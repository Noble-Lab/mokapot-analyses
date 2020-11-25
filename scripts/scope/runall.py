"""
This script reanlyzes the SCoPE2 dataset from the Slavov lab.

The idea is that training a joint model should work similarly to using
a static model when the number of datasets is large.
"""
import os
import sys
import logging
import subprocess

import mokapot
import numpy as np
import pandas as pd

sys.path.append(os.path.join("..", "..", "bin"))
import search
import download

# Setup -----------------------------------------------------------------------
np.random.seed(42)
FASTA = os.path.join("..", "..", "data", "fasta", "human_swissprot_2019-09.fasta")
MISSED_CLEAVAGES = 2


# Functions ------------------------------------------------------------------
def make_index(fasta_file, name="human.index"):
    """Create a tide index"""
    params = {
        "mods-spec": "C+0,K+229.162932,3M+15.99492,3N+0.9840155848",
        "nterm-peptide-mods-spec": "X+229.162932",
        "nterm-protein-mods-spec": "1X+42.01056",
        "enzyme": "trypsin/p",
        "missed-cleavages": str(MISSED_CLEAVAGES),
        "output-dir": "index-out",
    }

    if not os.path.isdir(name):
        search.tide_index(fasta_file, name, **params)

    return name


def tide2pin(target, name):
    """Convert tide results to a pin file."""
    out_file = f"pin-out/{name}.make-pin.pin"

    if not os.path.isfile(out_file):
        cmd = [
            "crux",
            "make-pin",
            "--max-charge-feature",
            "5",
            "--output-dir",
            "pin-out",
            "--top-match",
            "5",
            "--fileroot",
            name,
            target,
        ]

        subprocess.run(cmd, check=True)

    return out_file


def run_tide(mzml_files, name, index):
    """Perform a tide search"""
    if isinstance(mzml_files, str):
        mzml_files = [mzml_files]

    out_file = f"tide-out/{name}.tide-search.target.txt"

    params = {
        "remove-precursor-peak": "T",
        "auto-precursor-window": "warn",
        "isotope-error": "1",
        "precursor-window": "50",
        "precursor-window-type": "ppm",
        "fragment-tolerance": "0.02",
        "mz-bin-width": "1.0005079",
        "compute-sp": "T",
        "exact-p-value": "T",
        "score-function": "both",
        "fileroot": name,
        "output-dir": "tide-out",
        "concat": "F",
    }

    if not os.path.isfile(out_file):
        search.tide(mzml_files, index, **params)

    return out_file


def load_pins(pin_dir="pin-out", train=False, fasta=None):
    """Load the pin files."""
    if train:
        pin_files = ["qc.make-pin.pin"]
    else:
        pin_files = [
            p for p in os.listdir(pin_dir) if p.endswith("pin") and "qc" not in p
        ]

    small_files = [
        "190222S_LCA9_X_FP94BD.make-pin.pin",
        "190228S_LCA9_X_FP94BE.make-pin.pin",
        "190222S_LCA9_X_FP94AO.make-pin.pin",
    ]
    pin_files = [p for p in pin_files if p not in small_files]

    pin_files = [os.path.join(pin_dir, p) for p in pin_files]
    psms = [mokapot.read_pin(p) for p in pin_files]

    if fasta is not None:
        prots = mokapot.FastaProteins(fasta, missed_cleavages=MISSED_CLEAVAGES)
        for dat in psms:
            dat.add_proteins(prots)

    return psms, pin_files


def run_mokapot(model_type, fasta):
    """Run mokapot with a certain type of model"""
    out_dir = "mokapot-out"
    out_files = [
        os.path.join(out_dir, model_type + l + ".txt.gz")
        for l in (".psms", ".peptides", ".proteins")
    ]

    if all([os.path.isfile(f) for f in out_files]):
        return tuple(pd.read_csv(f, sep="\t") for f in out_files)

    os.makedirs(out_dir, exist_ok=True)
    if model_type == "static":
        train, _ = load_pins(train=True)
        model = mokapot.PercolatorModel()
        model.fit(train[0])
        del train
        test, pins = load_pins(fasta=fasta)
        results = [d.assign_confidence(model.predict(d)) for d in test]
        results = aggregate_results(results, pins)

    elif model_type == "independent":
        test, pins = load_pins(fasta=fasta)
        model = mokapot.PercolatorModel(override=True)
        results = []
        for d in test:
            try:
                res = mokapot.brew(d, model)[0]
            except (ValueError, RuntimeError) as e:
                res = None
                logging.warning("Brew failed for %s.", d)
                logging.warning("\t- Caught: %s", e)

            results.append(res)

        results = aggregate_results(results, pins)

    elif model_type == "joint":
        test, pins = load_pins(fasta=fasta)
        model = mokapot.PercolatorModel(override=True)
        results = aggregate_results(mokapot.brew(test, model)[0], pins)

    elif model_type == "tide":
        test, pins = load_pins(fasta=fasta)
        results = [d.assign_confidence() for d in test]
        results = aggregate_results(results, pins)

    else:
        raise ValueError("Unrecognized model_type")

    _ = [r.to_csv(f, sep="\t", index=False) for r, f in zip(results, out_files)]

    return results


def aggregate_results(results, pins):
    """Aggregate the results"""
    psms = []
    peps = []
    prots = []
    for res, pin in zip(results, pins):
        if res is None:
            continue

        res.psms["pin_file"] = pin
        res.peptides["pin_file"] = pin
        res.proteins["pin_file"] = pin
        psms.append(res.psms)
        peps.append(res.peptides)
        prots.append(res.proteins)

    return pd.concat(psms), pd.concat(peps), pd.concat(prots)


# MAIN ------------------------------------------------------------------------
def main():
    """Run the analyses"""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logging.info("##### Getting Files #####")
    mzml = download.scope2()

    logging.info("##### Making Index #####")
    index = make_index(FASTA)

    logging.info("##### Performing QC searches #####")
    qc_name = "qc"
    qc_files = [f for f in mzml if "_QC_" in f and f.endswith("mzML.gz")]
    qc_tide = run_tide(qc_files, qc_name, index)
    qc_pin = tide2pin(qc_tide, name=qc_name)

    logging.info("##### Performing main searches #####")
    x_files = [f for f in mzml if "_X_" in f and f.endswith("mzML.gz")]
    x_names = [os.path.split(f.replace(".mzML.gz", ""))[-1] for f in x_files]
    x_tide = [run_tide(f, n, index) for f, n in zip(x_files, x_names)]
    x_pins = [tide2pin(f, name=n) for f, n in zip(x_tide, x_names)]

    logging.info("##### Static Model #####")
    static = run_mokapot("static", FASTA)

    logging.info("##### Independent Models #####")
    independent = run_mokapot("independent", FASTA)

    logging.info("##### Joint Models #####")
    joint = run_mokapot("joint", FASTA)

    logging.info("##### No model #####")
    tide = run_mokapot("tide", FASTA)

    logging.info("##### DONE! #####")


if __name__ == "__main__":
    main()
