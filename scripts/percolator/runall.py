"""
Sample PSMs from the Kim et al then run Percolator and Mokapot
"""
import os
import logging
import subprocess

import mokapot
import pandas as pd

# Setup -----------------------------------------------------------------------
PIN = os.path.join(
    "..", "scope", "pin-out", "190222S_LCA9_X_FP94_col22.make-pin.pin"
)
FASTA = os.path.join(
    "..", "..", "data", "fasta", "human_swissprot_2019-09.fasta"
)


# Functions -------------------------------------------------------------------
def run_mokapot(pin, fasta):
    """Run mokapot"""
    out_dir = "mokapot-out"
    os.makedirs(out_dir, exist_ok=True)

    out_base = os.path.join(out_dir, "mokapot.{level}.txt")
    out_files = [out_base.format(level=l)
                 for l in ["psms", "peptides", "proteins"]]

    if not all([os.path.isfile(f) for f in out_files]):
        cmd = [
            "mokapot",
            "--dest_dir", out_dir,
            "--proteins", fasta,
            pin,
        ]

        subprocess.run(cmd, check=True)

    return out_files


def run_percolator(pin, fasta):
    """Run Percolator"""
    out_dir = "perc-out"
    os.makedirs(out_dir, exist_ok=True)

    out_base = os.path.join(out_dir, "percolator.{level}.txt")
    out_files = [out_base.format(level=l)
                 for l in ["psms", "peptides", "proteins"]]

    if not all([os.path.isfile(f) for f in out_files]):
        cmd = [
            "percolator",
            "--results-psms", out_files[0],
            "--results-peptides", out_files[1],
            "--results-proteins", out_files[2],
            "--picked-protein", fasta,
            "--protein-decoy-pattern", "decoy_",
            "--post-processing-tdc",
            "--seed", "42",
            pin
        ]

        subprocess.run(cmd, check=True)

    return out_files


def parse_results(res_files):
    """Parse and save the combined results"""
    out_dir = "combined-out"
    os.makedirs(out_dir, exist_ok=True)
    for level in ["psms", "peptides", "proteins"]:
        moka_file = [f for f in res_files if f"mokapot.{level}" in f][0]
        perc_file = [f for f in res_files if f"percolator.{level}" in f][0]

        perc_res = mokapot.read_percolator(perc_file)
        moka_res = pd.read_table(moka_file)

        if level == "proteins":
            moka_res["ProteinId"] = (moka_res["mokapot protein group"]
                                     .str.split(", ", expand=True)[0])
            perc_res = perc_res.rename(
                columns={
                    "q-value": "percolator q-value",
                    "posterior_error_prob": "percolator PEP",
                    "ProteinId": "mokapot protein group",
                }
            )
        else:
            perc_res = perc_res.rename(
                columns={
                    "PSMId": "SpecId",
                    "peptide": "Peptide",
                    "proteinIds": "Proteins",
                    "score": "percolator score",
                    "q-value": "percolator q-value",
                    "posterior_error_prob": "percolator PEP",
                }
            )

        merged = pd.merge(moka_res, perc_res)
        moka_sum = (moka_res["mokapot q-value"] <= 0.01).sum()
        perc_sum = (perc_res["percolator q-value"] <= 0.01).sum()

        logging.info("------------------------------")
        logging.info("%s", level)
        logging.info("  - Mokapot: %i (%i)", len(moka_res), moka_sum)
        logging.info("  - Percolator: %i (%i)", len(perc_res), perc_sum)
        logging.info("  - Merged: %i", len(merged))

        merged.to_csv(os.path.join(out_dir, f"{level}.txt"), sep="\t",
                      index=False)


# Main ------------------------------------------------------------------------
def main():
    """The main function"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    perc_res = run_percolator(PIN, FASTA)
    moka_res = run_mokapot(PIN, FASTA)
    parse_results(perc_res + moka_res)


if __name__ == "__main__":
    main()
