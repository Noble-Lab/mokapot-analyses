"""
Analyze the RNA data with a non-linear model
"""
import os
import sys
import copy
import shutil
import logging

import mokapot
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.inspection import permutation_importance

# local modules
sys.path.append(os.path.join("..", "..", "bin"))
import download
import search

# Constants and Setup ---------------------------------------------------------
MISSED_CLEAVAGES = 2
TOP_MATCH = 5
EXPERIMENT = "yeast"
FASTA = os.path.join("..", "..", "data", "fasta", "yeast_sp_2020-11.fasta")

logging.basicConfig(
    format=("[%(levelname)s]: %(message)s"),
    level=logging.INFO,
)


# Functions -------------------------------------------------------------------
def run_fragger(mzml_files, fasta, force_=False):
    """Run MSFragger"""
    os.makedirs("fragger-out", exist_ok=True)
    if isinstance(mzml_files, str):
        mzml_files = [mzml_files]

    base = [os.path.splitext(f)[0] for f in mzml_files]
    out_files = [
        os.path.join("fragger-out", os.path.split(f)[-1]) + ".pin" for f in base
    ]
    all_exist = all([os.path.isfile(f) for f in out_files])
    if all_exist and not force_:
        return out_files

    fragger_args = {
        "database_name": fasta,
        "precursor_mass_lower": "-150",
        "precursor_mass_upper": "1000",
        "precursor_mass_units": "0",
        "isotope_error": "1",
        "add_C_cysteine": "0",
        "decoy_prefix": "decoy_",
        "calibrate_mass": "2",
        "localize_delta_mass": "1",
        "clip_nTerm_M": "0",
        "allowed_missed_cleavage": MISSED_CLEAVAGES,
        "output_format": "tsv_pin",
        "output_report_topN": "5",
        "report_alternative_proteins": "1",
    }

    search.msfragger(mzml_files, **fragger_args)
    for frag in base:
        out = os.path.join("fragger-out", os.path.split(frag)[-1])
        for ext in (".tsv", ".pin"):
            shutil.move(frag + ext, out + ext)

    return out_files


def update_fragger(pin_files, top_match=5):
    """Update the pin files."""
    out_files = []
    os.makedirs("pin-out", exist_ok=True)
    for pin in pin_files:
        out_file = os.path.join("pin-out", os.path.basename(pin))
        out_files.append(out_file)
        if os.path.isfile(out_file):
            continue

        tsv = pin.replace(".pin", ".tsv")
        tsv_df = pd.read_csv(tsv, sep="\t")
        pin_df = mokapot.read_pin(pin, to_df=True).drop("ExpMass", axis=1)
        pin_df = pin_df.rename(columns={"ExpMass": "CalcMass"})
        cols = [
            "scannum",
            "hit_rank",
            "precursor_neutral_mass",
            "massdiff",
        ]
        tsv_df = tsv_df.loc[:, cols].rename(
            columns={
                "scannum": "ScanNr",
                "hit_rank": "rank",
                "precursor_neutral_mass": "ExpMass",
            }
        )

        pin_df = pd.merge(tsv_df, pin_df)
        pin_df["group"] = "unmodified"
        pin_df.loc[pin_df["abs_ppm"] > 50, "group"] = "modified"

        new_cols = ["group"] + list(pin_df.columns)[: pin_df.shape[1] - 1]
        pin_df = pin_df.loc[:, new_cols]
        pin_df["Peptide"] = (
            pin_df["Peptide"] + "[" + pin_df["massdiff"].round(2).astype(str) + "]"
        )

        pin_df = pin_df.loc[pin_df["rank"] <= top_match, :].drop(
            columns=["delta_hyperscore", "massdiff"]
        )

        pin_df.to_csv(out_file, index=False, sep="\t")

    return out_files


def run_mokapot(psms, model="linear", force_=False):
    """Run mokapot with a various or no model."""
    out_res = os.path.join("mokapot-out", f"{model}.modified.mokapot.psms.txt")
    if os.path.isfile(out_res) and not force_:
        return out_res

    if model == "fragger":
        np.random.seed(3)
        logging.info("======================")
        logging.info("====== Baseline ======")
        logging.info("======================")
        res = psms.assign_confidence()

    elif model == "linear":
        np.random.seed(3)
        logging.info("======================")
        logging.info("===== Linear SVM =====")
        logging.info("======================")
        mod = mokapot.PercolatorModel()
        res, mods = mokapot.brew(psms, mod)

    elif model == "xgb":
        np.random.seed(3)
        logging.info("======================")
        logging.info("====== XGBoost =======")
        logging.info("======================")
        grid = {
            "scale_pos_weight": np.logspace(0, 2, 3),
            "max_depth": [1, 3, 6],
            "min_child_weight": [1, 10, 100],
            "gamma": [0, 1, 10],
        }
        xgb_mod = GridSearchCV(
            XGBClassifier(),
            param_grid=grid,
            n_jobs=1,
            cv=3,
            scoring="roc_auc",
        )
        mod = mokapot.Model(xgb_mod)
        res, mods = mokapot.brew(psms, mod)

    else:
        raise ValueError("Must be a typo in model.")

    res.to_txt("mokapot-out", f"{model}")
    if model != "fragger":
        for i, mod in enumerate(mods):
            mod.save(os.path.join("mokapot-out", f"{model}.model{i}.pkl"))

    return out_res


def parse_results(res_files):
    """Parse the result files"""
    psms = {}
    peptides = {}
    proteins = {}
    models = {}
    for psm_file in res_files:
        label = os.path.split(psm_file)[-1].split(".")[0]
        psms[label] = pd.read_csv(psm_file, sep="\t")
        pep_file = psm_file.replace("psms", "peptides")
        prot_file = psm_file.replace("psms", "proteins")
        peptides[label] = pd.read_csv(pep_file, sep="\t")
        proteins[label] = pd.read_csv(prot_file, sep="\t")
        if label != "fragger" and label != "mod":
            model_file = psm_file.split(".")[0] + ".model1.pkl"
            models[label] = mokapot.load_model(model_file)

    results = (psms, peptides, proteins, models)
    logging.info("")
    logging.info("=== Results ===")
    logging.info("%s", res_files)
    for level, res in zip(("PSMs", "peptides", "proteins"), results[:3]):
        logging.info(level)
        for label, df in res.items():
            num_passing = (df["mokapot q-value"] <= 0.01).sum()
            logging.info("\t%s:  %i", label, num_passing)

    return results


def calc_importance(dset, models, force_=False):
    """Do all the feature importance calculations"""
    base_dset = copy.copy(dset)
    modified = base_dset._data["group"] == "modified"
    base_dset._data = base_dset._data.loc[modified, :]

    out_dir = "featimp-out"
    os.makedirs(out_dir, exist_ok=True)
    imp_out = os.path.join(out_dir, "importance.txt")

    if not os.path.isfile(imp_out) or force_:
        imp = pd.concat([feature_importance(dset, m, l) for l, m in models.items()])
        imp.to_csv(imp_out, sep="\t", index=False)

    return imp_out


def feature_importance(dset, model, label):
    """Get the permutation feature importance"""
    feat = model.scaler.transform(dset.features.values)
    labels = dset.targets
    logging.info("Calculating importance for %s...", label)
    imp = permutation_importance(
        model.estimator,
        feat,
        labels,
        scoring="roc_auc",
        n_jobs=-1,
    )

    imp_df = pd.DataFrame(imp.importances, index=dset.features.columns)
    imp_df = (
        imp_df.reset_index(drop=False)
        .melt("index", var_name="rep", value_name="importance")
        .rename(columns={"index": "feature"})
    )
    imp_df["model"] = label
    return imp_df


# MAIN ------------------------------------------------------------------------
def main():
    """The main function"""
    np.random.seed(1)
    out_dirs = ["figures", "mokapot-out"]
    [os.makedirs(d, exist_ok=True) for d in out_dirs]

    models = ["fragger", "linear", "xgb"]

    mzml_files = download.rnaxl(EXPERIMENT)
    td_fasta = "yeast_target-decoy.fasta"
    if not os.path.isfile(td_fasta):
        mokapot.make_decoys(FASTA, td_fasta)

    logging.info("Performing Searches...")
    search_res = run_fragger(mzml_files, td_fasta)
    pins = update_fragger(search_res, top_match=TOP_MATCH)

    logging.info("Reading Search Results...")
    psms = mokapot.read_pin(pins, group_column="group")
    logging.info("\n%s", psms._data.groupby("group")["Label"].value_counts())

    psms.add_proteins(td_fasta, missed_cleavages=MISSED_CLEAVAGES)

    logging.info("Training models...")
    res_files = [run_mokapot(psms, m) for m in models]
    _, _, _, trained_mods = parse_results(res_files)
    calc_importance(psms, trained_mods)


if __name__ == "__main__":
    main()
