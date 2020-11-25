"""
Download the required proteomics data and fasta files.
"""
import os
import sys
from typing import List
import ppx
import logging
import subprocess

DATA_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
)


def rnaxl(experiment="human") -> List[str]:
    """
    Download the mzML files from PXD000513

    Parameters
    ----------
    experiment : str
        The experiment to use.

    Returns
    -------
    A list of mzML files
    """
    out_dir = os.path.join(DATA_DIR, "rnaxl", "mzML")
    os.makedirs(out_dir, exist_ok=True)

    dataset = ppx.PXDataset("PXD000513")
    all_files = [f for f in dataset.list_files() if f.endswith(".mzML")]

    if experiment == "human":
        mzml_files = [f for f in all_files if "XL_human_RBPs" in f]
    elif experiment == "yeast":
        mzml_files = [f for f in all_files if "XL_yeast_RBPs_Ex" in f
                      and "UV" in f]
    elif experiment == "invivo":
        mzml_files = [f for f in all_files if "XL_yeast_RBPs_invivo_4SU_Ex" in f]

    logging.info("Downloading...")
    downloaded = dataset.download(mzml_files, dest_dir=out_dir)
    return downloaded


def scope2() -> List[str]:
    """
    Download and convert the raw files from PXD001468.

    If the files already exist they will not be re-downloaded and
    converted again.

    Returns
    -------
    A list of mzML files.
    """
    raw_dir = os.path.join(DATA_DIR, "scope2", "raw")
    mzml_dir = os.path.join(DATA_DIR, "scope2", "mzML")

    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mzml_dir, exist_ok=True)

    skip = ["190228S_LCA9_X_FP94BF.raw",
            "190321S_LCA10_X_FP97BE.raw"]

    dataset = ppx.MSVDataset("MSV000083945")
    raw_files = [f for f in dataset.list_files("raw/scope2_raw")
                 if f.endswith(".raw") and f not in skip]
    out_files = [os.path.join(mzml_dir, f.replace(".raw", ".mzML.gz"))
                 for f in raw_files]

    raw_files = ["raw/scope2_raw/" + f for f in raw_files]

    print([f for f in out_files if not os.path.isfile(f)])

    if not all([os.path.isfile(f) for f in out_files]):
        logging.info("Downloading...")
        downloaded = dataset.download(raw_files, dest_dir=raw_dir)
        logging.info("Converting...")
        for raw_file in downloaded:
            _ = convert(raw_file, mzml_dir, gzip=True)

    return out_files


def convert(raw, dest_dir=".", gzip=False):
    """
    Convert raw files to mzML format

    Parameters
    ----------
    raw : str
        The raw file to convert.

    dest_dir : str
        The destination directory for output.

    Returns
    -------
    str
        The mzML file.
    """
    ex_path = os.path.join(os.path.dirname(sys.executable), "ThermoRawFileParser.exe")

    ext = ".mzML"
    cmd = [
        "mono",
        ex_path,
        "-i=" + raw,
        "-o=" + dest_dir,
        "-f=2",
        "-m=1"
    ]

    if gzip:
        cmd += ["-g"]
        ext += ".gz"

    out_file = os.path.join(
        dest_dir, os.path.split(raw)[-1].replace(".raw", ext)
    )
    if not os.path.isfile(out_file):
        subprocess.run(cmd, check=True)

    assert os.path.isfile(out_file) # verify that it is actually there!
    return out_file


def _download_and_convert(dataset, out_dir, raw_files=None):
    """Do the download and conversion"""
    raw_dir = os.path.join(out_dir, "raw")
    mzml_dir = os.path.join(out_dir, "mzML")

    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mzml_dir, exist_ok=True)

    if raw_files is None:
        raw_files = [f for f in dataset.list_files() if f.endswith(".raw")]

    out_files = [os.path.join(mzml_dir, f.replace(".raw", ".mzML.gz"))
                 for f in raw_files]

    out_files = [os.path.join(mzml_dir, f.replace(".raw", ".mzML.gz"))
                 for f in raw_files]

    if not all([os.path.isfile(f) for f in out_files]):
        logging.info("Downloading...")
        downloaded = dataset.download(raw_files, dest_dir=raw_dir)
        logging.info("Converting...")
        for raw_file in downloaded:
            _ = convert(raw_file, mzml_dir, gzip=True)

    return out_files
