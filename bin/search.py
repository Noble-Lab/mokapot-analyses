"""
Python wrappers for various search engines and utilities.
"""
import os
import shutil
import subprocess
from tempfile import TemporaryDirectory


def msfragger(ms_files, max_mem="32G", jar_path=None, **kwargs):
    """
    Conduct a search using MSFragger.

    Parameters
    ----------
    ms_files : str or list of str
        The MS data files to search.
    max_mem : str
        Passed to java -Xmx
    jar_path : str
        The path to the MSFragger jar file.
    **kwargs : dict
        Arguments passed to MSFragger.
    """
    if jar_path is None:
        jar_path = os.getenv("MSFRAGGER_PATH")
        if not jar_path:
            jar_path = "~/bin/MSFragger-3.1.1/MSFragger-3.1.1.jar"

    jar_path = os.path.abspath(os.path.expanduser(jar_path))
    fragger_args = [f"--{k} {v}" for k, v in kwargs.items()]
    cmd = ["java", f"-Xmx{max_mem}", "-jar", jar_path] + fragger_args

    if isinstance(ms_files, str):
        ms_files = [ms_files]

    subprocess.run(" ".join(cmd + ms_files), check=True, shell=True)


def tide(ms_files, fasta, **kwargs):
    """
    Conduct a search using crux tide.

    Parameters
    ----------
    ms_files : str or list of str
        The MS data files to search.
    fasta : str
        The protein database to search.
    **kwargs : dict
        Arguments passed to MSFragger.
    """
    tide_args = [f"--{k} {v}" for k, v in kwargs.items()]
    cmd = ["crux", "tide-search"] + tide_args

    if isinstance(ms_files, str):
        ms_files = [ms_files]

    subprocess.run(" ".join(cmd + ms_files + [fasta]), check=True, shell=True)


def make_decoys(fasta, outfile, concat=True, **kwargs):
    """
    Run crux generate-peptides to create concatenated target-decoy database.

    Parameters
    ----------
    fasta : str
        The location of a fasta file for which to generate decoys.
    outfile : str
        The name of the resulting fasta file.
    concat : bool
        Return a concatenated database or just the decoys?
    **kwargs : dict
        Arguments passed to 'crux generate-peptides'
    """
    with TemporaryDirectory() as tmp:
        crux_args = [f"--{k} {v}" for k, v in kwargs.items()]
        cmd = ["crux", "generate-peptides", "--output-dir", tmp, fasta]
        cmd += crux_args

        subprocess.run(" ".join(cmd), check=True, shell=True)
        decoys = os.path.join(tmp, "generate-peptides.proteins.decoy.txt")

        with open(outfile, "wb") as fout:
            if concat:
                with open(fasta, "rb") as target_file:
                    shutil.copyfileobj(target_file, fout)

            with open(decoys, "rb") as decoy_file:
                shutil.copyfileobj(decoy_file, fout)

    return outfile


def tide_index(fasta_file, name, **kwargs):
    """
    Run crux tide-index.

    Parameters
    ----------
    fasta_file : str
        The protein database
    name : str
        The name of the index directory
    **kwargs : dict
        Arguments to pass to tide-index
    """
    index_args = [f"--{k} {v}" for k, v in kwargs.items()]
    cmd = ["crux", "tide-index"] + index_args
    subprocess.run(" ".join(cmd + [fasta_file, name]), check=True, shell=True)

    return name
