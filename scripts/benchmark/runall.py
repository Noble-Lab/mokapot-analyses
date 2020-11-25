"""
Benchmark Percolator and mokapot
"""
import os
import logging
import subprocess

import tqdm
import numpy as np

# Setup -----------------------------------------------------------------------
REPS = 3
PIN = os.path.join(os.getenv("TMPDIR"), "test.pin")
LENGTH = 23330312

# Functions -------------------------------------------------------------------
def get_results(log_file):
    """Extract the wall clock time and maximum RSS from a GNU time log file"""
    with open(log_file) as log:
        for line in log:
            if "Elapsed (wall clock) time" in line:
                time = line.split(" ")[-1].split(":")
                if len(time) == 3:
                    time = (
                        float(time[0]) * 60 ** 2 + float(time[1]) * 60 + float(time[2])
                    )
                else:
                    time = float(time[0]) * 60 + float(time[1])
            if "Maximum resident set size" in line:
                mem = float(line.split(" ")[-1]) / 1000

    return (time, mem)


def sample_psms(num, pin_file, out_file, total):
    """Sample PSMs"""
    if os.path.isfile(out_file):
        return out_file

    psm_set = set(np.random.choice(np.arange(total), num, replace=False))

    logging.info(f"Sampling {num} PSMs...")
    # Read through PIN, appending selected rows to the list.
    samp = []
    with open(pin_file, "r") as pin:
        samp.append(pin.readline())

        for idx, psm in tqdm.tqdm(enumerate(pin), total=total):
            if idx in psm_set:
                samp.append(psm)

    with open(out_file, "w+") as out_pin:
        out_pin.writelines(samp)

    return out_file


def benchmark(pin, mokapot=True, rep=None):
    """Benchmark a command"""
    prefix = ["/usr/bin/time", "-v"]
    suffix = [pin, "2>"]

    if rep is None:
        rep = ""
    else:
        rep = f"_{rep}"

    fileroot = os.path.split(pin)[-1].replace(".pin", f"{rep}.log.txt")
    out_dir = os.getenv("TMPDIR")

    if mokapot:
        out_file = [f"logs/mokapot_{fileroot}"]
        cmd = [f"mokapot -d {out_dir}"]
    else:
        out_file = [f"logs/percolator_{fileroot}"]
        cmd = [
            f"percolator -Y --results-psms {out_dir}/percolator.psms.txt "
            f"--results-peptides {out_dir}/percolator.peptides.txt"
        ]

    if not os.path.isfile(out_file[0]):
        logging.info(f"Running {cmd[0]} on {pin}")
        subprocess.run(
            " ".join(prefix + cmd + suffix + out_file), shell=True, check=True
        )
    else:
        logging.info(f"{out_file[0]} exist. Skipping...")

    return out_file[0]


# MAIN ------------------------------------------------------------------------
def main():
    """The main function"""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    np.random.seed(42)

    pin_dir = os.path.join(os.getenv("TMPDIR"), "pin-out")
    os.makedirs(pin_dir, exist_ok=True)
    nums = list(np.logspace(4, 7, 7)) + [LENGTH]
    pins = [
        sample_psms(int(n), PIN, f"{pin_dir}/sampled_{int(n)}.pin", LENGTH)
        for n in nums
    ]

    os.makedirs("logs", exist_ok=True)
    perc_benchmark = []
    mp_benchmark = []
    for r in range(REPS):
        mp_benchmark += [benchmark(p, True, r) for p in pins]
        perc_benchmark += [benchmark(p, False, r) for p in pins]

    logging.info("DONE!")


if __name__ == "__main__":
    main()
