# Code for evaluating mokapot
This repository contains the code for reproducing the results from "mokapot:
Fast and flexible semi-supervised learning for peptide detection**. 

## Reproducing the manuscript
The provided can fully reproduce the figures and analyses presented in the
manuscript, provided that the necessary software are installed and data are
present. Additionally, some analyses (such as the benchmarking experiments** will
provide different results depending on the hardware they are run on.

### Requirements  
**Operating System:** Our code was written for CentOS7 Linux machines, but
should be compatible with MacOS and other Linux distributions as well. The code
is unlikely to work on
Windows.  

**Hardware:** To most accurately reproduce our results, a 12-core machine with a
minimum of 32 Gb of memory should be used.

**Installed Software:** The analysis scripts are written for Python 3.7+. Many
of the software tools are installed automatically when executing the analysis,
however some need to be installed beforehand:  

- [Anaconda Python 3.7+](https://www.anaconda.com/products/individual) (Our
  analysis scripts install dependencies using the conda package manager)  

- [The Crux Mass Spectrometry Toolkit](http://crux.ms/) (Use the "I agree with
  to the licensing terms, download the most recent build of Crux" button)  
  
- [MSFragger 3.1.1](https://github.com/Nesvilab/MSFragger/wiki) (Note that
  MSFragger also requires the [Java 8 Runtime Environment](https://www.java.com))
  
You can then check that these have been successfully installed and configured
with the following commands (the example output is from my machine, but yours
may be slightly different).

Verify Python Version and configured:
```bash
$ python3 --version
Python 3.8.6
```

Verify conda is installed and configured:
```bash
$ conda --version
conda 4.9.2
```

Verify Crux is installed and configured:
```bash
$ crux version
INFO: Beginning version.
====================
Crux version 3.2-9d35092f
====================
Proteowizard version 3.0.20213
====================
Percolator version 3.05.nightly-1-e16f49a-dirty, Build Date Jul 30 2020 22:14:27
Copyright (c) 2006-9 University of Washington. All rights reserved.
Written by Lukas Käll (lukall@u.washington.edu) in the
Department of Genome Sciences at the University of Washington.
====================
Comet version 2019.0X rev. X
====================
Boost version 1_67
====================
INFO: Elapsed time: 0.000323 s
INFO: Finished crux version.
INFO: Return Code:0
```

Verify that MSFragger is installed and configured (your path may be different):
```
$java -jar ~/bin/MSFragger-3.1.1/MSFragger-3.1.1.jar --version
MSFragger version MSFragger-3.1.1
Batmass-IO version 1.19.5
timsdata library version timsdata-2-7-0
(c) University of Michigan
RawFileReader reading tool. Copyright (c) 2016 by Thermo Fisher Scientific, Inc. All rights reserved.
System OS: Linux, Architecture: amd64
Java Info: 14.0.2, OpenJDK 64-Bit Server VM, Red Hat, Inc.
```

Finally, you'll need to specify the path to the MSFragger jar file:
```
export MSFRAGGER=~/bin/MSFragger-3.1.1/MSFragger-3.1.1.jar
```

## Running the Analyses
Running the analyses is easy, but will potentially take days. First, we highly
recommend creating a new conda environment for the analyses. The scripts will
install packages into this environment:

```bash
$ conda create -n mokapot && conda activate mokapot
```

Then the analyses are run using the GNU make:
```bash
$ make
```

### Results
Once complete, all of the figures will be present in the `figures` directory.

## Questions?
If you have problems or questions, feel free to ask Will Fondrie (wfondrie@uw.edu).




