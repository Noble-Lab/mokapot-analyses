#!/usr/bin/bash
#$ -cwd
#$ -l h_rt=36:0:0
#$ -l mfree=8G
#$ -l disk_free=50G
#$ -l gpgpu=FALSE
#$ -pe serial 12
#$ -o benchmark.stdout.txt
#$ -e benchmark.stderr.txt
#$ -N benchmark
set -e
echo "Start - `date`"

export GOMP_CPU_AFFINITY="${SGE_BINDING}"

if [ ! -f "${TMPDIR}/test.pin" ] ; then
    cp -v ../../data/pin/kim.pin.gz ${TMPDIR}/test.pin.gz
    gunzip ${TMPDIR}/test.pin.gz
fi

python runall.py

echo "Finish - `date`"
