# Change this to point to your version of MSFragger
export MSFRAGGER_PATH ?= ~/bin/MSFragger-3.1.1/MSFragger-3.1.1.jar

# Make sure TMPDIR is set:
export TMPDIR ?= /tmp

kim = data/pin/kim.pin.gz
percolator = scripts/percolator/make_figures.html
scope = scripts/scope/make_figures.html
rna = scripts/rna-xl/make_figures.html
benchmark = scripts/benchmark/make_figures.html

all: install ${kim} ${scope} ${rna} ${percolator} ${benchmark} wrapup

install:
	conda install -y -c conda-forge \
		mkl \
		tqdm \
		numpy \
		pandas \
		matplotlib \
		seaborn \
		scikit-learn \
		numba \
		mono \
		nbconvert \
		xgboost \
		wget \
	&& conda install -y -c bioconda \
		thermorawfileparser \
		percolator \
		triqler \
	&& pip install \
		git+git://github.com/wfondrie/mokapot \
	 	git+git://github.com/wfondrie/wispy \
		ppx

${kim}:
	mkdir -p data/pin && \
	wget -N -O data/pin/kim.pin.gz https://ndownloader.figshare.com/files/19068101

${benchmark}: ${kim}
	cd scripts/benchmark && \
	./cluster.sh && \
	jupyter nbconvert --to html make_figures.ipynb

${percolator}: ${scope}
	cd scripts/percolator && \
	python3 runall.py && \
	jupyter nbconvert --to html make_figures.ipynb

${scope}:
	cd scripts/scope && \
	python3 runall.py && \
	jupyter nbconvert --to html make_figures.ipynb

${rna}:
	cd scripts/rna-xl && \
	python3 runall.py && \
	jupyter nbconvert --to html make_figures.ipynb

wrapup:
	mkdir -p figures && \
	cp scripts/*/figures/*.png figures
