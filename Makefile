# Change this to point to your version of MSFragger
export MSFRAGGER_PATH ?= ~/bin/MSFragger-3.1.1/MSFragger-3.1.1.jar

# Make sure TMPDIR is set:
export TMPDIR ?= /tmp

kim = data/pin/kim.pin.gz
percolator = scripts/percolator
scope = scripts/scope
rna = scripts/rna-xl
benchmark = scripts/benchmark

all: install ${kim} ${scope}/make_figures.html ${rna}/make_figures.html \
	${percolator}/make_figures.html ${benchmark}/make_figures.html wrapup

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


${benchmark}/make_figures.html: ${kim} ${benchmark}/cluster.sh \
	${benchmark}/runall.py ${benchmark}/make_figures.ipynb

	cd scripts/benchmark && \
	./cluster.sh && \
	jupyter nbconvert --to html --execute make_figures.ipynb


${percolator}/make_figures.html: ${percolator}/runall.py \
	${percolator}/make_figures.ipynb \
	${scope}/pin-out/190222S_LCA9_X_FP94_col22.make-pin.pin

	cd scripts/percolator && \
	python3 runall.py && \
	jupyter nbconvert --to html --execute make_figures.ipynb


${scope}/make_figures.html 	${scope}/pin-out/190222S_LCA9_X_FP94_col22.make-pin.pin: \
	${scope}/runall.py ${scope}/make_figures.ipynb

	cd scripts/scope && \
	python3 runall.py && \
	jupyter nbconvert --to html --execute make_figures.ipynb


${rna}/make_figures.html: ${rna}/runall.py ${rna}/make_figures.ipynb
	cd scripts/rna-xl && \
	python3 runall.py && \
	jupyter nbconvert --to html --execute make_figures.ipynb


wrapup:
	mkdir -p figures && \
	cp scripts/*/figures/*.png figures
