# quasar

quasar is a C++ software package for performing expression quantitative trait loci (eQTL) mapping.

Compared to other eQTL mapping software, quasar: 

* is ~6-400 times faster when run on CPUs,
* implements a much wider range of models, including both count distribution models and mixed models,
* implements the Cox-Reid adjusted profile likelihood for estimating the negative bimomial dispersion, and
* implements a novel trace-based approximation of the mixed-model score test variance.

## Citation

If you use quasar in your research, please cite

> TODO

## Installation

To install quasar run: 

```sh
git clone https://github.com/jeffreypullin/quasar.git
cd quasar
cmake 
make
```

## Quickstart

The following invocation runs quasar with our recommended settings for cis-eQTL mapping:

```
./quasar --plink plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --mode cis \
    --model nb_glm \
    --use-apl \
    --out quasar_out
```

For further information see the [Documenation](https://jeffreypullin.github.io/quasar/).

Please contact jp2045[at]cam.ac.uk for assistance running quasar for other questions.
