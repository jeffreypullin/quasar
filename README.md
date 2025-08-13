# quasar

quasar is a C++ software package for performing expression quantitative trait loci (eQTL) mapping.

Compared to other eQTL mapping software, quasar: 

* is ~5-45 times faster when run on CPUs,
* implements a much wider range of models, including both count distribution models and mixed models,
* implements the Cox-Reid adjusted profile likelihood for estimating the negative bimomial dispersion, and
* implements a novel trace-based approximation of the mixed-model score test variance.

## Citation

If you use quasar in your research, please cite

> Flexible and efficient count-distribution and mixed-model methods for eQTL mapping with quasar
> Jeffrey M. Pullin, Chris Wallace
> medRxiv 2025.07.17.25331702; doi: https://doi.org/10.1101/2025.07.17.25331702

## Installation

To install quasar run: 

```sh
git clone https://github.com/jeffreypullin/quasar.git
cd quasar
mkdir -p build
cd build
cmake ..
make
```

This process generates the binary in the `build/` subdirectory of quasar. To verify the installation has completed succesfully, run,

```
./quasar --version
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
    --out nb_fit
```

Here, `phenotype_data.bed` should contain bulk/pseudobulk gene counts.

For further information see the [Documentation](https://jeffreypullin.github.io/quasar/). 

Please contact jp2045[at]cam.ac.uk for assistance running quasar or other questions.
