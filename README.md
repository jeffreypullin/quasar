# quasar

quasar is a C++ software package for performing expression quantitative trait loci (eQTL) mapping.

Compared to other eQTL mapping software, quasar: 

* implements a much wider range of models, including both count distribution models and mixed models,
* implements the Cox-Reid adjusted profile likelihood for estimating the negative binomial dispersion,
* implements a novel trace-based approximation of the mixed-model score test variance, and 
* is generally faster.

As of version 2.0, quasar also:

* accepts single-cell resolution data as input,
* allows interaction-QTL testing,
* is generally faster than existing single-cell resolution approaches.

## Citation

If you use quasar in your research, please cite:

> Flexible and efficient count-distribution and mixed-model methods for eQTL mapping with quasar
> Jeffrey M. Pullin, Jarny Choi, Davis J. McCarthy, Chris Wallace  
> medRxiv 2025.07.17.25331702; doi: https://doi.org/10.1101/2025.07.17.25331702  

If you use quasar version 2.X+ in your research please also cite:

> TODO

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

The following invocation runs quasar with our recommended settings for cis-eQTL mapping of bulk or pseudobulk counts:

```
./quasar --plink plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --mode cis \
    --model p_glmm_id \
    --out p_glmm_id_fit
```

Here, `phenotype_data.bed` should contain bulk/pseudobulk gene expression counts. For single-cell resolution data, use:

```
./quasar --plink plink_prefix \
    --sc-pheno single-cell-pheno.tsv \
    --anno annotations.tsv \
    --cov single_cell_covariate_data.tsv \
    --mode cis \
    --model p_glmm_sc \
    --out single-cell-out
```

`--model p_glmm_sc` expects raw counts. For interaction-QTL mapping, pass one or more covariate names with `--interaction` (comma-separated for multiple covariates):

```
./quasar --plink plink_prefix \
    --sc-pheno single-cell-pheno.tsv \
    --anno annotations.tsv \
    --cov single_cell_covariate_data.tsv \
    --mode cis \
    --model p_glmm_sc \
    --interaction pseudotime \
    --out sc-int-out
```

For further information see the [Documentation](https://jeffreypullin.github.io/quasar/). 

## FAQs

quasar is named after the [astronomical object](https://en.wikipedia.org/wiki/Quasar): a short, memorable word that starts with Q, for QTL mapping. 

Please contact jp2045[at]cam.ac.uk for assistance running quasar or other questions.