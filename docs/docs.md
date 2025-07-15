# Usage

To run **quasar**, use the command `./quasar` on the command line after installation. Flags and options specify how quasar will run. 

To list all the possible options and see quasar's help you can run: 

```
./quasar --help
```

To see the version of quasar you are using run: 

```
./quasar --version
```

## Quickstart

We recommend running quasar to perform cis-eQTL mapping using the negative binomial model and adjusted profile likelihood estimation of the negative binomial dispersion parameter. To run quasar in this mode use the command: 

```
./quasar --plink plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --mode cis \
    --model nb_glm \
    --use-apl \
    --out quasar_out
```

## QTL mapping modes

The quasar software can be run in three modes: `cis`, `trans`, `gwas`. These modes specify which varaints are tested for asscociation with a particular feature. 

In mode `cis` variants within +- the window size of the gene (see phenotype data format for details). By default the window size is set to 1Mb but can be specified using the `--window_size` flag. We refer to this set of variants as the cis-window for that feature. 

In mode `trans` all variants except those in the cis window are tested for assoication. In mode `gwas` all variants are are tested for association. 

Note that during the development of quasar the `cis` mode was tested more extensively than the `trans` and `gwas` modes and that there are methodological issues with trans-eQTL mapping due to reads mapping to multiple locations causing false-positive associations. 

## Statistical models

The quasar software package supports a wide range of statistical models used to resiudalise the expression values. The supported models are:

* `lm`: linear model
* `nb_glm`: negative binomial GLM 
* `lmm`: linear mixed model
* `p_glmm`: Poisson generalised linear mixed model (GLMM)
* `p_glm`: Poisson GLM (**not recommended** due to producing a very high rate of false positives)
* `nb_glmm`: negative binomial GLMM (**not generally recommended** due to producing highly similar results to the Poisson GLMM while being slower, can be used if there is known to be high relatedness between samples)

When the model is a mixed model i.e. is specified to be any of `lmm`, `p_glmm`, `nb_glmm` the --grm flag (see below) must be used to specify a genetic relatedness matrix used in the covarariance matrix of the random effects. 

When the `nb_glm` or `nb_glmm` flags are specified the --use-apl flag can be specified to use the Cox-Reid adjusted profile likelihood (APL) when estimating the negative binomial dispersion paraemter. Use of the APL is recommended as it reduces the number of false-positives but is slightly slower than standard maximum likelihood estimation.

## Data formats

### Genotype data

--plink/-p

The genotype data should be in plink2 binary format with .bed/.bim/.fam files named plink_prefix.bed/.bim/.fam

The .bed/.bim/.fam files can be generated from vcf using the following command

```sh
plink2 \
    --output-chr chrX \
    --vcf ${plink_prefix}.vcf.gz \
    --out ${plink_prefix}
```

If using --make-bed with PLINK 1.9 or earlier, add the --keep-allele-order flag.

### Phenotype data

--bed/-p

The phenotype data is a tab-seperated file with where rows are features and the first four
columns give feature inforamtion and the rest are sample ids are the sample ids. For example, 

```
#chr      start         end      phenotype_id  sample_1   sample_2   sample_3 ...
   1  113871759   113813811   ENSG00000134242      39.1      435.4      435.8 ...
 ...
```

The start and end values are used to specify the centre of the cis-window. To specify the gene TSS as the centre of the window, set TSS = start, end = start + 1, so that the cis-window is [TSS - window, TSS + window + 1] or alternatively set the start and end values to the start and end of the gene so that the cis-window is [start - window, end + window].

### Covariate data

--cov/-c

The covariate data is a tab-separated file with rows as samples and first column `sample_id` and other columns containing the covariates. For example,

```
sample_id   covariate_1     covariate_2 ...
 sample_1             1             5.4 ...
 sample_2             1             3.1 ...
      ...
```

### Genetic relatedness matrix

--grm/-g

A tab separated text file contaning the genetic relatedness-matrix in matrix fomat. For example, 

```
sample_id sample_1 sample_2 sample_3 sample_4 ...
  sample_1       1	   0.18     0.03     -0.3 
  sample_2     0.1	      1      0.4      0.1
  sample_3    0.04	   0.45        1      0.1
  sample_4    -0.1	    0.4        0        1 ...
       ...
```

To construct the GRM we recommend using the plink2 --make-king command after pruning variants.

## Option list

* --plink/-p: Plink files prefix
* --cov/-c: Covariate data file
* --bed/-b: Phenotype bed file
* --grm/-g: A (dense) genetic relatedness matrix
* --out/-o: The output file prefix
* --mode: The mode used to run quasar in. One of: 
    - cis
    - trans
    - gwas
* --model: The model used to residualise phenotype data. One of:
    - lm: Linear model
    - lmm: Linear mixed model
    - p_glm: Poissom GLM
    - nb_glm: Negative binomial GLM
    - p_glmm: Poisson GLMM
    - nb_glmm: Negative binomial GLMM
* --window_size/-w: The size of the cis window in base pairs. Default: 1000000
* --use-apl: Use Cox-Reid adjusted profile likelihood when estimating negative binomial dispersion
* --verbose: Write additional information to the console
