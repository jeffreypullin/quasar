# Usage

To run **quasar**, use the command `./quasar` on the command line. Flags and options specify how quasar will run. 

To list all the possible options and see quasar's help you can run: 

```
./qusar --help
```

A normal use of quasar, providing the necessary data with all other options at their default values is:

```
./quasar --plink_prefix plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --grm grm.tsv
```

## Data formats

### Genotype data

The genotype data should be in binary plink format with .bed/.bim/.fam files named `plink_prefix`.

### Phenotype data

The phenotype data is a tab-seperated file with where rows are features and the first four
columns give feature inforamtion and the rest are sample ids are the sample ids. For example, 

```
#chr      start         end      phenotype_id  sample_1   sample_2   sample_3 ...
   1  113871759   113813811   ENSG00000134242      39.1      435.4      435.8 ...
 ...
```

### Covariate data

The covariate data is a tab-separated file with rows as samples and first column `sample_id` 
and other columns containing the covariates. For example,

```
sample_id   covariate_1     covariate_2 ...
 sample_1             1             5.4 ...
 sample_2             1             3.1 ...
      ...
```
### Genetic relatedness matrix

A tab separated text file contaning the genetic relatedness-matrix in matrix fomat. 
For example, 

```
sample_id sample_1 sample_2 sample_3 sample_4 ...
  sample_1       1	   0.18     0.03     -0.3 
  sample_2     0.1	      1      0.4      0.1
  sample_3    0.04	   0.45        1      0.1
  sample_4    -0.1	    0.4        0        1 ...
       ...
```

## Other options

* --model/-m: The model used to residualise phenotype data. Can be:
    - lm: Linear model
    - lmm: Linear mixed model
    - glm: Negative binomial generalised linear model
    - glmm: Poisson generalised linear mixed model
* --window/-w: The size of the cis window in base pairs. Default: 1000000
* --use-apl: Use Cox-Reid adjusted profile likelihood when estimating negative binomial dispersion