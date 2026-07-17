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

To run quasar with bulk or pseudobulk data use: 

### Bulk/pseduobulk data

```
./quasar \
    --plink plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --mode cis \
    --model bulk-out \
    --use-apl \
    --out nb_fit
```

### Single-cell resolution data

As of quasar 2.0, quasar can take single-cell data as input. To use this functionality run: 

```
./quasar \
    --plink plink_prefix \
    --sc-pheno single-cell-pheno.tsv \
    --anno annotations.tsv \
    --cov single_cellcovariate_data.tsv  \
    --mode cis \
    --model p_glmm_sc \
    --out single-cell-out
```

Single-cell data requires different formatted data to bulk/pseudobulk data, for more information see below.

### Interaction QTLs

As of quasar 2.0, quasar can compute interaction-QTLs. To use this functionality, simply specify which covariate the interaction should be tested for.

#### Bulk/pseudobulk interaction QTLs

```
./quasar \
    --plink plink_prefix \
    --bed phenotype_data.bed \
    --cov covariate_data.tsv \
    --mode cis \
    --interaction sex \
    --model lm \
    --use-apl \
    --out sex-int-out
```

#### Single-cell interaction QTLs

```
./quasar \
    --plink plink_prefix \
    --sc-pheno single-cell-pheno.tsv \
    --anno annotations.tsv \
    --cov single_cell_covariate_data.tsv \
    --mode cis \
    --model p_glmm_sc \
    --interaction pseudotime \
    --out sc-int-out
```

Interaction testing can be performed for bulk/pseudobulk data with the linear mixed model (`lmm`), linear model (`lm`) and Negative Binomial GLM (`nb_glm`) models, and for single-cell data with the Poisson GLMM (`p_glmm_sc`).

For single-cell data, interaction testing requires `--model p_glmm_sc`, which fits a random-slope null model (a per-donor random intercept and random slope on the interaction covariate). When covariates are provided in single-cell (cell-level) format, quasar splits the interaction covariate into between-donor (`{interaction_cov}_b`) and within-donor (`{interaction_cov}_w`) components and tests the interaction on the within-donor component (`{interaction_cov}_w`). `--interaction` cannot be combined with `--cell-groups.

When interaction testing is enabled, quasar inspects the interaction covariate values and automatically augments the nuisance covariates as follows:

* if the interaction covariate has `<=10` unique finite values, it is treated as categorical and no squared nuisance term is added;
* if it has `>10` unique finite values, it is treated as continuous and quasar automatically adds a nuisance covariate named `{interaction_cov}_sq` which contains the square of the covariate

## QTL mapping modes

The quasar software can be run in three modes: `cis`, `trans`, `gwas`. These modes specify which varaints are tested for asscociation with a particular feature. 

In mode `cis` variants within +- the window size of the gene (see phenotype data format for details). By default the window size is set to 1Mb but can be specified using the `-w/--window` flag. We refer to this set of variants as the cis-window for that feature. 

In mode `trans` all variants except those in the cis window are tested for assoication. 

In mode `gwas` all variants are are tested for association. 

Note that trans-eQTL mapping can be affected by cross-mapping reads and paralogous genes causing spurious associations.

## Statistical models

The quasar software package supports a wide range of statistical models used to resiudalise the expression values. The supported models are:

* `lm`: linear model
* `nb_glm`: negative binomial GLM 
* `lmm`: linear mixed model
* `p_glmm`: Poisson generalised linear mixed model (GLMM)
* `p_glm`: Poisson GLM (**not recommended** due to producing a very high rate of false positives)
* `nb_glmm`: negative binomial GLMM (**not generally recommended** due to producing highly similar results to the Poisson GLMM while being slower, can be used if there is known to be high relatedness between samples)
* `p_glmm_sc`: A Poisson GLMM accounting for repeated measures that can be used for single-cell level data. Without `--interaction` this is a random-intercept model; with `--interaction` it fits a per-donor random intercept and random slope on the interaction covariate.

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

#### Bulk/pseudobulk bed file 

--bed/-b

The phenotype data  a tab-seperated file with where rows are features and the first four
columns give feature information and the rest are sample ids are the sample ids. For example, 

```
#chr      start         end      phenotype_id  sample_1   sample_2   sample_3 ...
   1  113871759   113813811   ENSG00000134242      39           43         45 ...
 ...
```

The start and end values are used to specify the centre of the cis-window. To specify the gene TSS as the centre of the window, set TSS = start, end = start + 1, so that the cis-window is [TSS - window, TSS + window + 1] or alternatively set the start and end values to the start and end of the gene so that the cis-window is [start - window, end + window].

In mode `gwas`, genomic coordinates are optional. The bed file may use either the standard header above or a header with only `phenotype_id` followed by sample columns:

```
phenotype_id  sample_1   sample_2   sample_3 ...
ENSG00000134242      39         43         45 ...
```

When coordinates are omitted, all variants are tested and residual output (`{out-prefix}-resids.bed`) is written without `#chr`, `start`, or `end` columns. When coordinates are provided in `gwas` mode they are read and written back but are not used to define cis-windows. Modes `cis`, `trans`, and `residualise` require the four-column annotation format.

For the count based models (i.e. `nb_glm`, `p_glm`, `p_glmm` and `nb_glmm`) count data should be passed to quasar. This can be either RNA-seq counts or pseudobulk scRNA-seq counts (the sum of the counts over the inidivdual). For the linear models (i.e. `lm` and `lmm`) we recommend that when analysing scRNA-seq counts the mean over individuals is passed to quasar.

#### Single-cell phenotype data

--sc-pheno

The single-cell phenotype data is provided as tab-seperated files, where each row is one cell and all columns except the first two are different genes. The first two column hold the sample_id and (unique) cell_id of each cell.

```
sample_id     cell_id     gene_1     gene_2    gene_3    ...
 sample_1      cell_1          0          1         0    ... 
 sample_1      cell_2          0          0         0    ... 
 sample_2      cell_3          2          0         0    ... 
 sample_2      cell_4          0          0         1    ... 
     ...
```

### Covariate data

--cov/-c

The covariate data can be specified in bulk or single-cell formats. Both bulk and single-cell foramts can be used with single-cell data (although the bulk format cannot hold single-cell resolution covariates) but only the bulk format can be used with bulk/pseudobulk data.

In bulk format, the covariate data is a tab-separated file with rows as samples and first column `sample_id` and other columns containing the covariates. For example,

```
sample_id   covariate_1     covariate_2 ...
 sample_1             1             5.4 ...
 sample_2             1             3.1 ...
      ...
```

In single-cell format, the covariate data is a tab-separated file where rows are cells, the first column is `sample_id`, the second column is `cell_id` and other columns contain the covariates. For example, 

```
sample_id   cell_id    covariate_1     covariate_2 ...
 sample_1    cell_1              1             5.4 ...
 sample_1    cell_2              1             3.1 ...
 sample_2    cell_3              1             2.7 ...
 sample_2    cell_4              1             1.0 ...
      ...
```

If an intercept is not present in the covariate data it will be added automatically as of quasar 1.1

### Annotation data

--anno

In single-cell mode an annotation file containing information about features/genes must be passed to as this information is not stored in the phenotype file, as in bulk data. The annotation file should be a tab-seperated bed file with columns #chr, start, end and phenotype_id. For example,

```
#chr      start         end      phenotype_id      ...
   1  113871759   113813811   ENSG00000134242      ...
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

To construct the GRM we recommend using the plink2 --make-king command after pruning variants. The resulting matrix will then need to be multiplied by 2, and possibly slightly altered, for example by setting negative eigenvalues to 0, to ensure it is positive definite. Other methods for constructing the GRM should work but have not been evaluated.

## Output

In cis mode, quasar produces two files:

* {out-prefix}-quasar-variant.txt which contains variant information
* {out-prefix}-quasar-cis-region.txt which contains gene information

### Variant output

The variant level output of quasar has the following basic format:

```
     feature_id            snp_id     chrom       pos      alt     ref     maf     beta     se     pvalue
ENSG00000100181    22:16849971A-T        22  16849971        T      A     0.39    0.012  0.038     0.7385
           ...
```

In the output the `alt` allele is the effect allele. Other columns including `glm_converged`, `glmm_converged`, `phi`, `phi_converged` encode information about the gene-level models fit to the expression data and are included depending on the type of model used.

### Region level output

In cis mode, quasar produces a summary 

```
     feature_id     chrom       start       end     pvalue
ENSG00000100181        22    17082776  17082777   0.796123
ENSG00000069998        22    17646176  17646177  0.0388123
            ...
```

#### Interaction testing output

For variant level output of interaction testing quasar produces two effect sizes, standard errors and pvalues for each variant. These correspond to the estimates of $\beta$ and $\gamma$ in the model

$$
y = g $\beta + (g x x) \gamma
$$ 

where $g$ is the vector of genotypes and $x$ is the vector of interaction covariate. These two values for each variant are denoted as `snp_*` and `snp_x_{interaction covariate}` respectively. For example if the interaction covariate was `sex` the output would be of the form:

```
     feature_id            snp_id     chrom       pos      alt     ref     maf     snp_beta     snp_se     snp_pvalue    snp_x_sex_beta    snp_x_sex_se     snp_x_pvalue
ENSG00000100181    22:16849971A-T        22  16849971        T      A     0.39        0.012      0.038         0.7385              0.02            0.03              0.5
           ...
```

## Option list

| Option | Argument | Type | Description|
|--------|-------|------|----|
|`--plink` | FILE | Required | Plink files prefix, assumes that `{prefix}.bed`, `{prefix}.bim`, `{prefix}.fam` exist |
|`--cov` | FILE | Required | Covariate data file |
|`--bed` | FILE | Required (bulk/pseduobulk) | Phenotype bed file |
|`--sc-pheno` | FILE | Required (single-cell) | Single-cell phenotype file |
|`--anno` | FILE | Required (single-cell) | Annotation file |
|`--grm` | FILE | Optional | A (dense) genetic relatedness matrix |
|`--out` | STRING | Optional | The output file prefix |
|`--mode`  | STRING | Required | The mode used to run quasar in. One of: `cis`, `trans`, `gwas`. |
|`--interaction` | STRING | Optional | The covariate to perform interaction-QTL testing with. Must be a column header in the covariate data. |
|`--model` | STRING | Required | The model used to residualise phenotype data. One of: `lm`, `lmm`, `p_glm`, `nb_glm`, `p_glmm`, `nb_glmm`, `p_glmm_sc`. |
|`--window_size` | NUMBER | Optional | The size of the cis window in base pairs. Default: 1000000 |
|`--use-apl` | FLAG | Optional | Use Cox-Reid adjusted profile likelihood when estimating negative binomial dispersion |
|`--verbose` | FLAG | Optional | Write additional information to the console |
