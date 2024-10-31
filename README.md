# quasar

A software package for performing quantitative trait loci (QTL) mapping

## Installation

TODO

## Usage

```
./quasar --bed plink_name \\
    --cov covariate_data.tsv \\ 
    --feat_anno features.tsv \\
    --pheno phenotype_data.tsv \\
    --grm grm.txt
    --o quasar_output
```

### Genotype data

The genotype data should be in binary plink format with .bed/.bim/.fam files named `plink_name`.

### Covariate data

The covariate data is a tab-separated file with rows as samples and first column `sample_id` 
and other columns containing the covariates. For example,

```
sample_id   covariate_1     covariate_2 ...
 sample_1             1             5.4 ...
 sample_2             1             3.1 ...
      ...
```

### Phenotype data

The phenotype data is a tab-seperated file with where rows are features and first column 
`feature_id` and other columns are the sample ids. For example, 

```
feature_id  sample_1    sample_2    sample_3    sample_4    sample_5 ...
    PTPN22      39.2	   435.4       435.8        78.1      578.2 ...
       ...
```

### Feature annotation file

The feature annotation file contains information about the genetic locations of the 
quatative traits being mapped. It is tab separated file with columns: `feature_id`, 
`chrom`, `start`, `end`. Optional columns names include: `gene_name`. For example,

```
     feature_id   chrom        start         end   gene_name
ENSG00000134242       1    113871759   113813811      PTPN22
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

## Inspiration

Elements of the code and interface of quasar are based off

- (APEX)[https://corbinq.github.io/apex/doc/]
- regenie
- limix_qtl
- PQLseq2

