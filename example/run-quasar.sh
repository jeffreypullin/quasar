#!/bin/bash

# Linear model.
/home/jp2045/quasar/build/quasar \
  --plink chr22-n100 \
  --bed mean-pheno-n100.bed \
  --cov cov-n100.tsv \
  -o lm-example \
  --model lm \
  --mode cis

# Linear model with pgen.
/home/jp2045/quasar/build/quasar \
  --pgen chr22-n100 \
  --bed mean-pheno-n100.bed \
  --cov cov-n100.tsv \
  -o lm-pgen-example \
  --model lm \
  --mode cis

# Negative binomial GLM.
/home/jp2045/quasar/build/quasar \
  -p chr22-n100 \
  -b sum-pheno-n100.bed \
  -c cov-n100.tsv \
  -o nb_glm-example \
  --model nb_glm \
  --use-apl \
  --mode cis

# Linear mixed model.
/home/jp2045/quasar/build/quasar \
  -p chr22-n100 \
  -b mean-pheno-n100.bed \
  -c cov-n100.tsv \
  -o lmm-example \
  -g grm-n100.tsv \
  --model lmm \
  --mode cis

