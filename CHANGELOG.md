# Change log

## [1.2.0] - 24/8/2026

### Breaking changes

Due to the two fixes below the sign of effects in the negative binomial and Poisson GLMM models has reversed. The sign of effects in the linear model and linear mixed model is unchanged.

### Fixed

- Fixed direction of tails in qnorm and pnorm functions
- Fixed reading of alleles so behaviour now aligns with documentation

### Changed

- Relaxed reuquirement to always used residualised expression in trans/gwas mapping 

## [1.1.0] - 9/1/2026

### Added

- Add MAF in variant-level output
- Add check if phi <0
- Add better error message for non-autosomal chromosomes
- Add snp_id to variant-level output
- Add example data

### Changed

- Return both gene start and end in gene-level output
- Automatically add intercept if not present in the covariates
- Fix ref/alt order in variant-level output header

### Fixed

- Fix beta sizes when using LMM model
- Fix cis-window computation when genes are close to the end of the chromosome

## [1.0.0] - 29/7/2025

Initial release