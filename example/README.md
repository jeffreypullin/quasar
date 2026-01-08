# Example

The example data in this directory is a subset of the OneK1K data, processed as derscribed in

> Flexible and efficient count-distribution and mixed-model methods for eQTL mapping with quasar  
> Jeffrey M. Pullin, Jarny Choi, Davis J. McCarthy, Chris Wallace  
> medRxiv 2025.07.17.25331702; doi: https://doi.org/10.1101/2025.07.17.25331702  

The subset contains data for the B cell naive cluster filtered to the first 100 individuals. The variant data is restricted to chromosome 22 and the phenotype data is restrctied to the first 20 genes on chromosome 22. Phenotype data is provided for both mean and sum processing/aggregation approaches.

Examples of running quasar on this data are contained in `run-quasar.sh`.
