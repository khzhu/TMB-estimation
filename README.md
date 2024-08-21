# Tumor Mutational Burden (TMB) Project

## Overview
Tumor Mutational Burden (TMB) - the total number of mutations (changes) found in the DNA of cancer cells. TMB varies from different sequencing platforms to different type of cancer. Hence, it is important to harmonize variant selection methods and choose reasonable thresholds for the alternative allele count and the read depth. Accurate estimate of the total number of mutations in the targeted regions may help doctors plan the best treatment for each individual patient.

We present a nextflow workflow that calculates tumour mutational burden based on the criteria suggested by the [TMB Harmonization Consortium](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7174078/).
