# Msc-GWAMA_Project-code-logs
# GWAMA Project

A pipeline for Genome-Wide Association Meta-Analysis (GWAMA) and downstream genetic analyses.

## Pipeline Steps

1. **Original GWAS Files**
   
https://www.ebi.ac.uk/gwas/studies/GCST004773
https://www.ebi.ac.uk/gwas/studies/GCST006801
https://www.ebi.ac.uk/gwas/studies/GCST90029024
https://www.ebi.ac.uk/gwas/studies/GCST90086068
https://www.ebi.ac.uk/gwas/studies/GCST90026417
https://www.ebi.ac.uk/gwas/studies/GCST90077725
https://www.ebi.ac.uk/gwas/studies/GCST90668070
https://www.ebi.ac.uk/gwas/studies/GCST90435704
https://www.ebi.ac.uk/gwas/studies/GCST007517
https://www.ebi.ac.uk/gwas/studies/GCST005413
https://www.ebi.ac.uk/gwas/studies/GCST005898
https://www.ebi.ac.uk/gwas/studies/GCST005047

2. **Meta-analysis using GWAMA tool**
   - Combines summary statistics across multiple GWAS datasets.

3. **LD Clumping (PLINK command line)**
   - Removes variants in high linkage disequilibrium to identify independent association signals.

4. **Sentinel Selection (R script)**
   - Selects the lead (sentinel) variants for each independent signal.

5. **Loci Merging (R script)**
   - Combines overlapping or neighboring genomic regions into distinct loci.

6. **TwoSampleMR (R script)**
   - Performs Mendelian Randomization analyses using the selected instruments.

7. **VEP (Variant Effect Predictor) online tool**
   - Annotates the functional consequences of the identified variants.

8. **LDSC (Linkage Disequilibrium Score Regression)**
   - Evaluates heritability and genetic correlation.
