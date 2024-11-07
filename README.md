# Pcitrea_fm_microbes 

This repo includes all the code necessary to reproduce results presented in the associated manuscript. Below, you will find an overview of the repo's contents, along with guidance on individual README files that provide detailed instructions and descriptions for specific sections.

**Associated manuscript:** *Matthews, A.E., B.K. Trevelline, A.J. Wijeratne, and T.J. Boves. 2024. Picky eaters: selective
microbial diet of avian ectosymbionts. Journal of Animal Ecology. DOI: 10.1111/1365-2656.14215*

The purpose of this study was to understand the dietary selection of vane-dwelling feather mites found on Prothonotary Warblers (*Protonotaria citrea*) and how their dietary selection impacts their functional nature with hosts (i.e., are they parasites, mutualists, or commensals?). We found that feather mites exhibit a selective diet of feather microbes, many of which have known feather-degrading functions. Our findings support the idea that mites act as mutualistic feather “cleaners” and highlight the synergistic interactions between co-occurring symbionts in animal ecosystems.

---
## Contents overview
This repo is organized into two main sections:

### 16S_qiime_20221004
- Contains analyses for 16S (bacteria) sequences of feathers and mite data
- This directory is divided into several subdirectories that are numbered in order of operation:
  - 01_batch_effects: test for sequencing batch effects
  - 02_Run1Run2_combo: combine data from two sequencing runs, run initial analyses (e.g., trim primers, DADA2, taxonomy/SILVA)
  - 03_decontam: decontaminate reads (i.e., leave only bacteria)
  - 04_SeparateSamples: separate biological samples from control samples
  - 05_AlphaBetaDiversity: analyze alpha and beta diversity stats
  - 06_RPCA: run RPCA analyses
  - 07_TaxaBarplots: create taxa barplots
  - 08_R: plotting and analyses run in R (e.g., 'pretty plots,' dietR)
  - 09_Maaslin: differential abundance analyses
  - 11_picrust: predictive metagenomic analyses (removed, but available upon request)

### ITS_qiime_20230220
- Contains analyses for ITS (fungi) sequences of feathers and mite data
- This directory is divided into several subdirectories that are numbered in order of operation:
  - 01_ForwardData: trim adapters
  - 02_dada2: run DADA2
  - 03_taxonomy: run taxonomy/UNITE
  - 04_decontam: decontaminate reads (i.e., leave only fungi)
  - 05_SeparateSamples: separate biological samples from control samples
  - 06_AlphaBetaDiversity: analyze alpha and beta diversity stats
  - 07_TaxaBarplots: create taxa barplots
  - 08_RPCA: run RPCA analyses
  - 09_R: plotting and analyses run in R (e.g., 'pretty plots,' dietR)
  - 10_Maaslin: differential abundance analyses
  - 11_FUNguild: predictive functional analyses (removed, but available upon request)

## Individual subdirectory READMEs
Each subdirectory has a dedicated README file providing more detailed information, such as the scripts used to generate results/output files, analysis guidelines, data formats, as well as variable definitions. Please consult these individual READMEs for specific instructions related to each component of this repo.

---

**NOTE:** This repo has been "cleaned up" from a private repo for the public.
