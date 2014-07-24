Ancestry pipeline
=================
Overall description: Run PCA on a VCF, phase data, convert phased data to RFMix format, run local ancestry with RFMix, collapse RFMix output into bed files, alter bed files, plot karyograms, estimate global ancestry proportions, run TRACTS, generate PCAMask input, run PCA on PCAMask output (ASPCA).

This repo gives information about how to run through phasing, local ancestry inference, generate collapsed bed files, plot karyograms, estimate global ancestry proportion from local ancestry proportions, generate and run ASPCA, and run TRACTS to model migration events, proportions, and timings.

## Pipeline Map ##
#### 0.) Phase #####
* Overview
* SHAPEIT2 check
* SHAPEIT2 phasing
Make RFMix input
##### 1.) Infer local ancestry #####
  * Run RFMix
##### 2.1) Collapse inferred data #####
* Collapse RFMix output into TRACTS-compatible bed files
* Posthoc bed file filter (OPTIONAL)
* Plot ancestry karyograms
* Estimate global ancestry proportions from local ancestry inference
##### 2.2) Generate ASPCA input #####
* Make ASPCA input per chromosome
* Combine ASPCA input across chromosomes
##### 3.1) Model migration timings with TRACTS #####
* Model migration timings with TRACTS
* Helpful tips and tricks
##### 3.2) Run ASPCA #####
* Run PCAMask
Orthogonalize PCs
