Ancestry pipeline
=================
Overall description: Run PCA on a VCF, phase data, convert phased data to RFMix format, run local ancestry with RFMix, collapse RFMix output into bed files, alter bed files, plot karyograms, estimate global ancestry proportions, run TRACTS, generate PCAMask input, run PCA on PCAMask output (ASPCA).

This repo gives information about how to run through phasing, local ancestry inference, generate collapsed bed files, plot karyograms, estimate global ancestry proportion from local ancestry proportions, generate and run ASPCA, and run TRACTS to model migration events, proportions, and timings.

## Pipeline Map ##
#### 0.) Phase #####
* Overview
* SHAPEIT2 check
* SHAPEIT2 phasing
* Make RFMix input

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
* Orthogonalize PCs

## 0.) Phase ###
##### Overview #####
Any phasing tool will do, such as SHAPEIT2, BEAGLE, or MaCH. I usually run SHAPEIT2 (documentation here: http://www.shapeit.fr/) in two steps, first in check mode to identify SNPs that are incompatible and need to be excluded (also useful for summary info, which tells you about mendelian inconsistencies). 

#### 1) SHAPEIT2 check ####
Run according to the documentation. Examples as follows:
```
for i in {1..22}; 
do shapeit \
-check \
--input-ref output.SHAPEIT.20140226.chr${i}.hap.gz \
output.SHAPEIT.20140226.chr${i}.legend.gz \
output.SHAPEIT.20140226.chr${i}.sample \
-B nam_fwd_cleaned_hg19_ref_chr${i} \
--input-map genetic_map_chr${i}_combined_b37.txt \
--output-log nam_fwd_cleaned_hg19_ref_chr${i}.mendel; 
done
```

#### 2) SHAPEIT2 phasing ####
Second, I run the phasing algorithm itself (usually with a reference panel like phase 1 1000 Genomes), which takes some time and memory, for example as follows:
```
for i in {1..22}; 
do shapeit \
--input-ref output.SHAPEIT.20140226.chr${i}.hap.gz \
output.SHAPEIT.20140226.chr${i}.legend.gz \
output.SHAPEIT.20140226.chr${i}.sample \
-B nam_fwd_cleaned_hg19_ref_chr${i} \
--input-map genetic_map_chr${i}_combined_b37.txt \
--exclude-snp nam_fwd_cleaned_hg19_ref_chr${i}.mendel.snp.strand.exclude \
--output-max nam_fwd_cleaned_hg19_ref_chr${i}.haps \
nam_fwd_cleaned_hg19_ref_chr${i}.sample; done
```
Phasing should be parallelized across chromosomes and can be run with plink files or VCF files. The output logs from SHAPEIT2 indicate whether Medelian errors are present in family designations in plink .fam files. I have not yet found a way to run SHAPEIT2 incorporating family information with VCF files. The output files are *.haps and *.sample files. 

#### 3) Make RFMix input ####

I wrote a script to convert SHAPEIT2 output to RFMix input (as well as a .map file with 3 columns: physical position, cM position, rsID, which will be useful downstream), which you can run as follows:

```
for i in {1..22}; do python shapeit2rfmix.py \
--shapeit_hap1 CEU_chr${i}.haps \
--shapeit_hap2 LWK_chr${i}.haps \
--shapeit_hap3 SAN_all_chr${i}.haps \
--shapeit_hap_admixed SAN_all_chr${i}.haps \
--shapeit_sample1 CEU_chr${i}.sample \
--shapeit_sample2 LWK_chr${i}.sample \
--shapeit_sample3 SAN_all_chr${i}.sample \
--shapeit_sample_admixed SAN_all_chr${i}.sample \
--ref_keep CEU+LWK+SA_90.ref \
--admixed_keep CEU+LWK+SA_90.inf \
--chr ${i} \
--out CEU_LWK_SA; done
```

This script assumes that different reference and inference panels are phased separately and performs strand flip checking across files. It's totally fine (and actually better) if you phase them all together. You just need to beware that the output *.classes file may not be defined as you want, and I wrote another script that allows you to fix posthoc (currently only supports 3 reference panels, but should be easily made more flexible). To run the above script when all individuals are phased together, you can run as follows:

```
for i in {1..22}; do python shapeit2rfmix.py \
--shapeit_hap1 SA_550_OMNI_CEU_LWK_phase3_chr${i}.haps \
--shapeit_hap_admixed SA_550_OMNI_CEU_LWK_phase3_chr${i}.haps \
--shapeit_sample1 SA_550_OMNI_CEU_LWK_phase3_chr${i}.sample \
--shapeit_sample_admixed SA_550_OMNI_CEU_LWK_phase3_chr${i}.sample \
--ref_keep CEU_LWK_SAN.ref \
--admixed_keep CEU_LWK_SAN.notref \
--chr ${i} \
--out CEU_LWK_SA_phase3; done
```

After the fact, you can fix the classes file as follows:

```
python classes.py \
--ref pop1.ref,pop2.ref,pop3.ref \
--sample pops.sample \
--out pops.classes
```

## 1.) Infer local ancestry ##
#### Run RFMix ####
After RFMix input is generated, run RFMix. There is nice documentation on what all of the options mean here: https://sites.google.com/site/rfmixlocalancestryinference/. Here is an example run that I went through for 1000 Genomes:

```
for i in {1..22}; do for POP in ACB ASW CLM MXL PEL PUR; 
do python RunRFMix.py \
-e 5 \
-w 0.5 \
--num-threads 4 \
--use-reference-panels-in-EM \
--forward-backward \
PopPhased \
${POP}_chr${i}.alleles \
${POP}.classes \
${POP}_chr${i}.snp_locations \
${POP}_3_chr${i}.rfmix; done; done
```

## 2.1) Collapse inferred data ##
#### Collapse RFMix output into TRACTS-compatible bed files ####
After running RFMix, I always collapse the output into bed files and generate karyogram plots to ensure that there weren't upstream issues, such as class file errors, phasing technical artifacts, sample mixups, etc. I wrote a script to do this, which can be run for example as follows:
```
for POP in ACB ASW CLM MXL PEL PUR; do sed '1,143d' ${POP}.sample | while read line; do python collapse_ancestry.py \
--rfmix ${POP}_3_chrX.rfmix.5.Viterbi.txt \
--snp_locations ${POP}_chrX.snp_locations \
--fbk ${POP}_3_chrX.rfmix.5.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind ${line} \
--ind_info ${POP}.sample \
--pop_labels AFR,EUR,NAT \
--chrX \
--out ${line}"; done; done
```

#### Posthoc bed file filter (OPTIONAL) ####
After the first bed file is created, it might be desirable to mask certain regions, for example if a particular region is shown to be frequently misspecified empirically in the reference panel. I have only used this script once, so it almost assuredly has some bugs, but I have provided it here as a starting point in case posthoc masking is a desirable feature. An example run is as follows:

```
cat bed_files.txt | while read line; do \
python mask_lcr_bed.py \
--bed ${line} \
--out ${line}2;
done
```

#### Plot ancestry karyograms ####
I have also written a visualization script to plot karyograms, which can be run for example as follows:
```
line=bed_a.bed; line2=bed_b.bed; IND='NA12878'; python plot_karyogram.py \
--bed_a ${line} \
--bed_b ${line2} \
--ind ${IND} \
--out ${IND}.png"; done
```
Example output looks like the following:

![alt tag](https://aliciarmartindotcom.files.wordpress.com/2012/02/hg01892.png?w=800)

This script accepts a centromere bed file.

To do:
* Fix plot_karyogram.py so that the rounding at the ends of chromosomes occurs because the first and last chromosome tracts have been identified in the script, rather than required in the centromere bed file

#### Estimate global ancestry proportions from local ancestry inference ####

The last step is to calculate global ancestry proportions from the tracts. This can be useful to compare to orthogonal methods, i.e. ADMIXTURE, to see how well the ancestry estimates agree. This step can be run as follows:
```
for POP in ACB ASW CLM MXL PEL PUR; do python lai_global.py \
--bed_list bed_list_${POP}.txt \
--ind_list ${POP}.inds \
--pops AFR,EUR,NAT \
--out lai_global_${POP}.txt; done
```
The bed_list input here is a text file with a list of bed files, two per line and separated by whitespace, where each row corresponds to a single individual. For example:
```
ind1_a.bed    ind1_b.bed
ind2_a.bed    ind2_b.bed...
```
The ind_list input has individual IDs that will be used to summarize the output. The pops option specifies all of the populations to estimate global proportion ancestry for. I created this option so that UNK tracts could be easily dropped from global proportion ancestry estimated. An example txt output file is attached.

## 2.2) Generate ASPCA input ##
#### Make ASPCA input per chromosome ####

Once the local ancestry output has been generated from RFMix, ASPCA input can be generated for example as follows:
```
for i in {1..22}; do python make_aspca_inputs.py \
--alleles SA_550_Omni_CEU_LWK_SA_phase3_chr${i}.rfmix.allelesRephased5.txt \
--vit SA_550_Omni_CEU_LWK_SA_phase3_chr${i}.rfmix.5.Viterbi.txt \
--haps SA_550_Omni_CEU_LWK_SA_phase3_chr${i}.rfmix.5.map \
--sample SA_550_Omni_CEU_LWK_SA_phase3.sample \
--classes SA_550_Omni_CEU_LWK_SA_phase3.classes \
--out SA_550_Omni_CEU_LWK_SA_phase3_san_multi_chr${i}; done
```
There is currently a hard posterior threshold set. We have seen that this posterior thresholding can generate tighter clusters assuming a reasonably large number of markers are available.

To do: 
* This script could probably be made more flexible.

#### Combine ASPCA input across chromosomes ####

After this script is run, the input files need to be combined across chromosomes, for example as follows:
```
python combine_aspca_chrs.py \
--aspca_prefix SA_550_Omni_CEU_LWK_SA_phase3_san_multi_chr \
--keep_anc HGDP_Schuster_ref.inds \
--anc san \
--out SA_550_Omni_CEU_LWK_SA_phase3_san
```

## 3.1) Model migration timings with TRACTS ##
#### Model migration timings with TRACTS ####

After a cleaned set of bed files has been generated, migration timings can be modeled using TRACTS across all bed files in a given admixed population. I have attached an iPython notebook I generated that runs through one type of model. Typically multiple models will need to be tested to find the most likely fit (maximum likelihood). The resulting parameters from a single model (not necessarily the best one) have been plotted and are attached.

#### Helpful tips and tricks ####

Andres Moreno-Estrada has provided the following notes to run TRACTS:

if same models want to be applied to a new population = just adapt main script to indicate location of files, etc ==>  "MXL_samples_ppx_xxp.py" or "complete_CLM_ccp.py"

if want to test different models, then need to modify "cont_pulse_model_unif_params.py"  = this has the demographic models to be tested 
*contact Simon if want to modify...

also, no need to modify "admix_recomb.py" = it is the former name of tracts and before was imported from the main script, while now is still evoked but as "import tracts as" instead of "import" only… 

admix_recomb was the fomer name of tracts! 
You can do:

import tracts as admix_recomb instead of 
import admix_recomb,


Nomenclature:

p=pulse
c=continuous
x=nothing

number of parameters:
ancestral proportion pop 1
time pop 1
ancestral proportion pop 2
time pop 2

This model has for example 4 parameters

Increasing parameters will generally improve the likelihood, so choosing the best model is not simply taking that with the highest likelihood…
unless it is at least 5 log-likelihood units above the simple pulse model according to Simon's calculations…

input needed:
.bed files from PCAdmix
 

once output has been generated:

can get simple Tract length plot by typing "coin.mig" in python
currently, plots are generated in Mathematica by Simon



-- Python scripts by simon:

New ones have the "fix tot" suffix meaning (fix total ancestry numerically)

which means that these models fix the total amount of ancestry in the population
eg:

miami_scripts_ppx_xxp_fixtot_num


fpop function:
 
python is zero-based and "fpop" function is specified per population

fpop(0) means fpop applied to the first population


to run:

python miami_scripts_ppx_xxp_fixtot_num.py

no flags needed, everything is specified inside the script


contains all the models that we want to test:
cont_pulse_model_unif_params.py


whats the expected tract length distribution given the parameters
do they look like the observed data?
if not, optimize the parameters
find the highest likelihood

always try the simplest model as the base line…

keep the same architecture of where the data is in the computer

can be run locally….

maybe take overnight…

does not take much more memory

not optimized for use in the qsub…


-- scripts need to import packages that need to be previously installed:

import os,sys
sys.path.append("/Users/simongravel/Documents/tractlength/python/")
import tracts as admix_recomb 
import cont_pulse_model_unif_params
import numpy,pylab

ONLY one that pends to be installed:
numpy


Python needs to import all these from the computer into the memory of python

keep the name of admix_recomb and use "import tracts as admix_recomb"
or 
change all instances of "admix_recomb" by tracts

## 3.2) Run ASPCA
#### Run PCAMask

After ASPCA input is generated, run PCAMask. Reading the manual and considering the question you're trying to ask are important here. You can then run PCAMask, for example as follows:
```
i=550_Omni; PCAmask_linux \
-anc SA_${i}_CEU_LWK_SA_phase3_san_san.beagle \
-adm SA_${i}_CEU_LWK_SA_phase3_san_adm.beagle \
-vit SA_${i}_CEU_LWK_SA_phase3_san.vit \
-m SA_${i}_CEU_LWK_SA_phase3_san.markers \
-o SA_${i}_CEU_LWK_SA_phase3_san_posterior99.aspca \
-mask 3
```
#### Orthogonalize PCs

After running PCAMask, run PCA on the output to orthogonalize the results. You can use the simple_pca.py script (written by Chris Gignoux) to quickly run PCA on the output. For example:

```
python simple_pca.py SA_${i}_CEU_LWK_SA_phase3_san_posterior99.aspca_3.pca.txt
```
