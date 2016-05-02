Ancestry pipeline
=================
Overall description: Run PCA on a VCF, phase data, convert phased data to RFMix format, run local ancestry with RFMix, collapse RFMix output into bed files, alter bed files, plot karyograms, estimate global ancestry proportions, run TRACTS, generate PCAMask input, run PCA on PCAMask output (ASPCA).

This repo gives information about how to run through phasing, local ancestry inference, generate collapsed bed files, plot karyograms, estimate global ancestry proportion from local ancestry proportions, generate and run ASPCA, and run TRACTS to model migration events, proportions, and timings.

**Each step can be followed with the 1000 Genomes data (ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/hd_genotype_chip/) and systematically through this pipeline by downloading a toy dataset here: https://www.dropbox.com/sh/zbwka9u09f73gwo/AABc6FNl9fVBPjby8VQWzyeXa?dl=0. Slides from a tutorial I gave from ASHG 2015 using 1000 Genomes data are available here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20151008_ASHG15_tutorial/20151008_ASHG15_admixture.pdf.**

## Pipeline Map ##
#### 0.) Phase #####
* Overview
  * HAPI-UR
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
Any phasing tool will do, such as SHAPEIT2, BEAGLE, or HAPI-UR. I usually run SHAPEIT2 (documentation here: http://www.shapeit.fr/) for ~100s of samples in two steps, first in check mode to identify SNPs that are incompatible and need to be excluded (also useful for summary info, which tells you about mendelian inconsistencies). For larger datasets (~1000s+), HAPI-UR is more accurate and orders of magnitude faster.

```
for i in {1..22};
do plink \
--bfile ACB_example \
--chr ${i} \
--make-bed \
--out ACB_example_chr${i};
done
```

#### 1) SHAPEIT2 check ####
Run according to the documentation. Examples as follows:

```
for i in {1..22}; 
do shapeit \
-check \
--input-ref 1000GP_Phase3_chr${i}.hap.gz \
1000GP_Phase3_chr${i}.legend.gz \
1000GP_Phase3.sample \
-B ACB_example_chr${i} \
--input-map genetic_map_chr${i}_combined_b37.txt \
--output-log ACB_example_chr${i}.mendel; 
done
```

Note: 1000 Genomes reference panels and genetic maps can be downloaded here: https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html.
Be sure to remove fully missing SNPs.

#### 2) SHAPEIT2 phasing ####
Second, I run the phasing algorithm itself (usually with a reference panel like phase 1 1000 Genomes), which takes some time and memory, for example as follows:

```
for i in {1..22}; 
do shapeit \
--input-ref 1000GP_Phase3_chr${i}.hap.gz \
1000GP_Phase3_chr${i}.legend.gz \
1000GP_Phase3.sample \
-B ACB_example_chr${i} \
--duohmm \
--input-map genetic_map_chr${i}_combined_b37.txt \
--exclude-snp ACB_example_chr${i}.mendel.snp.strand.exclude \
--output-max ACB_example_chr${i}.haps.gz \
ACB_example_chr${i}.sample; done
```

Phasing should be parallelized across chromosomes and can be run with plink files or VCF files. The output logs from SHAPEIT2 indicate whether Medelian errors are present in family designations in plink .fam files. I have not yet found a way to run SHAPEIT2 incorporating family information with VCF files. The output files are *.haps and *.sample files. 

#### 3) Make RFMix input ####

I have written quick converters for two phasing programs: SHAPEIT2 and HAPI-UR.

##### HAPI-UR #####
The output from HAPI-UR is almost exactly what you need to run RFMix. Before running 

##### SHAPEIT2 #####
I also wrote a script to convert SHAPEIT2 output to RFMix input (as well as a .map file with 3 columns: physical position, cM position, rsID, which will be useful downstream), which you can run as follows:

```
for i in {1..22}; do python shapeit2rfmix.py \
--shapeit_hap_ref ACB_example_chr${i}.haps.gz \
--shapeit_hap_admixed ACB_example_chr${i}.haps.gz \
--shapeit_sample_ref ACB_example_chr${i}.sample \
--shapeit_sample_admixed ACB_example_chr${i}.sample \
--ref_keep ACB_example.ref \
--admixed_keep ACB_example.notref \
--chr ${i} \
--genetic_map genetic_map_chr${i}_combined_b37.txt \
--out ACB_example; done
```

However, if you phased them separately, the script will check for strand flip errors across phased files and correct them where possible when run e.g. as follows:

```
for i in {1..22}; do python shapeit2rfmix.py \
--shapeit_hap_ref CEU_example_chr${i}.haps.gz,YRI_example_chr${i}.haps.gz \
--shapeit_hap_admixed ACB_only_chr${i}.haps \
--shapeit_sample_ref CEU_example_chr${i}.sample,CEU_example_chr${i}.sample \
--shapeit_sample_admixed ACB_only_chr${i}.sample \
--ref_keep ACB.ref \
--admixed_keep ACB.notref \
--chr ${i} \
--genetic_map genetic_map_chr${i}_combined_b37.txt \
--out CEU_YRI_ACB; done
```

If you run the script as in the first example, you will need to fix the classes file , e.g. as follows:

```
python classes.py \
--ref CEU_example_chr22.keep,YRI_example_chr22.keep \
--sample CEU_YRI_ACB.sample \
--out CEU_YRI_ACB.classes
```

Note that the sample file option here denotes a list of individuals as output by the previous script (consistent with individual order in alleles file), and NOT the shapeit sample file. Each ref file has at least two columns (similar to a plink keep file), with the relevant 2nd column corresponding with the individual ID.

## 1.) Infer local ancestry ##
#### Run RFMix ####
After RFMix input is generated, run RFMix. There is nice documentation on what all of the options mean here: https://sites.google.com/site/rfmixlocalancestryinference/. Here is an example run that I went through for 1000 Genomes:

```
for i in {1..22}; do \
python RunRFMix.py \
-e 2 \
-w 0.2 \
--num-threads 4 \
--use-reference-panels-in-EM \
--forward-backward \
PopPhased \
CEU_YRI_ACB_chr${i}.alleles \
CEU_YRI_ACB.classes \
CEU_YRI_ACB_chr${i}.snp_locations \
-o CEU_YRI_ACB_chr${i}.rfmix; done
```

## 2.1) Collapse inferred data ##
#### Collapse RFMix output into TRACTS-compatible bed files ####
After running RFMix, I always collapse the output into bed files and generate karyogram plots to ensure that there weren't upstream issues, such as class file errors, phasing technical artifacts, sample mixups, etc. I wrote a script to do this, which can be run for example as follows:
```
python collapse_ancestry.py \
--rfmix CEU_YRI_ACB_chr1.rfmix.2.Viterbi.txt \
--snp_locations CEU_YRI_ACB_chr1.snp_locations \
--fbk CEU_YRI_ACB_chr1.rfmix.5.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind HG02481 \
--ind_info CEU_YRI_ACB.sample \
--pop_labels EUR,AFR \
--chrX \
--out HG02481; done; done
```

Note: all autosomes must have successfully completed, and including chromosome X is optional with the flag. The order of the population labels correspond with the order of labels in the classes file.

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
IND='HG02481'; python plot_karyogram.py \
--bed_a ${IND}_A.bed \
--bed_b ${IND}_B.bed \
--ind ${IND} \
--out ${IND}.png
```
Example output looks like the following:

![alt tag](https://aliciarmartindotcom.files.wordpress.com/2012/02/hg02481.png?w=800)

This script accepts a centromere bed file (see Dropbox data).

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
--alleles ACB_chr${i}.rfmix.allelesRephased1.txt \
--vit ACB_chr${i}.rfmix.1.Viterbi.txt \
--markers ACB_chr${i}.map \
--inds ACB.sample \
--classes ACB.classes \
--out ACB_chr${i}; done
```
There is currently a hard posterior threshold set. We have seen that this posterior thresholding can generate tighter clusters assuming a reasonably large number of markers are available.

#### Combine ASPCA input across chromosomes ####

After this script is run, the input files need to be combined across chromosomes, for example as follows:
```
python combine_aspca_chrs.py \
--aspca_prefix ACB_chr \
--anc EUR \
--out ACB
```

## 2.3) Run ASPCA ##
#### Run PCAMask ####

After ASPCA input is generated, run PCAMask (https://sites.google.com/site/pcamask/dowload). Reading the manual and considering the question you're trying to ask are important here. You can then run PCAMask, for example as follows:
```
PCAmask_linux \
-anc ACB_anc.beagle \
-adm ACB_adm.beagle \
-vit ACB.vit \
-m ACB.markers \
-o ACB.aspca \
-mask 2
```

#### Orthogonalize PCs ####

After running PCAMask, run PCA on the output to orthogonalize the results. You can use the simple_pca.py script (written by Chris Gignoux) to quickly run PCA on the output. For example:

```
python simple_pca.py ACB.aspca_2.pca.txt
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
