#takes in a sample file and writes classes
from optparse import  OptionParser
import os

USAGE = """
intersect_snps.py   --anc1
                    --anc2
                    --anc3
                    --sample
                    --out
"""
parser = OptionParser(USAGE)

parser.add_option('--anc1', default='/home/armartin/rare/chip_collab/admixed/affy6/shapeit_in/YRI.inds')
parser.add_option('--anc2', default='/home/armartin/rare/chip_collab/admixed/affy6/shapeit_in/CEU.inds')
parser.add_option('--anc3', default='/home/armartin/rare/chip_collab/admixed/NATAM/nam_fwd_cleaned_hg19.fam')
parser.add_option('--sample', default='/home/armartin/rare/chip_collab/admixed/affy6/rfmix_input/ACB/ACB.sample')
parser.add_option('--out', default='/home/armartin/rare/chip_collab/admixed/affy6/rfmix_input/ACB/ACB.classes')

(options, args) = parser.parse_args()

anc1_set = set()
anc1 = open(options.anc1)
for line in anc1:
    anc1_set.add(line.strip())

anc2_set = set()
anc2 = open(options.anc2)
for line in anc2:
    anc2_set.add(line.strip())

anc3_set = set()
anc3 = open(options.anc3)
for line in anc3:
    anc3_set.add(line.strip().split()[1])
    #anc3_set.add(line.strip())
    
ind_order = []
sample = open(options.sample)
sample.readline()
sample.readline()
for line in sample:
    line = line.strip().split()
    ind_order.append(line[0])

out = open(options.out, 'w')
sample = open(options.sample)
for line in sample:
    line = line.strip().split()
    if line[0] in anc1_set:
        out.write('1 1 ')
    elif line[0] in anc2_set:
        out.write('2 2 ')
    elif line[0] in anc3_set:
        out.write('3 3 ')
    else:
        out.write('0 0 ')
out.write('\n')
out.close()
