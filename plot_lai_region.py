from optparse import  OptionParser
import gzip
from datetime import datetime
import time
from string import maketrans
import bedparser
import operator
import brewer2mpl
import pylab
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol

def split_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

USAGE = """
plot_lai_region.py  --bed_files
                    --pheno
                    --pheno_col
                    --position
                    --out
"""

parser = OptionParser(USAGE)

parser.add_option('--bed_files', default='/vault/henn/people/alicia/lai_files/bed_files.txt')
parser.add_option('--pheno', default='/vault/henn/people/alicia/phenos/PhenotypeMatrix_Avg_Oct2013.csv')
parser.add_option('--pheno_col', default='M.Arm.Avg')
parser.add_option('--position', help='chr,start,stop', action='callback', type='str', callback=split_callback)
parser.add_option('--title', default='')
parser.add_option('--out', default='/home/armartin/sa_analysis/ihs/SA_Omni_only_phase3_chr22')

(options,args)=parser.parse_args()

print options
#read all phenotype files and store phenotype of interest into a dict
pheno = open(options.pheno)
pheno_column = pheno.readline().strip().split(',').index(options.pheno_col)
pheno_dict = {}
for line in pheno:
    line = line.strip().split(',')
    if line[pheno_column] != 'NA':
        pheno_dict[line[0]] = float(line[pheno_column])

#store paths to all bed_files
bed_files = open(options.bed_files)
ind_paths = {}
for line in bed_files:
    line = line.strip()
    ind = line.split('/')[-1].split('_')[0]
    ind_path = '/'.join(line.split('/')[0:-1])
    ind_paths[ind] = ind_path + '/' + ind

dicbedsA={}
dicbedsB={}

for indiv in pheno_dict:
    if indiv in ind_paths:
        dicbedsA[indiv]=bedparser.bed(ind_paths[indiv]+'_A.bed')
        dicbedsB[indiv]=bedparser.bed(ind_paths[indiv]+'_B.bed')

#get tract indices for most 5' positions. for each individual, store a list of start index and end indices
pos = options.position
tract_A = {}
tract_B = {}
for ind in pheno_dict:
    if ind in ind_paths:
        tract_A[ind] = [dicbedsA[ind].loc(int(pos[1]), dicbedsA[ind].chrstarts[pos[0]], dicbedsA[ind].chrends[pos[0]]),
                        dicbedsA[ind].loc(int(pos[2]), dicbedsA[ind].chrstarts[pos[0]], dicbedsA[ind].chrends[pos[0]])]
        tract_B[ind] = [dicbedsB[ind].loc(int(pos[1]), dicbedsB[ind].chrstarts[pos[0]], dicbedsB[ind].chrends[pos[0]]),
                        dicbedsB[ind].loc(int(pos[2]), dicbedsB[ind].chrstarts[pos[0]], dicbedsB[ind].chrends[pos[0]])]

#sort individuals by phenotype
sorted_phenos = sorted(pheno_dict.iteritems(), key=operator.itemgetter(1))
print sorted_phenos

#define plotting space
brewer_vec = brewer2mpl.get_map('Set1', 'qualitative', 3).hex_colors
brewer_vec = {'CEU': brewer_vec[0], 'LWK': brewer_vec[1], 'SAN': brewer_vec[2], 'UNK': 'black'}

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(float(pos[1])/1e6,float(pos[2])/1e6)
ax.set_ylim(len(pheno_dict.keys())*2,0)
plt.xlabel('Physical position (Mb)')
plt.ylabel('Haplotype (ordered by increasing ' + options.pheno_col + ')')
plt.title(options.title)
#plt.yticks(range(1,23))

##plot rectangles
#for line in bed_a:
#    line = line.strip().split()
#    try:
#      plot_rects(line[3], int(line[0]), line[4], line[5], 'A', pop_order, colors, chrX)
#    except ValueError: #flexibility for chrX
#      plot_rects(line[3], 23, line[4], line[5], 'A', pop_order, colors, chrX)
#for line in bed_b:
#    line = line.strip().split()
#    try:
#      plot_rects(line[3], int(line[0]), line[4], line[5], 'B', pop_order, colors, chrX)
#    except ValueError: #flexibility for chrX
#      plot_rects(line[3], 23, line[4], line[5], 'B', pop_order, colors, chrX)

fig.savefig(options.out)
#dicbedsA['SA029'].chrdict['1'][0].name #this gives the name of the tract for this position
#dicbedsA['SA029'].loc(70513248, dicbedsA['SA029'].chrstarts['3'], dicbedsA['SA029'].chrends['3']) #finds the index of the tract that this position resides in
