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

#gets individuals with both genotype and phenotype data
pheno_geno = set(pheno_dict.keys()).intersection(set(dicbedsA.keys()))

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
brewer_vec = {'CEU': brewer_vec[0], 'LWK': brewer_vec[1], 'SAN': brewer_vec[2], 'UNK': 'k'}

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(float(pos[1])/1e6,float(pos[2])/1e6)
ax.set_ylim(len(pheno_dict.keys())*2,0)
plt.xlabel('Physical position (Mb)')
plt.ylabel('Haplotype (ordered by increasing ' + options.pheno_col + ')')
plt.title(options.title)
plt.yticks(range(len(pheno_dict.keys())*2,0))

##plot rectangles
def plot_rects(bottom, left, top, right, color):
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    verts = [
            (left, bottom), #left, bottom
            (left, top), #left, top
            (right, top), #right, top
            (right, bottom), #right, bottom
            (0, 0), #ignored
    ]
    
    clip_path = Path(verts, codes)
    col=mcol.PathCollection([clip_path],facecolor=color, linewidths=0)
    ax.add_collection(col)

#num_inds =
counter = len(pheno_geno)*2
for ind in sorted_phenos:
    this_ind = ind[0]
    #print this_ind
    if this_ind in pheno_geno:
        #if this_ind == 'SA1000':
        print this_ind
        for tract in range(len(dicbedsA[this_ind].chrstarts[pos[0]])):
            print [counter, dicbedsA[this_ind].chrstarts[pos[0]][tract]/1e6, counter+1, dicbedsA[this_ind].chrends[pos[0]][tract]/1e6,
                   brewer_vec[dicbedsA[this_ind].chrdict[pos[0]][tract].name]]
            plot_rects(counter, dicbedsA[this_ind].chrstarts[pos[0]][tract]/1e6, counter+1, dicbedsA[this_ind].chrends[pos[0]][tract]/1e6,
                       brewer_vec[dicbedsA[this_ind].chrdict[pos[0]][tract].name])
        print counter
        counter = counter - 1
        print counter
        for tract in range(len(dicbedsB[this_ind].chrstarts[pos[0]])):
            print [counter, dicbedsB[this_ind].chrstarts[pos[0]][tract]/1e6, counter+1, dicbedsB[this_ind].chrends[pos[0]][tract]/1e6,
                       brewer_vec[dicbedsB[this_ind].chrdict[pos[0]][tract].name]]
            plot_rects(counter, dicbedsB[this_ind].chrstarts[pos[0]][tract]/1e6, counter+1, dicbedsB[this_ind].chrends[pos[0]][tract]/1e6,
                       brewer_vec[dicbedsB[this_ind].chrdict[pos[0]][tract].name])
        counter = counter - 1

p = []
pops = ['CEU', 'LWK', 'SAN', 'UNK']
for i in range(len(pops)):
    p.append(plt.Rectangle((0, 0), 1, 1, color=brewer_vec[pops[i]]))
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

ax.legend(p, pops, loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=5)

#leg = ax.legend(p, pops, loc=4, fancybox=True, bbox_to_anchor=(0.5, -0.05), shadow=True)
#leg.get_frame().set_alpha(0)

fig.savefig(options.out)
