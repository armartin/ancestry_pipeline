__author__ = 'armartin'
from optparse import  OptionParser
import matplotlib.pyplot as plt
import pylab
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol
import brewer2mpl
import os

def splitstr(option, opt, value, parser):
  return(setattr(parser.values, option.dest, value.split(',')))

USAGE = """
plot_karyogram.py   --bed_a
                    --bed_b
                    --pop_order
                    --out
"""
parser = OptionParser(USAGE)

parser.add_option('--bed_a', default='/Users/alicia/Dropbox/Shared/SA_Phenotypes/Ancestry/RFMix_LocalAncestry/CEU_LWK_SAN/rfmix1.5.4/550/bed_files/SA006_A.bed')
parser.add_option('--bed_b', default='/Users/alicia/Dropbox/Shared/SA_Phenotypes/Ancestry/RFMix_LocalAncestry/CEU_LWK_SAN/rfmix1.5.4/550/bed_files/SA006_B.bed')
parser.add_option('--ind', default=None)
parser.add_option('--chrX', help='include chrX?', default=False, action="store_true")
parser.add_option('--centromeres', default='/home/armartin/rare/chip_collab/admixed/affy6/lai_output/centromeres.bed')
parser.add_option('--pop_order', default=['AFR','EUR','NAT'], type='string', action='callback', callback=splitstr,
                  help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
parser.add_option('--out')

(options, args) = parser.parse_args()

def plot_rects(anc, chr, start, stop, hap, pop_order, colors, chrX):    
    centro_coords = map(float, centromeres[str(chr)])
    if len(centro_coords) == 3: #acrocentric chromosome
        mask = [
        (centro_coords[1]+2,chr-0.4), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chr-0.4),
        (centro_coords[2]+2,chr),
        (centro_coords[2]-2,chr+0.4),
        (centro_coords[1]+2,chr+0.4),
        (centro_coords[1]-2,chr),
        (centro_coords[1]+2,chr-0.4)
        ]
        
        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)
    
    else: #need to write more complicated clipping mask with centromere masked out
        mask = [
        (centro_coords[1]+2,chr-0.4), #add +/- 2 at the end of either end
        (centro_coords[2]-2,chr-0.4),
        (centro_coords[2]+2,chr+0.4),
        (centro_coords[3]-2,chr+0.4),
        (centro_coords[3]+2,chr),
        (centro_coords[3]-2,chr-0.4),
        (centro_coords[2]+2,chr-0.4),
        (centro_coords[2]-2,chr+0.4),
        (centro_coords[1]+2,chr+0.4),
        (centro_coords[1]-2,chr),
        (centro_coords[1]+2,chr-0.4)
        ]
        
        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)
        
    if hap == 'A': #bed_a ancestry goes on top
        verts = [
            (float(start), chr), #left, bottom
            (float(start), chr + 0.4), #left, top
            (float(stop), chr + 0.4), #right, top
            (float(stop), chr), #right, bottom
            (0, 0), #ignored
        ]
    else: #bed_b ancestry goes on bottom
        verts = [
            (float(start), chr - 0.4), #left, bottom
            (float(start), chr), #left, top
            (float(stop), chr), #right, top
            (float(stop), chr - 0.4), #right, bottom
            (0, 0), #ignored
        ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    clip_path = Path(verts, codes)
    if anc in pop_order:
        col=mcol.PathCollection([clip_path],facecolor=colors[pop_order.index(anc)], linewidths=0)
    else:
        col=mcol.PathCollection([clip_path],facecolor=colors[-1], linewidths=0)
    if 'clip_mask' in locals():
        col.set_clip_path(clip_mask, ax.transData)
    ax.add_collection(col)


#read in bed files and get individual name
bed_a = open(options.bed_a)
bed_b = open(options.bed_b)
pop_order = options.pop_order
if options.ind is None:
    ind = '_'.join(options.bed_a.split('/')[-1].split('.')[0].split('_')[0:-1])
else:
    ind = options.ind

#define plotting space
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,300)
chrX = options.chrX
if chrX:
  ax.set_ylim(24,0)
else:
  ax.set_ylim(23,0)
plt.xlabel('Genetic position (cM)')
plt.ylabel('Chromosome')
plt.title(ind)
if options.chrX:
  plt.yticks(range(1,24))
  yticks = range(1,23)
  yticks.append('X')
  ax.set_yticklabels(yticks)
else:
  plt.yticks(range(1,23))

#define colors
bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
colors=bmap.mpl_colors
colors.append((0,0,0))

#define centromeres
centro = open(options.centromeres)
centromeres = {}
for line in centro:
    line = line.strip().split()
    if chrX and line[0] == 'X':
      line[0] = '23'
    centromeres[line[0]] = line

#plot rectangles
for line in bed_a:
    line = line.strip().split()
    try:
      plot_rects(line[3], int(line[0]), line[4], line[5], 'A', pop_order, colors, chrX)
    except ValueError: #flexibility for chrX
      plot_rects(line[3], 23, line[4], line[5], 'A', pop_order, colors, chrX)
for line in bed_b:
    line = line.strip().split()
    try:
      plot_rects(line[3], int(line[0]), line[4], line[5], 'B', pop_order, colors, chrX)
    except ValueError: #flexibility for chrX
      plot_rects(line[3], 23, line[4], line[5], 'B', pop_order, colors, chrX)

#write a legend
p = []
for i in range(len(pop_order)):
    p.append(plt.Rectangle((0, 0), 1, 1, color=colors[i]))
p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
labs = list(pop_order)
labs.append('UNK')
leg = ax.legend(p, labs, loc=4, fancybox=True)
leg.get_frame().set_alpha(0)

#get rid of annoying plot features
spines_to_remove = ['top', 'right']
for spine in spines_to_remove:
    ax.spines[spine].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

fig.savefig(options.out)