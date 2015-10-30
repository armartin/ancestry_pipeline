"""
Plots the output of smartpca
"""

__author__ = 'armartin'
import matplotlib.pyplot as plt
import argparse
import brewer2mpl
import pandas as pd
import numpy as np
from matplotlib.font_manager import FontProperties
import matplotlib.cm as cm
from datetime import datetime
import time

parser = argparse.ArgumentParser(description='Parse some args')
parser.add_argument('--evec') #assumes eval in same place
parser.add_argument('--out', help='uses png or pdf suffix to determine type of file to plot')
parser.add_argument('--projected')
parser.add_argument('--title')
parser.add_argument('--noGrid')
parser.add_argument('--which_pcs', default='1,2', help='comma-separated list of PCs')

args = parser.parse_args()

evec = open(args.evec)
my_eval = open(args.evec.replace('evec', 'eval'))
which_pcs = map(int, args.which_pcs.split(','))
print ['PC' + str(p) + ' ' for p in which_pcs]

if args.out is None:
    out = args.evec.replace('.evec', '_pca' + str(which_pcs[0]) + '_' + str(which_pcs[1]) + '.png')
else:
    out = args.out
if not (args.out.endswith('png') or args.out.endswith('pdf')):
    raise ValueError('--out does not end in png or pdf')

eval_per = []
for line in my_eval:
    eval_per.append(float(line.strip()))
eval_per = [x/sum(eval_per)*100 for x in eval_per]

eigs = {}
evec.readline()
evec_all = pd.read_csv(evec, header=None, sep='\s+')
pcs = ['ID']
pcs.extend(['PC' + str(x) for x in range(1, evec_all.shape[1]-1)])
pcs.append('case')
evec_all.columns = pcs
iid_fid = pd.DataFrame(evec_all['ID'].str.split(':').tolist())
df = pd.concat([iid_fid, evec_all], axis=1)
df.rename(columns={0: 'FID', 1: 'IID'}, inplace=True)
projected = args.projected
#fix this
#if projected is not None:
#    new_df = df.loc[df.loc[:,0]==projected,:]
#    new_df2 = df.loc[df.loc[:,0]!=projected,:]
    #new_df3 = Frame.append(pandas.DataFrame(data=SomeNewLineOfData),ignore_index=True)

pops = list(pd.Series.unique(df.iloc[:,0]))
print pops

try:
    if projected is not None:
        colors = [(0, 0, 0)]
        bmap = (brewer2mpl.get_map('Set1', 'qualitative', len(pops)-1))
        colors.extend(bmap.mpl_colors)
    else:
        bmap = brewer2mpl.get_map('Set1', 'qualitative', len(pops))
        colors=bmap.mpl_colors
except ValueError: #too many colors for color brewer
    if projected is not None:
        colors = np.vstack((np.array([0, 0, 0, 1]), cm.rainbow(np.linspace(0, 1, len(pops)-1))))
    else:
        colors = cm.rainbow(np.linspace(0, 1, len(pops)))

fig = plt.figure()
ax = plt.subplot(111)
#pop_col = [colors[pops.index(x)] for x in df.iloc[:,0]]
markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'] #probably need to cycle through
while len(markers) < len(pops):
    markers.extend(markers)
#print markers
#pop_mark = [markers[pops.index(x)] for x in df.iloc[:,0]]

#for _pop_mark, _pop_col, _x, _y in zip(pop_mark, pop_col, df.iloc[:,which_pcs[0]+2], df.iloc[:,which_pcs[1]+2]):
#    ax.scatter(_x, _y, s=12, c=_pop_col, marker=_pop_mark, alpha=0.5, edgecolors='none', label=df.iloc[:,0])

for pop in range(len(pops)):
    print pops[pop] + ' [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    ax.scatter(df.ix[df.index[df['FID'] == pops[pop]],which_pcs[0]+2],
               df.ix[df.index[df['FID'] == pops[pop]],which_pcs[1]+2],
               label=df.ix[df.index[df['FID'] == pops[pop]],0],
               s=20, c=colors[pop], marker=markers[pop], alpha=0.5, edgecolors=colors[pop])
    
#ax.set_xlabel(list(df.columns.values)[which_pcs[0]+2] + ' (' + "%.1f" % eval_per[which_pcs[0]+2] + '% var explained)', fontsize=14)
ax.set_xlabel(list(df.columns.values)[which_pcs[0]+2], fontsize=14)
#ax.set_ylabel(list(df.columns.values)[which_pcs[1]+2] + ' (' + "%.1f" % eval_per[which_pcs[1]+2] + '% var explained)', fontsize=14)
ax.set_ylabel(list(df.columns.values)[which_pcs[1]+2], fontsize=14)
ax.set_title(args.title, fontsize=20)

if not args.noGrid:
    ax.grid(True, color='0.2', linestyle='--')
    ax.set_axisbelow(True)
box = ax.get_position()
p = []
for i in range(len(pops)):
    p.append(plt.Line2D(range(1), range(1), color='white', markerfacecolor=colors[i], marker=markers[i], markeredgecolor=colors[i]))
#if len(pops) < 7:
#    leg = ax.legend(p, pops, loc='best', fancybox=True, ncol=len(pops)/6, prop={'size':9}) #bbox_to_anchor=(0.5, 1),
#    leg.get_frame().set_alpha(0.5)
#else:
leg = ax.legend(p, pops, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':9}, numpoints=1) #, fancybox=True, ncol=len(pops)
#leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':9}, numpoints=1) #, fancybox=True, ncol=len(pops)
#fig.tight_layout()
plt.savefig(out, bbox_extra_artists=(leg,), bbox_inches='tight')