__author__ = 'armartin'
import argparse

USAGE = """
lai_global.py   --bed_list
                --ind_list
                --pops
                --out
"""
parser = argparse.ArgumentParser()

parser.add_argument('--bed_list')
parser.add_argument('--ind_list')
parser.add_argument('--pops', default='AFR,EUR,NAT,UNK',
                  help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
parser.add_argument('--out')

args = parser.parse_args()

bed_list = open(args.bed_list)
ind_list = open(args.ind_list)
out = open(args.out, 'w')
pops = args.pops
out.write('ID\t' + '\t'.join(pops) + '\n')
lai_props = [0]*len(pops)

for line in bed_list:
    line = line.strip().split()
    ind = ind_list.readline().strip()
    bed_a = open(line[0])
    bed_b = open(line[1])
    for tract in bed_a:
      tract = tract.strip().split()
      if tract[3] in pops: #this excludes stuff not listed in pops
        lai_props[pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))
    for tract in bed_b:
      tract = tract.strip().split()
      if tract[3] in pops: #this excludes stuff not listed in pops
        lai_props[pops.index(tract[3])] += (float(tract[5]) - float(tract[4]))
    
    out.write(ind + '\t' + '\t'.join(map(str, [round(i/sum(lai_props), 4) for i in lai_props])) + '\n')
    lai_props = [0]*len(pops)
        
out.close()