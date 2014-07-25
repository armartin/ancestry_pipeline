__author__ = 'armartin'
from optparse import  OptionParser

def splitstr(option, opt, value, parser):
  return(setattr(parser.values, option.dest, value.split(',')))

USAGE = """
lai_global.py   --bed_list
                --ind_list
                --pops
                --out
"""
parser = OptionParser(USAGE)

parser.add_option('--bed_list')
parser.add_option('--ind_list')
parser.add_option('--pops', default=['AFR','EUR','NAT','UNK'], type='string', action='callback', callback=splitstr,
                  help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
parser.add_option('--out')

(options, args) = parser.parse_args()

bed_list = open(options.bed_list)
ind_list = open(options.ind_list)
out = open(options.out, 'w')
pops = options.pops
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