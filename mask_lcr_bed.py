from optparse import  OptionParser
from collections import defaultdict

USAGE = """
mask_lcr_bed.py --bed
                --out
"""
parser = OptionParser(USAGE)

parser.add_option('--bed', default='/home/armartin/rare/chip_collab/admixed/affy6/lai_output/ACB/seq/BF001_A.bed')
parser.add_option('--lcr', default='/home/armartin/rare/chip_collab/admixed/affy6/lai_output/lcr_250kb.bed')
parser.add_option('--out', default='/home/armartin/rare/chip_collab/admixed/affy6/lai_output/ACB/seq/BF001_A_2.bed')

(options, args) = parser.parse_args()

lcr = open(options.lcr)
bed = open(options.bed)
out = open(options.out, 'w')

chrs = map(str, range(1,23))

lcr_dict = defaultdict(dict)
for line in lcr:
    line = line.strip().split()
    if line[4] != line[5]:
        lcr_dict[line[0]][int(line[1])] = line

print lcr_dict

out = open(options.out, 'w')
while True:
    line = bed.readline().strip().split()
    if not line: break
    if line[0] in lcr_dict:
        #print lcr_dict[line[0]]
        for lcr_start in lcr_dict[line[0]]:
            #lcr is the same tract
            if int(line[1]) == lcr_start and int(line[2]) == int(lcr_dict[line[0]][lcr_start][2]):
                out.write('\t'.join(lcr_dict[line[0]][lcr_start]) + '\n')
            #lcr is within a single tract
            elif int(line[1]) <= lcr_start and int(line[2]) >= int(lcr_dict[line[0]][lcr_start][2]):
                tract1 = list(line)
                tract1[2] = str(lcr_start)
                tract1[5] = lcr_dict[line[0]][lcr_start][4]
                tract2 = list(line)
                tract2[1] = lcr_dict[line[0]][lcr_start][2]
                tract2[4] = lcr_dict[line[0]][lcr_start][5]
                tract2[5] = line[5]
                #out.write('\t'.join(tract1) + '\n')
                out.write('\t'.join(lcr_dict[line[0]][lcr_start]) + '\n')
                out.write('\t'.join(tract2) + '\n')
            #lcr is within two tracts
            elif int(line[1]) <= lcr_start and int(line[2]) <= int(lcr_dict[line[0]][lcr_start][2]) and int(line[2]) >= lcr_start:
                tract1 = list(line)
                tract1[2] = str(lcr_start)
                tract1[5] = lcr_dict[line[0]][lcr_start][4]
                out.write('\t'.join(tract1) + '\n')
                out.write('\t'.join(lcr_dict[line[0]][lcr_start]) + '\n')
                #print tract1
                #print lcr_dict[line[0]][lcr_start]
                #while int(line[2]) < int(lcr_dict[line[0]][lcr_start][2]):
                line = bed.readline().strip().split()
                tract2 = list(line)
                tract2[1] = lcr_dict[line[0]][lcr_start][2]
                tract2[4] = lcr_dict[line[0]][lcr_start][5]
                tract2[5] = line[5]
                #print tract2
                out.write('\t'.join(tract1) + '\n')
                out.write('\t'.join(lcr_dict[line[0]][lcr_start]) + '\n')
                out.write('\t'.join(tract2) + '\n')
            else:
                out.write('\t'.join(line) + '\n')
            #tract ends in middle of lcr
    else:
        out.write('\t'.join(line) + '\n')