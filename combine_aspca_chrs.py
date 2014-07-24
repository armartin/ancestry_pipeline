#takes in vit, ancestral, marker, and admixed aspca files split by chromosome and aggregates them
from optparse import  OptionParser
import os

USAGE = """
combine_aspca_chrs.py   --aspca_prefix
                        --keep_anc
                        --out
"""
parser = OptionParser(USAGE)

parser.add_option('--aspca_prefix', default='/home/armartin/gsfs0/sa_analysis/aspca/SA_Omni_CEU_LWK_SA_phase3_chr')
parser.add_option('--keep_anc', default='/home/armartin/sa_analysis/local_ancestry/hmp3_ref/CEU_all.inds')
parser.add_option('--anc', default='eur')
parser.add_option('--out', default='/home/armartin/gsfs0/sa_analysis/aspca/SA_Omni_CEU_LWK_SA_phase3')

(options, args) = parser.parse_args()

chrs = map(str, range(1,23))

keep_inds = []
if options.keep_anc != 'None':
    keep = open(options.keep_anc)
    for line in keep:
        line = line.strip()
        keep_inds.append(line)

print keep_inds

anc = open(options.aspca_prefix + '1_anc.beagle')
out_anc = open(options.out + '_' + options.anc + '.beagle', 'w')
anc_header = anc.readline().strip().split()
print anc_header
keep_index = [0, 1]
if options.keep_anc != 'None':
    out_anc.write('I\tid\t')
    for i in range(len(anc_header)):
        ind = '_'.join(anc_header[i].split('_')[:-1])
        if ind in keep_inds:
            keep_index.append(i)
            out_anc.write(anc_header[i] + '\t')
    out_anc.write('\n')
else:
    out_anc.write('\t'.join(anc_header) + '\n')

keep_index = set(keep_index)



adm = open(options.aspca_prefix + '1_adm.beagle')
out_adm = open(options.out + '_adm.beagle', 'w')
adm_header = adm.readline().strip().split()
out_adm.write('\t'.join(adm_header) + '\n')

vit_dict = {}
for ind in anc_header[2:len(anc_header)]:
    vit_dict[ind] = []
for ind in adm_header[2:len(adm_header)]:
    vit_dict[ind] = []

out_markers = open(options.out + '.markers', 'w')
out_vit = open(options.out + '.vit', 'w')

window_count = 0
for chr in chrs:
    current_anc = open(options.aspca_prefix + chr + '_anc.beagle')
    current_adm = open(options.aspca_prefix + chr + '_adm.beagle')
    current_anc.readline()
    current_adm.readline()
    for line in current_anc:
        if len(keep_inds) > 1:
            i = 0
            line = line.strip().split()
            for hap in line:
                if i in keep_index:
                    out_anc.write(hap + '\t')
                i += 1
            out_anc.write('\n')
        else:
            line = line.strip().split()
            out_anc.write('\t'.join(line) + '\n')
    for line in current_adm:
        line = line.strip().split()
        out_adm.write('\t'.join(line) + '\n')
    current_vit = open(options.aspca_prefix + chr + '.vit')
    for line in current_vit:
        line = line.strip().split()
        vit_dict[line[0]].extend(line[1:len(line)])
    current_markers = open(options.aspca_prefix + chr + '.markers')
    for line in current_markers: #might need to recheck if monomorphic
        line = line.strip().split()
        out_markers.write('window' + str(window_count) + '\t' + line[1] + '\n')
        window_count += 1

out_anc.close()
out_adm.close()
out_markers.close()

#for ind in anc_header[2::2]:
#    ind = '_'.join(ind.split('_')[:-1])
#    out_vit.write(ind + '\t' + '\t'.join(vit_dict[ind]) + '\n')
for ind in adm_header[2:]:
    out_vit.write(ind + '\t' + '\t'.join(vit_dict[ind]) + '\n')

out_vit.close()
