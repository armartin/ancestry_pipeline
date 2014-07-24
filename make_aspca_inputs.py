#takes in rephased alleles rfmix Viterbi file and outputs ASPCA inputs files
from optparse import  OptionParser
import os

USAGE = """
intersect_snps.py   --alleles
                    --vit
                    --haps
                    --sample
                    --classes
                    --out
"""
parser = OptionParser(USAGE)

parser.add_option('--alleles', default='/home/armartin/gsfs0/sa_analysis/local_ancestry/rfmix_output/SA_550_CEU_LWK_SA_phase3_chr22.rfmix.allelesRephased5.txt')
parser.add_option('--vit', default='/home/armartin/gsfs0/sa_analysis/local_ancestry/rfmix_output/SA_550_CEU_LWK_SA_phase3_chr22.rfmix.5.Viterbi.txt')
parser.add_option('--fbk_threshold', default=0.99, type='float')
parser.add_option('--mono_class', default='3')
parser.add_option('--keep_unrel', default='/home/armartin/sa_analysis/pedigree/SA.max_unrel')
parser.add_option('--haps', default='/home/armartin/sa_analysis/phasing/SA_550_sites_CEU_LWK_phase3_chr22.haps')
parser.add_option('--sample', default='/home/armartin/sa_analysis/local_ancestry/hmp3_ref/rfmix_input/SA_550_CEU_LWK_SA_phase3.sample')
parser.add_option('--classes', default='/home/armartin/sa_analysis/local_ancestry/hmp3_ref/rfmix_input/SA_550_CEU_LWK_SA_phase3.classes')
parser.add_option('--out', default='/home/armartin/gsfs0/sa_analysis/aspca/SA_550_CEU_LWK_SA_phase3_chr22')

(options, args) = parser.parse_args()

#find admixed versus reference panel individuals
ind_order = []
all_inds = []
sample = open(options.sample)
classes = open(options.classes)
classes = classes.readline().strip().split()
classes = classes[0::2] #get odds
print len(classes)
print classes
i = 0
ind_class = {}

keep_unrel = set()#
unrel = open(options.keep_unrel)#
for line in unrel:#
    keep_unrel.add(line.strip())#

for line in sample:
    line = line.strip()
    all_inds.append(line)
    if (line.startswith('SA') and line in keep_unrel) or not line.startswith('SA'):#this might mess up a lot of things
        ind_order.append(line)
    if classes[i] == '3' and line.startswith('SA'):
        ind_class[line] = '0'
    else:
        ind_class[line] = classes[i]
    i += 1

print keep_unrel#
print ind_order#
print ind_class
#write beagle header files for admixed and ancestral beagle files
out_anc = open(options.out + '_anc.beagle', 'w')
out_adm = open(options.out + '_adm.beagle', 'w')
out_anc.write('I\tid\t')
out_adm.write('I\tid\t')
vit_all = {}
vit_prob = {}
for ind in ind_order:
    if (ind.startswith('SA') and ind in ind_order) or not ind.startswith('SA'):#
        for hap in ['_A', '_B']:
            vit_all[ind + hap] = []
            vit_prob[ind + hap] = {}
        if ind_class[ind] == '0':
            out_adm.write(ind + '_A\t' + ind + '_B\t')
        else:
            out_anc.write(ind + '_A\t' + ind + '_B\t')
out_adm.write('\n')
out_anc.write('\n')

def check_anc(anc): #TO DO: take in a number corresponding to the ancestral class of interest
    sum_anc = 0
    for i in range(2,len(anc)):
        sum_anc += int(anc[i])
    if sum_anc == (len(anc) - 2) or sum_anc == 0:
        #print 'monomorphic'
        #print anc
        return False
    else:
        return True

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

#write all markers that are not monomorphic in reference panel
haps = open(options.haps)
alleles = open(options.alleles)
vit = open(options.vit)
fbk = open(options.vit.replace('Viterbi.txt', 'ForwardBackward.txt'))
out_vit = open(options.out + '.vit', 'w')
out_markers = open(options.out + '.markers', 'w')
marker_count = 0
for line in haps:
    line = line.strip().split()
    markers = alleles.readline().strip().split()
    vit_line = vit.readline().strip().split()
    fbk_line = fbk.readline().strip().split()
    fbk_chunked = chunker(fbk_line, 3)
    max_chunks = []
    for chunk in fbk_chunked:
        max_chunks.append(max(map(float, chunk)))
    #print line
    #print markers
    anc = ['M', line[2]]
    adm = ['M', line[2]]
    i=0
    #print markers[0]
    #print vit_line
    anc_mono = []
    for ind in all_inds: #this would be a problem because indexing is off
        if ind in ind_order:
            if ind_class[ind] == '0':
            #if ind_class[ind] == '0' and ind in keep_unrel:
                adm.append(markers[0][i])
                adm.append(markers[0][i+1])
            else:
            #elif ind.startswith('SA') and ind in keep_unrel:
                anc.append(markers[0][i])
                anc.append(markers[0][i+1])
                if options.mono_class is not None and ind_class[ind] == options.mono_class:
                    anc_mono.append(markers[0][i])
                    anc_mono.append(markers[0][i+1])
        i += 2
    i=0
    multimorphic = (check_anc(anc_mono) if anc_mono != [] else check_anc(anc))
    if multimorphic: #by default, check all for monomorphic, alternative, check subset for monomorphic
        for ind in all_inds:
            if ind in ind_order:
                #print ind + ' ' + vit_line[i]
                if max_chunks[i] > options.fbk_threshold:
                    vit_all[ind + '_A'].append(vit_line[i])
                else:
                    vit_all[ind + '_A'].append('-9')
                if max_chunks[i+1] > options.fbk_threshold:
                    vit_all[ind + '_B'].append(vit_line[i+1])
                else:
                    vit_all[ind + '_B'].append('-9')
            i += 2
        out_adm.write('\t'.join(adm) + '\n')
        out_anc.write('\t'.join(anc) + '\n')
        out_markers.write('window' + str(marker_count) + '\t' + line[2] + '\n')
        marker_count += 1
    
for ind in ind_order:
    out_vit.write(ind + '_A ' + ' '.join(vit_all[ind + '_A']) + '\n')
    out_vit.write(ind + '_B ' + ' '.join(vit_all[ind + '_B']) + '\n')

out_adm.close()
out_anc.close()
out_vit.close()