#sets up RFMix (alleles, snp_locations, and classes) files from shapeit output
from optparse import OptionParser
from datetime import datetime
import time
import gzip

parser = OptionParser()
parser.add_option('--shapeit_hap1', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/shapeit_out/YRI_chr20.haps')
parser.add_option('--shapeit_hap2', default=None)
parser.add_option('--shapeit_hap3', default=None)
parser.add_option('--shapeit_hap4', default=None)
parser.add_option('--shapeit_hap_admixed', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/shapeit_out/ASW_chr20.haps')

parser.add_option('--shapeit_sample1', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/shapeit_out/YRI_chr20.sample')
parser.add_option('--shapeit_sample2', default=None)
parser.add_option('--shapeit_sample3', default=None)
parser.add_option('--shapeit_sample4', default=None)
parser.add_option('--shapeit_sample_admixed', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/shapeit_out/ASW_chr20.sample')

parser.add_option('--ref_keep', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/shapeit_in/CEU+YRI.trios') #a list of individual IDs
parser.add_option('--admixed_keep', default=None) #a list of individual IDs
#could probably write an option in here to keep admixed individuals for ease of use with other datasets

parser.add_option('--chr', default='20')

parser.add_option('--genetic_map', default='/home/armartin/bustamante/reference_panels/recombination_rates_hapmap_b37/genetic_map_chr')
parser.add_option('--out', default='/home/armartin/rare/chip_collab/admixed/ACB+ASW/rfmix_input/ASW')

(options,args)=parser.parse_args()

#open haps files
def open_shapeit():
    global shapeit_hap1
    if options.shapeit_hap1.endswith('gz'):
        shapeit_hap1 = gzip.open(options.shapeit_hap1)
    else:
        shapeit_hap1 = open(options.shapeit_hap1)
    global shapeit_hap2
    if options.shapeit_hap2 is not None:
        if options.shapeit_hap2.endswith('gz'):
            shapeit_hap2 = gzip.open(options.shapeit_hap2)
        else:
            shapeit_hap2 = open(options.shapeit_hap2)
    global shapeit_hap3
    if options.shapeit_hap3 is not None:
        if options.shapeit_hap3.endswith('gz'):
            shapeit_hap3 = gzip.open(options.shapeit_hap3)
        else:
            shapeit_hap3 = open(options.shapeit_hap3)
    global shapeit_hap4
    if options.shapeit_hap4 is not None:
        if options.shapeit_hap4.endswith('gz'):
            shapeit_hap4 = gzip.open(options.shapeit_hap4)
        else:
            shapeit_hap4 = open(options.shapeit_hap4)
    global shapeit_hap_admixed
    if options.shapeit_hap_admixed.endswith('gz'):
        shapeit_hap_admixed = gzip.open(options.shapeit_hap_admixed)
    else:
        shapeit_hap_admixed = open(options.shapeit_hap_admixed)

#fill sites and genos variables with info from each shapeit hap files
def fill_sites_genos(sites, genos, hap, pos_id):
    for line in hap:
        myLine = line.strip().split()
        pos_id[int(myLine[2])] = myLine[1]
        sites.add(int(myLine[2]))
        genos[int(myLine[2])] = [myLine[3], myLine[4]]

#get indices corresponding with reference individuals to keep
#note: index order corresponds with the input keep file. i.e. they are output in same order as ref_keep and admixed_keep files
def find_indices(ref_set, shapeit_sample):
    shapeit_sample.readline()
    shapeit_sample.readline() #get rid of header
    indices = []
    ind_order = []
    for line in shapeit_sample:
        myLine = line.strip().split()
        ind_order.append(myLine[1])
    for ind in ref_set:
        #print ind_order
        if ind in ind_order:
            indices.extend([5+2*ind_order.index(ind), 6+2*ind_order.index(ind)])
            out_sample.write(ind + '\n')
    return (indices, ind_order)

#uses sample indices to get only the haplotypes of interest
def get_haps(hap_dict, shapeit_hap, indices):
    for line in shapeit_hap:
        myLine = line.strip().split()
        if int(myLine[2]) in full_intersection:
            hap_dict[int(myLine[2])] = [myLine[i] for i in indices] #this is the problematic line for the admixed groups
    return hap_dict

#need to go through the intersection ordered sites and make sure the ref/alt alleles are the exact same between all groups!
def write_or_flip(snp_alleles, snp_haps, allele_order, file):
    if len(allele_order) == 1:
        file.write(''.join(snp_haps))
    else:
        if '0' in snp_alleles: #0 is always the 2nd SNP
            if snp_alleles[0] == '0':
                snp_alleles[0] = allele_order[abs(allele_order.index(snp_alleles[1])-1)]
            else:
                snp_alleles[1] = allele_order[abs(allele_order.index(snp_alleles[0])-1)]
                
            which_allele = allele_order.index(snp_alleles[0])
            if which_allele == 0:
                file.write(''.join(snp_haps))
            else: #need to flip
                flipped = map(lambda x: '0' if x == '1' else '1', snp_haps)
                file.write(''.join(flipped))
        elif snp_alleles == allele_order:
            file.write(''.join(snp_haps))
        else:
        #need to flip alleles before writing
            flipped = map(lambda x: '0' if x == '1' else '1', snp_haps)
            file.write(''.join(flipped))

#set up variables
hap1_sites = set()
hap1_genos = {}
hap2_sites = set()
hap2_genos = {}
hap3_sites = set()
hap3_genos = {}
hap4_sites = set()
hap4_genos = {}
hap_admixed_sites = set()
hap_admixed_genos = {}

#open shapeit haps and sample files as well as reference inds to keep
open_shapeit()
shapeit_sample1 = open(options.shapeit_sample1)
if options.shapeit_hap2 is not None:
    shapeit_sample2 = open(options.shapeit_sample2)
if options.shapeit_hap3 is not None:
    shapeit_sample3 = open(options.shapeit_sample3)
if options.shapeit_hap4 is not None:
    shapeit_sample4 = open(options.shapeit_sample4)
if options.admixed_keep is not None:
    shapeit_sample_admixed = open(options.shapeit_sample_admixed)
    admixed_keep = open(options.admixed_keep)
ref_keep = open(options.ref_keep)

#fill sites and genos variables with info from each shapeit file
print 'Filling sites and genos from shapeit haps files [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
global pos_id
pos_id = {}
fill_sites_genos(hap1_sites, hap1_genos, shapeit_hap1, pos_id)
if options.shapeit_hap2 is not None:
    fill_sites_genos(hap2_sites, hap2_genos, shapeit_hap2, pos_id)
if options.shapeit_hap3 is not None:
    fill_sites_genos(hap3_sites, hap3_genos, shapeit_hap3, pos_id)
if options.shapeit_hap4 is not None:
    fill_sites_genos(hap4_sites, hap4_genos, shapeit_hap4, pos_id)
fill_sites_genos(hap_admixed_sites, hap_admixed_genos, shapeit_hap_admixed, pos_id)

#open the genetic map
global genetic_map
genetic_map = open(options.genetic_map + options.chr + '_combined_b37.txt')

#find inds to keep in reference set
ref_set = []
for line in ref_keep:
    ref_set.append(line.strip())
#print ref_set
admixed_set = []
if options.admixed_keep is not None:
    for line in admixed_keep:
        admixed_set.append(line.strip())
#print admixed_set

out_sample = open(options.out + '.sample', 'w')

(sample1_indices,sample1_order) = find_indices(ref_set, shapeit_sample1)

if options.shapeit_hap2 is not None:
    (sample2_indices,sample2_order) = find_indices(ref_set, shapeit_sample2)
if options.shapeit_hap3 is not None:
    (sample3_indices,sample3_order) = find_indices(ref_set, shapeit_sample3)
if options.shapeit_hap4 is not None:
    (sample4_indices,sample4_order) = find_indices(ref_set, shapeit_sample4)
if admixed_keep is not None:
    (sample_admixed_indices,sample_admixed_order) = find_indices(admixed_set, shapeit_sample_admixed)
open_shapeit()

#get intersection between all reference panel and admixed sites
print 'Getting intersection between reference and admixed sites [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
full_intersection = hap1_sites & hap_admixed_sites
if options.shapeit_hap2 is not None:
    full_intersection = full_intersection & hap2_sites
if options.shapeit_hap3 is not None:
    full_intersection = full_intersection & hap3_sites
if options.shapeit_hap4 is not None:
    full_intersection = full_intersection & hap4_sites
    
for snp in list(full_intersection):
    all_genos = []
    all_genos.extend(hap1_genos[snp])
    if options.shapeit_hap2 is not None:
        all_genos.extend(hap2_genos[snp])
    if options.shapeit_hap3 is not None:
        all_genos.extend(hap3_genos[snp])
    if options.shapeit_hap4 is not None:
        all_genos.extend(hap4_genos[snp])
    all_genos.extend(hap_admixed_genos[snp])
    #get rid of sites that are not biallelic (but not monomorphic sites)
    if not (len(set(all_genos)) == 2 or (len(set(all_genos)) == 3 and '0' in set(all_genos))): #filter to biallelic sites
        full_intersection.remove(snp)
intersection_ordered = sorted(list(full_intersection))

#write all genetic positions that are known in a recombination map in a dict
pos_gen_map = {}
gen_pos = set()
genetic_map.readline()
for line in genetic_map:
    myLine = line.strip().split()
    pos_gen_map[int(myLine[0])] = float(myLine[2])
    gen_pos.add(int(myLine[0]))

#write physical positions with unknown genetic positions in a dict mapping to None to interpolate
union_pos = sorted(list(full_intersection | gen_pos))
for line in (full_intersection - gen_pos):
    pos_gen_map[line] = None

#get genetic positions for all sites using known genetic positions and interpolating when not possible
print 'Getting genetic positions through known genetic map and interpolation [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
current_position = 0
for position in range(len(union_pos)):
    if pos_gen_map[union_pos[position]] is not None: #we already have the genetic position
        continue
    elif position == 0 and pos_gen_map[union_pos[position]] is None: #we're at the beginning of the file, which doesn't have a genetic position
        next_position = position
        while(pos_gen_map[union_pos[next_position]]) is None:
            next_position += 1
        ####running into None issues at the next line. probably need a try except here, but first figure out why it's not getting the right index
        interpolated = float(pos_gen_map[union_pos[next_position]]*union_pos[position])/union_pos[next_position]
        pos_gen_map[union_pos[position]] = interpolated
    elif pos_gen_map[union_pos[position]] is None: #we're between known genetic positions or at the end
        next_position = position
        try:
            while(pos_gen_map[union_pos[next_position]]) is None:
                next_position += 1
            interpolated = (float(union_pos[position] - union_pos[position-1]) * (pos_gen_map[union_pos[next_position]] - pos_gen_map[union_pos[position-1]]) /
                float(union_pos[next_position] - union_pos[position-1])) + pos_gen_map[union_pos[position-1]]
            pos_gen_map[union_pos[position]] = interpolated
        except IndexError: #we're near the end, need to perform linear interpolation like we did at the beginning of the fil
            last_position = position
            while(pos_gen_map[union_pos[last_position]]) is None:
                last_position = last_position - 1
            interpolated = float(pos_gen_map[union_pos[last_position]]*union_pos[position])/union_pos[last_position]
            pos_gen_map[union_pos[position]] = interpolated

alleles = open(options.out + '_chr' + options.chr + '.alleles', 'w')
snps = open(options.out + '_chr' + options.chr + '.snp_locations', 'w')
my_map = open(options.out + '_chr' + options.chr + '.map', 'w') #maps genetic to physical position
classes = open(options.out + '.classes', 'w')

print 'Writing snp_locations file [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
for snp in intersection_ordered:
    snps.write(str(pos_gen_map[snp]) + '\n')
    my_map.write(str(snp) + '\t' + str(pos_gen_map[snp]) + '\t' + pos_id[snp] + '\n')

print 'Writing classes file [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
for hap in range(len(sample1_indices)):
    classes.write('1 ')
if options.shapeit_hap2 is not None:
    for hap in range(len(sample2_indices)):
        classes.write('2 ')
if options.shapeit_hap3 is not None:
    for hap in range(len(sample3_indices)):
        classes.write('3 ')
if options.shapeit_hap4 is not None:
    for hap in range(len(sample4_indices)):
        classes.write('4 ')
if options.admixed_keep is None:
    shapeit_sample_admixed = open(options.shapeit_sample_admixed)
    shapeit_sample_admixed.readline()
    shapeit_sample_admixed.readline()
    sample_admixed_indices = []
    current_line = 5
    for line in shapeit_sample_admixed:
        classes.write('0 0 ')
        myLine = line.strip().split()
        sample_admixed_indices.extend([current_line, current_line+1])
        current_line+=2
else:
    for hap in range(len(sample_admixed_indices)):
        classes.write('0 ')
classes.write('\n')

#dicts mapping position -> haplotype
print 'Getting haplotypes of interest [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
hap1_ref = {}
hap2_ref = {}
hap3_ref = {}
hap4_ref = {}
hap_ref_admixed = {}
hap1_ref = get_haps(hap1_ref, shapeit_hap1, sample1_indices)
if options.shapeit_hap2 is not None:
    hap2_ref = get_haps(hap2_ref, shapeit_hap2, sample2_indices)
if options.shapeit_hap3 is not None:
    hap3_ref = get_haps(hap3_ref, shapeit_hap3, sample3_indices)
if options.shapeit_hap4 is not None:
    hap4_ref = get_haps(hap4_ref, shapeit_hap4, sample4_indices)
hap_ref_admixed = get_haps(hap_ref_admixed, shapeit_hap_admixed, sample_admixed_indices)

print 'Writing alleles file (total sites=' + str(len(intersection_ordered)) + ') [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
for snp in intersection_ordered:
    if options.shapeit_hap3 is not None and options.shapeit_hap4 is None:
        if '0' not in hap1_genos[snp] and '0' not in hap2_genos[snp] and '0' not in hap3_genos[snp] and '0' not in hap_admixed_genos[snp]:#all pops have variable snps, so just pick an order.
            allele_order = hap1_genos[snp]
        else: #figure out what to do with 0's
            allele_order = set(hap1_genos[snp]) | set(hap2_genos[snp]) | set(hap3_genos[snp]) | set(hap_admixed_genos[snp])
            allele_order.remove('0')
            allele_order = list(allele_order)
        write_or_flip(hap1_genos[snp], hap1_ref[snp], allele_order, alleles)            
        write_or_flip(hap2_genos[snp], hap2_ref[snp], allele_order, alleles)
        write_or_flip(hap3_genos[snp], hap3_ref[snp], allele_order, alleles)
        write_or_flip(hap_admixed_genos[snp], hap_ref_admixed[snp], allele_order, alleles) #hap_ref_admixed is empty when it shouldn't be
    elif options.shapeit_hap3 is not None and options.shapeit_hap4 is not None:
        if '0' not in hap1_genos[snp] and '0' not in hap2_genos[snp] and '0' not in hap3_genos[snp] and '0' not in hap4_genos[snp] and '0' not in hap_admixed_genos[snp]:#all pops have variable snps, so just pick an order.
            allele_order = hap1_genos[snp]
        else: #figure out what to do with 0's
            allele_order = set(hap1_genos[snp]) | set(hap2_genos[snp]) | set(hap3_genos[snp]) | set(hap4_genos[snp]) | set(hap_admixed_genos[snp])
            allele_order.remove('0')
            allele_order = list(allele_order)
        write_or_flip(hap1_genos[snp], hap1_ref[snp], allele_order, alleles)
        write_or_flip(hap2_genos[snp], hap2_ref[snp], allele_order, alleles)
        write_or_flip(hap3_genos[snp], hap3_ref[snp], allele_order, alleles)
        write_or_flip(hap4_genos[snp], hap4_ref[snp], allele_order, alleles)
        write_or_flip(hap_admixed_genos[snp], hap_ref_admixed[snp], allele_order, alleles)
    else: #only 2 ref panel
        if options.shapeit_hap2 is None:
            if '0' not in hap1_genos[snp] and '0' not in hap_admixed_genos[snp]:#all pops have variable snps, so just pick an order.
                allele_order = hap1_genos[snp]
            else: #figure out what to do with 0's
                allele_order = set(hap1_genos[snp]) | set(hap_admixed_genos[snp])
                allele_order.remove('0')
                allele_order = list(allele_order)
            write_or_flip(hap1_genos[snp], hap1_ref[snp], allele_order, alleles)            
            write_or_flip(hap_admixed_genos[snp], hap_ref_admixed[snp], allele_order, alleles) #hap_ref_admixed is empty when it shouldn't be
        else:
            if '0' not in hap1_genos[snp] and '0' not in hap2_genos[snp] and '0' not in hap_admixed_genos[snp]:#all pops have variable snps, so just pick an order.
                allele_order = hap1_genos[snp]
            else: #figure out what to do with 0's
                allele_order = set(hap1_genos[snp]) | set(hap2_genos[snp]) | set(hap_admixed_genos[snp])
                allele_order.remove('0')
                allele_order = list(allele_order)
            write_or_flip(hap1_genos[snp], hap1_ref[snp], allele_order, alleles)
            write_or_flip(hap2_genos[snp], hap2_ref[snp], allele_order, alleles)
            write_or_flip(hap_admixed_genos[snp], hap_ref_admixed[snp], allele_order, alleles)
    alleles.write('\n')

print 'Done! [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
alleles.close()
snps.close()
classes.close()