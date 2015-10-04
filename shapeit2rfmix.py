"""
sets up RFMix (alleles, snp_locations, and classes) files from shapeit output
to do: test that triallelics and strand flips are dealt with properly
admixed_keep can be None if keeping all individuals
"""
__author__ = 'armartin'

import argparse
from datetime import datetime
import time
import gzip

def open_shapeit(filename):
    """
    open haps and sample files
    """
    if filename.endswith('gz'):
        shapeit_file = gzip.open(filename)
    else:
        shapeit_file = open(filename)
    return(shapeit_file)


def fill_genos(hap, pos_id):
    """
    fill genos variables with info from each shapeit hap files
    assumes only one SNP per position
    """
    genos = {}
    for line in hap:
        myLine = line.strip().split()
        pos_id[int(myLine[2])] = myLine[1]
        genos[int(myLine[2])] = [myLine[3], myLine[4]]
    return(genos, pos_id)


def find_indices(keep_set, shapeit_sample, out_sample):
    """
    get indices corresponding with reference individuals to keep
    note: index order corresponds with the input keep file.
    i.e. they are output in same order as ref_keep and admixed_keep files
    """
    header0 = shapeit_sample.readline().strip().split()
    header1 = shapeit_sample.readline()
    if header0 != ['ID_1', 'ID_2', 'missing', 'father', 'mother', 'sex', 'plink_pheno']:
        raise RuntimeError('Shapeit sample file appears to be incorrect')
    indices = []
    ind_order = []
    for line in shapeit_sample:
        myLine = line.strip().split()
        ind_order.append(myLine[1])
    for ind in keep_set:
        if ind in ind_order:
            indices.extend([5+2*ind_order.index(ind), 6+2*ind_order.index(ind)])
            out_sample.write(ind + '\n')
    return indices

def get_haps(shapeit_hap, indices, full_intersection):
    """
    uses sample indices to get only the haplotypes of interest
    """
    hap_dict = {}
    for line in shapeit_hap:
        myLine = line.strip().split()
        if int(myLine[2]) in full_intersection:
            hap_dict[int(myLine[2])] = [myLine[i] for i in indices] #this is the problematic line for the admixed groups
    return hap_dict

def write_or_flip(snp_alleles, snp_haps, allele_order, file):
    """
    need to go through the intersection ordered sites and make sure the ref/alt alleles are the exact same between all groups!
    this probably breaks for fixed differences
    TO DO: get rid of triallelics but allow strand flips and don't consider 0's in triallelics
    """
    if len(allele_order) == 1: #should only happen if monomorphic across all files
        file.write(''.join(snp_haps))
    else: #TO DO change this to elif len()  = 2 else exception
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

def main(args):
    ref_keep = open(args.ref_keep)
    #find inds to keep in reference set
    ref_set = []
    for line in ref_keep:
        ref_set.append(line.strip())
    
    #require admixed keep file
    admixed_set = []
    if args.admixed_keep is not None:
        admixed_keep = open(args.admixed_keep)
        for line in admixed_keep:
            admixed_set.append(line.strip())
    
    ref_haps = args.shapeit_hap_ref.split(',')
    ref_samples = args.shapeit_sample_ref.split(',')
    
    #fill sites and genos variables with info from each shapeit file
    print 'Filling sites and genos from shapeit haps files [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    pos_id = {}
    genos = []
    for hap_file in ref_haps:
        sites_vars = fill_genos(open_shapeit(hap_file), pos_id)
        genos.append(sites_vars[0])
        pos_id = sites_vars[1]
    admixed_vars = fill_genos(open_shapeit(args.shapeit_hap_admixed), pos_id)
    genos.append(admixed_vars[0])
    pos_id = admixed_vars[1]
    
    out_sample = open(args.out + '.sample', 'w')
    
    sample_indices = []
    sample_admixed_indices = []
    for sample in ref_samples:
        sample_indices.append(find_indices(ref_set, open_shapeit(sample), out_sample))
    if args.admixed_keep is not None:
        sample_admixed_indices = find_indices(admixed_set, open_shapeit(args.shapeit_sample_admixed), out_sample)
    
    #get intersection between all reference panel and admixed sites
    print 'Getting intersection between reference and admixed sites [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    full_intersection = set(genos[0].keys())
    for site in range(1,len(genos)):
        full_intersection = set(genos[site].keys()) & full_intersection
        
    intersection_ordered = sorted(list(full_intersection))
    
    #save all genetic positions that are known in a recombination map in a dict
    genetic_map = open(args.genetic_map)
    pos_gen_map = {}
    header = genetic_map.readline().strip().split()
    if len(header) != 3:
        raise RuntimeError('Genetic map file does not have expected columns')
    for line in genetic_map:
        myLine = line.strip().split()
        pos_gen_map[int(myLine[0])] = float(myLine[2])

    #write physical positions with unknown genetic positions in a dict mapping to None to interpolate
    union_pos = sorted(list(full_intersection | set(pos_gen_map.keys())))
    for snp in (full_intersection - set(pos_gen_map.keys())):
        pos_gen_map[snp] = None
    
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
    
    alleles = open(args.out + '_chr' + args.chr + '.alleles', 'w')
    snps = open(args.out + '_chr' + args.chr + '.snp_locations', 'w')
    my_map = open(args.out + '_chr' + args.chr + '.map', 'w') #maps genetic to physical position
    classes = open(args.out + '.classes', 'w')
    
    print 'Writing snp_locations file [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    for snp in intersection_ordered:
        snps.write(str(pos_gen_map[snp]) + '\n')
        my_map.write(str(snp) + '\t' + str(pos_gen_map[snp]) + '\t' + pos_id[snp] + '\n')
    
    print 'Writing classes file [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    for index in range(len(sample_indices)):
        for i in sample_indices[index]:
            classes.write(str(index + 1) + ' ')

    if args.admixed_keep is None:
        for ind in admixed_set: ###changed this quite a bit, make sure it's ok
            classes.write('0 0 ')
    else:
        for hap in range(len(sample_admixed_indices)):
            classes.write('0 ')
    classes.write('\n')
    
    
    #dicts mapping position -> haplotype
    print 'Getting haplotypes of interest [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    all_haps = [] ######check that this is ok
    for i,hap in enumerate(ref_haps): # open all haps files
        all_haps.append(get_haps(open_shapeit(hap), sample_indices[i], full_intersection))
    all_haps.append(get_haps(open_shapeit(args.shapeit_hap_admixed), sample_admixed_indices, full_intersection))
    
    #checks that the intersection of ref_keep and ref_sample + admixed_keep and admixed_sample is equal to the number of haplotypes I will print
    assert sum([len(x[intersection_ordered[0]]) for x in all_haps]) == sum([len(x) for x in sample_indices]) + len(sample_admixed_indices)
    
    
    print 'Writing alleles file (total sites=' + str(len(intersection_ordered)) + ') [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    for snp in intersection_ordered:
        allele_pairs = [tuple(hap[snp]) for hap in genos]
        if len(set(allele_pairs)) == 1: #preferentially choose 2 alleles if bialleleic, if only monomorphic, choose that
            if '0' in allele_pairs[0]: #if all monomorphic, remove 2nd allele
                allele_pairs[0].remove('0')
            allele_order = list(allele_pairs[0])
        else:
            ap_boolean = ['0' in x for x in allele_pairs]
            if all(ap_boolean): #all have monomorphic snps
                allele_order = allele_pairs[0]
            else: #choose first non-monomorphic snps
                allele_order = allele_pairs[ap_boolean.index(False)]
        
        for i in range(len(genos)):
            write_or_flip(genos[i][snp], all_haps[i][snp], allele_order, alleles)
            
        alleles.write('\n')
    
    print 'Done! [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    alleles.close()
    snps.close()
    classes.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--shapeit_hap_ref', help='can be a single filename or comma-separated list if using multiple phased haplotypes', required=True)
    parser.add_argument('--shapeit_hap_admixed', required=True)
    
    parser.add_argument('--shapeit_sample_ref', required=True)
    parser.add_argument('--shapeit_sample_admixed', required=True)
    
    parser.add_argument('--ref_keep', help='a list of individual IDs', required=True)
    parser.add_argument('--admixed_keep', help='a list of individual IDs', required=True)
    
    parser.add_argument('--chr', required=True)
    
    parser.add_argument('--genetic_map', help='3 column file: physical position, recombination rate, genetic position', required=True)
    parser.add_argument('--out', required=True)

    args = parser.parse_args()
    main(args)
    

