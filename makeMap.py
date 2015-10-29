"""
Interpolates genetic positions from a genetic map given a map or bim plink file
"""

__author__ = 'armartin'
import argparse
import sys

#sys.setrecursionlimit(1500) #this can really crash things. careful. default on broad node was 1000

def main():
    parser = argparse.ArgumentParser(description='Parse some args')
    parser.add_argument('--chr', default='22') #will replace chr[\d] with chr[chr]
    parser.add_argument('--genmap', default='/home/unix/armartin/atgu/shared_resources/recombination_maps/recombination_rates_hapmap_b37/genetic_map_chr22_combined_b37.txt')
    parser.add_argument('--bim', default='/psych/genetics_data/armartin/dbs/haplotypes/plink.dbs_db10_eur_rw-qc.hg19.ch.fl.chr22_.bim')
    parser.add_argument('--map_bim', default='bim')
    parser.add_argument('--out', default='/psych/genetics_data/armartin/dbs/haplotypes/plink.dbs_db10_eur_rw-qc.hg19.ch.fl.chr22_.map')
    
    args = parser.parse_args()
    
    genmap = open(args.genmap)
    genmap.readline()
    #position COMBINED_rate(cM/Mb) Genetic_Map(cM)
    #72765 0.1245577896 0
    start = genmap.readline().strip().split()
    (start_bp, start_cM) = (int(start[0]), float(start[2]))
    end = genmap.readline().strip().split()
    (end_bp, end_cM) = (int(end[0]), float(end[2]))
    
    bim = open(args.bim)
    if args.out is not None:
        my_map = open(args.out, 'w')
    else:
        my_map = open(args.bim.replace('bim', 'map'), 'w')
    
    bim_line = bim.readline().strip().split()
    chr = args.chr
    print chr
    
    (rsid, phys_pos, a0, a1) = (bim_line[1], int(bim_line[3]), bim_line[4], bim_line[5])
    while phys_pos < start_bp:
        proportion = (float(phys_pos) * float(start_cM)) / float(start_bp)
        write_map(my_map, [chr, rsid, str(proportion), phys_pos, a0, a1])
        bim_line = bim.readline().strip().split()
        (rsid, phys_pos, a0, a1) = (bim_line[1], bim_line[3], bim_line[4], bim_line[5])

    current_args = [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, args.chr, a0, a1]
    current_args = check_conditions(*current_args)
    
    for bim_line in bim:
        bim_line = bim_line.strip().split()
        (rsid, phys_pos, a0, a1) = (bim_line[1], bim_line[3], bim_line[4], bim_line[5])
        try:
            (current_args[0], current_args[3], current_args[10], current_args[11]) = (phys_pos, rsid, a0, a1)
        except TypeError:
            print current_args
            print [phys_pos, rsid, a0, a1]
        current_args = check_conditions(*current_args)
    
    my_map.close()

#don't start genetic positions quite at 0 because this throws off program (e.g. hapi-ur) assumptions
def write_map(my_map, write_vars):
    write_vars[2] = str(write_vars[2])
    if str(write_vars[2]) == '0.0':
        write_vars[2] = 1e-4
    my_map.write('\t'.join(map(str, write_vars)) + '\n')

def check_conditions(phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1):
    my_var = True
    my_map.flush()
    while my_var:
        if int(phys_pos) > int(end_bp):
            #print 'Criteria 1 - genotypes ahead of genetic map' #chr19 getting to this point but not subloops... ?
            while int(phys_pos) > int(end_bp):
                (start_bp, start_cM) = (end_bp, end_cM)
                end = genmap.readline().strip().split() #make sure this can happen
                if end == []: #might there be many snps meeting this criteria?
                    #print 'Criteria 1a'
                    ###Note: this one hasn't been checked yet
                    proportion = (float(phys_pos) * float(end_cM)) / float(end_bp)
                    write_map(my_map, [chr, rsid, str(proportion), phys_pos, a0, a1])
                    my_var = False
                    return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
                    #break
                else:
                    #print 'Criteria 1b'
                    (end_bp, end_cM) = (end[0], end[2])
                    #return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1] #this breaks things
        else:
            #print 'Criteria 2 - genotypes not ahead of genetic map'
            if phys_pos == start_bp:
                #print 'Criteria 2a'
                write_map(my_map, [chr, rsid, str(start_cM), phys_pos, a0, a1])
                return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
            elif phys_pos == end_bp:
                #print 'Criteria 2b'
                #print [chr, rsid, end_cM, phys_pos, a0, a1]
                write_map(my_map, [chr, rsid, str(end_cM), phys_pos, a0, a1])
                return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
            elif int(phys_pos) > int(start_bp) and int(phys_pos) < int(end_bp):
                #print 'Criteria 2c'
                interpolate = float((int(phys_pos) - int(start_bp))*(float(end_cM) - float(start_cM)))/(int(end_bp) - int(start_bp)) + float(start_cM)
                write_map(my_map, [chr, rsid, str(interpolate), phys_pos, a0, a1])
                return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
            else:
                print 'Criteria 2d'
                #this will happen if the first genotype in the bim file is before the first position in the genetic map
            #    print [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
            ##########!!!!!!!! Something wrong when genetic data is before physical positions
            #elif phys_pos < start_bp:
            #    proportion = (float(phys_pos) * float(start_cM)) / float(start_bp)
            #    my_map.write('\t'.join([chr, rsid, str(proportion), phys_pos, a0, a1]) + '\n')
            #    return [phys_pos, start_bp, end_bp, rsid, bim, start_cM, end_cM, genmap, my_map, chr, a0, a1]
            break
            
if __name__ == '__main__':
    main()
