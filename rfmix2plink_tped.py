#takes in rfmix output and writes ancestry plink tped files for admixture mapping, etc
__author__='armartin'
import argparse
import re

def main(args):
    chrs = map(str, range(1,23))
    inds = open(args.fam)
    ind_order = []
    for ind in inds:
        ind_order.append(ind.strip())
        
    pop_labels = args.pop_labels.split(',')
    out_fam = []
    for pop in pop_labels:
        out_fam.append(open(args.out + '_' + pop + '.tfam', 'w'))
    
    out_tped = []
    for pop in pop_labels:
        out_tped.append(open(args.out + '_' + pop + '.tped', 'w'))
    
    #go through every chromosome and population, printing tped along the way. counting number of tracts of each ancestry: G = absent, A = present
    for i in chrs:
        rfmix = open(re.sub(r'chr[X0-9]+', 'chr' + str(i), args.rfmix))
        snp_map = open(re.sub(r'chr[X0-9]+', 'chr' + str(i), args.snp_map))
        for line in snp_map:
            rf_line = rfmix.readline().strip().split()
            snp_line = line.strip().split()
            for out in out_tped:
                out.write(' '.join([i, snp_line[2], snp_line[1], snp_line[0]]) + ' ')
            for j in range(len(rf_line)/2):
                current_anc = [rf_line[2*j], rf_line[2*j+1]]
                for pop in range(len(pop_labels)):
                    pop_count = current_anc.count(str(pop+1))
                    if pop_count == 0:
                        out_tped[pop].write('G G ')
                    elif pop_count == 1:
                        out_tped[pop].write('G A ')
                    else:
                        out_tped[pop].write('A A ')
            for pop in range(len(pop_labels)):
                out_tped[pop].write('\n')

    for pop in range(len(pop_labels)):
        out_tped[pop].close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rfmix')
    parser.add_argument('--fam')
    parser.add_argument('--pop_labels')
    parser.add_argument('--snp_map')
    parser.add_argument('--out')
    args = parser.parse_args()
    main(args)