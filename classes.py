#takes in a sample file output by shapeit2rfmix and writes classes
import argparse
import os

def read_ref(anc):
    anc_set = set()
    anc = open(anc)
    for line in anc:
        try:
            anc_set.add(line.strip().split()[1]) #split?
        except IndexError:
            raise IOError('Input misspecified. 2nd column needs to correspond with individual ID.')
    return(anc_set)

def main(args):
    ref = args.ref.strip().split(',')
    
    ancs = []
    for anc in ref:
        ancs.append(read_ref(anc))

    ind_order = []
    out = open(args.out, 'w')
    sample = open(args.sample)
    for line in sample:
        line = line.strip()
        if line == 'ID_1 ID_2 missing father mother sex plink_pheno':
            raise IOError('sample file must be list of individual IDs output by shapeit2rfmix with order of inds in alleles file, not shapeit sample file.')
        ind_order.append(line)
        in_ref = 0
        for anc in range(len(ancs)):
            in_ref = in_ref + int(line in ancs[anc])
            if line in ancs[anc]:
                out.write(str(anc + 1) + ' ' + str(anc + 1) + ' ')

        if in_ref == 0:
            out.write('0 0 ')
    out.write('\n')
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--ref', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--out', required=True)
    
    args = parser.parse_args()

    main(args)