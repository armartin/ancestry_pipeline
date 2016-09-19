#takes in vit, ancestral, marker, and admixed aspca files split by chromosome and aggregates them
import argparse
import os

def main(args):
    chrs = map(str, range(1,23))
    
    keep_inds = []
    if args.keep_anc is not None:
        keep = open(args.keep_anc)
        for line in keep:
            line = line.strip()
            keep_inds.append(line)
    
    if args.extract is not None:
        keep_snps = set()
        extract = open(args.extract)
        for line in extract:
            keep_snps.add(line.strip())
            
    anc = open(args.aspca_prefix + '1_anc.beagle')
    out_anc = open(args.out + '_' + args.anc + '.beagle', 'w')
    anc_header = anc.readline().strip().split()
    keep_index = [0, 1]
    if args.keep_anc is not None:
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
    
    
    
    adm = open(args.aspca_prefix + '1_adm.beagle')
    out_adm = open(args.out + '_adm.beagle', 'w')
    adm_header = adm.readline().strip().split()
    out_adm.write('\t'.join(adm_header) + '\n')
    
    vit_dict = {}
    for ind in anc_header[2:len(anc_header)]:
        vit_dict[ind] = []
    for ind in adm_header[2:len(adm_header)]:
        vit_dict[ind] = []
    
    out_markers = open(args.out + '.markers', 'w')
    out_vit = open(args.out + '.vit', 'w')
    
    window_count = 0
    for chr in chrs:
        current_anc = open(args.aspca_prefix + chr + '_anc.beagle')
        current_adm = open(args.aspca_prefix + chr + '_adm.beagle')
        current_anc.readline()
        current_adm.readline()
        count_marker = 1
        marker_indices = []
        for line in current_anc:
            line = line.strip().split()
            if args.extract is not None and line[1] in keep_snps or args.extract is None:
                marker_indices.extend([count_marker, count_marker + 1])
                if len(keep_inds) > 1:
                    i = 0
                    for hap in line:
                        if i in keep_index:
                            out_anc.write(hap + ' ')
                        i += 1
                    out_anc.write('\n')
                else:
                    out_anc.write(' '.join(line) + '\n')
            count_marker += 2
        for line in current_adm:
            line = line.strip().split()
            if args.extract is not None and line[1] in keep_snps or args.extract is None:
                out_adm.write(' '.join(line) + '\n')
        current_vit = open(args.aspca_prefix + chr + '.vit')
        for line in current_vit: #check if there are markers to extract
            line = line.strip().split()
            if args.extract is not None:
                vit_dict[line[0]].extend([line[x] for x in marker_indices])
            else:
                vit_dict[line[0]].extend(line[1:len(line)])
        current_markers = open(args.aspca_prefix + chr + '.markers')
        for line in current_markers: #might need to recheck if monomorphic
            line = line.strip().split()
            if args.extract is not None and line[1] in keep_snps or args.extract is None:
                out_markers.write('window' + str(window_count) + '\t' + line[1] + '\n')
                window_count += 1
    
    out_anc.close()
    out_adm.close()
    out_markers.close()
    
    for ind in adm_header[2:]:
        out_vit.write(ind + ' ' + ' '.join(vit_dict[ind]) + '\n')
    
    out_vit.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--aspca_prefix')
    parser.add_argument('--keep_anc')
    parser.add_argument('--anc')
    parser.add_argument('--extract')
    parser.add_argument('--out')
    args = parser.parse_args()
    main(args)