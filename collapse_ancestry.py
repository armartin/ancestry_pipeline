__author__ = 'armartin'
#takes in RFMix output and creates collapsed 2 haploid bed files per individual

from datetime import datetime
import time
import argparse
from itertools import izip_longest
import gzip
import re

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('--rfmix', help='path to RFMix Viterbi output, chr expected in filename', required=True)
  parser.add_argument('--snp_locations', help='path to snp_locations file required by RFMix, chr expected in filename', required=True)
  parser.add_argument('--fbk', default=None)
  parser.add_argument('--fbk_threshold', type=float, default = 0.9)
  parser.add_argument('--ind', help='Individual ID, must match a line in --ind_info option', required=True)
  parser.add_argument('--ind_info', help='Individual IDs listed in the order they appear in RFMix Viterbi output', required=True)
  parser.add_argument('--pop_labels', default='AFR,EUR,NAT',
                    help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
  parser.add_argument('--chrX', help='include chrX?', default=False, action="store_true")
  parser.add_argument('--out', help='prefix to bed file, _A.bed and _B.bed will be appended', required=True)
  
  args = parser.parse_args()
  return(args)
  
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def check_gt_posterior(fbk_max, fbk_threshold, index, add, hap_anc, line, current_anc):
  if fbk_max[index*2+add] >= fbk_threshold:
      hap_anc = myLine[index*2+add]
  else:
      hap_anc = -9
      current_anc[add] = -9
  return (current_anc, hap_anc)
  
def find_haplotype_bounds(index, add, pop_order, hap,npop):
  for chr in chrs:
    print str(chr) + ' [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    rfmix = open(re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.rfmix))
    if args.fbk is not None:
      fbk_file = re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.fbk)
      if fbk_file.endswith('gz'):
        fbk = gzip.open(fbk_file)
      else:
        fbk = open(fbk_file)
    snp_file = re.sub(r'chr[X0-9]+', 'chr' + str(chr), args.snp_locations)
    snp_locations = open(snp_file)
    snp_map = open(snp_file.replace('snp_locations', 'map')) #map of physical position -> genetic position
    
    last_anc_pos_cm = [None, None, 0]
    
    counter = 0
    for line in rfmix:
      counter += 1
      myLine = line.strip().split()
      my_pos = snp_locations.readline().strip()
      my_map = snp_map.readline().strip().split()
      if args.fbk is not None:
        fbk_line = fbk.readline().strip().split()
        fbk_line = map(float, fbk_line)
        fbk_max = []
        for i in grouper(npop, fbk_line):
          fbk_max.append(max(i))
        if fbk_max[index*2+add] < args.fbk_threshold:
          myLine[index*2+add] = -9
      #fencepost for start of the chromosome
      if counter == 1:
        last_anc_pos_cm = [myLine[2*index + add], my_map[0], my_pos]
        post_anc_pos_cm = [myLine[2*index + add], my_map[0], my_pos]
        continue
      
      #start regular iterations
      current_anc_pos_cm = [myLine[2*index + add], my_map[0], my_pos]
      if current_anc_pos_cm[0] == last_anc_pos_cm[0]:
        last_anc_pos_cm = current_anc_pos_cm
        continue
      else:
        #we've reached the end of a region. Need to print.
        if last_anc_pos_cm[0] == -9:
          hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] +
                    '\tUNK\t' + post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
        else:
          hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + last_anc_pos_cm[1] + '\t' +
                    pop_order[int(last_anc_pos_cm[0])-1] + '\t' +
                    post_anc_pos_cm[2] + '\t' + last_anc_pos_cm[2] + '\n')
        post_anc_pos_cm = current_anc_pos_cm
      
      last_anc_pos_cm = current_anc_pos_cm
    
    #last iteration, still need to print
    if last_anc_pos_cm[0] == -9:
      hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
                '\tUNK\t' + post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')
    else:
      hap.write(str(chr) + '\t' + post_anc_pos_cm[1] + '\t' + current_anc_pos_cm[1] +
                '\t' + pop_order[int(current_anc_pos_cm[0])-1] + '\t' +
                post_anc_pos_cm[2] + '\t' + current_anc_pos_cm[2] + '\n')

def main(current_ind, index, pop_order):
    
  #open bed files (2 haplotypes per individual)
  hap_a = open(args.out + '_A.bed', 'w')
  hap_b = open(args.out + '_B.bed', 'w')
  
  find_haplotype_bounds(index, 0, pop_order, hap_a,npop)
  hap_a.close()
  
  find_haplotype_bounds(index, 1, pop_order, hap_b,npop)
  hap_b.close()    

if __name__ == '__main__':
  #load parameters and files
  args = parse_args()
  pop_labels = args.pop_labels.split(',')
  npop = len(pop_labels)
  ind_info = open(args.ind_info)
  if args.ind is None:
    raise Exception('individual not set')
  current_ind = args.ind
  
  #reading order of individuals
  ind_list = []
  for line in ind_info:
      myLine = line.strip().split()
      ind_list.append(myLine[0])
  try:
    ind_index = ind_list.index(current_ind)
  except ValueError:
    raise Exception('Individual is not in the list provided')
  
  #set up chromosome variables
  chrs = range(1,23)
  if args.chrX:
    chrs.append('X')
  
  print 'Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
  main(current_ind, ind_index, pop_labels)
          
  print 'Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
