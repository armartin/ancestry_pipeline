__author__ = 'armartin'
#takes in RFMix output and creates collapsed 2 haploid bed files per individual

from datetime import datetime
import time
from optparse import  OptionParser
from itertools import izip_longest

def parse_options():
  USAGE = """
  collapse_ancestry.py    --rfmix
                          --snp_locations
                          --fbk
                          --fbk_threshold
                          --ind
                          --ind_info
                          --pop_labels
                          --chrX
                          --out
  """
  
  parser = OptionParser(USAGE)
  
  parser.add_option('--rfmix', help='path to RFMix Viterbi output, chr expected in filename')
  parser.add_option('--snp_locations', help='path to snp_locations file required by RFMix, chr expected in filename')
  parser.add_option('--fbk', default=None)
  parser.add_option('--fbk_threshold', type='float', default = 0.9)
  parser.add_option('--ind', help='Individual ID, must match a line in --ind_info option')
  parser.add_option('--ind_info', help='Individual IDs listed in the order they appear in RFMix Viterbi output')
  parser.add_option('--pop_labels', type='string', action='callback', callback=splitstr, default=['AFR','EUR','NAT'],
                    help='comma-separated list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels')
  parser.add_option('--chrX', help='include chrX?', default=False, action="store_true")
  parser.add_option('--out', help='prefix to bed file, _A.bed and _B.bed will be appended')
  
  (options,args)=parser.parse_args()
  return(options)
  
def splitstr(option, opt, value, parser):
  return(setattr(parser.values, option.dest, value.split(',')))

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

#def track_lai_proportions(last_bound, current_bound, current_info, lai_proportions):
#    current_bound = float(current_bound)
#    last_bound = float(last_bound)
#    for hap in range(1,3):
#        if current_info[hap] == '1':
#            lai_proportions[0] = lai_proportions[0] + current_bound - last_bound
#        elif current_info[hap] == '2':
#            lai_proportions[1] = lai_proportions[1] + current_bound - last_bound
#        elif current_info[hap] == '3':
#            lai_proportions[2] = lai_proportions[2] + current_bound - last_bound
#        elif current_info[hap] == '4':
#            lai_proportions[3] = lai_proportions[3] + current_bound - last_bound
#    return lai_proportions
    

#def plot_chromosomes(last_anc, chr, last_plot_bound, current_plot_bound, current_info, ax):
#    if last_anc is None:
#        pass
#    else:
#        #first ancestry goes on the bottom
#        verts1 = [
#            (float(last_plot_bound), chr - 0.4), #left, bottom
#            (float(last_plot_bound), chr), #left, top
#            (float(current_plot_bound), chr), #right, top
#            (float(current_plot_bound), chr - 0.4), #right, bottom
#            (0, 0), #ignored
#        ]
#        #second ancestry goes on the top
#        verts2 = [
#            (float(last_plot_bound), chr), #left, bottom
#            (float(last_plot_bound), chr + 0.4), #left, top
#            (float(current_plot_bound), chr + 0.4), #right, top
#            (float(current_plot_bound), chr), #right, bottom
#            (0, 0), #ignored
#        ]
#        codes = [
#            Path.MOVETO,
#            Path.LINETO,
#            Path.LINETO,
#            Path.LINETO,
#            Path.CLOSEPOLY,
#        ]
#        
#        bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
#        colors=bmap.mpl_colors
#        colors.append((0,0,0))
#        
#        path = Path(verts1, codes)
#        path2 = Path(verts2, codes)
#        if last_anc[0] != -9:
#            patch = patches.PathPatch(path, color=colors[int(last_anc[0])-1], lw=0)
#        else:
#            patch = patches.PathPatch(path, color=colors[-1], lw=0)
#        ax.add_patch(patch)
#        if last_anc[1] != -9:
#            patch = patches.PathPatch(path2, color=colors[int(last_anc[1])-1], lw=0)
#        else:
#            patch = patches.PathPatch(path2, color=colors[-1], lw=0)
#        ax.add_patch(patch)

def main(current_ind, ind):
    
    ##generate figure settings
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.set_xlim(-5,300)
    #ax.set_ylim(23,0)
    #plt.xlabel('Genetic position (cM)')
    #plt.ylabel('Chromosome')
    #plt.title(current_ind)
    #plt.yticks(range(1,23))
    
    #open bed files (2 haplotypes per individual)
    hap_a = open(options.out + current_ind + '_A.bed', 'w')
    hap_b = open(options.out + current_ind + '_B.bed', 'w')
    
    ##admixed_pop
    #lai_proportions = [0,0,0] #changed this to only have 3-way ancestry (admixed).
    
    counter = 0
    for chr in chrs:
        print chr
        rfmix = open(options.rfmix) #1.rfmix.0.Viterbi.txt
        if options.fbk is not None: #will need to change global proportion calculator, make new anc possibility (-9?), and color black
            fbk = open(options.fbk)
        snp_locations = open(options.snp_locations + str(chr) + '.snp_locations') ###fix this filename
        snp_map = open(options.snp_locations + str(chr) + '.map') #map of physical position -> genetic position
        ind = ind_list.index(current_ind)

        last_anc = None
        last_plot_bound = None
        last_hapa_anc = 0
        last_hapa_pos = None
        last_hapb_anc = 0
        last_hapb_pos = None
        
        for line in rfmix:
            myLine = line.strip().split()
            my_pos = snp_locations.readline().strip()
            my_map = snp_map.readline().strip().split()
            #finds max of 3 ancestry posterior probabilities
            if options.fbk is not None:
                fbk_line = fbk.readline().strip().split()
                fbk_line = map(float, fbk_line)
                fbk_max = []
                for i in grouper(3, fbk_line):
                    fbk_max.append(max(i))
            #print fbk_max
            current_anc = [myLine[ind*2], myLine[ind*2+1]]
            current_info = [my_pos, myLine[ind*2], myLine[ind*2+1]]
            if fbk_max[ind*2] >= fbk_threshold:
                current_hapa_anc = myLine[ind*2]
            else:
                current_hapa_anc = -9
                current_anc[0] = -9
            if fbk_max[ind*2+1] >= fbk_threshold:
                current_hapb_anc = myLine[ind*2+1]
            else:
                current_hapb_anc = -9
                current_anc[1] = -9
            #current_hapb_anc = myLine[ind*2+1]
            current_hapa_pos = my_map[0]
            current_hapb_pos = my_map[0]
            current_hapa_cm = my_map[1]
            current_hapb_cm = my_map[1]
            current_plot_bound = my_pos
            if last_hapa_anc == 0:
                last_hapa_anc = current_hapa_anc
                last_hapa_pos = current_hapa_pos
                last_hapa_cm = current_hapa_cm
            if last_hapb_anc == 0:
                last_hapb_anc = current_hapb_anc
                last_hapb_pos = current_hapb_pos
                last_hapb_cm = current_hapb_cm
            if current_anc == last_anc:
                pass
            else:
                print [last_anc, chr, last_plot_bound, current_plot_bound, current_info, [current_hapa_anc, current_hapb_anc]]
                #plot_chromosomes(last_anc, chr, last_plot_bound, current_plot_bound, current_info, ax)
                #print chr
                #print last_anc
                if last_plot_bound is not None:
                    #lai_proportions = track_lai_proportions(last_plot_bound, current_plot_bound, [current_plot_bound, last_anc[0], last_anc[1]], lai_proportions)
                    counter = counter + 2*(float(current_plot_bound) - float(last_plot_bound))
                    if last_hapa_anc is not None and last_hapa_anc != current_hapa_anc:
                        print [chr, last_hapa_pos, current_hapa_pos, last_hapa_anc]
                        if last_hapa_anc != -9:
                            hap_a.write(str(chr) + '\t' + last_hapa_pos + '\t' + current_hapa_pos + '\t' + pop_labels[int(last_hapa_anc)-1] + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
                        else:
                            hap_a.write(str(chr) + '\t' + last_hapa_pos + '\t' + current_hapa_pos + '\tUNK\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
                        last_hapa_anc = current_hapa_anc
                        last_hapa_pos = current_hapa_pos
                        last_hapa_cm = current_hapa_cm
                    if last_hapb_anc is not None and last_hapb_anc != current_hapb_anc:
                        print [chr, last_hapb_pos, current_hapb_pos, last_hapb_anc]
                        if last_hapb_anc != -9:
                            hap_b.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\t' + pop_labels[int(last_hapb_anc)-1] + '\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
                        else:
                            hap_b.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\tUNK\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
                        last_hapb_anc = current_hapb_anc
                        last_hapb_pos = current_hapb_pos
                        last_hapb_cm = current_hapb_cm
                last_plot_bound = current_plot_bound
                last_anc = current_anc #might need to change this for plotting purposes
        if current_hapa_anc == -9:
            hap_a.write(str(chr) + '\t' + last_hapa_pos + '\t' + current_hapa_pos + '\tUNK\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
        else:
            hap_a.write(str(chr) + '\t' + last_hapa_pos + '\t' + current_hapa_pos + '\t' + pop_labels[int(current_hapa_anc)-1] + '\t' + last_hapa_cm + '\t' + current_hapa_cm + '\n')
        if current_hapb_anc == -9:
            hap_b.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\tUNK\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
        else:
            hap_b.write(str(chr) + '\t' + last_hapb_pos + '\t' + current_hapb_pos + '\t' + pop_labels[int(current_hapb_anc)-1] + '\t' + last_hapb_cm + '\t' + current_hapb_cm + '\n')
        #lai_proportions = track_lai_proportions(last_plot_bound, current_plot_bound, [current_plot_bound, last_anc[0], last_anc[1]], lai_proportions)
        print [last_anc, chr, last_plot_bound, current_plot_bound, current_info]
        #plot_chromosomes(last_anc, chr, last_plot_bound, current_plot_bound, current_info, ax)
    
    hap_a.close()
    hap_b.close()
    
    #bmap = brewer2mpl.get_map('Set1', 'qualitative', 4)
    #colors=bmap.mpl_colors
    #
    #p1 = plt.Rectangle((0, 0), 1, 1, color=colors[0])
    #p2 = plt.Rectangle((0, 0), 1, 1, color=colors[1])
    #p3 = plt.Rectangle((0, 0), 1, 1, color='k')
    #if len(options.pop_labels) > 2:
    #    p3 = plt.Rectangle((0, 0), 1, 1, color=colors[2])
    #    p4 = plt.Rectangle((0, 0), 1, 1, color='k')
    #if len(options.pop_labels) > 3:
    #    p4 = plt.Rectangle((0, 0), 1, 1, color=colors[3])
    #    p5 = plt.Rectangle((0, 0), 1, 1, color='k')
    #labs = list(options.pop_labels)
    #labs.append('UNK (Prob < 0.9)')
    #if len(options.pop_labels) ==2:
    #    ax.legend([p1, p2, p3], labs, loc=4, fancybox=True)
    #elif len(options.pop_labels) ==3:
    #    ax.legend([p1, p2, p3, p4], labs, loc=4, fancybox=True)
    #elif len(options.pop_labels) ==4:
    #    ax.legend([p1, p2, p3, p4, p5], labs, loc=4, fancybox=True)
    #
    ## Remove top and right axes lines ("spines")
    #spines_to_remove = ['top', 'right']
    #for spine in spines_to_remove:
    #    ax.spines[spine].set_visible(False)
    #
    ## Get rid of ticks. The position of the numbers is informative enough of
    ## the position of the value.
    #ax.xaxis.set_ticks_position('none')
    #ax.yaxis.set_ticks_position('none')
    #
    #fig.savefig(options.out + current_ind + '.png', transparent=True)
    #return(lai_proportions)

if __name__ == '__main__':
  #load parameters and files
  options = parse_options()
  pop_labels = options.pop_labels
  ind_info = open(options.ind_info)
  if options.ind is None:
    raise Exception('individual not set')
  current = options.ind
  
  #reading order of individuals
  ind_list = []
  for line in ind_info:
      myLine = line.strip().split()
      ind_list.append(myLine[0])
  print ind_list
  try:
    ind_index = ind_list.index(ind)
  except ValueError:
    raise Exception('Individual is not in the list provided')
  
  #set up chromosome variables
  chrs = range(1,23)
  chrs.append('X') #remove this if only autosomes were run
  
  print 'Starting [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
  out = open(options.out, 'w')
  for pop in options.pop_labels:
      out.write('\tRF_' + pop)
  out.write('\n')
  
  for ind in current_ind:
      if ind == '0':
          pass
      else:
          print ind
          lai_proportions = main(ind, ind_list.index(ind))
          out.write(ind + '\t')
          out.write(str(lai_proportions[0]/sum(lai_proportions)) + '\t' + str(lai_proportions[1]/sum(lai_proportions)) + '\t')
          if len(options.pop_labels) > 2:
              out.write(str(lai_proportions[2]/sum(lai_proportions)) + '\t')
          if len(options.pop_labels) > 3:
              out.write(str(lai_proportions[3]/sum(lai_proportions)))
          out.write('\n')
          
  print 'Finished [' + datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S') + ']'
    
import matplotlib.pyplot as plt
import pylab
from matplotlib.path import Path
import matplotlib.patches as patches
import brewer2mpl
