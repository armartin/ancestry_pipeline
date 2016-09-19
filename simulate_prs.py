from __future__ import print_function
import sys
import msprime
import numpy as np
import math
import argparse
from datetime import datetime
import random
import gzip
import os
from scipy import stats
from collections import defaultdict
from collections import Counter
from itertools import izip
from tqdm import tqdm

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def current_time():
    return(' [' + datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + ']')

def out_of_africa(nhaps):
    """
    Specify the demographic model used in these simulations (Gravel et al, 2011 PNAS)
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Gravel et al, 2011 Table 2.
    N_A = 7300
    N_B = 1861
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time
    T_B = 51e3 / generation_time
    T_EU_AS = 23e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.0038
    r_AS = 0.0048
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=nhaps[0], initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=nhaps[1], initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=nhaps[2], initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dp = msprime.DemographyDebugger(
        Ne=N_A,
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    dp.print_history()
    
    return(population_configurations, migration_matrix, demographic_events)

def simulate_ooa(population_configurations, migration_matrix, demographic_events, recomb):
    """
    simulate according to the specified demographic model with recombination
    """
    eprint('Starting simulations' + current_time())
    simulation = msprime.simulate(
        population_configurations = population_configurations,
        migration_matrix=migration_matrix,
        mutation_rate=2e-8,
        recombination_map = msprime.RecombinationMap.read_hapmap(recomb)
    )
    eprint('Ending simulations' + current_time())
    return(simulation)

def compute_prs(simulation, ncausal, h2, nhaps, out):
    """
    choose some number of causal alleles
    assign these alleles effect sizes
    from these effect sizes, compute polygenic risk scores for everyone
    """
    eprint('Reading all site info' + current_time())
    #my_sim2 = msprime.simulate(sample_size=60, Ne=1000, length=1e5, recombination_rate=2e-8, mutation_rate=2e-8)
    causal_mut_index = np.linspace(0, simulation.get_num_mutations()-1, ncausal, dtype=int)
    causal_mut_pos = set()
    causal_vars = set(causal_mut_index)
    
    # compute means and standard deviations for every variant
    means = np.zeros(simulation.get_num_mutations()) # allele frequency
    sds = np.zeros(simulation.get_num_mutations()) # pq/2N - from the binomial, assumes diploid, just N because haploid seqs
    
    out_sites = gzip.open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal) + '.sites.gz', 'w')
    out_sites.write('\t'.join(['Index', 'Pos', 'AFR_count', 'EUR_count', 'EAS_count', 'Total', 'beta']) + '\n')
    mut_info = {}
    pop_count = 0
    for pop_leaves in [simulation.get_samples(population_id=0), simulation.get_samples(population_id=1), simulation.get_samples(population_id=2)]:
        mut_num = 0
        for tree in simulation.trees(tracked_leaves=pop_leaves):
            #if mut_num == 0:
                #eprint(pop_leaves[0:10])
            for pos, node in tree.mutations():
                
                if mut_num in causal_vars:
                    if pop_count == 0:
                        mut_info[mut_num] = [pos, tree.get_num_tracked_leaves(node)]
                    else:
                        mut_info[mut_num].append(tree.get_num_tracked_leaves(node))
                mut_num += 1
        pop_count += 1
    
    mut_num = 0               
    for tree in simulation.trees():
        #print(tree.get_interval(), str(tree))
        for pos, node in tree.mutations():
            
            if mut_num in causal_vars:
                mut_info[mut_num].append(tree.get_sample_size())
                causal_mut_pos.add(pos)
            means[mut_num] = float(tree.get_num_leaves(node)) / tree.get_sample_size()
            sds[mut_num] = np.sqrt(means[mut_num] * (1 - means[mut_num]) / (tree.get_sample_size()))
            mut_num += 1
            if not mut_num % 1000000:
                print("Mutation {} @ position {} has frequency {} / {}".format(
                mut_num, pos, tree.get_num_leaves(node), tree.get_sample_size()))
    
    causal_means = [means[x] for x in causal_mut_index]
    causal_sds = [sds[x] for x in causal_mut_index]
    causal_effects = np.random.normal(loc=0,scale=h2/ncausal, size=ncausal)# / causal_sds
    
    eprint('Writing all site info' + current_time())
    for var in range(len(causal_mut_index)):
        out_sites.write(str(causal_mut_index[var]) + '\t' + '\t'.join(map(str, mut_info[causal_mut_index[var]])) + '\t')
        out_sites.write(str(causal_effects[var]) + '\n')
        #
    out_sites.close()
    prs_haps = np.zeros(sum(nhaps)) #score for each individual
    
    eprint('Computing true PRS' + current_time())
    mut_num = 0
    for pos, variant in tqdm(simulation.variants(), total=simulation.get_num_mutations()):
        if pos in causal_mut_pos:
            for ind in range(len(variant)):
                prs_haps[ind] += int(variant[ind]) * causal_effects[mut_num] # mean center variants w global frequencies?
            mut_num += 1
        
    prs_true = prs_haps[0::2] + prs_haps[1::2]
    return(prs_true)
    
def case_control(prs_true, h2, nhaps, prevalence, ncontrols, out):
    """
    get cases assuming liability threshold model
    get controls from non-cases in same ancestry
    write file containing cases, controls, polygenic scores for every individual
    """
    eprint('Defining cases/controls' + current_time())
    env_effect = np.random.normal(loc=0,scale=1-h2, size=sum(nhaps)/2)
    prs_norm = (prs_true - np.mean(prs_true)) / np.std(prs_true)
    env_norm = (env_effect - np.mean(env_effect)) / np.std(env_effect)
    total_liability = math.sqrt(h2) * prs_norm + math.sqrt(1 - h2) * env_norm
    eur_liability = total_liability[nhaps[0]/2:(nhaps[0]+nhaps[1])/2]
    sorted_liability = sorted(eur_liability)
    cases = [i for (i, x) in enumerate(eur_liability) if x >= sorted_liability[int((1-prevalence)*len(sorted_liability))]]
    controls = set(range(nhaps[1]/2))
    for case in cases:
        controls.remove(case)
    controls = random.sample(controls, ncontrols)
    
    case_ids = map(lambda(x): x+nhaps[0]/2, cases)
    control_ids = sorted(map(lambda(x): x+nhaps[0]/2, controls))
    
    #print 'Writing output' + current_time()
    #out_pheno = open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal) + '.pheno', 'w')
    #out_prs = open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal)  + '.prs', 'w')
    #
    #for case in case_ids:
    #    current_id = 'msp_' + str(case)
    #    out_pheno.write('\t'.join([current_id, current_id, '2']) + '\n')
    #
    #for control in control_ids:
    #    current_id = 'msp_' + str(control)
    #    out_pheno.write('\t'.join([current_id, current_id, '1']) + '\n')
    #
    #for score in range(len(prs_true)):
    #    out_prs.write('msp_' + str(score) + '\t' + str(prs_true[score]) + '\t' + str(total_liability[score]) + '\n')
    #
    #out_pheno.close()
    #out_prs.close()
    return(case_ids, control_ids, prs_norm, env_norm)

def run_gwas(simulation, diploid_cases, diploid_controls, p_threshold):
    """
    use cases and controls to compute OR, log(OR), and p-value for every variant
    """
    eprint('Running GWAS (' + str(len(diploid_cases)) + ' cases, ' + str(len(diploid_controls)) + ' controls)' + current_time())
    summary_stats = {} # pos -> OR, p-value
    pos_case_control = {} # pos -> ncases w mut, ncontrols w mut
    
    cases = [2*x for x in diploid_cases] + [2*x+1 for x in diploid_cases]
    controls = [2*x for x in diploid_controls] + [2*x+1 for x in diploid_controls]
    
    eprint('Counting case mutations' + current_time())
    
    for tree in simulation.trees(tracked_leaves=cases):
        for pos, node in tree.mutations():
            pos_case_control[pos] = [tree.get_num_tracked_leaves(node)]
            
    eprint('Counting control mutations' + current_time())
    for tree in simulation.trees(tracked_leaves=controls):
        for pos, node in tree.mutations():
            pos_case_control[pos].append(tree.get_num_tracked_leaves(node))
    
    # only keep sites with non-infinite or nan effect size
    num_var = 0
    eprint('Computing fisher\'s exact test' + current_time())
    num_controls = float(len(controls))
    num_cases = float(len(cases))
    for pos in tqdm(pos_case_control):
        case_maf = min(pos_case_control[pos][0]/num_cases, (num_cases - pos_case_control[pos][0])/num_cases)
        control_maf = min(pos_case_control[pos][1]/num_controls, (num_controls - pos_case_control[pos][1])/num_controls)
        if case_maf > 0.01 and control_maf > 0.01:
            contingency = [[pos_case_control[pos][0], num_controls - pos_case_control[pos][0]],
                [pos_case_control[pos][1], num_controls - pos_case_control[pos][1]]]
            fisher = stats.fisher_exact(contingency)
            if not np.isnan(fisher[0]) and not np.isinf(fisher[0]) and fisher[1] <= p_threshold:
                summary_stats[pos] = fisher
                num_var += 1
                #if not num_var % 100: #remove this eventually
                #    break
            
    eprint('Done with GWAS! (' + str(len(summary_stats)) + ' amenable sites)' + current_time())
        
    return(summary_stats, cases, controls)

def clump_variants(simulation, summary_stats, nhaps, r2_threshold, window_size):
    """
    perform variant clumping in a greedy fasion with p-value and r2 threshold in windows
    return only those variants meeting some nominal threshold
    
    1: make a dict of pos -> variant for subset of sites meeting criteria
    2: make an r2 dict of all pairs of snps meeting p-value threshold and in same window
    """
    # make a list of SNPs ordered by p-value
    eprint('Subsetting variants to usable list' + current_time())
    usable_subset = {}
    usable_variants = {}
    max_per_pop = min(min(nhaps), 1000) #only compute r2 on a max of 1000 haps per population
    for pos, variant in tqdm(simulation.variants(), total=simulation.get_num_mutations()):
        if pos in summary_stats:
            usable_variants[pos] = variant
            usable_subset[pos] = variant[nhaps[0]:(nhaps[0]+max_per_pop)]
    
    r2 = defaultdict(dict)
    eprint('Computing r2 between all variants within window meeting significance threshold (' + str(max_per_pop) + ' European haplotypes)' + current_time())
    #num_vars = some constant so you don't 
    for posA in tqdm(usable_subset):
        variantA = usable_subset[posA]
        variantA_len = float(len(usable_subset[posA]))
        for posB in usable_subset:
            if posA > posB and posA - posB <= window_size:
                variantB = usable_subset[posB]
                variantB_len = float(len(usable_subset[posB]))
                hap = Counter(izip(variantA, variantB)) # get the counts of covariats
                D = hap[('0', '0')]/variantA_len * hap[('1', '1')]/variantA_len - hap[('0', '1')]/variantA_len * hap[('1', '0')]/variantA_len
                countA = Counter(variantA)
                countB = Counter(variantB)
                r2_denom = countA['0']/variantA_len*countA['1']/variantA_len*countB['0']/variantB_len*countB['1']/variantA_len
                if r2_denom > 0:
                    r2[posA][posB] = float(D**2)/r2_denom
                else:
                    r2[posA][posB] = 1 #make sure there are no divide by zero errors                    
                #if float(countA['1'])/(countA['0'] + countA['1']) > 0.05:
                #    print [posA, posB, countA, countB, D, my_r2]
                #print [posA, posB, countA, countB, hap, summary_stats[posA][-1], summary_stats[posB][-1], D, r2[posA][posB]]
    
    eprint('Greedily go through variants in p-value order and remove those in LD' + current_time())
    # get a list of positions, sorted by p-value
    ordered_pos = sorted(summary_stats.keys(), key=lambda x: summary_stats[x][-1])
    clumped_snps = set(ordered_pos)
    for posA in ordered_pos:
        for posB in ordered_pos:
            if abs(posA - posB) < window_size:
                p_A = summary_stats[posA][1]
                p_B = summary_stats[posB][1]
                if posA > posB:
                    current_r2 = r2[posA][posB]
                elif posB > posA:
                    current_r2 = r2[posB][posA]
                else:
                    continue
                if current_r2 > r2_threshold:
                    if p_A < p_B:
                        try:
                            clumped_snps.remove(posB)
                        except KeyError: #already removed by a previous variant
                            pass
                    elif p_B < p_A:
                        try:
                            clumped_snps.remove(posA)
                        except KeyError: #already removed by a previous variant
                            pass
    eprint('Starting SNPs: ' + str(len(ordered_pos)) + '; SNPs after clumping: ' + str(len(clumped_snps)) + current_time())
    
    return(clumped_snps, usable_variants)
    

def infer_prs(simulation, nhaps, clumped_snps, summary_stats, usable_variants):
    """
    use clumped variants from biased gwas to compute inferred prs for everyone
    """
    eprint('Computing inferred PRS' + current_time())
    prs_haps = np.zeros(sum(nhaps))
    for pos in usable_variants:
        variant = usable_variants[pos]
        if pos in clumped_snps:
            if summary_stats[pos][0] > 0:
                for ind in range(len(variant)):
                    prs_haps[ind] += int(variant[ind]) * math.log(summary_stats[pos][0])
    
    prs_infer = prs_haps[0::2] + prs_haps[1::2]
    return(prs_infer)
    
def write_summaries(out, prs_true, prs_infer, nhaps, cases, controls, h2, ncausal, environment):
    eprint('Writing output!' + current_time())
    scaled_prs = math.sqrt(h2) * prs_true
    scaled_env = math.sqrt(1 - h2) * environment
    out_prs = gzip.open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal)  + '.prs.gz', 'w')
    out_prs.write('\t'.join(['Ind', 'Pop', 'PRS_true', 'PRS_infer', 'Pheno', 'Environment']) + '\n')
    for ind in range(len(prs_true)):
        if ind in cases:
            pheno = 1
        elif ind in controls:
            pheno = 0
        else:
            pheno = 'NA'
        if ind in range(nhaps[0]/2):
            pop = 'AFR'
        elif ind in range(nhaps[0]/2, nhaps[0]/2+nhaps[1]/2):
            pop = 'EUR'
        else:
            pop = 'EAS'
        out_prs.write('\t'.join(map(str, [ind+1, pop, scaled_prs[ind], prs_infer[ind], pheno, scaled_env[ind]])) + '\n')
    out_prs.close()
        
def main(args):
    nhaps = map(int, args.nhaps.split(','))
    recomb = args.recomb_map
    ncausal = args.ncausal
    
    # generate/load coalescent simulations
    if args.tree is None:
        (pop_config, mig_mat, demog) = out_of_africa(nhaps)
        simulation = simulate_ooa(pop_config, mig_mat, demog, recomb)
        simulation.dump(args.out+ '_nhaps_' + '_'.join(map(str, nhaps)) + '.hdf5', True) # UPDATE FILENAME FROM FULL RUN (not out + '.dump')
    else:
        simulation = msprime.load(args.tree)
    
    eprint(simulation)
    eprint('Number of haplotypes: ' + ','.join(map(str, nhaps)))
    eprint('Number of trees: ' + str(simulation.get_num_trees()))
    eprint('Number of mutations: ' + str(simulation.get_num_mutations()))
    eprint('Sequence length: ' + str(simulation.get_sequence_length()))
    
    
    prs_true = compute_prs(simulation, args.ncausal, args.h2, nhaps, args.out)
    cases_diploid, controls_diploid, prs_norm, environment = case_control(prs_true, args.h2, nhaps, args.prevalence, args.ncontrols, args.out)
    summary_stats, cases_haploid, controls_haploid = run_gwas(simulation, cases_diploid, controls_diploid, args.p_threshold)
    clumped_snps, usable_variants = clump_variants(simulation, summary_stats, nhaps, args.r2, args.window_size)
    prs_infer = infer_prs(simulation, nhaps, clumped_snps, summary_stats, usable_variants)
    write_summaries(args.out, prs_norm, prs_infer, nhaps, cases_diploid, controls_diploid, args.h2, args.ncausal, environment)
    
    #print summary_stats


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree')
    parser.add_argument('--nhaps', help='AFR,EUR,EAS', default='400000,400000,400000')
    parser.add_argument('--recomb_map', default='/Users/alicia/Documents/Grad_School/Rotations/Bustamante_Rotation/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr20.txt')
    parser.add_argument('--ncausal', type=int, default=200)
    parser.add_argument('--ncontrols', type=int, default=10000)
    parser.add_argument('--h2', type=float, default=float(2)/3)
    parser.add_argument('--prevalence', type=float, default=0.05)
    parser.add_argument('--p_threshold', type=float, default=0.01)
    parser.add_argument('--r2', type=float, default=0.5)
    parser.add_argument('--window_size', type=int, default=250e3)
    parser.add_argument('--out', default='/Users/alicia/rare/chip_design/prs/simulations/sim0')
    
    args = parser.parse_args()
    main(args)
