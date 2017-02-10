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

def true_prs(simulation, ncausal, h2, nhaps, out):
    """
    choose some number of causal alleles
    assign these alleles effect sizes
    from these effect sizes, compute polygenic risk scores for everyone
    """
    eprint('Reading all site info' + current_time())
    #my_sim2 = msprime.simulate(sample_size=60, Ne=1000, length=1e5, recombination_rate=2e-8, mutation_rate=2e-8)
    causal_mut_index = np.linspace(0, simulation.get_num_mutations()-1, ncausal, dtype=int)
    causal_mutations = set()
    
    # go through each population's trees
    out_sites = gzip.open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal) + '.sites.gz', 'w')
    out_sites.write('\t'.join(['Index', 'Pos', 'AFR_count', 'EUR_count', 'EAS_count', 'Total', 'beta']) + '\n')
    mut_info = {} # index -> position, afr count, eur count, eas count
    pop_count = 0
    for pop_leaves in [simulation.get_samples(population_id=0), simulation.get_samples(population_id=1), simulation.get_samples(population_id=2)]:
        for tree in simulation.trees(tracked_leaves=pop_leaves):
            for mutation in tree.mutations():
                if mutation.index in causal_mut_index:
                    causal_mutations.add(mutation)
                    if pop_count == 0:
                        mut_info[mutation.index] = [mutation.position, tree.get_num_tracked_leaves(mutation.node)]
                    else:
                        mut_info[mutation.index].append(tree.get_num_tracked_leaves(mutation.node))
        pop_count += 1
    
    causal_effects = {mutation.index:np.random.normal(loc=0,scale=h2/ncausal) for mutation in causal_mutations}
    for mutation in causal_mutations:
        causal_effects[mutation.index] = np.random.normal(loc=0,scale=h2/ncausal)
    
    eprint('Writing all site info' + current_time())
    for mutation in causal_mutations:
        out_sites.write(str(mutation.index) + '\t' + '\t'.join(map(str, mut_info[mutation.index])) + '\t' + str(simulation.get_sample_size()) + '\t')
        out_sites.write(str(causal_effects[mutation.index]) + '\n')
    out_sites.close()
    
    prs_haps = np.zeros(sum(nhaps)) #score for each haplotype
    eprint('Computing true PRS' + current_time())
    for variant in tqdm(simulation.variants(), total=simulation.get_num_mutations()):
        if variant.index in causal_mut_index:
            prs_haps += variant.genotypes * causal_effects[variant.index] # multiply vector of genotypes by beta for given variant
    prs_true = prs_haps[0::2] + prs_haps[1::2] #add to get individuals
    return(prs_true)
    
def case_control(prs_true, h2, nhaps, prevalence, ncontrols, out):
    """
    get cases assuming liability threshold model
    get controls from non-cases in same ancestry
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
    
    return(case_ids, control_ids, prs_norm, env_norm)

def run_gwas(simulation, diploid_cases, diploid_controls, p_threshold, cc_maf):
    """
    use cases and controls to compute OR, log(OR), and p-value for every variant
    """
    eprint('Running GWAS (' + str(len(diploid_cases)) + ' cases, ' + str(len(diploid_controls)) + ' controls)' + current_time())
    summary_stats = {} # position -> OR, p-value
    case_control = {} # position -> ncases w mut, ncontrols w mut
    
    cases = [2*x for x in diploid_cases] + [2*x+1 for x in diploid_cases]
    controls = [2*x for x in diploid_controls] + [2*x+1 for x in diploid_controls]
    
    eprint('Counting case mutations' + current_time())
    for tree in simulation.trees(tracked_leaves=cases):
        for mutation in tree.mutations():
            case_control[mutation.position] = [tree.get_num_tracked_leaves(mutation.node)]
            
    eprint('Counting control mutations' + current_time())
    for tree in simulation.trees(tracked_leaves=controls):
        for mutation in tree.mutations():
            case_control[mutation.position].append(tree.get_num_tracked_leaves(mutation.node))
            #case_control[mutation.position].append(mutation.position)
    
    # only keep sites with non-infinite or nan effect size with case and control maf > .01
    num_var = 0
    eprint('Computing fisher\'s exact test' + current_time())
    num_controls = float(len(controls))
    num_cases = float(len(cases))
    for position in tqdm(case_control):
        case_maf = min(case_control[position][0]/num_cases, (num_cases - case_control[position][0])/num_cases)
        control_maf = min(case_control[position][1]/num_controls, (num_controls - case_control[position][1])/num_controls)
        case_control_maf = min((case_control[position][0]+case_control[position][1])/(num_cases+num_controls), (num_cases + num_controls - case_control[position][0] - case_control[position][1])/(num_cases + num_controls))
        if case_control_maf > cc_maf:
            contingency = [[case_control[position][0], num_cases - case_control[position][0]],
                [case_control[position][1], num_controls - case_control[position][1]]]
            (OR, p) = stats.fisher_exact(contingency) #OR, p-value
            if not np.isnan(OR) and not np.isinf(OR) and OR != 0 and p <= p_threshold:
                summary_stats[position] = [OR, p]
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
    usable_positions = {} # position -> variant (simulation indices)
    
    sim_pos_index = {}
    for variant in tqdm(simulation.variants(), total=simulation.get_num_mutations()):
        if variant.position in summary_stats:
            usable_positions[variant.position] = variant
            sim_pos_index[variant.position] = variant.index
    
    # order all snps by p-value
    ordered_positions = sorted(summary_stats.keys(), key=lambda x: summary_stats[x][-1])
    #[(x, (x in usable_positions.keys())) for x in ordered_positions]
    
    eur_subset = simulation.subset(range(nhaps[0], (nhaps[0]+nhaps[1])))
    eur_index_pos = {}
    eur_pos_index = {}
    for mutation in tqdm(eur_subset.mutations(), total=eur_subset.get_num_mutations()):
        eur_index_pos[mutation.index] = mutation.position
        eur_pos_index[mutation.position] = mutation.index
    ordered_eur_index = sorted(eur_index_pos.keys())
    ld_calc = msprime.LdCalculator(eur_subset)
    #ld_calc = msprime.LdCalculator(simulation)
    
    # compute LD and prune in order of significance (popping index of SNPs)
    for position in ordered_positions:
        if position in usable_positions:
            r2_forward = ld_calc.get_r2_array(eur_pos_index[position], direction=msprime.FORWARD, max_distance=125e3)
            #print([position, np.where(r2_forward > r2_threshold)[0], np.where(r2_reverse > r2_threshold)[0]])
            for i in np.where(r2_forward > r2_threshold)[0]:
                usable_positions.pop(eur_index_pos[eur_pos_index[position]+i+1], None) #identify next position in eur space
            r2_reverse = ld_calc.get_r2_array(eur_pos_index[position], direction=msprime.REVERSE, max_distance=125e3)
            for i in np.where(r2_reverse > r2_threshold)[0]:
                usable_positions.pop(eur_index_pos[eur_pos_index[position]-i-1], None)
    
    clumped_snps = set(usable_positions.keys())
    
    eprint('Starting SNPs: ' + str(len(ordered_positions)) + '; SNPs after clumping: ' + str(len(clumped_snps)) + current_time())

    return(clumped_snps, usable_positions)
    

def infer_prs(simulation, nhaps, clumped_snps, summary_stats, usable_positions, h2, ncausal, out):
    """
    use clumped variants from biased gwas to compute inferred prs for everyone
    """
    eprint('Computing inferred PRS' + current_time())
    prs_haps = np.zeros(sum(nhaps))
    for variant in tqdm(simulation.variants(), total=simulation.get_num_mutations()):
        if variant.position in usable_positions:
            for ind in range(simulation.get_sample_size()):
                prs_haps[ind] += int(variant.genotypes[ind]) * math.log(summary_stats[variant.position][0])
    
    prs_infer = prs_haps[0::2] + prs_haps[1::2]
    
    # go through each population's trees
    out_sites = gzip.open(out + '_nhaps_' + '_'.join(map(str, nhaps)) + '_h2_' + str(round(h2, 2)) + '_m_' + str(ncausal) + '.infer_sites.gz', 'w')
    out_sites.write('\t'.join(['Index', 'Pos', 'AFR_count', 'EUR_count', 'EAS_count', 'Total', 'beta']) + '\n')
    mut_info = {}
    causal_mutations = set()
    pop_count = 0
    for pop_leaves in [simulation.get_samples(population_id=0), simulation.get_samples(population_id=1), simulation.get_samples(population_id=2)]:
        for tree in simulation.trees(tracked_leaves=pop_leaves):
            for mutation in tree.mutations():
                if mutation.position in usable_positions:
                    causal_mutations.add(mutation)
                    if pop_count == 0:
                        mut_info[mutation.index] = [mutation.position, tree.get_num_tracked_leaves(mutation.node)]
                    else:
                        mut_info[mutation.index].append(tree.get_num_tracked_leaves(mutation.node))
        pop_count += 1
    
    eprint('Writing all site info' + current_time())
    for mutation in causal_mutations:
        out_sites.write(str(mutation.index) + '\t' + '\t'.join(map(str, mut_info[mutation.index])) + '\t' + str(simulation.get_sample_size()) + '\t')
        out_sites.write(str(summary_stats[mutation.position][0]) + '\n')
        #
    out_sites.close()
    
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
        out_prs.write('\t'.join(map(str, [ind+1, pop, prs_true[ind], prs_infer[ind], pheno, scaled_env[ind]])) + '\n')
    out_prs.close()
        
def main(args):
    nhaps = map(int, args.nhaps.split(','))
    recomb = args.recomb_map
    ncausal = args.ncausal
    
    # generate/load coalescent simulations
    if args.tree is None:
        (pop_config, mig_mat, demog) = out_of_africa(nhaps)
        simulation = simulate_ooa(pop_config, mig_mat, demog, recomb)
        simulation.dump(args.out+ '_nhaps_' + '_'.join(map(str, nhaps)) + '.hdf5', True)
    else:
        simulation = msprime.load(args.tree)
    
    eprint(simulation)
    eprint('Number of haplotypes: ' + ','.join(map(str, nhaps)))
    eprint('Number of trees: ' + str(simulation.get_num_trees()))
    eprint('Number of mutations: ' + str(simulation.get_num_mutations()))
    eprint('Sequence length: ' + str(simulation.get_sequence_length()))

    
    prs_true = true_prs(simulation, args.ncausal, args.h2, nhaps, args.out)
    cases_diploid, controls_diploid, prs_norm, environment = case_control(prs_true, args.h2, nhaps, args.prevalence, args.ncontrols, args.out)
    summary_stats, cases_haploid, controls_haploid = run_gwas(simulation, cases_diploid, controls_diploid, args.p_threshold, args.cc_maf)
    clumped_snps, usable_positions = clump_variants(simulation, summary_stats, nhaps, args.r2, args.window_size)
    prs_infer = infer_prs(simulation, nhaps, clumped_snps, summary_stats, usable_positions, args.h2, args.ncausal, args.out)
    write_summaries(args.out, prs_true, prs_infer, nhaps, cases_diploid, controls_diploid, args.h2, args.ncausal, environment)
    
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
    parser.add_argument('--cc_maf', type=float, default=0.01)
    parser.add_argument('--r2', type=float, default=0.5)
    parser.add_argument('--window_size', type=int, default=250e3)
    parser.add_argument('--out', default='/Users/alicia/rare/chip_design/prs/simulations/sim0')
    
    args = parser.parse_args()
    main(args)
