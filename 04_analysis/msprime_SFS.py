#!/usr/bin/env python3
"""
Coalescent simulations for WBDC site frequency spectrum analysis.

This script runs msprime simulations to generate coalescent simulations
for comparison to the empirical SFS.

NOTE: We seem to have lost the original script but it looked similar to this.
"""

import msprime
import numpy as np
from collections import Counter


def run_msprime_simulation(num_replicates=1000, num_sites=100000, num_individuals=281, theta=0.008):
    """
    Run msprime coalescent simulations.
    
    Parameters
    ----------
    num_replicates : int
        Number of simulation replicates to run
    num_sites : int
        Number of sites to simulate (similar to cutdown VCF)
    num_individuals : int
        Number of haploid individuals
    theta : float
        Mutation parameter (4*N_e*mu). Default: 0.008
    
    Returns
    -------
    sfs_cumulative : dict
        Cumulative site frequency spectrum
    all_variants : list
        All variants across replicates
    """
    
    all_variants = []
    
    # theta = 0.008 per base pair (given parameter)
    # theta = 4*N_e*u, so N_e = theta/(4*u)
    # Mutation rate per site per generation (typical for plants: 5e-9)
    u = 5e-9
    Ne = theta / (4 * u)
    
    # Recombination rate: rho = 1.5 * theta
    rho = 1.5 * theta
    r = rho / (4 * Ne)
    
    print(f"Running {num_replicates} msprime replicates...")
    print(f"Simulating {num_sites} sites with {num_individuals} haploid individuals")
    print(f"Parameters: theta={theta}, rho={rho:.4f}")
    print(f"Mutation rate (u): {u:.2e}, Recombination rate (r): {r:.2e}")
    
    for rep in range(num_replicates):
        if (rep + 1) % 100 == 0:
            print(f"  Completed replicate {rep + 1}/{num_replicates}")
        
        # Run coalescent simulation
        ts = msprime.simulate(
            sample_size=num_individuals,  # haploid
            length=num_sites,
            recombination_rate=r,
            mutation_rate=u
        )
        
        # Extract variants
        for variant in ts.variants():
            # Count alternate alleles
            alt_count = np.sum(variant.genotypes)
            all_variants.append(alt_count)
    
    # Calculate cumulative SFS
    sfs_cumulative = calculate_cumulative_sfs(all_variants, num_individuals)
    
    return sfs_cumulative, all_variants


def calculate_cumulative_sfs(variants, num_chromosomes):
    """
    Calculate cumulative site frequency spectrum.
    
    Parameters
    ----------
    variants : list
        List of alternate allele counts
    num_chromosomes : int
        Total number of haploid individuals
    
    Returns
    -------
    sfs_cumulative : dict
        Dictionary mapping frequency class to cumulative count
    """
    
    sfs = Counter(variants)
    
    # Calculate frequency classes (as a fraction of n)
    sfs_cumulative = {}
    
    for count in sorted(sfs.keys()):
        if 0 < count < num_chromosomes:  # Skip monomorphic sites
            frequency = count / num_chromosomes
            cumulative_count = sum(sfs[c] for c in sfs if c <= count)
            sfs_cumulative[count] = cumulative_count
    
    return sfs_cumulative


def print_sfs_table(sfs_cumulative, num_individuals):
    """
    Print the SFS in a readable table format.
    
    Parameters
    ----------
    sfs_cumulative : dict
        Cumulative SFS from simulations
    num_individuals : int
        Total number of haploid individuals
    """
    
    print("\n" + "="*60)
    print("CUMULATIVE SITE FREQUENCY SPECTRUM")
    print("="*60)
    print(f"{'Allele Count':<15} {'Frequency':<15} {'Cumulative SNPs':<20}")
    print("-"*60)
    
    for count in sorted(sfs_cumulative.keys()):
        frequency = count / num_individuals
        cum_count = sfs_cumulative[count]
        print(f"{count:<15} {frequency:<15.6f} {cum_count:<20}")
    
    print("="*60)


if __name__ == "__main__":
    # Run simulations
    num_replicates = 1000
    num_sites = 100000
    num_individuals = 281  # haploid individuals
    theta = 0.008
    
    sfs_cumulative, all_variants = run_msprime_simulation(
        num_replicates=num_replicates,
        num_sites=num_sites,
        num_individuals=num_individuals,
        theta=theta
    )
    
    # Print results
    print_sfs_table(sfs_cumulative, num_individuals)
    
    # Additional summary statistics
    print("\nSUMMARY STATISTICS")
    print("-"*60)
    print(f"Total variants simulated: {len(all_variants)}")
    print(f"Segregating sites: {len([v for v in all_variants if 0 < v < num_individuals])}")
    print(f"Mean allele frequency: {np.mean(all_variants):.4f}")
    print(f"Median allele frequency: {np.median(all_variants):.4f}")
