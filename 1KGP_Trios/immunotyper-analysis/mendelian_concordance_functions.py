import pandas as pd
import os
from collections import defaultdict
import itertools
from glob import glob
import seaborn as sns
import matplotlib.pyplot as plt
from statistics import mean
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path
from Bio import SeqIO
from typing import Tuple, List

def load_genotype_data(genotype_dir, gene_type):
    """
    Load genotype data for each individual, structure as Gene -> Alleles list.
    
    Args:
        genotype_dir (str): Directory containing genotype files.
        gene_type (str): Type of gene to load data for (e.g., 'TRBV').
        
    Returns:
        dict: Dictionary where keys are sample IDs and values are dictionaries 
              where 'Gene' is the key and 'Alleles' is a list of alleles.
              
    Raises:
        FileNotFoundError: If no genotype files are found matching the pattern.
    """
    glob_str = f"{genotype_dir}*{gene_type.upper()}_functional_allele_calls.txt"
    genotype_files = glob(glob_str)
    if not genotype_files:
        raise FileNotFoundError(
            f"No genotype files found matching pattern '{glob_str}' "
        )
    
    genotype_dfs = {}
    for file_path in genotype_files:
        sample_id = os.path.split(file_path)[-1].split('.')[0]  # Extract sample ID from the filename

        # Create a dictionary to store alleles for each gene
        gene_alleles = defaultdict(list)
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()  # Remove newline characters
                if '*' in line:
                    gene, allele = line.split('*', 1)  # Split on the first '*' character
                    gene_alleles[gene].append(allele)
        
        # Convert dictionary into a DataFrame with Gene as the index and Alleles as the column
        gene_alleles_df = pd.DataFrame(list(gene_alleles.items()), columns=['Gene', 'Allele'])
        gene_alleles_df.set_index('Gene', inplace=True)
        genotype_dfs[sample_id] = gene_alleles_df
            
    return genotype_dfs

def calculate_mendelian_concordance(child_id, father_id, mother_id, genotype_dfs):
    """
    Calculate Mendelian concordance for a given trio, handling cases where child or parent alleles may be missing.
    
    Args:
        child_id (str): The child's sample ID.
        father_id (str): The father's sample ID.
        mother_id (str): The mother's sample ID.
        genotype_dfs (dict): Dictionary with genotype dataframes for each individual.

    Returns:
        concordance_dict (dict): Dictionary with gene as key and concordance (0/1) as value.
    """
    concordance_dict = {}
    
    if child_id not in genotype_dfs or father_id not in genotype_dfs or mother_id not in genotype_dfs:
        return concordance_dict  # Return empty if any of the trio is missing
    
    child_genotype = genotype_dfs[child_id]
    father_genotype = genotype_dfs[father_id]
    mother_genotype = genotype_dfs[mother_id]
    
    # Iterate through each gene in the child's genotype
    for gene, child_alleles in child_genotype['Allele'].items():
        if gene in father_genotype.index and gene in mother_genotype.index:
            father_alleles = father_genotype.loc[gene, 'Allele']
            mother_alleles = mother_genotype.loc[gene, 'Allele']
            
            # Handle cases where the child has two or more alleles
            if len(child_alleles) >= 2:
                # Count how many alleles match parents
                matching_alleles = 0
                for pair in itertools.combinations(child_alleles, 2):
                    if (pair[0] in father_alleles and pair[1] in mother_alleles) or \
                       (pair[1] in father_alleles and pair[0] in mother_alleles):
                        matching_alleles = 2
                        break
                    elif (pair[0] in father_alleles or pair[0] in mother_alleles) or \
                         (pair[1] in father_alleles or pair[1] in mother_alleles):
                        matching_alleles = 1
                
                concordance_dict[gene] = matching_alleles / 2  # 0, 0.5, or 1
            elif len(child_alleles) == 1:
                # Check if the single allele is from either parent
                if child_alleles[0] in father_alleles or child_alleles[0] in mother_alleles:
                    concordance_dict[gene] = 1  # Concordant (partial match)
                else:
                    concordance_dict[gene] = 0  # Non-concordant
        else:
            concordance_dict[gene] = 0  # Non-concordant if parent's data is missing
    
    return concordance_dict

def calculate_trio_concordance(trio_df, genotype_dfs):
    """
    Calculate concordance for all trios.
    
    Args:
        trio_df (pd.DataFrame): DataFrame containing trio information.
        genotype_dfs (dict): Dictionary with genotype dataframes for each individual.
        
    Returns:
        trio_concordance (defaultdict): Dictionary with gene as key and list of concordance values as value.
        trio_concordance_rates (dict): Dictionary with trio ID as key and concordance rate as value.
    """
    trio_concordance = defaultdict(list)
    trio_concordance_rates = {}  # Will store concordance per trio

    for _, trio in trio_df.iterrows():
        child_id = trio['sampleID']
        father_id = trio['fatherID']
        mother_id = trio['motherID']
        trio_id = "-".join([child_id, father_id, mother_id])
        
        if father_id != '0' and mother_id != '0':  # Only process trios with non-zero parents
            concordance = calculate_mendelian_concordance(child_id, father_id, mother_id, genotype_dfs)
            for gene, is_concordant in concordance.items():
                trio_concordance[gene].append(is_concordant)
            
            # Calculate concordance for the trio (across all genes)
            trio_level_concordance = sum(concordance.values()) / len(concordance) if concordance else 0
            if trio_level_concordance > 0:
                trio_concordance_rates[trio_id] = trio_level_concordance
                
    return trio_concordance, trio_concordance_rates


def load_novel_variants(variants_dir):
    """
    Load novel variants TSV files into a dictionary of dataframes.
    
    Args:
        variants_dir (str): Directory containing the .final-novel-variants.tsv files
        
    Returns:
        dict: Dictionary with sample IDs as keys and variant dataframes as values
    """
    variant_dfs = {}
    
    # Get all TSV files in the directory
    for file in Path(variants_dir).glob("*.final-novel-variants.tsv"):
        sample_id = file.stem.split('.')[0]  # Extract sample ID from filename
        
        # Read TSV file
        df = pd.read_csv(file, sep='\t')
        
        # Group variants by gene to create a more useful structure
        variant_dict = {}
        for gene, group in df.groupby('gene'):
            # Create variant identifiers combining position and change
            variants = []
            for _, row in group.iterrows():
                if len(row['ref']) == len(row['alt']):  # SNP
                    variant = f"{row['position']}_{row['ref']}>{row['alt']}"
                else:  # INDEL
                    variant = f"{row['position']}_{row['ref']}->{row['alt']}"
                variants.append(variant)
            variant_dict[gene] = variants
            
        variant_dfs[sample_id] = pd.DataFrame({'Variants': variant_dict}).T

    return variant_dfs

def calculate_variant_concordance(child_id, father_id, mother_id, variant_dfs):
    """
    Calculate Mendelian concordance for novel variants in a trio.
    
    Args:
        child_id (str): The child's sample ID
        father_id (str): The father's sample ID
        mother_id (str): The mother's sample ID
        variant_dfs (dict): Dictionary with variant dataframes for each individual
        
    Returns:
        dict: Dictionary with gene as key and concordance (0/1) as value.
               Returns empty dict if child is missing.
               Treats missing parents as having no variants.
    """
    concordance_dict = {}
    
    # If child is missing, we can't calculate concordance
    if child_id not in variant_dfs:
        return concordance_dict
    
    child_variants = variant_dfs[child_id]
    
    # Get parent variants if available, otherwise use empty DataFrames
    father_variants = variant_dfs.get(father_id, pd.DataFrame({'Variants': []}).T)
    mother_variants = variant_dfs.get(mother_id, pd.DataFrame({'Variants': []}).T)
    
    # Iterate through each gene in the child's variants
    for gene in child_variants.columns:
        try:
            child_vars = set(child_variants[gene].loc['Variants'])
        except KeyError:
            print(child_id)
            print(gene)
            raise
        
        # Get parent variants for this gene, or empty set if gene not present
        father_vars = set(father_variants[gene].loc['Variants']) if gene in father_variants.columns else set()
        mother_vars = set(mother_variants[gene].loc['Variants']) if gene in mother_variants.columns else set()
        
        # Check if each child variant is present in either parent
        concordant = True
        for variant in child_vars:
            if variant not in father_vars and variant not in mother_vars:
                concordant = False
                break
        
        concordance_dict[gene] = 1 if concordant else 0
            
    return concordance_dict

def calculate_trio_variant_concordance(trio_df, variant_dfs):
    """
    Calculate variant concordance across all trios and generate summary statistics.
    
    Args:
        trio_df (pd.DataFrame): DataFrame containing trio information
        variant_dfs (dict): Dictionary of variant dataframes for each sample
        
    Returns:
        tuple: (concordance_df, trio_concordance_rates, overall_concordance_rate)
            - concordance_df: DataFrame with gene-level concordance rates
            - trio_concordance_rates: Dictionary of concordance rates per trio
            - overall_concordance_rate: Float representing overall concordance
    """
    trio_concordance = defaultdict(list)
    trio_concordance_rates = {}  # Will store concordance per trio
    
    complete_trio_count = 0
    for _, trio in trio_df.iterrows():
        child_id = trio['sampleID']
        father_id = trio['fatherID']
        mother_id = trio['motherID']
        trio_id = "-".join([child_id, father_id, mother_id])
        
        if father_id != '0' and mother_id != '0':  # Only process trios with non-zero parents
            concordance = calculate_variant_concordance(child_id, father_id, mother_id, variant_dfs)
            for gene, is_concordant in concordance.items():
                trio_concordance[gene].append(is_concordant)
            
            # Calculate concordance for the trio (across all genes)
            trio_level_concordance = sum(concordance.values()) / len(concordance) if concordance else 0
            if concordance:
                complete_trio_count += 1
                trio_concordance_rates[trio_id] = trio_level_concordance
    
    # Summarize concordance per gene
    concordance_summary = {gene: sum(values)/len(values) 
                         for gene, values in trio_concordance.items()}
    concordance_df = pd.DataFrame.from_dict(concordance_summary, 
                                          orient='index', 
                                          columns=['Variant_Concordance'])
    
    # Calculate overall concordance rate
    overall_concordance_rate = sum(trio_concordance_rates.values()) / len(trio_concordance_rates)
    
    print(f"Processed concordance from {complete_trio_count} trios")
    return concordance_df, trio_concordance_rates, overall_concordance_rate

def calculate_total_gene_copies(genotype_dfs):
    """
    Calculate the total number of distinct gene copies (alleles) for each sample.
    
    Args:
        genotype_dfs (dict): Dictionary of genotype DataFrames for each sample.
    
    Returns:
        total_gene_copies (dict): Dictionary with sample IDs as keys and total number of gene copies as values.
    """
    total_gene_copies = {}
    
    for sample_id, df in genotype_dfs.items():
        total_copies = df['Allele'].apply(lambda alleles: len(set(alleles))).sum()
        total_gene_copies[sample_id] = total_copies
    
    return total_gene_copies

def map_samples_to_trios(trio_df):
    """
    Create a mapping from each sample ID to the corresponding trio (child-father-mother).
    
    Args:
        trio_df (pd.DataFrame): DataFrame containing the trio information (sampleID, fatherID, motherID).
    
    Returns:
        sample_to_trio (dict): Dictionary mapping sample IDs to trio strings (e.g., "child-father-mother").
    """
    sample_to_trio = {}
    
    for _, row in trio_df.iterrows():
        child_id = row['sampleID']
        father_id = row['fatherID']
        mother_id = row['motherID']
        
        if father_id != '0' and mother_id != '0':  # Only process valid trios
            trio_string = f"{child_id}-{father_id}-{mother_id}"
            sample_to_trio[child_id] = trio_string
            sample_to_trio[father_id] = trio_string
            sample_to_trio[mother_id] = trio_string
    
    return sample_to_trio





def analyze_copies_vs_concordance(total_gene_copies, trio_concordance_rates, sample_to_trio):
    """
    Analyze the relationship between the total gene copies per sample and the trio concordance rate.
    
    Args:
        total_gene_copies (dict): Dictionary with sample IDs as keys and total gene copies as values.
        trio_concordance_rates (dict): Dictionary with trio strings as keys and trio concordance rates as values.
        sample_to_trio (dict): Dictionary mapping sample IDs to trio strings.
    
    Returns:
        analysis_df (pd.DataFrame): DataFrame containing the total gene copies and trio concordance rate per sample.
    """
    analysis_data = []
    
    for sample_id, total_copies in total_gene_copies.items():
        trio_string = sample_to_trio.get(sample_id)
        if trio_string in trio_concordance_rates:
            trio_concordance_rate = trio_concordance_rates[trio_string]
            analysis_data.append({
                'SampleID': sample_id,
                'Trio': trio_string,
                'Total_Copies': total_copies,
                'Trio_Concordance_Rate': trio_concordance_rate
            })
    
    analysis_df = pd.DataFrame(analysis_data)
    return analysis_df

def identify_outliers_copies_vs_concordance(analysis_df):
    """
    Identify outliers where the number of gene copies or concordance rate is unusually high or low.
    
    Args:
        analysis_df (pd.DataFrame): DataFrame containing the total gene copies and trio concordance rate per sample.

    Returns:
        outliers_df (pd.DataFrame): DataFrame containing outliers based on total gene copies and concordance rate.
    """
    # Define thresholds for outliers based on quantiles
    high_copy_threshold = analysis_df['Total_Copies'].quantile(0.95)  # Top 5% highest gene copies
    low_copy_threshold = analysis_df['Total_Copies'].quantile(0.05)   # Bottom 5% lowest gene copies
    low_concordance_threshold = analysis_df['Trio_Concordance_Rate'].quantile(0.05)  # Lowest 5% concordance
    
    # Identify outliers based on copies and concordance rate
    outliers_df = analysis_df[
        ((analysis_df['Total_Copies'] > high_copy_threshold) | 
         (analysis_df['Total_Copies'] < low_copy_threshold)) |
        (analysis_df['Trio_Concordance_Rate'] < low_concordance_threshold)
    ]
    
    return outliers_df

def get_gene_concordance_data(gene_types, genotype_dirs, trio_df):
    """
    Collect concordance and incidence data for different gene types.
    
    Args:
        gene_types (list): List of gene types (e.g., ['TRBV', 'IGHV'])
        genotype_dirs (list): List of directories containing genotype files, matching gene_types order
        trio_df (pd.DataFrame): DataFrame containing trio information
    
    Returns:
        pd.DataFrame: DataFrame containing concordance and incidence data with columns:
            - Gene Type: The type of gene
            - Concordance: Mean concordance value
            - Gene: Gene name
            - Incidence_Count: Number of samples containing the gene
            - Incidence_Percentage: Percentage of samples containing the gene
            - Total_Samples: Total number of samples analyzed
    """
    plot_data = []
    
    for gene_type, genotype_dir in zip(gene_types, genotype_dirs):
        # Load genotype data for this gene type
        genotype_dfs = load_genotype_data(genotype_dir, gene_type)
        
        # Calculate concordance
        trio_concordance, _ = calculate_trio_concordance(trio_df, genotype_dfs)
        
        # Calculate incidence
        total_samples = len(genotype_dfs)
        gene_incidence = defaultdict(int)
        
        for sample_df in genotype_dfs.values():
            for gene in sample_df.index:
                gene_incidence[gene] += 1
        
        # For each gene, get its concordance values and incidence data
        for gene, concordance_values in trio_concordance.items():
            incidence = gene_incidence.get(gene, 0)
            plot_data.append({
                'Gene Type': gene_type,
                'Concordance': mean(concordance_values),
                'Gene': gene,
                'Incidence_Count': incidence,
                'Incidence_Percentage': (incidence / total_samples) * 100,
                'Total_Samples': total_samples
            })
    
    return pd.DataFrame(plot_data)

def plot_gene_concordance_distribution(concordance_df):
    """
    Generate a box plot comparing per-gene Mendelian concordance across different gene types.
    
    Args:
        concordance_df (pd.DataFrame): DataFrame containing concordance data with columns:
            - Gene Type: The type of gene
            - Concordance: Mean concordance value
            - Gene: Gene name
    
    Returns:
        fig: matplotlib Figure object
        ax: matplotlib Axes object
    """
    # Create figure and axis objects
    fig, ax = plt.subplots(figsize=(20, 6))
    
    # Create the violin plot with different colors for each gene type
    sns.boxplot(data=concordance_df, x='Gene Type', y='Concordance', 
                palette='deep', ax=ax)  # Added palette parameter
    
    # Add individual points and store the scatter plot object
    scatter = sns.stripplot(data=concordance_df, x='Gene Type', y='Concordance', 
                          color='red', alpha=0.5, size=4, jitter=0.2, ax=ax)
    
    # Calculate outlier threshold for each gene type
    for gene_type_idx, gene_type in enumerate(concordance_df['Gene Type'].unique()):
        gene_type_data = concordance_df[concordance_df['Gene Type'] == gene_type]
        mean_concordance = gene_type_data['Concordance'].mean()
        std_concordance = gene_type_data['Concordance'].std()
        low_threshold = mean_concordance - std_concordance
        
        # Get outlier points
        outliers = gene_type_data[gene_type_data['Concordance'] <= low_threshold]
        
        # Find the collection for this gene type by matching number of points
        expected_points = len(gene_type_data)
        for collection in scatter.collections:
            if isinstance(collection, plt.matplotlib.collections.PathCollection):
                points = collection.get_offsets()
                if len(points) == expected_points:
                    # Label outliers using points from this collection
                    for _, row in outliers.iterrows():
                        # Find matching points by y-coordinate
                        matching_points = points[np.abs(points[:, 1] - row['Concordance']) < 0.001]
                        if len(matching_points) > 0:
                            point_pos = matching_points[0]
                            
                            # Determine if point is left or right of midline
                            relative_pos = point_pos[0] - gene_type_idx
                            
                            # Remove first 4 characters from gene name
                            short_label = row['Gene'][4:]
                            
                            if relative_pos < 0:
                                offset = 0.001 * len(short_label)  # Scale offset by shortened text length
                                ax.text(point_pos[0] - offset, point_pos[1], 
                                      short_label,
                                      horizontalalignment='right',
                                      verticalalignment='center',
                                      fontsize=13)
                            else:
                                ax.text(point_pos[0] + 0.01, point_pos[1], 
                                      short_label,
                                      horizontalalignment='left',
                                      verticalalignment='center',
                                      fontsize=13)

    # Customize the plot
    ax.set_xlabel('Gene Type', fontsize=17)
    ax.tick_params(axis='x', labelsize=17)
    
    # Fix the warning by getting current tick positions first
    current_labels = [label.get_text().upper() for label in ax.get_xticklabels()]
    ax.set_xticks(ax.get_xticks())  # Set tick positions explicitly
    ax.set_xticklabels(current_labels)  # Now set the labels
    
    ax.set_ylabel('Concordance Rate', fontsize=17)
    ax.set_ylim(0, 1)
    
    # Adjust layout
    fig.tight_layout()
    
    return fig, ax

def analyze_gene_incidence(concordance_df):
    """
    Analyze and plot the relationship between gene incidence and concordance for all genes.
    
    Args:
        concordance_df (pd.DataFrame): DataFrame containing concordance and incidence data
        
    Returns:
        tuple: (fig, ax, analysis_df) where:
            - fig: matplotlib Figure object
            - ax: matplotlib Axes object
            - analysis_df: DataFrame containing analysis results
    """
    # Create scatter plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create scatter plot with different colors for each gene type
    for gene_type in concordance_df['Gene Type'].unique():
        mask = concordance_df['Gene Type'] == gene_type
        ax.scatter(concordance_df[mask]['Incidence_Percentage'], 
                  concordance_df[mask]['Concordance'],
                  label=gene_type.upper(),
                  alpha=0.7)
    
    # Add labels for points with low concordance or unusual incidence
    concordance_threshold = concordance_df['Concordance'].mean() - concordance_df['Concordance'].std()
    incidence_threshold = concordance_df['Incidence_Percentage'].mean() - concordance_df['Incidence_Percentage'].std()
    
    for _, row in concordance_df.iterrows():
        if row['Concordance'] < concordance_threshold or row['Incidence_Percentage'] < incidence_threshold:
            # Remove gene type prefix for cleaner labels
            short_label = row['Gene'][4:]
            ax.annotate(short_label,
                       (row['Incidence_Percentage'], row['Concordance']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=13)  # Increased from 10 to 13
    
    # Calculate linear regression and add trend line with confidence interval
    slope, intercept, corr, p_value, std_err = stats.linregress(
        concordance_df['Incidence_Percentage'],
        concordance_df['Concordance']
    )
    x_trend = np.linspace(concordance_df['Incidence_Percentage'].min(),
                         concordance_df['Incidence_Percentage'].max(), 100)
    y_trend = slope * x_trend + intercept
    
    # Calculate 95% confidence interval
    n = len(concordance_df)
    x_mean = concordance_df['Incidence_Percentage'].mean()
    
    # Sum of squares of x deviations
    x_dev = concordance_df['Incidence_Percentage'] - x_mean
    ss_xx = np.sum(x_dev ** 2)
    
    # Calculate standard error of regression line
    y_hat = slope * concordance_df['Incidence_Percentage'] + intercept
    mse = np.sum((concordance_df['Concordance'] - y_hat) ** 2) / (n - 2)
    
    # Calculate confidence intervals
    ci = stats.t.ppf(0.975, n - 2) * np.sqrt(mse * (1/n + (x_trend - x_mean)**2 / ss_xx))
    
    # Plot regression line and confidence interval
    ax.plot(x_trend, y_trend, "r-", alpha=0.8, label='Regression line')
    ax.fill_between(x_trend, y_trend - ci, y_trend + ci, color='red', alpha=0.1, label='95% CI')
    
    # Add correlation information (bottom right)
    ax.text(0.95, 0.05,
            f'Correlation: {corr:.3f}\n$p$-value: {p_value:.3e}',
            transform=ax.transAxes,
            horizontalalignment='right',
            verticalalignment='bottom',
            fontsize=13,  # Added fontsize (increased by 3 from default ~10)
            bbox=dict(facecolor='white', alpha=0.8))
    
    # Customize plot
    ax.set_xlabel('Gene Incidence (%)', fontsize=17)  # Increased from 14 to 17
    ax.set_ylabel('Concordance Rate', fontsize=17)    # Increased from 14 to 17
    # ax.set_title('Relationship between Gene Incidence and Concordance')
    ax.legend(fontsize=13)  # Added fontsize (increased by 3 from default ~10)
    ax.grid(True, alpha=0.3)
    
    # Increase tick label font size
    ax.tick_params(axis='both', which='major', labelsize=13)  # Added for tick labels
    
    # Adjust layout
    fig.tight_layout()
    
    return fig, ax, concordance_df
def get_trio_concordance_data(gene_types, genotype_dirs, trio_df):
    """
    Collect trio-level concordance data for different gene types.
    
    Args:
        gene_types (list): List of gene types (e.g., ['TRBV', 'IGHV'])
        genotype_dirs (list): List of directories containing genotype files, matching gene_types order
        trio_df (pd.DataFrame): DataFrame containing trio information
    
    Returns:
        pd.DataFrame: DataFrame containing concordance data with columns:
            - Gene Type: The type of gene
            - Concordance: Trio-level concordance rate
            - Trio: Trio identifier
    """
    plot_data = []
    
    for gene_type, genotype_dir in zip(gene_types, genotype_dirs):
        # Load genotype data for this gene type
        genotype_dfs = load_genotype_data(genotype_dir, gene_type)
        
        # Calculate concordance
        _, trio_concordance_rates = calculate_trio_concordance(trio_df, genotype_dfs)
        
        # Add each trio's concordance rate to plot data
        for trio_id, concordance_rate in trio_concordance_rates.items():
            plot_data.append({
                'Gene Type': gene_type,
                'Concordance': concordance_rate,
                'Trio': trio_id
            })
    
    return pd.DataFrame(plot_data)

def plot_trio_concordance_distribution(concordance_df):
    """
    Generate a box plot with strip points comparing trio-level concordance across different gene types.
    
    Args:
        concordance_df (pd.DataFrame): DataFrame containing concordance data with columns:
            - Gene Type: The type of gene
            - Concordance: Trio-level concordance rate
            - Trio: Trio identifier
    
    Returns:
        fig: matplotlib Figure object
        ax: matplotlib Axes object
    """
    # Create figure and axis objects
    fig, ax = plt.subplots(figsize=(20, 6))
    
    # Create the box plot
    sns.boxplot(data=concordance_df, x='Gene Type', y='Concordance', ax=ax, palette='deep')
    
    # Add individual points
    sns.stripplot(data=concordance_df, x='Gene Type', y='Concordance', 
                 color='red', alpha=0.5, size=4, jitter=0.2, ax=ax)
    
    # Customize the plot
    # ax.set_title('Distribution of Trio-Level Concordance by Gene Type')
    ax.set_xlabel('Gene Type', fontsize=17)
    ax.set_ylabel('Concordance Rate', fontsize=17)
    ax.tick_params(axis='x', labelsize=17)
    ax.set_xticklabels([label.get_text().upper() for label in ax.get_xticklabels()])
    ax.set_ylim(0, 1)
    
    # Adjust layout
    fig.tight_layout()
    
    return fig, ax

def load_sample_metadata(genotype_dir, gene_type, high_coverage_index, related_index):
    """
    Load read statistics and 1000G metadata for samples.
    
    Args:
        genotype_dir (str): Directory containing genotype and debug log files
        gene_type (str): Type of gene (e.g., 'IGHV', 'TRBV')
        high_coverage_index (str): Path to 1000G_2504_high_coverage.sequence.index
        related_index (str): Path to 1000G_698_related_high_coverage.sequence.index
        
    Returns:
        pd.DataFrame: DataFrame with columns:
            - sample_id: Sample identifier
            - recruited_reads: Number of reads that passed filters
            - total_reads: Total number of reads processed
            - read_depth: Estimated single copy read depth
            - population: Sample population/ethnicity
            - sequencing_reads: Total sequencing read count
            - sequencing_bases: Total sequencing base count
            - source_files: Which 1000G index files contained this sample
    """
    # Load read stats from debug logs
    stats_data = []
    debug_files = glob(f"{genotype_dir}*-{gene_type.lower()}-immunotyper-debug.log")
    
    for file_path in debug_files:
        sample_id = os.path.basename(file_path).split('-')[0].split('.')[0]
        recruited_reads = None
        total_reads = None
        read_depth = None
        
        with open(file_path, 'r') as f:
            for line in f:
                # Parse read recruitment info
                if "passed and" in line and "failed filter" in line:
                    parts = line.split()
                    total_reads = int(parts[parts.index("Of") + 1].rstrip(','))
                    recruited_reads = int(parts[parts.index("Of") + 2])
                
                # Parse read depth info
                if "Estimated single copy read depth:" in line:
                    read_depth = int(line.split()[-1])
                
                if recruited_reads is not None and read_depth is not None:
                    break
        
        if recruited_reads is not None and read_depth is not None:
            stats_data.append({
                'sample_id': sample_id,
                'recruited_reads': recruited_reads,
                'total_reads': total_reads,
                'read_depth': read_depth
            })
    
    # Create DataFrame from read stats
    stats_df = pd.DataFrame(stats_data)
    
    # Load and process 1000G metadata
    metadata_df = load_1000g_metadata(high_coverage_index, related_index)
    
    # Extract relevant columns and rename for clarity
    metadata_subset = metadata_df[[
        'SAMPLE_NAME',
        'POPULATION',
        'READ_COUNT',
        'BASE_COUNT',
        'source_files'
    ]].rename(columns={
        'SAMPLE_NAME': 'sample_id',
        'POPULATION': 'population',
        'READ_COUNT': 'sequencing_reads',
        'BASE_COUNT': 'sequencing_bases'
    })
    
    # Merge read stats with metadata
    merged_df = pd.merge(
        stats_df,
        metadata_subset,
        on='sample_id',
        how='left'  # Keep all samples from stats_df even if not in metadata
    )
    
    # Add warning for samples without metadata
    missing_metadata = merged_df[merged_df['population'].isna()]['sample_id'].tolist()
    if missing_metadata:
        print(f"\nWarning: No 1000G metadata found for {len(missing_metadata)} samples: "
              f"{', '.join(missing_metadata[:5])}{'...' if len(missing_metadata) > 5 else ''}")
    
    # Report samples with metadata but no read stats
    extra_metadata = set(metadata_subset['sample_id']) - set(stats_df['sample_id'])
    if extra_metadata:
        print(f"\nNote: {len(extra_metadata)} samples in 1000G metadata but not in read stats: "
              f"{', '.join(sorted(list(extra_metadata))[:5])}{'...' if len(extra_metadata) > 5 else ''}")
    
    return merged_df

def load_1000g_metadata(high_coverage_index: str, related_index: str) -> pd.DataFrame:
    """
    Load and combine metadata from 1000 Genomes Project index files.
    
    Args:
        high_coverage_index: Path to 1000G_2504_high_coverage.sequence.index
        related_index: Path to 1000G_698_related_high_coverage.sequence.index
        
    Returns:
        pd.DataFrame: Combined metadata with columns:
            - SAMPLE_NAME: Sample identifier
            - POPULATION: Sample population/ethnicity
            - READ_COUNT: Total sequencing read count
            - BASE_COUNT: Total sequencing base count
            - source_files: Which index files contained this sample
    """
    def read_index_file(file_path):
        # Read comments to get column definitions
        comments = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    comments.append(line.strip())
                elif line.startswith('#'):
                    # This is the header line - remove the # and split
                    header = line[1:].strip().split('\t')
                    break
        
        # Read the data using pandas, skipping comment lines
        df = pd.read_csv(file_path, 
                        sep='\t',
                        comment='#',
                        names=header,
                        dtype={'READ_COUNT': 'Int64', 'BASE_COUNT': 'Int64'},
                        index_col=False)
        
        # Select only the columns we need
        return df
    
    # Read both index files
    high_coverage_df = read_index_file(high_coverage_index)
    high_coverage_df['source_files'] = 'high_coverage'
    
    related_df = read_index_file(related_index)
    related_df['source_files'] = 'related'
    
    # Merge dataframes on SAMPLE_NAME
    combined_df = pd.concat([
        high_coverage_df,
        related_df[~related_df['SAMPLE_NAME'].isin(high_coverage_df['SAMPLE_NAME'])]
    ], ignore_index=True)
    
    return combined_df

def plot_concordance_vs_read_variation(trio_concordance_rates, read_stats_df, method='cv', ax=None):
    """
    Plot relationship between trio concordance and read recruitment variation using different metrics.
    
    Args:
        trio_concordance_rates (dict): Dictionary with trio ID as key and concordance rate as value
        read_stats_df (pd.DataFrame): DataFrame containing read statistics
        method (str): One of 'cv', 'pairwise', or 'spread'
    
    Returns:
        tuple: (fig, ax) matplotlib figure and axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure
        
    def get_max_pairwise_difference(trio_samples):
        """Calculate maximum pairwise relative difference between samples"""
        ratios = trio_samples['recruitment_ratio'].values
        max_diff = 0
        for i, j in itertools.combinations(range(len(ratios)), 2):
            diff = abs(ratios[i] - ratios[j]) / min(ratios[i], ratios[j])
            max_diff = max(max_diff, diff)
        return max_diff
        
    read_stats_df['recruitment_ratio'] = read_stats_df['recruited_reads'] / read_stats_df['read_depth']
    
    variation_metrics = {
        'cv': lambda x: x['recruitment_ratio'].std() / x['recruitment_ratio'].mean(),
        'pairwise': get_max_pairwise_difference,
        'spread': lambda x: (
            x['recruitment_ratio'].quantile(0.75) - x['recruitment_ratio'].quantile(0.25)
        ) / x['recruitment_ratio'].quantile(0.5)
    }
    
    metric_labels = {
        'cv': 'Coefficient of Variation',
        'pairwise': 'Maximum Pairwise Relative Difference',
        'spread': 'Normalized IQR Spread'
    }
    
    plot_data = []
    for trio_id, concordance_rate in trio_concordance_rates.items():
        child_id, father_id, mother_id = trio_id.split('-')
        trio_samples = read_stats_df[read_stats_df['sample_id'].isin([child_id, father_id, mother_id])]
        
        if not trio_samples.empty:
            variation = variation_metrics[method](trio_samples)
            plot_data.append({
                'trio_id': trio_id,
                'concordance': concordance_rate,
                'variation': variation * 100  # Convert to percentage
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create plot with regression line and confidence interval
    sns.scatterplot(data=plot_df, x='variation', y='concordance', ax=ax)
    
    # Add regression analysis
    x = plot_df['variation']
    y = plot_df['concordance']
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    
    # Plot regression line
    ax.plot(x, line, 'r--', label='Regression line')
    
    # Calculate and plot confidence interval
    def confidence_interval(x, y, slope, intercept, p=0.95):
        n = len(x)
        mean_x = np.mean(x)
        ss_x = np.sum((x - mean_x)**2)
        y_hat = slope * x + intercept
        sr = np.sqrt(np.sum((y - y_hat)**2) / (n-2))
        se_line = sr * np.sqrt(1/n + (x - mean_x)**2/ss_x)
        t = stats.t.ppf((1 + p)/2, n-2)
        ci = t * se_line
        return y_hat - ci, y_hat + ci

    lower_ci, upper_ci = confidence_interval(x, y, slope, intercept)
    ax.fill_between(x, lower_ci, upper_ci, color='red', alpha=0.1, label='95% CI')
    
    # Customize plot with increased font sizes
    ax.set_xlabel(f'{metric_labels[method]} (%)', fontsize=17)
    ax.set_ylabel('Trio Concordance Rate', fontsize=17)
    ax.set_title(f'Trio Concordance vs {metric_labels[method]}\n' +
                 f'Correlation: {r_value:.3f} (p={p_value:.2e})', fontsize=15)
    ax.legend(fontsize=12)
    
    # Increase tick label font sizes
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    return fig, ax

def plot_concordance_vs_child_parent_reduction(trio_concordance_rates, read_stats_df, ax=None):
    """
    Plot relationship between trio concordance and maximum reduction in read recruitment 
    ratio between child and either parent, as a percentage of child's ratio.
    
    Args:
        trio_concordance_rates (dict): Dictionary with trio ID as key and concordance rate as value
        read_stats_df (pd.DataFrame): DataFrame containing read statistics
    
    Returns:
        tuple: (fig, ax) matplotlib figure and axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    else:
        fig = ax.figure
        
    def get_max_child_parent_reduction(child_ratio, parent1_ratio, parent2_ratio):
        """Calculate maximum reduction in ratio between child and either parent,
        as a percentage of child's ratio"""
        # Calculate reduction relative to child for each parent
        # Positive values mean parent has higher ratio than child
        reduction1 = -(child_ratio - parent1_ratio) / child_ratio
        reduction2 = -(child_ratio - parent2_ratio) / child_ratio
        # Return the larger reduction
        return reduction1 if abs(reduction1) > abs(reduction2) else reduction2
    
    read_stats_df['recruitment_ratio'] = read_stats_df['recruited_reads'] / read_stats_df['read_depth']
    
    plot_data = []
    for trio_id, concordance_rate in trio_concordance_rates.items():
        child_id, father_id, mother_id = trio_id.split('-')
        trio_samples = read_stats_df[read_stats_df['sample_id'].isin([child_id, father_id, mother_id])]
        
        if not trio_samples.empty and len(trio_samples) == 3:  # Ensure we have all three samples
            child_ratio = trio_samples[trio_samples['sample_id'] == child_id]['recruitment_ratio'].iloc[0]
            father_ratio = trio_samples[trio_samples['sample_id'] == father_id]['recruitment_ratio'].iloc[0]
            mother_ratio = trio_samples[trio_samples['sample_id'] == mother_id]['recruitment_ratio'].iloc[0]
            
            max_reduction = get_max_child_parent_reduction(child_ratio, father_ratio, mother_ratio)
            plot_data.append({
                'trio_id': trio_id,
                'concordance': concordance_rate,
                'max_reduction': max_reduction * 100  # Convert to percentage
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create plot with regression line and confidence interval
    sns.scatterplot(data=plot_df, x='max_reduction', y='concordance', ax=ax)
    
    # Add regression analysis
    x = plot_df['max_reduction']
    y = plot_df['concordance']
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    
    # Plot regression line
    ax.plot(x, line, 'r--', label='Regression line')
    
    # Calculate and plot confidence interval
    def confidence_interval(x, y, slope, intercept, p=0.95):
        n = len(x)
        mean_x = np.mean(x)
        ss_x = np.sum((x - mean_x)**2)
        y_hat = slope * x + intercept
        sr = np.sqrt(np.sum((y - y_hat)**2) / (n-2))
        se_line = sr * np.sqrt(1/n + (x - mean_x)**2/ss_x)
        t = stats.t.ppf((1 + p)/2, n-2)
        ci = t * se_line
        return y_hat - ci, y_hat + ci

    lower_ci, upper_ci = confidence_interval(x, y, slope, intercept)
    ax.fill_between(x, lower_ci, upper_ci, color='red', alpha=0.1, label='95% CI')
    
    # Customize plot with increased font sizes (+3 from original)
    ax.set_xlabel('Maximum Child-Parent Difference in Read Recruitment Ratio (%)', fontsize=17)
    ax.set_ylabel('Trio Concordance Rate', fontsize=17)
    
    # Set tick label font sizes (+3 from default ~10)
    ax.tick_params(axis='both', which='major', labelsize=13)
    
    print('Trio Concordance vs Maximum Parent-to-Child Read Recruitment Ratio Increase\n' +
                 f'Correlation: {r_value:.3f} (p={p_value:.2e})')
    
    # Set legend font size (+3 from default ~10)
    ax.legend(fontsize=13)
    
    # Add vertical line at x=0 to help visualize positive vs negative differences
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    return fig, ax

POPULATION_GROUPS = {
    'African': ['YRI', 'ESN', 'GWD', 'MSL', 'ACB', 'ASW'],  # West African and African American
    'European': ['CEU', 'IBS'],  # European
    'South Asian': ['PJL', 'BEB', 'STU', 'ITU'],  # South Asian
    'East Asian': ['CHS', 'KHV'],  # East Asian
    'American': ['PUR', 'CLM', 'PEL', 'MXL']  # Latino/Admixed American
}

def get_population_group(population):
    """Helper function to get the group for a given population"""
    if pd.isna(population):
        return 'Unknown'
    for group, populations in POPULATION_GROUPS.items():
        if population in populations:
            return group
    return 'Unknown'


def plot_concordance_by_population(trio_concordance_rates, metadata_df):
    """
    Create boxplot visualization of concordance rates by population, grouped by continental ancestry.
    """
    # Prepare data for plotting with group information
    plot_data = []
    for trio_id, concordance in trio_concordance_rates.items():
        child_id = trio_id.split('-')[0]
        population = metadata_df[metadata_df['sample_id'] == child_id]['population'].iloc[0] \
            if child_id in metadata_df['sample_id'].values else None
        group = get_population_group(population)
        
        plot_data.append({
            'trio_id': trio_id,
            'population': population,
            'group': group,
            'concordance': concordance
        })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Sort populations by group for adjacent plotting
    population_order = []
    group_colors = {
        'African': '#FFA07A',    # Light salmon
        'European': '#98FB98',   # Pale green
        'South Asian': '#87CEFA', # Light sky blue
        'East Asian': '#DDA0DD',  # Plum
        'American': '#F0E68C',    # Khaki
        'Unknown': '#D3D3D3'     # Light gray
    }
    
    for group in group_colors.keys():
        group_pops = [pop for pop in plot_df['population'].unique() 
                     if get_population_group(pop) == group]
        population_order.extend(sorted(group_pops))
    
    # Calculate statistics with grouped ordering
    stats_df = plot_df.groupby('population').agg({
        'concordance': ['count', 'mean', 'std', 'median', 'min', 'max'],
        'trio_id': 'count'
    }).round(3)
    stats_df.columns = ['concordance_count', 'concordance_mean', 'concordance_std', 
                       'concordance_median', 'concordance_min', 'concordance_max', 
                       'n_trios']
    stats_df = stats_df.reindex(population_order)
    
    # Create figure
    fig, (ax_top, ax_main) = plt.subplots(
        2, 1, 
        gridspec_kw={'height_ratios': [1, 3]},
        figsize=(15, 10)  # Made wider to accommodate group labels
    )
    fig.subplots_adjust(hspace=0.1)
    
    # Plot sample sizes with group colors
    x_positions = np.arange(len(population_order))
    sample_sizes = stats_df['n_trios']
    
    for i, pop in enumerate(population_order):
        group = get_population_group(pop)
        ax_top.bar([i], sample_sizes[pop], width=0.8, 
                  color=group_colors[group], alpha=0.6)
    
    ax_top.set_xlim(-0.5, len(population_order) - 0.5)
    ax_top.set_xticks([])
    ax_top.set_ylabel('Number of Trios', fontsize=15)  # Increased from 12 to 15
    ax_top.tick_params(axis='y', labelsize=12)  # Added tick label font size
    
    # Create color palette mapping each population to its group color
    pop_colors = {pop: group_colors[get_population_group(pop)] 
                 for pop in population_order}
    
    # Create boxplot with explicit color palette
    bp = sns.boxplot(data=plot_df, x='population', y='concordance', ax=ax_main, 
                    order=population_order, whis=1.5,
                    palette=pop_colors)
    
    # Add horizontal mean lines for each group
    for group in group_colors.keys():
        if group != 'Unknown':  # Skip Unknown group
            group_pops = [pop for pop in population_order if get_population_group(pop) == group]
            if group_pops:
                # Calculate group mean
                group_data = plot_df[plot_df['population'].isin(group_pops)]['concordance']
                mean_value = group_data.mean()
                
                # Find x-coordinates for line
                start_idx = population_order.index(group_pops[0])
                end_idx = population_order.index(group_pops[-1])
                
                # Add horizontal line with consistent, dark color
                ax_main.hlines(mean_value, start_idx - 0.4, end_idx + 0.4,
                             colors='#404040',  # Dark gray
                             linestyles='dashed',
                             linewidth=2.5)  # Removed label here
    
    # Add swarmplot
    sns.swarmplot(data=plot_df, x='population', y='concordance', ax=ax_main,
                  order=population_order, size=4, alpha=0.5, color='0.3')
    
    # Customize main plot
    ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=45, ha='right')
    ax_main.set_xlabel('Population', fontsize=15)  # Increased from 12 to 15
    ax_main.set_ylabel('Trio Concordance Rate', fontsize=15)  # Increased from 12 to 15
    ax_main.tick_params(axis='both', labelsize=12)  # Added tick label font size
    
    # Add group separators and labels
    current_x = -0.5
    for group in group_colors.keys():
        group_pops = [pop for pop in population_order 
                     if get_population_group(pop) == group]
        if group_pops:
            group_width = len(group_pops)
            # Add separator line
            ax_main.axvline(x=current_x, color='black', linestyle='--', alpha=0.3)
            # Add group label
            ax_main.text(current_x + group_width/2, ax_main.get_ylim()[1],
                        group, ha='center', va='bottom', fontsize=13)  # Added fontsize
            current_x += group_width
    
    # Create legend elements
    legend_elements = [
        plt.Rectangle((0,0),1,1, facecolor=color, alpha=0.6, label=group)
        for group, color in group_colors.items()
        if group != 'Unknown'  # Skip Unknown group
    ]
    # Add single mean line entry
    legend_elements = [(plt.Line2D([0], [0], color='#404040', linestyle='--', 
                                    label='Mean Continental Concordance', linewidth=2.5))]
    
    # Add single combined legend in bottom right
    ax_main.legend(handles=legend_elements, 
                #   title='Continental Groups',
                  loc='lower right',
                  bbox_to_anchor=(0.98, 0.02),
                  fontsize=13)  # Added legend font size
    
    # Add statistical test results to figure title
    title = f'Trio Concordance by Population (Grouped by Continental Ancestry)\n'
    if len(plot_df['population'].unique()) > 1:
        populations = [group for _, group in plot_df.groupby('population')['concordance']]
        h_stat, p_value = stats.kruskal(*populations)
        title += f'Kruskal-Wallis H-test: p={p_value:.2e}'
    # fig.suptitle(title, y=0.93, fontsize=16)  # If you uncomment, increased font size
    print(title)
    
    return fig, ax_main, stats_df

def calculate_aligned_distances(fasta_file: str) -> Tuple[np.ndarray, List[str]]:
    """
    Calculate pairwise distances between pre-aligned sequences in a FASTA file.
    
    Args:
        fasta_file: Path to aligned FASTA file
        
    Returns:
        Tuple containing:
        - np.ndarray: Square matrix of pairwise distances
        - List[str]: List of sequence IDs in same order as matrix
    """
    # Read sequences and store in dictionary
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract allele name from header (e.g. "M83134|IGHV3-30*01" -> "IGHV3-30*01")
        allele = record.description.split("|")[1]
        sequences[allele] = str(record.seq)
    
    # Get sorted list of allele names
    alleles = sorted(sequences.keys())
    n = len(alleles)
    
    # Initialize distance matrix
    distances = np.zeros((n, n))
    
    # Calculate pairwise distances
    for i in range(n):
        for j in range(i+1, n):
            seq1 = sequences[alleles[i]]
            seq2 = sequences[alleles[j]]
            
            # Count positions where sequences differ
            dist = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
            
            # Store distance in both halves of symmetric matrix
            distances[i,j] = dist
            distances[j,i] = dist
            
    return distances, alleles

def calculate_distribution_distances(distance_matrix: np.ndarray, allele_list: List[str]):
    """Calculate distribution-based distance metrics between genes.
    
    Args:
        distance_matrix: Square numpy array of pairwise allele distances
        allele_list: List of allele IDs in same order as distance matrix
        
    Returns:
        dict: Dictionary with gene pairs as keys and dict of distribution metrics as values
    """
    import numpy as np
    from scipy import stats
    
    # Create mapping of alleles to their genes
    gene_to_alleles = defaultdict(list)
    for allele in allele_list:
        gene = allele.split('*')[0]  # Extract gene name from allele ID
        gene_to_alleles[gene].append(allele)
    
    gene_distances = {}
    genes = sorted(gene_to_alleles.keys())
    
    # Get all gene pairs
    for i, gene1 in enumerate(genes):
        for gene2 in genes[i+1:]:
            distances = []
            
            # Get indices for alleles of both genes
            gene1_indices = [allele_list.index(allele) for allele in gene_to_alleles[gene1]]
            gene2_indices = [allele_list.index(allele) for allele in gene_to_alleles[gene2]]
            
            # Extract all pairwise distances between alleles of these genes
            for idx1 in gene1_indices:
                for idx2 in gene2_indices:
                    distances.append(distance_matrix[idx1, idx2])
            
            if distances:
                # Convert to numpy array for easier calculation
                distances = np.array(distances)
                
                # Calculate distribution metrics
                metrics = {
                    'num_comparisons': len(distances),
                    'min': np.min(distances),
                    'max': np.max(distances),
                    'mean': np.mean(distances),
                    'median': np.median(distances),
                    'std': np.std(distances),
                    'q1': np.percentile(distances, 25),
                    'q3': np.percentile(distances, 75),
                    'iqr': np.percentile(distances, 75) - np.percentile(distances, 25),
                    'skewness': stats.skew(distances),
                    'kurtosis': stats.kurtosis(distances)
                }
                
                # Calculate mode(s)
                mode_result = stats.mode(distances)
                metrics['mode'] = mode_result.mode[0]
                metrics['mode_count'] = mode_result.count[0]
                
                # Store raw distances for potential visualization
                metrics['distances'] = distances
                
                gene_distances[(gene1, gene2)] = metrics
                
    return gene_distances

def calculate_gene_distance_metrics(distance_matrix: np.ndarray, allele_list: List[str]):
    """Calculate intra-gene compactness and inter-gene distinctiveness metrics.
    
    Args:
        distance_matrix: Square numpy array of pairwise allele distances
        allele_list: List of allele IDs in same order as distance matrix
        
    Returns:
        dict: Dictionary with genes as keys and metrics as values:
            - intra_distances: List of distances between alleles within the gene
            - mean_intra_distance: Average distance between alleles within the gene
            - compactness: Inverse of mean intra-distance (normalized)
            - nearest_gene: Name of closest gene by mean distance
            - min_inter_distance: Minimum distance to any other gene's alleles
            - mean_inter_distance: Mean distance to other genes' alleles
            - distinctiveness: Ratio of min_inter_distance to mean_intra_distance
            - num_alleles: Number of alleles for this gene
    """
    # Create mapping of alleles to their genes
    gene_to_alleles = defaultdict(list)
    for allele in allele_list:
        gene = allele.split('*')[0]  # Extract gene name from allele ID
        gene_to_alleles[gene].append(allele)
    
    genes = sorted(gene_to_alleles.keys())
    gene_metrics = {}
    
    # Calculate metrics for each gene
    for gene in genes:
        # Get indices for this gene's alleles
        gene_indices = [allele_list.index(allele) for allele in gene_to_alleles[gene]]
        
        # Calculate intra-gene distances
        intra_distances = []
        for i, idx1 in enumerate(gene_indices):
            for idx2 in gene_indices[i+1:]:
                intra_distances.append(distance_matrix[idx1, idx2])
        
        # Calculate inter-gene distances and find closest gene
        inter_gene_distances = {}  # mean distance to each other gene
        all_inter_distances = []
        
        for other_gene in genes:
            if other_gene != gene:
                other_indices = [allele_list.index(allele) 
                               for allele in gene_to_alleles[other_gene]]
                distances = []
                
                for idx1 in gene_indices:
                    for idx2 in other_indices:
                        dist = distance_matrix[idx1, idx2]
                        distances.append(dist)
                        all_inter_distances.append(dist)
                
                inter_gene_distances[other_gene] = np.mean(distances)
        
        # Calculate metrics
        mean_intra = np.mean(intra_distances) if intra_distances else 0
        nearest_gene = min(inter_gene_distances.items(), key=lambda x: x[1])[0]
        min_inter = min(all_inter_distances) if all_inter_distances else float('inf')
        mean_inter = np.mean(all_inter_distances) if all_inter_distances else float('inf')
        
        # Store metrics
        gene_metrics[gene] = {
            'intra_distances': intra_distances,
            'mean_intra_distance': mean_intra,
            'compactness': 1 / (mean_intra + 1e-10),  # Add small constant to avoid division by zero
            'nearest_gene': nearest_gene,
            'min_inter_distance': min_inter,
            'mean_inter_distance': mean_inter,
            'distinctiveness': min_inter / (mean_intra + 1e-10),
            'num_alleles': len(gene_indices)
        }
    
    # Normalize compactness scores
    max_compactness = max(m['compactness'] for m in gene_metrics.values())
    for metrics in gene_metrics.values():
        metrics['compactness'] = metrics['compactness'] / max_compactness
    
    return gene_metrics

def analyze_concordance_vs_distances(gene_metrics: dict, trio_concordance: dict):
    """Analyze relationship between gene distance metrics and concordance.
    
    Args:
        gene_metrics: Output from calculate_gene_distance_metrics
        trio_concordance: Dictionary with gene as key and list of concordance values
        
    Returns:
        tuple: (DataFrame with analysis results, correlation statistics)
    """
    # Calculate mean concordance per gene
    gene_concordance = {
        gene: np.mean(concordance_values)
        for gene, concordance_values in trio_concordance.items()
    }
    
    # Prepare data for analysis
    analysis_data = []
    for gene in gene_metrics:
        if gene in gene_concordance:
            metrics = gene_metrics[gene]
            analysis_data.append({
                'gene': gene,
                'concordance': gene_concordance[gene],
                'compactness': metrics['compactness'],
                'distinctiveness': metrics['distinctiveness'],
                'num_alleles': metrics['num_alleles'],
                'mean_intra_distance': metrics['mean_intra_distance'],
                'min_inter_distance': metrics['min_inter_distance'],
                'mean_inter_distance': metrics['mean_inter_distance']
            })
    
    analysis_df = pd.DataFrame(analysis_data)
    
    # Calculate correlations with concordance
    metrics_to_correlate = ['compactness', 'distinctiveness', 'num_alleles',
                           'mean_intra_distance', 'min_inter_distance', 
                           'mean_inter_distance']
    
    correlations = []
    for metric in metrics_to_correlate:
        corr, p_value = stats.pearsonr(analysis_df[metric], analysis_df['concordance'])
        correlations.append({
            'metric': metric,
            'correlation': corr,
            'p_value': p_value
        })
    
    # Adjust p-values for multiple testing
    corr_df = pd.DataFrame(correlations)
    corr_df['p_value_adjusted'] = multipletests(corr_df['p_value'], 
                                              method='fdr_bh')[1]
    
    return analysis_df, corr_df

def plot_distance_concordance_analysis(analysis_df: pd.DataFrame, corr_df: pd.DataFrame,
                                     figsize: tuple = (15, 10)):
    """Create visualization of distance metrics vs concordance relationships.
    
    Args:
        analysis_df: DataFrame from analyze_concordance_vs_distances
        corr_df: Correlation statistics from analyze_concordance_vs_distances
        figsize: Figure size tuple
        
    Returns:
        matplotlib Figure object
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 3)
    
    # Plot relationships for each metric
    metrics_to_plot = ['compactness', 'distinctiveness', 'num_alleles',
                      'mean_intra_distance', 'min_inter_distance', 
                      'mean_inter_distance']
    
    for i, metric in enumerate(metrics_to_plot):
        ax = fig.add_subplot(gs[i//3, i%3])
        
        # Create scatter plot
        sns.scatterplot(data=analysis_df, x=metric, y='concordance', ax=ax)
        
        # Add trend line
        sns.regplot(data=analysis_df, x=metric, y='concordance', 
                   scatter=False, color='red', ax=ax)
        
        # Get correlation info
        corr_info = corr_df[corr_df['metric'] == metric].iloc[0]
        
        # Add correlation annotation
        ax.text(0.05, 0.95, 
                f"r = {corr_info['correlation']:.3f}\n"
                f"p = {corr_info['p_value_adjusted']:.2e}",
                transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8))
        
        ax.set_title(f'{metric.replace("_", " ").title()} vs Concordance')
        
    plt.tight_layout()
    return fig