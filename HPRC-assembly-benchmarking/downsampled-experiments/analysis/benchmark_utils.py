"""Utility functions for genotyping accuracy benchmarking."""

from collections import OrderedDict, Counter
from pathlib import Path
import pandas as pd
import logging
from Bio import SeqIO
from typing import Tuple
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)


def load_genotype_calls(sample_id: str, gene_type: str, call_base_dir: str="../immunotyper-output/") -> list:
    """
    Load ImmunoTyper genotype calls for a given sample and gene type.
    
    Args:
        sample_id: Sample identifier
        gene_type: Type of gene (e.g., 'iglv', 'trav')
        
    Returns:
        List of called alleles
        
    Raises:
        FileNotFoundError: If no matching files are found
        RuntimeError: If multiple matching files are found
    """
    # Create glob pattern for the calls file
    pattern = os.path.join(call_base_dir, f"{gene_type.lower()}/{sample_id}.final*{gene_type.upper()}_functional_allele_calls.txt")
    matching_files = list(Path().glob(pattern))
    
    if not matching_files:
        log.error(f"No genotype calls file found for {sample_id}. Pattern: {pattern}")
        raise FileNotFoundError(f"No genotype calls file found matching pattern: {pattern}")
    
    if len(matching_files) > 1:
        log.error(f"Multiple genotype calls files found for {sample_id}: {matching_files}")
        raise RuntimeError(f"Expected exactly one genotype calls file, found {len(matching_files)}")
    
    calls_path = matching_files[0]
    with open(calls_path, 'r') as f:
        calls = [line.strip() for line in f if line.strip()]
    
    return calls



def load_assembly_annotations(gene_type, path="../digger-functional-annotations-top-contig/"):
    """
    Load annotations for a specific gene_type across all samples and merge
    maternal/paternal haplotypes.
    
    Args:
        gene_type (str): The gene type to load (e.g., 'IGH', 'TRB')
        path (str): Base path where sample directories are located
        
    Returns:
        dict: Dictionary with merged annotations keyed by sample ID
    """
    import os
    import re
    
    # Dictionary to store results
    results = {}
    
    # Dictionary to track haplotype data before merging
    haplotype_data = {}
    
    # Regular expression to extract sample name from directory
    pattern = r"^(.*?)\.(maternal|paternal)\.f1_assembly_v2_genbank$"
    print(f"Loading annotations for {gene_type} from {path}")
    # Walk through directories
    for sample_dir in os.listdir(path):
        # Skip non-directories or directories that don't match pattern
        sample_dir_path = os.path.join(path, sample_dir)

        if not os.path.isdir(sample_dir_path):
            log.error(f"Skipping non-directory: {sample_dir_path}")
            continue
            
        # Match the directory name against our pattern
        match = re.match(pattern, sample_dir)

        if not match:
            log.error(f"Skipping directory that doesn't match pattern: {sample_dir}")
            continue
            
        # Extract sample ID and haplotype
        sample_id = match.group(1)
        haplotype = match.group(2)
        
        # Path to the gene type annotations file
        gene_file = os.path.join(sample_dir_path, f"{gene_type[:3]}_functional_alleles.txt")
        
        # Skip if file doesn't exist (gene type might be missing for this sample)
        if not os.path.exists(gene_file):
            log.error(f"Skipping missing gene file: {gene_file}")
            continue
            
        # Load alleles from file
        alleles = []
        with open(gene_file, 'r') as f:
            for line in f:
                # Skip header lines and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                alleles.append(line.strip())
        
        # Store in haplotype data dictionary
        if sample_id not in haplotype_data:
            haplotype_data[sample_id] = {}
        haplotype_data[sample_id][haplotype] = alleles
    
    # Merge maternal and paternal haplotypes
    for sample_id, haplotypes in haplotype_data.items():
        all_alleles = set()
        
        # Add maternal alleles if available
        if 'maternal' in haplotypes:
            all_alleles.update(haplotypes['maternal'])
        
        # Add paternal alleles if available
        if 'paternal' in haplotypes:
            all_alleles.update(haplotypes['paternal'])
        
        # Store merged results
        results[sample_id] = sorted(list(all_alleles))
    
    return results


def calculate_accuracy_metrics(results: list, truth: list, cnv: bool = False) -> tuple:
    """
    Calculate true positives, false positives, and false negatives.
    
    Args:
        results: List of predicted alleles (may contain duplicates)
        truth: List of true alleles (may contain duplicates)
        cnv: If True, consider copy number variations (count duplicates)
        
    Returns:
        Tuple of (true_positives, false_positives, false_negatives)
    """
    if not cnv:
        # Original set-based logic for presence/absence
        results_set = set(results)
        truth_set = set(truth)
        
        tp = list(results_set & truth_set)  # Intersection
        fp = list(results_set - truth_set)  # In results but not in truth
        fn = list(truth_set - results_set)  # In truth but not in results
    else:
        # Use Counter objects to handle multiple copies
        results_counter = Counter(results)
        truth_counter = Counter(truth)
        
        tp = []
        fp = []
        fn = []
        
        # Process all unique alleles from both results and truth
        all_alleles = set(results) | set(truth)
        for allele in all_alleles:
            result_count = results_counter[allele]
            truth_count = truth_counter[allele]
            
            if truth_count > 0:
                # Add true positives (minimum of predicted and true counts)
                tp.extend([allele] * min(result_count, truth_count))
                
                # Add false negatives (missing copies)
                if truth_count > result_count:
                    fn.extend([allele] * (truth_count - result_count))
            
            # Add false positives (extra copies)
            if result_count > truth_count:
                fp.extend([allele] * (result_count - truth_count))
    
    return tp, fp, fn


def calculate_precision_recall(tp: list, fp: list, fn: list) -> tuple:
    """
    Calculate precision and recall metrics.
    
    Args:
        tp: List of true positives
        fp: List of false positives
        fn: List of false negatives
        
    Returns:
        Tuple of (precision, recall)
    """
    tp_count = float(len(tp))
    fp_count = float(len(fp))
    fn_count = float(len(fn))
    
    precision = round(tp_count / (tp_count + fp_count), 3) if (tp_count + fp_count) > 0 else 0
    recall = round(tp_count / (tp_count + fn_count), 3) if (tp_count + fn_count) > 0 else 0
    
    return precision, recall
