"""Utility functions for genotyping accuracy benchmarking."""

from collections import OrderedDict, Counter
from pathlib import Path
import pandas as pd
import logging
from Bio import SeqIO
from typing import Tuple

# Configure logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

def load_assembly_annotations(sample_id: str, gene_type: str) -> Tuple[pd.DataFrame, list]:
    """
    Load ground truth annotations for a given sample and gene type from TSV files
    into a single DataFrame, and add functionality status from FASTA reference.
    
    Args:
        sample_id: Sample identifier
        gene_type: Type of gene (e.g., 'IGLV', 'TRAV')
        
    Returns:
        Tuple containing:
        - DataFrame with columns: [START, END, ALLELE, EDIT DISTANCE, IS_REVERSE, 
                                 CLOSEST_ALLELE, source_file, is_functional]
        - List of unique functional alleles
    """
    # Use glob to find all TSV files for this sample
    annotation_dir = Path(f"data/assembly-annotations/{gene_type}-reference-annotations/{sample_id}")
    tsv_files = list(annotation_dir.glob("*.tsv"))
    
    if not tsv_files:
        log.warning(f"No TSV annotation files found for {sample_id} in {annotation_dir}")
        return pd.DataFrame(), []
    
    # Read and combine all TSV files
    dfs = []
    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep='\t')
        df['source_file'] = tsv_file.name
        dfs.append(df)
    
    # Concatenate all DataFrames
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Load FASTA file and create functionality lookup dictionary
    fasta_path = Path(f"data/allele-databases/{gene_type}-IMGT-allele-db.fa")
    functionality_dict = {}
    
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            allele_info = record.id.split('|')
            allele_name = f"{allele_info[1]}"  # e.g., "IGLV(I)-11-1*01"
            is_functional = allele_info[3] in ["F", "(F)", "[F]"]
            functionality_dict[allele_name] = is_functional
            
        # Add functionality status to DataFrame
        combined_df['is_functional'] = combined_df['ALLELE'].map(functionality_dict)
        
        # Report alleles with missing functionality information
        missing_alleles = combined_df[combined_df['is_functional'].isna()]['ALLELE'].unique()
        if len(missing_alleles) > 0:
            log.warning(f"Alleles not found in FASTA reference for {sample_id}: {missing_alleles}")
        
        # Get list of unique functional alleles, handling NaN values
        functional_alleles = combined_df[
            combined_df['is_functional'].fillna(False)
        ]['ALLELE'].tolist()
        
    except Exception as e:
        log.error(f"Error processing FASTA file {fasta_path}: {e}")
        combined_df['is_functional'] = None
        functional_alleles = []
    
    return combined_df, functional_alleles

def load_genotype_calls(sample_id: str, gene_type: str) -> list:
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
    pattern = f"data/genotypes/{gene_type.upper()}/{sample_id}*/{sample_id}*{gene_type.upper()}_functional_allele_calls.txt"
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
