#!/usr/bin/env python3
"""Script to calculate allele frequencies per chromosome from a VCF file.

Reads a VCF file, parses sample genotypes (GT) to calculate allele frequency
for each variant, and outputs both individual variant frequencies and average
allele frequencies grouped by chromosome.
"""

import sys
import gzip
from pathlib import Path
from collections import defaultdict


def parse_vcf(vcf_path):
    """Parses a VCF file (gzip supported) and computes allele frequencies.
    
    Returns:
        tuple: (variant_records, chrom_summaries)
            - variant_records: list of dicts with variant-level information
            - chrom_summaries: dict mapping chromosome -> statistics
    """
    variant_records = []
    chrom_allele_freqs = defaultdict(list)
    
    # Select appropriate file opener based on file extension
    open_func = gzip.open if str(vcf_path).endswith(".gz") else open
    mode = "rt" if str(vcf_path).endswith(".gz") else "r"
    
    with open_func(vcf_path, mode, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            fields = line.strip().split("\t")
            if len(fields) < 10:
                # VCF genotype columns start at 9 (0-indexed)
                continue
                
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # Identify where the GT (Genotype) subfield is in the FORMAT field
            format_col = fields[8].split(":")
            try:
                gt_idx = format_col.index("GT")
            except ValueError:
                # GT not found in format column
                continue
                
            total_alleles = 0
            alt_alleles = 0
            
            # Process sample genotype columns (indices 9 to end)
            for sample_col in fields[9:]:
                sample_format = sample_col.split(":")
                if len(sample_format) <= gt_idx:
                    continue
                
                gt = sample_format[gt_idx]
                # Normalize genotype separator (could be / or |)
                alleles = gt.replace("|", "/").split("/")
                
                for allele in alleles:
                    if allele == "." or allele == "":
                        # Skip missing alleles
                        continue
                    try:
                        allele_int = int(allele)
                        total_alleles += 1
                        # 0 is REF, >0 is ALT
                        if allele_int > 0:
                            alt_alleles += 1
                    except ValueError:
                        # Non-integer allele (malformed GT)
                        continue
            
            if total_alleles > 0:
                af = alt_alleles / total_alleles
                variant_records.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "af": af
                })
                chrom_allele_freqs[chrom].append(af)
                
    # Aggregate stats per chromosome
    chrom_summaries = {}
    for chrom, freqs in chrom_allele_freqs.items():
        chrom_summaries[chrom] = {
            "count": len(freqs),
            "avg_af": sum(freqs) / len(freqs) if freqs else 0.0
        }
        
    return variant_records, chrom_summaries


def main():
    if len(sys.argv) > 1:
        vcf_path = Path(sys.argv[1])
    else:
        vcf_path = Path(__file__).parent / "sample.vcf"
        
    if not vcf_path.exists():
        print(f"Error: File '{vcf_path}' not found.", file=sys.stderr)
        print("Usage: python vcf_allele_freq.py [path_to_vcf_file]", file=sys.stderr)
        sys.exit(1)
        
    print(f"Parsing VCF file: {vcf_path}\n")
    try:
        variants, summaries = parse_vcf(vcf_path)
    except Exception as e:
        print(f"Error parsing VCF: {e}", file=sys.stderr)
        sys.exit(1)
        
    # Print variant-level allele frequencies
    print("--- Variant Allele Frequencies ---")
    print(f"{'CHROM':<10}\t{'POS':<10}\t{'REF':<5}\t{'ALT':<5}\t{'Allele Frequency':<16}")
    print("-" * 60)
    for var in variants:
        print(f"{var['chrom']:<10}\t{var['pos']:<10}\t{var['ref']:<5}\t{var['alt']:<5}\t{var['af']:.4f}")
    print()
    
    # Print chromosome-level summaries
    print("--- Chromosome Summary ---")
    print(f"{'CHROM':<10}\t{'Variants Count':<15}\t{'Average Allele Frequency':<25}")
    print("-" * 60)
    for chrom, summary in summaries.items():
        print(f"{chrom:<10}\t{summary['count']:<15}\t{summary['avg_af']:.4f}")


if __name__ == "__main__":
    main()
