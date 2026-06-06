#!/usr/bin/env python3
import sys
import os

# Standard Genetic Code mapping
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

COMPLEMENT_TABLE = str.maketrans('ATCGatcgNn', 'TAGCtagcNn')

def reverse_complement(seq):
    return seq.translate(COMPLEMENT_TABLE)[::-1]

def translate(seq):
    seq = seq.upper()
    protein = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i+3]
        protein.append(GENETIC_CODE.get(codon, 'X'))
    return "".join(protein)

def read_fasta(fasta_path):
    sequences = {}
    current_id = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = "".join(current_seq)
    return sequences

def parse_gff3(gff3_path):
    cds_features = []
    with open(gff3_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes_str = parts
            if feature_type != 'CDS':
                continue
                
            start = int(start_str)
            end = int(end_str)
            
            # Parse attributes
            attrs = {}
            for item in attributes_str.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k.strip()] = v.strip()
            
            cds_features.append({
                'seqid': seqid,
                'start': start,
                'end': end,
                'strand': strand,
                'phase': phase,
                'attrs': attrs
            })
    return cds_features

def main():
    fasta_path = 'genome.fa'
    gff3_path = 'annotations.gff3'
    
    if len(sys.argv) > 1:
        fasta_path = sys.argv[1]
    if len(sys.argv) > 2:
        gff3_path = sys.argv[2]
        
    if not os.path.exists(fasta_path):
        print(f"Error: FASTA file not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(gff3_path):
        print(f"Error: GFF3 file not found: {gff3_path}", file=sys.stderr)
        sys.exit(1)
        
    # Read sequences
    sequences = read_fasta(fasta_path)
    
    # Parse GFF3 CDS features
    cds_features = parse_gff3(gff3_path)
    
    # Group CDS segments.
    # Group by Parent first, then by ID, then fallback to index.
    groups = {}
    for idx, cds in enumerate(cds_features):
        attrs = cds['attrs']
        group_key = attrs.get('Parent') or attrs.get('ID') or f"cds_unnamed_{idx}"
        if group_key not in groups:
            groups[group_key] = []
        groups[group_key].append(cds)
        
    # Process and print each CDS group
    for group_key, segments in groups.items():
        first_seg = segments[0]
        seqid = first_seg['seqid']
        strand = first_seg['strand']
        attrs = first_seg['attrs']
        
        # Get the gene name using fallback options
        gene_name = attrs.get('Name') or attrs.get('gene') or attrs.get('gene_name') or attrs.get('ID') or group_key
        
        if seqid not in sequences:
            print(f"Warning: Sequence ID '{seqid}' not found in FASTA file. Skipping.", file=sys.stderr)
            continue
            
        seq = sequences[seqid]
        
        # Sort segments by coordinate (ascending)
        segments_sorted = sorted(segments, key=lambda x: x['start'])
        
        # Concatenate exons
        nt_parts = []
        for seg in segments_sorted:
            start = seg['start']
            end = seg['end']
            if start < 1 or end > len(seq):
                print(f"Warning: Coordinates {start}-{end} out of bounds for sequence {seqid} (length {len(seq)}). Skipping.", file=sys.stderr)
                continue
            part = seq[start - 1 : end]
            nt_parts.append(part)
            
        nt_seq = "".join(nt_parts)
        
        # Handle negative strand
        if strand == '-':
            nt_seq = reverse_complement(nt_seq)
            
        # Translate to protein
        protein_seq = translate(nt_seq)
        
        # Print tab-separated: gene_name <TAB> nt_sequence <TAB> protein_sequence
        print(f"{gene_name}\t{nt_seq}\t{protein_seq}")

if __name__ == '__main__':
    main()
