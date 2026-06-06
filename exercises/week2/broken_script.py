"""Small intentionally broken FASTA summarizer for Week 2.

Use an agentic IDE to fix this script one issue at a time. The intended
behavior is: read a FASTA file, print one row per sequence, and report sequence
length plus GC percentage.
"""

from pathlib import Path


def read_fasta(path):
    content = Path(path).read_text()
    records = {}
    
    for entry in content.split(">")[1:]:
        lines = entry.splitlines()
        if lines:
            name = lines[0].strip()
            sequence = "".join(line.strip() for line in lines[1:])
            records[name] = sequence
            
    return records


def gc_percent(sequence):
    if not sequence:
        return 0.0
    sequence = sequence.upper()
    gc = sequence.count("G") + sequence.count("C")
    return gc / len(sequence)


def main():
    fasta_path = Path(__file__).parent / "example.fa"
    records = read_fasta(fasta_path)

    for name, sequence in records.items():
        print(name, len(sequence), gc_percent(sequence))


if __name__ == "__main__":
    main()
