import argparse
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def main():
    parser = argparse.ArgumentParser(description="Filter sequences based on length threshold and GC content")
    parser.add_argument("--input-sequences", required=True, help="Path to fasta file containing sequences")
    parser.add_argument("--length-threshold", type=int, default=50, help="Minimum length of sequences to retain")
    parser.add_argument("--gc-threshold", type=float, default=0.4, help="Minimum GC content of sequences to retain")
    parser.add_argument("--output", required=True, help="Path to output file for filtered sequences")

    args = parser.parse_args()

    retained = []
    filtered = 0
    
    with open(args.input_sequences, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if len(record.seq) >= args.length_threshold and gc_fraction(record.seq) >= args.gc_threshold:
                retained.append(record)
            else:
                filtered += 1
    
    # Write filtered sequences
    with open(args.output, 'w') as out_handle:
        SeqIO.write(retained, out_handle, 'fasta')

    print(f"Total sequences processed: {len(retained) + filtered}")
    print(f"Sequences retained: {len(retained)}")
    print(f"Sequences filtered out: {filtered}")

if __name__ == "__main__":
    main()