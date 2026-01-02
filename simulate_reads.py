import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def sample_short_reads(fasta_path, output_path, read_length, num_reads):
    long_reads = list(SeqIO.parse(fasta_path, "fasta"))
    total_sampled = 0
    simulated_reads = []

    while total_sampled < num_reads:
        record = random.choice(long_reads)
        seq = str(record.seq)

        if len(seq) < read_length:
            continue

        start = random.randint(0, len(seq) - read_length)
        short_seq = seq[start:start + read_length]

        # Create synthetic quality scores (Phred 40 = ASCII 'I')
        qualities = [40] * read_length

        read = SeqRecord(
            Seq(short_seq),
            id=f"sim_read_{total_sampled}",
            description="",
            letter_annotations={"phred_quality": qualities}
        )

        simulated_reads.append(read)
        total_sampled += 1

    # Write all reads to FASTQ file
    with open(output_path, "w") as out_f:
        SeqIO.write(simulated_reads, out_f, "fastq")

    print(f"Finished writing {total_sampled} FASTQ reads of length {read_length} with max quality.")

def main():
    parser = argparse.ArgumentParser(description="Simulate FASTQ short reads from long reads (FASTA).")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file with long reads")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ file for simulated short reads")
    parser.add_argument("-l", "--length", type=int, default=150, help="Length of each short read (default: 150)")
    parser.add_argument("-n", "--num_reads", type=int, default=100000, help="Number of short reads to simulate (default: 100000)")
    args = parser.parse_args()

    sample_short_reads(args.input, args.output, args.length, args.num_reads)

if __name__ == "__main__":
    main()