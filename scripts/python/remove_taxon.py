#!/usr/bin/env python3


import os
import argparse

def remove_sequence_from_fasta(input_file, output_file, seq_id):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        write_seq = True
        for line in infile:
            if line.startswith(">"):
                # Check if this is the sequence to remove
                header = line.strip()
                if header[1:].split()[0] == seq_id:
                    write_seq = False
                else:
                    write_seq = True

            if write_seq:
                outfile.write(line)


def process_directory(directory, seq_id, output_directory=None):
    if output_directory is None:
        output_directory = directory

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(directory):
        if filename.endswith((".fasta", ".fa", ".aln")):
            input_path = os.path.join(directory, filename)
            output_path = os.path.join(output_directory, filename)

            print(f"Processing {filename}...")
            remove_sequence_from_fasta(input_path, output_path, seq_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Remove a FASTA entry from all alignment files in a directory"
    )
    parser.add_argument("directory", help="Input directory with alignment files")
    parser.add_argument("seq_id", help="Sequence ID to remove (header without '>')")
    parser.add_argument(
        "-o", "--output", help="Output directory (default: overwrite input directory)"
    )

    args = parser.parse_args()

    process_directory(args.directory, args.seq_id, args.output)