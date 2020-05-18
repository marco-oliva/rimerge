#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path
from Bio import SeqIO

Description = """
Tesging script
"""


#------------------------------------------------------------
# Profiling options
time = "/usr/bin/time --verbose"

#------------------------------------------------------------
# Get first n sequences from the input fasta file
def split(input_file, seqs, blocks=1):
    file_paths = []
    filename, file_extension = os.path.splitext(input_file)
    input_file_iterator = SeqIO.parse(input_file, "fasta")

    seqs_per_block = int(seqs / blocks)
    for b in range (0, blocks):
        output_path = filename + "." + str(seqs) + "." + str(b) + ".seqs"
        file_paths.append(output_path)
        output = open(output_path, "a")
        for s in range(0, seqs_per_block - 1):
            record = next(input_file_iterator)
            output.write(str(record.seq))
            output.write('$')
        record = next(input_file_iterator)
        output.write(str(record.seq))
        output.close()

    return file_paths

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

# execute command: return True is everything OK, False otherwise
def execute_command(command):
    try:
        subprocess.check_call(command.split())
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    return True

#------------------------------------------------------------
# Run BigBWT
def run_big_bwt(file_path, n_of_sequences):
    #base_path="/s/apu/a/homes/marco.oliva/projects/Big-BWT/"
    base_path=""
    command = "{profiler} {bp}bigbwt -t 8 -N {ns} -s -e {file}".format(profiler=time, ns=n_of_sequences, file=file_path, bp=base_path)
    execute_command(command)


#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path):
    #base_path="/s/apu/a/homes/marco.oliva/projects/r-merge/build/bin/"
    base_path = ""
    command = "{profiler} {bp}main -a {left} -b {right} -o {out}".format(profiler=time, bp=base_path, left=left_file_path, right=right_file_path, out = out_file_path)
    execute_command(command)

#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file name', type=str, dest="input")
    parser.add_argument('-o', '--output', help='output file prefix', type=str, dest="output")
    parser.add_argument('-n', '--seqs', help='number of sequences', default=10, type=int, dest="seqs")
    parser.add_argument('-b', '--blocks', help='number of blocks', default=1, type=int, dest="blocks")
    parser.add_argument('-B',  help='Run Big-BWT ONLY',action='store_true', dest="bonly")

    args = parser.parse_args()

    # Split the input file
    seqs_per_block = int(args.seqs / args.blocks)
    print("Splitting the input file, {} sequences per block".format(seqs_per_block))
    file_paths = split(args.input, args.seqs, args.blocks)

    # Run Big-BWT on each block
    print("Running Big-BWT")
    for i in range(0, len(file_paths)):
        run_big_bwt(file_paths[i], seqs_per_block)

    # Clean up seqs files
    for path in file_paths:
        remove_file(path)

    # Iterative Merge
    if (not args.bonly):
        print("Merging")
        print("Output file prefix: {}".format(args.output))
        run_merge(file_paths[0], file_paths[1], args.output)
        for i in range(2, len(file_paths)):
            run_merge(args.output, file_paths[i], args.output)


if __name__ == '__main__':
    main()