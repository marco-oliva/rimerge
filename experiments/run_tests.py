#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, random, string, errno
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

Description = """
Testing script
"""

#------------------------------------------------------------
# Binaries
dirname = os.path.dirname(os.path.abspath(__file__))
rimerge_exe     =  os.path.join(dirname, "rimerge.x")
check_exe       =  os.path.join(dirname, "check.x")
estw_exe        =  os.path.join(dirname, "estw.x")
rle_exe         =  os.path.join(dirname, "rle.x")
bigbwt_exe      =  os.path.join(dirname, "bigbwt")

#------------------------------------------------------------
# Profiling options
#time = "/usr/bin/time --verbose"
#time = "valgrind --tool=callgrind"
time = ""
#time = "valgrind --tool=massif"
#time = "valgrind  --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose"

#------------------------------------------------------------
# Get n sequences from the input fasta file
def split(input_file, output_dir, seqs, blocks=1, ignore=0):
    print("Splitting input file in {} parts".format(blocks))
    file_paths = []
    filename, file_extension = os.path.splitext(input_file)
    input_file_iterator = SeqIO.parse(input_file, "fasta")

    for i in range(0, ignore):
        next(input_file_iterator)

    seqs_per_block = int(seqs / blocks)
    for b in range (0, blocks):
        output_path = output_dir + "/" + os.path.basename(filename) + "." + str(seqs) + "." + str(b) + ".seqs"
        file_paths.append(output_path)
        output = open(output_path, "w")
        for s in range(0, seqs_per_block - 1):
            record = next(input_file_iterator)
            output.write(str(record.seq))
            output.write('$')
        record = next(input_file_iterator)
        output.write(str(record.seq))
        output.close()

    return file_paths

#------------------------------------------------------------
# Generate random fasta file
def generate_random_fasta(num_of_sequences, length, error_rate = 0.15):
    sequences=[]
    origin_seq = ''.join(random.choices("ACGTN", k = length))
    for i in range(0, num_of_sequences):
        tmp_seq = MutableSeq(origin_seq)
        snps = random.randrange(error_rate * length)
        for r in range(snps):
            snp_pos = random.randrange(length)
            snp_char = random.choice("ACGTN")
            tmp_seq[snp_pos] = snp_char
        seq = Seq(str(tmp_seq))
        record = SeqRecord(seq, "{id}".format(id=i), '', '')
        sequences.append(record)
    cwd = os.getcwd()
    file_name = "{wd}/random_{n}.fasta".format(wd=cwd, n=num_of_sequences)
    SeqIO.write(sequences, file_name, "fasta")
    return file_name

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds):
    try:
        print("Command: {}".format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid)
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        #process.terminate()
        print("Command exceeded timeout:")
        print("\t"+ command)
        return False
    return True

#------------------------------------------------------------
# Run BigBWT
def run_big_bwt(file_path, n_of_sequences, seconds, window, module, check, bonly):
    command = "{profiler} {bigbwt} -t 32 -N {ns} -s -e {file} -p {mod} -w {win}".format(profiler=time, ns=n_of_sequences, file=file_path, bigbwt=bigbwt_exe, mod=module, win=window)
    execute_command(command, seconds)
    command = "{profiler} {estw} -i {i_prefix}".format(profiler=time, estw=estw_exe, i_prefix=file_path)
    execute_command(command, seconds)
    command = "{profiler} {rle} -i {i_prefix}".format(profiler=time, rle=rle_exe, i_prefix=file_path)
    execute_command(command, seconds)
    remove_file(file_path + ".ssa")
    remove_file(file_path + ".esa")
    remove_file(file_path + ".nsa")
    remove_file(file_path + ".bwt")
    remove_file(file_path)
    os.replace(file_path + ".rle", os.path.dirname(file_path) + "/bwt.rle")
    os.replace(file_path + ".rle.meta", os.path.dirname(file_path) + "/bwt.rle.meta")
    os.replace(file_path + ".saes", os.path.dirname(file_path) + "/samples.saes")
    if (check):
        command = "{profiler} {check} -i {i_prefix} -o {i_prefix}".format(profiler=time, i_prefix=file_path, check=check_exe)
        execute_command(command, seconds)

#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path, seconds, check, merge_jobs):
    base_path = ""
    command = "{profiler} {rimerge} -a {left} -b {right} -o {out} -j {mj}".format(profiler=time, rimerge=rimerge_exe, left=left_file_path, right=right_file_path, out=out_file_path, mj=merge_jobs)
    if (check):
        command = command + " -c"
    execute_command(command, seconds)


#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file name', type=str, dest="input")
    parser.add_argument('-w', '--wsize', help='sliding window size (def. 10)', default=10, type=int, dest="window")
    parser.add_argument('-p', '--mod', help='hash modulus (def. 100)', default=100, type=int, dest="module")
    parser.add_argument('-I', '--ignore', help='Ignore first I sequences', type=int, dest="ignore", default=0)
    parser.add_argument('-o', '--output', help='output file prefix', type=str, dest="output")
    parser.add_argument('-j', '--merge-jobs', help='number of merge jobs', type=int, default=4, dest="merge_jobs")
    parser.add_argument('-n', '--seqs', help='number of sequences as percentage', default=1, type=int, dest="seqs")
    parser.add_argument('-b', '--blocks', help='number of blocks', default=1, type=int, dest="blocks")
    parser.add_argument('-B',  help='run Big-BWT ONLY',action='store_true', dest="bonly")
    parser.add_argument('-C',  help='check index structure after each merge',action='store_true', dest="check")
    parser.add_argument('-T', '--timeout', help='command timeouts', default=18000, type=int, dest="seconds")
    parser.add_argument('-R',  help='Run on random sequences',action='store_true', dest="random")
    parser.add_argument('-S',  help='Skip Big-BWT',action='store_true', dest="skip")

    args = parser.parse_args()

    if (args.skip and args.bonly):
        print("Skip and Big-BWT only called toghether, error")

    print("Version: ")
    execute_command("{} --version".format(rimerge_exe), 10)

    # Generate the input file
    if (args.random and not args.skip):
        args.input = generate_random_fasta(args.seqs, 600000)
    elif (args.random and args.skip):
        cwd = os.getcwd()
        file_name = "{wd}/random_{n}.fasta".format(wd=cwd, n=args.seqs)
        args.input = file_name

    # Split the input file
    if (not args.skip):
        seqs_per_block = int(args.seqs / args.blocks)
        print("Splitting the input file, {} sequences per block".format(seqs_per_block))
        file_paths = split(args.input, args.output, args.seqs, args.blocks, args.ignore)
    else:
        filename, file_extension = os.path.splitext(args.input)
        for b in range (0, args.blocks):
            output_path = filename + "." + str(args.seqs) + "." + str(b) + ".seqs"
            file_paths.append(output_path)

    # Run Big-BWT on each block
    if (not args.skip):
        print("Running Big-BWT")
        for i in range(0, len(file_paths)):
            run_big_bwt(file_paths[i], seqs_per_block + 1, args.seconds, args.window, args.module, args.check, args.bonly)

    # Iterative Merge
    if (not args.bonly):
        print("Merging")
        print("Output file prefix: {}".format(args.output))
        run_merge(file_paths[0], file_paths[1], args.output, args.seconds, args.check, args.merge_jobs)
        for i in range(2, len(file_paths)):
            run_merge(args.output, file_paths[i], args.output, args.seconds, args.check, args.merge_jobs)


if __name__ == '__main__':
    main()