#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, random, string, shutil
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

Description = """
Rimerge
"""

#------------------------------------------------------------
# Absolute path of auxiliary executables
dirname = os.path.dirname(os.path.abspath(__file__))
rimerge_exe = os.path.join(dirname, "rimerge.x")
estw_exe    = os.path.join(dirname, "estw.x")
check_exe   = os.path.join(dirname, "check.x")
rle_exe     = os.path.join(dirname, "rle.x")
cfa_exe     = os.path.join(dirname, "cfa.x")
bigbwt_exe  = os.path.join(dirname, "bigbwt")


#------------------------------------------------------------
# Profiling options
#time = "/usr/bin/time --verbose"
#time = "valgrind --tool=callgrind"
time = ""
#time = "valgrind --tool=massif"
#time = "valgrind  --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose"

#------------------------------------------------------------
# Split fasta file
def split(input_file, output_dir, seqs, blocks=1, ignore=0):
    print("Splitting input file in {} parts".format(blocks))
    file_paths = []
    filename, file_extension = os.path.splitext(input_file)
    input_file_iterator = SeqIO.parse(input_file, "fasta")

    for i in range(0, ignore):
        next(input_file_iterator)

    seqs_per_block = 0
    remainder = 0
    if (seqs % blocks != 0):
        seqs_per_block = int(seqs / (blocks - 1))
        remainder = int(seqs % blocks)
        blocks = blocks - 1
    else :
        seqs_per_block = int(seqs / blocks)

    for b in range (0, blocks):
        output_path = output_dir + "/" + os.path.basename(filename) + "." + str(seqs) + "." + str(b) + ".fa"
        file_paths.append(output_path)
        output = open(output_path, "w")
        for s in range(0, seqs_per_block):
            record = next(input_file_iterator)
            output.write(">" + str(record.id) + "\n")
            output.write(str(record.seq) + "\n")
        output.close()

    # reamainder
    if (remainder != 0):
        output_path = output_dir + "/" + os.path.basename(filename) + "." + str(seqs) + "." + str(blocks) + ".fa"
        file_paths.append(output_path)
        output = open(output_path, "w")
        for s in range(0, remainder):
            record = next(input_file_iterator)
            output.write(">" + str(record.id) + "\n")
            output.write(str(record.seq) + "\n")
        output.close()

    return file_paths

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

def remove_dir(dir_path):
    shutil.rmtree(dir_path)

def move_dir_content(dir_source, dir_dest):
    files = os.listdir(dir_source)
    for f in files:
        shutil.move(dir_source + "/" + f, dir_dest)

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds):
    try:
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid)
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Command exceeded timeout:")
        print("\t"+ command)
        return False
    return True

#------------------------------------------------------------
# Run BigBWT
def run_big_bwt(file_path, n_of_sequences, seconds, window, module):
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

#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path, seconds, merge_jobs):
    command = "{profiler} {rimerge} -a {left} -b {right} -o {out} -j {mj}".format(profiler=time, rimerge=rimerge_exe, left=left_file_path, right=right_file_path, out=out_file_path, mj=merge_jobs)
    execute_command(command, seconds)


#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--A-directory-path', help='input A file name', type=str, dest="input_a", required=True)
    parser.add_argument('-b', '--B-directory-path', help='input B file name', type=str, dest="input_b")
    parser.add_argument('-w', '--wsize', help='sliding window size (def. 10)', default=10, type=int, dest="window")
    parser.add_argument('-p', '--mod', help='hash modulus (def. 100)', default=100, type=int, dest="module")
    parser.add_argument('-o', '--output', help='output index prefix', type=str, dest="output", required=True)
    parser.add_argument('-j', '--merge-jobs', help='number of merge jobs', type=int, default=4, dest="merge_jobs")
    parser.add_argument('-n', '--blocks', help='number of blocks', default=1, type=int, dest="blocks")
    parser.add_argument('-T', '--timeout', help='command timeouts', default=1800000, type=int, dest="seconds")
    parser.add_argument('-C', '--check', help='check output index structure', action='store_true',default=False, dest="check")
    args = parser.parse_args()

    # Check arguements
    if (not os.path.isdir(args.input_a) and not args.input_a.endswith(('.fa', '.fasta', '.FA', '.FASTA'))):
        sys.exit('Input A must be a directory containing the index or a fasta file')
    if ((args.input_b is not None) and (not os.path.isdir(args.input_b) and not args.input_b.endswith(('.fa', '.fasta', '.FA', '.FASTA')))):
        sys.exit('Input B must be a directory containing the index or a fasta file')

    print("Version: ")
    execute_command("{rimerge} --version".format(rimerge=rimerge_exe), 10)

    # 1 fasta file as input
    if ( (not os.path.isdir(args.input_a)) and (args.input_b is None) ):
        if (not os.path.isdir(args.output)):
            os.mkdir(args.output)
        sequences = 0
        with open(args.input_a, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences = sequences + 1
        tmp_dir_name = os.path.splitext(args.input_a)[0] + "_tmp"
        if (not os.path.isdir(tmp_dir_name)):
            os.mkdir(tmp_dir_name)
        file_paths = split(args.input_a, tmp_dir_name, sequences, args.blocks)

        file_paths_moved = []
        for i in range(0, len(file_paths)):
            split_dir_name = tmp_dir_name + "/split_{}".format(i)
            if (not os.path.isdir(split_dir_name)):
                os.mkdir(split_dir_name)
            os.replace(file_paths[i], split_dir_name + "/input.fa")
            file_paths_moved.append(split_dir_name + "/input.fa")

        # run cfa
        for i in range(len(file_paths_moved)):
            out_file = os.path.splitext(file_paths_moved[i])[0] + ".seqs"
            command  = "{profiler} {cfa} -i {f} -o {of}".format(profiler=time, cfa=cfa_exe, f=file_paths_moved[i], of=out_file)
            execute_command(command, args.seconds)
            file_paths_moved[i] = out_file

        for file in file_paths_moved:
            run_big_bwt(file, int((sequences / args.blocks) + 1), args.seconds, args.window, args.module)

        if (len(file_paths_moved) > 1):
            run_merge(os.path.dirname(file_paths_moved[0]), os.path.dirname(file_paths_moved[1]), args.output, args.seconds, args.merge_jobs)
            for i in range(2, len(file_paths_moved)):
                run_merge(args.output, os.path.dirname(file_paths_moved[i]), args.output, args.seconds, args.merge_jobs)
        else:
            move_dir_content(os.path.dirname(file_paths_moved[0]), args.output)

        remove_dir(tmp_dir_name)
    # 2 Indexes as input
    elif (os.path.isdir(args.input_a) and os.path.isdir(args.input_b)):
        if ( (not os.path.isfile(args.input_a + "/bwt.rle")) or (not os.path.isfile(args.input_a + "/samples.saes")) ):
            sys.exit('Input A wrong format')
        if ( (not os.path.isfile(args.input_b + "/bwt.rle")) or (not os.path.isfile(args.input_b + "/samples.saes")) ):
            sys.exit('Input B wrong format')
        if (not os.path.isdir(args.output)):
            os.mkdir(args.output)
        run_merge(args.input_a, args.input_b, args.output, args.seconds, args.merge_jobs)
    # Second index as fasta
    elif (os.path.isdir(args.input_a) and (not os.path.isdir(args.input_b))):
        sequences = 0
        with open(args.input_b, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences = sequences + 1
        tmp_dir_name = os.path.splitext(args.input_b)[0] + "_tmp"
        if (not os.path.isdir(tmp_dir_name)):
            os.mkdir(tmp_dir_name)
        file_paths = split(args.input_b, tmp_dir_name, sequences, args.blocks)

        file_paths_moved = []
        for i in range(0, len(file_paths)):
            split_dir_name = tmp_dir_name + "/split_{}".format(i)
            if (not os.path.isdir(split_dir_name)):
                os.mkdir(split_dir_name)
            os.replace(file_paths[i], split_dir_name + "/input.fa")
            file_paths_moved.append(split_dir_name + "/input.fa")

        # run cfa
        for i in range(len(file_paths_moved)):
            out_file = os.path.splitext(file_paths_moved[i])[0] + ".seqs"
            command  = "{profiler} {cfa} -i {f} -o {of}".format(profiler=time, cfa=cfa_exe, f=file_paths_moved[i], of=out_file)
            execute_command(command, args.seconds)
            file_paths_moved[i] = out_file

        for file in file_paths_moved:
            run_big_bwt(file, int((sequences / args.blocks) + 1), args.seconds, args.window, args.module)

        if (not os.path.isdir(args.output)):
            os.mkdir(args.output)
        run_merge(args.input_a, os.path.dirname(file_paths_moved[0]), args.output, args.seconds, args.merge_jobs)
        for i in range(1, len(file_paths_moved)):
            run_merge(args.output, os.path.dirname(file_paths_moved[i]), args.output, args.seconds, args.merge_jobs)

        remove_dir(tmp_dir_name)

    if (args.check):
        command = "{} -i {} -o {}".format(check_exe, args.output, args.output)
        execute_command(command, args.seconds)

if __name__ == '__main__':
    main()