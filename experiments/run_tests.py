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
check_sa_exe    = os.path.join(dirname, "check_sa.x")
estw_exe        =  os.path.join(dirname, "estw.x")
rle_exe         =  os.path.join(dirname, "rle.x")
bigbwt_exe      =  os.path.join(dirname, "bigbwt")
pfp_exe         =  os.path.join(dirname, "pfp++")
aupair_exe         =  os.path.join(dirname, "AuPair")

parsebwt_exe    = os.path.join(dirname, "bwtparse")
parsebwt_exe64    = os.path.join(dirname, "bwtparse64")

pfbwtNT_exe       = os.path.join(dirname, "pfbwtNT.x")
pfbwtNT_exe64       = os.path.join(dirname, "pfbwtNT64.x")

pfbwtSANT_exe       = os.path.join(dirname, "pfbwtSANT.x")
pfbwtSANT_exe64       = os.path.join(dirname, "pfbwtSANT64.x")

#------------------------------------------------------------
# Profiling options
#time = "/usr/bin/time --verbose"
#time = "valgrind --tool=callgrind"
profiler = ""
#profiler = "valgrind --tool=massif"
#profiler = "valgrind  --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose"

#------------------------------------------------------------
# Get n sequences from the input fasta file
def split(input_file, output_dir, seqs=1, blocks=1, ignore=0):
    print("Splitting input file in {} parts".format(blocks))
    file_paths = []
    filename, file_extension = os.path.splitext(input_file)
    input_file_iterator = SeqIO.parse(input_file, "fasta")

    for i in range(0, ignore):
        next(input_file_iterator)

    seqs_per_block = int(seqs / blocks)
    for b in range (0, blocks):
        output_path = output_dir + "/" + os.path.basename(filename) + "." + str(seqs) + "." + str(b) + ".fa"
        file_paths.append(output_path)
        output = open(output_path, "w")
        for s in range(0, seqs_per_block - 1):
            record = next(input_file_iterator)
            SeqIO.write(record, output, "fasta")
        record = next(input_file_iterator)
        SeqIO.write(record, output, "fasta")
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

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python ≥ 2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds = 10000000):
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
def build_rindex(file_path, num_of_sequences, window_length, modulo, check):
    input_prefix = os.path.basename(os.path.splitext(file_path)[0])
    command = "{profiler} {pfp} -f {file} -p {mod} -w {win} -o {out}".format(out=input_prefix, profiler=profiler, file=file_path, pfp=pfp_exe, mod=modulo, win=window_length)
    execute_command(command)
    command = "{profiler} {aupair} -w {win} -i {input} -o {out}".format(profiler=profiler, input=input_prefix, out=input_prefix + ".aup.deleted_ts", aupair=aupair_exe, mod=modulo, win=window_length)
    execute_command(command)
    input_prefix = input_prefix + ".aup"

    # ----------- computation of the BWT of the parsing
    start = time.time()
    parse_size = os.path.getsize(input_prefix + ".parse")/4
    if(parse_size >=  (2**32-1) ):
        print("Sorry, the parse contains %d words" %  parse_size )
        print("which is more than my current limit 2^32-2")
        print("Please re-run the program with a larger modulus")
        sys.exit(1)
    elif(parse_size >=  (2**31-1) ):
        command = "{exe} {file}".format(exe = parsebwt_exe64, file=input_prefix)
    else:
        command = "{exe} {file}".format(exe = parsebwt_exe, file=input_prefix)
    command += " -s"
    print("==== Computing BWT of parsing. Command:", command)
    if(execute_command(command)!=True):
        sys.exit(1)
    print("Elapsed time: {0:.4f}".format(time.time()-start));

    # ----------- compute final BWT using dictionary and BWT of parse
    start = time.time()
    if(os.path.getsize(input_prefix + ".dict") >=  (2**31-4) ):
        # 64 bit version with and without threads
        command = "{exe} -w {wsize} {file} -s -e".format(
            exe = pfbwtNT_exe64, wsize=window_length, file=input_prefix)
    else:  # 32 bit version
        command = "{exe} -w {wsize} {file} -s -e".format(
            exe = pfbwtNT_exe, wsize=window_length, file=input_prefix)

    print("==== Computing final BWT. Command:", command)
    if(execute_command(command)!=True):
        sys.exit(1)
    print("Elapsed time: {0:.4f}".format(time.time()-start))

    # ----------- compute first N samples using dictionary and BWT of parse
    if num_of_sequences == 0:
        print("The number of sequences needs to be specified for now")
        sys.exit(1)

    start = time.time()
    if(os.path.getsize(input_prefix + ".dict") >=  (2**31-4) ):
        # 64 bit version with and without threads
        command = "{exe} -N {N} -w {wsize} {file}".format(
            exe = pfbwtSANT_exe64, wsize=window_length, file=input_prefix, N=num_of_sequences)
    else:  # 32 bit version
        command = "{exe} -N {N} -w {wsize} {file}".format(
            exe = pfbwtSANT_exe, wsize=window_length, file=input_prefix, N=num_of_sequences)

    print("==== Computing first N SA samples. Command:", command)
    if(execute_command(command)!=True):
        sys.exit(1)
    print("Elapsed time: {0:.4f}".format(time.time()-start))
    print("Total construction time: {0:.4f}".format(time.time()-start))

    # ----------- compute
    command = "{estw} -i {i_prefix}".format(estw=estw_exe, i_prefix=input_prefix)
    execute_command(command)

    #remove_file(input_prefix + ".ssa")
    #remove_file(input_prefix + ".esa")
    #remove_file(input_prefix + ".nsa")
    #remove_file(input_prefix + ".ilist")
    #remove_file(input_prefix + ".bwlast")
    #remove_file(input_prefix + ".bwsai")

    mkdir_p(input_prefix)
    os.replace(input_prefix + ".rlebwt", input_prefix+ "/bwt.rle")
    os.replace(input_prefix + ".rlebwt.meta", input_prefix + "/bwt.rle.meta")
    os.replace(input_prefix + ".saes", input_prefix + "/samples.saes")

    if (check):
        command = "{profiler} {check} -i {i_prefix} -o {i_prefix}".format(profiler=profiler, i_prefix=input_prefix, check=check_exe)
        execute_command(command)
        command = "{} -i {}".format(check_sa_exe, input_prefix)
        execute_command(command)

    return input_prefix

#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path, check, merge_jobs, search_jobs):
    base_path = ""
    command = "{profiler} {rimerge} -a {left} -b {right} -o {out} -m {mj} -t {sj}".format(profiler=profiler, rimerge=rimerge_exe, left=left_file_path, right=right_file_path, out=out_file_path, mj=merge_jobs, sj=search_jobs)
    execute_command(command)

#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file name', type=str, dest="input")
    parser.add_argument('-w', '--wsize', help='sliding window size (def. 10)', default=10, type=int, dest="window")
    parser.add_argument('-p', '--mod', help='hash modulus (def. 100)', default=100, type=int, dest="module")
    parser.add_argument('-I', '--ignore', help='Ignore first I sequences', type=int, dest="ignore", default=0)
    parser.add_argument('-o', '--output', help='output file prefix', type=str, dest="output")
    parser.add_argument('-m', '--merge-jobs', help='number of merge jobs', type=int, default=4, dest="merge_jobs")
    parser.add_argument('-t', '--search-jobs', help='number of search jobs', type=int, default=1, dest="search_jobs")
    parser.add_argument('-n', '--seqs', help='number of sequences as percentage', default=1, type=int, dest="seqs")
    parser.add_argument('-b', '--blocks', help='number of blocks', default=1, type=int, dest="blocks")
    parser.add_argument('-B',  help='run Big-BWT ONLY',action='store_true', dest="bonly")
    parser.add_argument('-C',  help='check index structure after each merge',action='store_true', dest="check")
    parser.add_argument('-R',  help='Run on random sequences',action='store_true', dest="random")
    parser.add_argument('-S',  help='Skip Big-BWT',action='store_true', dest="skip")

    args = parser.parse_args()

    if (args.skip and args.bonly):
        print("Skip and Big-BWT only called toghether, error")

    print("Version: ")
    execute_command("{} --version".format(rimerge_exe), 10)

    # Generate the input file
    if (args.random and not args.skip):
        args.input = generate_random_fasta(args.seqs, 100000)
    elif (args.random and args.skip):
        cwd = os.getcwd()
        file_name = "{wd}/random_{n}.fasta".format(wd=cwd, n=args.seqs)
        args.input = file_name

    # Split the input file
    if (not args.skip):
        seqs_per_block = int(args.seqs / args.blocks)
        print("Splitting the input file, {} sequences per block".format(seqs_per_block))
        file_paths = split(args.input, os.getcwd(), args.seqs, args.blocks, args.ignore)
    else:
        file_paths = list()
        filename, file_extension = os.path.splitext(args.input)
        for b in range (0, args.blocks):
            output_path = filename + "." + str(args.seqs) + "." + str(b) + ".fa"
            file_paths.append(output_path)

    # Run Big-BWT on each block
    indexes = list()
    if (not args.skip):
        print("Running PFP++")
        for i in range(0, len(file_paths)):
            index_name = build_rindex(file_paths[i], seqs_per_block + 1, args.window, args.module, args.check)
            indexes.append(index_name)
    else:
        for file_name in file_paths:
            indexes.append(os.path.basename(os.path.splitext(file_name)[0]))

    # Iterative Merge
    if (not args.bonly):
        print("Merging")
        print("Output file prefix: {}".format(args.output))
        mkdir_p(args.output)
        run_merge(indexes[0], indexes[1], args.output, args.check, args.merge_jobs, args.search_jobs)
        if (args.check):
            command = "{} -i {} -o {}".format(check_exe, args.output, args.output)
            execute_command(command)
            if (len(file_paths) == 2):
                command = "{} -i {}".format(check_sa_exe, args.output)
                execute_command(command)
        for i in range(2, len(file_paths)):
            run_merge(args.output, indexes[i], args.output, args.check, args.merge_jobs, args.search_jobs)
            if (args.check):
                command = "{} -i {} -o {}".format(check_exe, args.output, args.output)
                execute_command(command)
                command = "{} -i {}".format(check_sa_exe, args.output)
                execute_command(command)

if __name__ == '__main__':
    main()