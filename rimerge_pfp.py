#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, random, string, shutil, errno
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq

Description = """
Rimerge from PFP++
"""

#------------------------------------------------------------
# Absolute path of auxiliary executables
dirname         = os.path.dirname(os.path.abspath(__file__))
rimerge_exe     = os.path.join(dirname, "rimerge.x")
estw_exe        = os.path.join(dirname, "estw.x")
check_exe       = os.path.join(dirname, "check.x")
check_sa_exe    = os.path.join(dirname, "check_sa.x")
rle_exe         = os.path.join(dirname, "rle.x")
cfa_exe         = os.path.join(dirname, "cfa.x")

parsebwt_exe    = os.path.join(dirname, "bwtparse")
parsebwt_exe64    = os.path.join(dirname, "bwtparse64")

pfbwtNT_exe       = os.path.join(dirname, "pfbwtNT.x")
pfbwtNT_exe64       = os.path.join(dirname, "pfbwtNT64.x")

pfbwtSANT_exe       = os.path.join(dirname, "pfbwtSANT.x")
pfbwtSANT_exe64       = os.path.join(dirname, "pfbwtSANT64.x")

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

def remove_dir(dir_path):
    shutil.rmtree(dir_path)

def move_dir_content(dir_source, dir_dest):
    files = os.listdir(dir_source)
    for f in files:
        shutil.move(dir_source + "/" + f, dir_dest)

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
def execute_command(command, seconds = 100000000):
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
def build_rindex(input_prefix, window_length, modulo, num_of_sequences, check):
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

    remove_file(input_prefix + ".ssa")
    remove_file(input_prefix + ".esa")
    remove_file(input_prefix + ".nsa")
    remove_file(input_prefix + ".ilist")
    remove_file(input_prefix + ".bwlast")
    remove_file(input_prefix + ".bwsai")

    mkdir_p(input_prefix + "_idx")
    os.replace(input_prefix + ".rlebwt", input_prefix + "_idx"+ "/bwt.rle")
    os.replace(input_prefix + ".rlebwt.meta", input_prefix + "_idx" + "/bwt.rle.meta")
    os.replace(input_prefix + ".saes", input_prefix + "_idx" + "/samples.saes")

    if (check):
        command = "{check} -i {i_prefix}_idx -o {i_prefix}_errors".format(i_prefix=input_prefix, check=check_exe)
        execute_command(command)
        command = "{} -i {}".format(check_sa_exe, input_prefix + "_idx")
        execute_command(command)

#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path, merge_jobs, search_jobs):
    command = "{rimerge} -a {left} -b {right} -o {out} -m {mj} -t {sj}".format(rimerge=rimerge_exe, left=left_file_path, right=right_file_path, out=out_file_path, mj=merge_jobs, sj=search_jobs)
    execute_command(command)

#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--A-prefix', help='input A prefix', type=str, dest="input_a", required=True)
    parser.add_argument('-b', '--B-prefix', help='input prefiix', type=str, dest="input_b")
    parser.add_argument('-w', '--wsize', help='sliding window size (def. 10)', default=10, type=int, dest="window")
    parser.add_argument('-p', '--modulo', help='modulo (def. 100)', default=100, type=int, dest="modulo")
    parser.add_argument('-o', '--output', help='output index prefix', type=str, dest="output", required=True)
    parser.add_argument('-m', '--merge-jobs', help='number of merge jobs', type=int, default=1, dest="merge_jobs")
    parser.add_argument('-t', '--search-jobs', help='number of search jobs', type=int, default=1, dest="search_jobs")
    parser.add_argument('-n', '--sequences', help='number of sequences', type=int, default=0, dest="num_of_sequences")
    parser.add_argument('-C', '--check', help='check output index structure', action='store_true',default=False, dest="check")
    args = parser.parse_args()

    print("Version: ")
    execute_command("{rimerge} --version".format(rimerge=rimerge_exe), 10)

    # build index from pfp, A given as an index B given as pfp
    if (args.input_b is None):
        print("Build index")
        build_rindex(args.input_a, args.window, args.modulo, args.num_of_sequences, args.check)
    else:
        if (not os.path.isdir(args.input_b + "_idx")):
            print("Build index of {}".format(args.input_b))
            build_rindex(args.input_b, args.window, args.modulo, args.num_of_sequences, args.check)
        else:
            print("input_b is already an index, not re-computing")
        print("Start merging")
        mkdir_p(args.output + "_idx")
        run_merge(args.input_a, args.input_b + "_idx", args.output + "_idx", args.merge_jobs, args.search_jobs)

    if (args.check):
        command = "{} -i {} -o {}".format(check_exe, args.output + "_idx", args.output + "_idx")
        execute_command(command)
        command = "{} -i {}".format(check_sa_exe, args.output + "_idx")
        execute_command(command)

if __name__ == '__main__':
    main()