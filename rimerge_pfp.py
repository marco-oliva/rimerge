#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, random, string, shutil, errno, logging, struct
from psutil import virtual_memory
import numpy as np

Description = '''
Rimerge from PFP++
'''

#------------------------------------------------------------
# Absolute path of auxiliary executables
dirname                 = os.path.dirname(os.path.abspath(__file__))
rimerge_exe             = os.path.join(dirname, 'rimerge.x')
estw_exe                = os.path.join(dirname, 'estw.x')
check_exe               = os.path.join(dirname, 'check.x')
check_sa_exe            = os.path.join(dirname, 'check_sa.x')
rle_exe                 = os.path.join(dirname, 'rle.x')
cfa_exe                 = os.path.join(dirname, 'cfa.x')

parsebwt_exe            = os.path.join(dirname, 'bwtparse')
parsebwt_exe64          = os.path.join(dirname, 'bwtparse64')

pfbwtNT_exe             = os.path.join(dirname, 'pfbwtNT.x')
pfbwtNT_exe64           = os.path.join(dirname, 'pfbwtNT64.x')

pfbwtSANT_exe           = os.path.join(dirname, 'pfbwtSANT.x')
pfbwtSANT_exe64         = os.path.join(dirname, 'pfbwtSANT64.x')

pfp_thresholds          = os.path.join(dirname, 'pfp_thresholds')
pfp_thresholds64        = os.path.join(dirname, 'pfp_thresholds64')

repair_exe              = os.path.join(dirname,'irepair')
largerepair_exe         = os.path.join(dirname,'largeb_irepair')
despair_exe             = os.path.join(dirname,'despair')
integer_despair_exe     = os.path.join(dirname,'idespair')
preprocess_exe          = os.path.join(dirname,'procdic')
integer_preprocess_exe  = os.path.join(dirname,'iprocdic')
postprocess_exe         = os.path.join(dirname,'postproc')
integer_postprocess_exe = os.path.join(dirname,'ipostproc')

shaped_slp              = os.path.join(dirname, 'SlpEncBuild')

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

def remove_dir(dir_path):
    shutil.rmtree(dir_path)

def copy_file(src_path, dest_path):
    shutil.copyfile(src_path, dest_path)

def move_file(src_path, dest_path):
    os.replace(src_path, dest_path)

def move_dir_content(dir_source, dir_dest):
    files = os.listdir(dir_source)
    for f in files:
        shutil.move(dir_source + '/' + f, dir_dest)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python ≥ 2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise # nop

# ------------------------------------------------------------
# execute command: return command's stdout if everything OK, None otherwise
def execute_command(command, out_file_path='', time_it=False, seconds=10000000):
    rootLogger = logging.getLogger()
    try:
        if time_it:
            command = '/usr/bin/time --verbose {}'.format(command)
        rootLogger.info('Executing: {}'.format(command))
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, err = process.communicate()
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        rootLogger.info('Error executing command line')
        return None
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        rootLogger.info('Command exceeded timeout')
        return None
    if output and out_file_path != '':
        rootLogger.info('Writing output to {}'.format(out_file_path))
        with open(out_file_path, 'wb') as out_file:
            out_file.write(output)
    elif output:
        output = output.decode('utf-8')
        rootLogger.info('CMD stoud:\n' + output)
    if err:
        err = err.decode('utf-8')
        rootLogger.info('CMD stderr:\n' + err)
    return output


#------------------------------------------------------------
# Run BigBWT
def build_rindex(input_prefix, output_prefix, window_length, num_of_sequences, check, text_start, text_end, build_moni=False):
    root_logger = logging.getLogger()

    # ----------- computation of the BWT of the parsing
    start = time.time()
    parse_size = os.path.getsize(input_prefix + '.parse')/4
    if parse_size >= (2 ** 32 - 1):
        root_logger.info('Sorry, the parse contains {} words'.format(parse_size))
        root_logger.info('which is more than my current limit 2^32-2')
        root_logger.info('Please re-run the program with a larger modulus')
        sys.exit(1)
    elif parse_size >= (2 ** 31 - 1):
        command = '{exe} {file}'.format(exe = parsebwt_exe64, file=input_prefix)
    else:
        command = '{exe} {file}'.format(exe = parsebwt_exe, file=input_prefix)
    command += ' -s'
    root_logger.info('==== Computing BWT of parsing.')
    execute_command(command)
    root_logger.info('==== Done')

    # ----------- compute first N samples using dictionary and BWT of parse
    if num_of_sequences == 0:
        root_logger.info('The number of sequences needs to be specified for now')
        sys.exit(1)

    start = time.time()
    if os.path.getsize(input_prefix + '.dict') >= (2 ** 31 - 4):
        # 64 bit version with and without threads
        command = '{exe} -N {N} -w {wsize} {file}'.format(
            exe = pfbwtSANT_exe64, wsize=window_length, file=input_prefix, N=num_of_sequences)
    else:  # 32 bit version
        command = '{exe} -N {N} -w {wsize} {file}'.format(
            exe = pfbwtSANT_exe, wsize=window_length, file=input_prefix, N=num_of_sequences)

    root_logger.info('==== Computing first N SA samples')
    execute_command(command)
    root_logger.info('==== Done')

    if build_moni:
        # ----------- compute final BWT, thresholds, ssa and esa
        start = time.time()
        parse_size = os.path.getsize(input_prefix + '.parse')/4
        dictionary_size = os.path.getsize(input_prefix + '.dict')

        if parse_size >= (2 ** 31 - 1) or dictionary_size >= (2 ** 31 - 4):
            command = '{exe} {file} -w {wsize}'.format(exe=pfp_thresholds64,wsize=window_length, file=input_prefix)
        else:
            command = '{exe} {file} -w {wsize}'.format(exe=pfp_thresholds,wsize=window_length, file=input_prefix)

        root_logger.info('==== Computing Thresholds. ')
        execute_command(command)
        root_logger.info('==== Done')

        # ----------- add random access

        # ---- preprocess the dictionary
        start = time.time()
        command = '{exe} {file}.dicz'.format(exe=preprocess_exe, file=input_prefix)
        root_logger.info('==== Preprocessing the dictionary')
        execute_command(command)
        preprocess_time = time.time()-start
        root_logger.info('==== Done')

        # ---- apply repair to the modified dictionary
        start = time.time()

        mem = virtual_memory()
        repair_mem  = round(mem.total / 1024 / 1024) # total physical memory available in MB
        root_logger.info('RePair maximum memory: {}'.format(repair_mem))

        command = '{exe} {file}.dicz.int {mb}'.format(mb=repair_mem, exe=largerepair_exe, file=input_prefix)
        root_logger.info('==== Repair dictionary.')
        execute_command(command)
        repair_time = time.time()-start
        root_logger.info('==== Done')

        # ---- apply repair to the parse
        start = time.time()
        command = '{exe} {file}.parse {mb}'.format(mb=repair_mem,exe=largerepair_exe, file=input_prefix)
        root_logger.info('==== Repair parse.')
        execute_command(command)
        repair_time = time.time()-start
        root_logger.info('==== Done')

        # ---- postprocess
        start = time.time()
        command = '{exe} {file}'.format(exe=postprocess_exe, file=input_prefix)
        root_logger.info('==== Postprocessing the dictionary.')
        execute_command(command)
        postprocess_time = time.time()-start
        root_logger.info('==== Done')

        # ---- ShapeSLP
        start = time.time()
        command = '{exe} -i {file} -o {outfile}.{ext} -e {grammar} -f Bigrepair'.format(
            exe=shaped_slp, file=input_prefix,ext='slp', outfile=input_prefix,
            grammar='SelfShapedSlp_SdSd_Sd')
        root_logger.info('==== ShapedSLP construction.')
        execute_command(command)
        preprocess_time = time.time()-start
        root_logger.info('==== Done')

        exts_to_remove = ['.ssa', '.esa', '.nsa', '.ilist', '.bwlast', '.bwsai', '.C', '.R',
         '.parse.C', '.parse.R', '.dicz.int', '.dicz.int.C', '.dicz.int.R']
        remove_file('rs_temp_output')

    else:
        # ----------- compute final BWT using dictionary and BWT of parse
        start = time.time()
        if(os.path.getsize(input_prefix + '.dict') >=  (2**31-4) ):
            # 64 bit version with and without threads
            command = '{exe} -w {wsize} {file} -s -e'.format(exe=pfbwtNT_exe64, wsize=window_length, file=input_prefix)
        else:  # 32 bit version
            command = '{exe} -w {wsize} {file} -s -e'.format(exe=pfbwtNT_exe, wsize=window_length, file=input_prefix)

        root_logger.info('==== Computing final BWT.')
        execute_command(command)
        root_logger.info('==== Done')

        exts_to_remove = ['.ssa', '.esa', '.nsa', '.ilist', '.bwlast', '.bwsai']

    # ----------- convert samples in rimerge format
    command = '{estw} -i {i_prefix}'.format(estw=estw_exe, i_prefix=input_prefix)
    execute_command(command)

    root_logger.info('==== Done')

    for ext_to_remove in exts_to_remove:
        remove_file(input_prefix + ext_to_remove)

    mkdir_p(output_prefix + '_idx')
    move_file(input_prefix + '.rlebwt', output_prefix + '_idx'+ '/bwt.rle')
    move_file(input_prefix + '.rlebwt.meta', output_prefix + '_idx' + '/bwt.rle.meta')
    move_file(input_prefix + '.saes', output_prefix + '_idx' + '/samples.saes')
    if build_moni:
        move_file(input_prefix + '.slp', output_prefix + '_idx' + '/grammar_0.slp')
        move_file(input_prefix + '.thr', output_prefix + '_idx' + '/thresholds.thr')
        move_file(input_prefix + '.thr_pos', output_prefix + '_idx' + '/thresholds.thr_pos')

        with open(output_prefix + '_idx/random_access.meta', 'wb') as metadata_file:
            metadata_file.write(int(1).to_bytes(8, 'little'))
            metadata_file.write(text_start.to_bytes(8, 'little'))
            metadata_file.write(text_end.to_bytes(8, 'little'))
            metadata_file.write(len('grammar_0.slp').to_bytes(8, 'little'))
            metadata_file.write('grammar_0.slp'.encode('UTF-8'))

    if (check):
        command = '{check} -i {i_prefix}_idx -o {i_prefix}_errors'.format(i_prefix=input_prefix, check=check_exe)
        execute_command(command)
        command = '{} -i {}'.format(check_sa_exe, input_prefix + '_idx')
        execute_command(command)

#------------------------------------------------------------
# Run merge
def run_merge(left_file_path, right_file_path, out_file_path, merge_jobs, search_jobs, build_moni=False):
    # merge indexes
    command = '{rimerge} -a {left} -b {right} -o {out} -m {mj} -t {sj}'.format(rimerge=rimerge_exe, left=left_file_path, right=right_file_path, out=out_file_path, mj=merge_jobs, sj=search_jobs)
    execute_command(command)

    if (build_moni):
        # merge grammars for random access
        grammars_left  = list()
        grammars_right = list()

        # from left
        max_from_left = np.uint64(0)
        with open(left_file_path + '/random_access.meta', 'rb') as ra_meta:
            n_of_grammar = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
            i = np.uint64(0)
            while i < n_of_grammar:
                i += 1
                text_start  = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
                text_end    = np.uint64(struct.unpack('<Q', ra_meta.read(8)))

                if (text_end > max_from_left):
                    max_from_left = text_end

                name_length = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
                name        = struct.unpack('{}s'.format(name_length.item()), ra_meta.read(name_length.item()))
                name        = name[0].decode('UTF-8')
                grammar = [name, text_start, text_end]
                grammars_left.append(grammar)

        # from right
        with open(right_file_path + '/random_access.meta', 'rb') as ra_meta:
            n_of_grammar = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
            i = np.uint64(0)
            while i < n_of_grammar:
                i += 1
                text_start  = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
                text_end    = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
                name_length = np.uint64(struct.unpack('<Q', ra_meta.read(8)))
                name        = struct.unpack('{}s'.format(name_length.item()), ra_meta.read(name_length.item()))
                name        = name[0].decode('UTF-8')
                grammar = [name, text_start, text_end]
                grammars_right.append(grammar)

        # merge grammars into the merged index folder
        id = 0
        with open(out_file_path + '/random_access.meta', 'wb') as metadata_file:
            metadata_file.write(int(len(grammars_left) + len(grammars_right)).to_bytes(8, 'little'))

            for grammar in grammars_left:
                grammar_new_name = 'grammar_{}.slp'.format(id)
                id += 1

                metadata_file.write(grammar[1].item().to_bytes(8, 'little'))
                metadata_file.write(grammar[2].item().to_bytes(8, 'little'))
                metadata_file.write(len(grammar_new_name).to_bytes(8, 'little'))
                metadata_file.write(grammar_new_name.encode('UTF-8'))

                copy_file(left_file_path + '/' + grammar[0], out_file_path + '/' + grammar_new_name)

            for grammar in grammars_right:
                grammar_new_name = 'grammar_{}.slp'.format(id)
                id += 1

                metadata_file.write(((grammar[1] + max_from_left).item()).to_bytes(8, 'little'))
                metadata_file.write(((grammar[2] + max_from_left).item()).to_bytes(8, 'little'))
                metadata_file.write(len(grammar_new_name).to_bytes(8, 'little'))
                metadata_file.write(grammar_new_name.encode('UTF-8'))

                copy_file(right_file_path + '/' + grammar[0], out_file_path + '/' + grammar_new_name)

#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-a', '--A-prefix', help='input A prefix', type=str, dest='input_a', required=True)
    parser.add_argument('-b', '--B-prefix', help='input prefiix', type=str, dest='input_b')
    parser.add_argument('-w', '--wsize', help='sliding window size (def. 10)', default=10, type=int, dest='window')
    parser.add_argument('-p', '--modulo', help='modulo (def. 100)', default=100, type=int, dest='modulo')
    parser.add_argument('-o', '--output', help='output index prefix', type=str, dest='output', required=True)
    parser.add_argument('-m', '--merge-jobs', help='number of merge jobs', type=int, default=4, dest='merge_jobs')
    parser.add_argument('-t', '--search-jobs', help='number of search jobs', type=int, default=1, dest='search_jobs')
    parser.add_argument('-n', '--sequences', help='number of sequences', type=int, default=0, dest='num_of_sequences')
    parser.add_argument('-C', '--check', help='check output index structure', action='store_true',default=False, dest='check')
    parser.add_argument('--text-start-0-based', type=int, dest='text_start', default=0)
    parser.add_argument('--text-end-0-based', type=int, dest='text_end', default=0)
    parser.add_argument('--moni', help='build other MONI structures', action='store_true',default=False, dest='build_moni')

    args = parser.parse_args()

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)

    root_logger.info('rimerge Version: ')
    execute_command('{rimerge} --version'.format(rimerge=rimerge_exe))

    # build index from pfp, A given as an index B given as pfp
    if (args.input_b is None):
        root_logger.info('Build index')
        build_rindex(args.input_a, args.output, args.window, args.num_of_sequences, args.check, args.text_start, args.text_end, args.build_moni)
    else:
        if (not os.path.isdir(args.input_b)):
            root_logger.info('Build index of {}'.format(args.input_b))
            build_rindex(args.input_b, args.output, args.window, args.num_of_sequences, args.check, args.text_start, args.text_end, args.build_moni)
        else:
            root_logger.info('input_b is already an index, not re-computing')
        root_logger.info('Start merging')
        mkdir_p(args.output + '_idx')
        run_merge(args.input_a, args.input_b, args.output + '_idx', args.merge_jobs, args.search_jobs, args.build_moni)

    if (args.check):
        command = '{} -i {} -o {}'.format(check_exe, args.output + '_idx', args.output + '_idx')
        execute_command(command)
        command = '{} -i {}'.format(check_sa_exe, args.output + '_idx')
        execute_command(command)

if __name__ == '__main__':
    main()