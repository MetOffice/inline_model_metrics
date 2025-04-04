import pandas as pd
import glob, os, subprocess, numpy

suites = ['u-cx129']
moo_path = 'moose:/crum/{suite}/ady.file/'
fnames_old = '{runid}_tracknodeedit_*_6h_tc_psl.csv'
fnames_new = '{runid}_track_{date}_6h_tc_psl.csv'
temp_dir = '/data/scratch/malcolm.roberts/nodeedit_files'
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)

extra_columns = ['radius_8', 'ace', 'pace', 'ike', 'pdi', 'max_core']

def run_cmd(
        cmd,
    check=False
):
    sts = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=check,
    )
    print('cmd output ', cmd, sts.stderr)
    if sts.stderr:
        if 'Warning' in sts.stderr or 'ALREADY_EXISTS' in sts.stderr:            
            msg = (
                "Warning found in cat output ", sts.stderr
            )
            print(msg)
        else:
            msg = (
                "Error found in cat output ", sts.stderr
            )
            raise RuntimeError(msg)
    return sts

def get_nodeedit_files(usite):
    cmd = 'moo get -i '+moo_path.format(suite=suite)+fnames_old.format(runid=suite[2:])+' '+temp_dir
    print(cmd)
    run_cmd(cmd)

    files = glob.glob(os.path.join(temp_dir, fnames_old.format(runid=suite[2:])))
    print(files)
    if len(files) > 0:
        return files
    else:
        raise Exception('No files to convert')

def put_files(fnames):
    for f in fnames:
        cmd = 'moo put '+f+' '+moo_path.format(suite=suite)
        print(cmd)
        run_cmd(cmd)

def convert_files(files, suite):
    files_to_archive = []
    for f in files:
        fbase = os.path.basename(f)
        time_period = f.split('_')[3]
        fname_new = fnames_new.format(runid=suite[2:], date=time_period)
        df = pd.read_csv(f, engine='python', sep=', ')
        df = df.drop(columns=extra_columns)
        #df.to_csv(os.path.join(temp_dir, fname_new), sep=',', index=False)
        csv_file_delimiter = ', '
        numpy.savetxt(os.path.join(temp_dir, fname_new), df, delimiter=csv_file_delimiter, header=csv_file_delimiter.join(df.columns.values), fmt=['%i', '%i', '%i', '%i', '%i', '%i', '%i', '%s', '%s', '%.6e', '%.6e', '%.6e', '%.6e', '%.6e', '%.6e'], comments='', encoding=None)
        files_to_archive.append(os.path.join(temp_dir, fname_new))
        return files_to_archive

def work(suite):
    files = get_nodeedit_files(suite)
    files_to_archive = convert_files(files, suite)
    put_files(files_to_archive)

if __name__ == '__main__':
    for suite in suites:
        work(suite)

    
