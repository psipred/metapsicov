import sys
from shutil import copyfile
import os
import subprocess
#  Test if hhaln has enough sequences (3000).
#  If not run jack_hhblits
#  either way output .jackaln
# sys.argv[1] - the hhaln file_len
# sys.argv[2] - tmp location
# sys.argv[3] - db, /scratch1/NOT_BACKED_UP/dbuchan/uniref/unirefmain.fasta


def file_len(fname):
    with open(fname) as f:
        i = None
        for i, l in enumerate(f):
            pass
        if i:
            return i + 1
        else:
            return 0


def run_exe(args, name):
    """
        Function takes a list of command line args and executes the subprocess.
        Sensibly tries to catch some errors too!
    """
    code = 0
    print("Running "+name)
    try:
        print(' '.join(args))
        code = subprocess.call(' '.join(args), shell=True,
                               env={"PATH": os.environ.get('PATH')+":"+args[2],
                                    "HHLIB": os.environ.get('HHLIB'),
                                    }
                               )
    except Exception as e:
        print('type is:', e.__class__.__name__)
        print("Call Error "+str(e))
        sys.exit(1)
    if code != 0:
        print(name+" Non Zero Exit status: "+str(code))
        sys.exit(code)


hhaln_length = file_len(sys.argv[1])
file_id = sys.argv[1][:-6]
if hhaln_length < 3000:
    bin_dir = os.path.dirname(os.path.realpath(sys.argv[0]))[:-3] + \
              "bin/"
    jackhhblits = bin_dir+"jack_hhblits"
    tmpdir = sys.argv[2]+"/"
    processJackhhblits_args = [jackhhblits,
                               file_id,
                               bin_dir,
                               tmpdir,
                               sys.argv[3],
                               '4',
                               # '12'
                               ]
    # print(processJackhhblits_args)
    run_exe(processJackhhblits_args, "jackhhblits")
else:
    print(file_id)
    f = open(file_id+".jackaln", "w")
    f.close()

jackaln_length = file_len(file_id+".jackaln")

if jackaln_length > hhaln_length:
    copyfile(file_id+".jackaln", file_id+".aln")
else:
    copyfile(sys.argv[1], file_id+".aln")
