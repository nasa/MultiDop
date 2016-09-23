from __future__ import print_function
import subprocess
import shlex
import os


def run_command(command):
    """
    This function captures text output from a command-line program and
    displays it in the Python shell. Currently does not print the output
    unitl after the command has completed.
    """
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE,
                               bufsize=0)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc


def do_analysis(param_file, cmd_path='../src/DDA'):
    """
    Execute a 3DVAR multi-Doppler analysis. The analysis will be output to
    the output file specified within the param_file.

    Parameters
    ----------
    param_file : str
        Name of file containing the parameters needed by the DDA engine.

    Other Parameters
    ----------------
    cmd_path : str
        Full path to the DDA executable, including the executable itself.
    """
    os.system('cp ' + cmd_path + ' .')
    rc = run_command('./DDA ' + param_file)
    os.remove('./DDA')
