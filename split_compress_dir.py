#!/usr/bin/env python3

from doctest import OutputChecker
import os
import sys
import re
import argparse
import shutil
import subprocess

"""Split a large directory tree into sub tar.gz files

This script is intended for splitting up a large minION output dir
into equal sized, compressed files that can then be uploaded to a cloud
storage in parts

This script requires `fpart` to be installed in the operating system
"""

def get_args() -> argparse.ArgumentParser:
    """ Get commandline arguments """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input directory", required=True)
    parser.add_argument("-o", "--output", help="output path to store all tar.gz compressed subdirs", default="out")
    parser.add_argument("-s", "--partition-size", 
                help="approximate size of each subfile before compression", 
                default="20Gb",
                dest="size")
    parser.add_argument("-t", "--threads", default=1, help="number of threads for xargs")
    parser.add_argument("--overwrite", action="store_true")

    return parser

def check_inputs(
    args: argparse.ArgumentParser,
    inputs: dict
    ) -> dict:
    """ Checks if the inputs are supplied """

    #check input
    if not os.path.exists(args.input):
        sys.exit("Input folder path does not exist")

    #check output
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    else:
        if args.overwrite:
            shutil.rmtree(args.output)
            os.makedirs(args.output)
        else:
            print(f"Output folder {os.path.abspath(args.output)} exists and '--overwrite' not provided")
            sys.exit(1)
    size = args.size
    size_by = parse_size(size)
    #set inputs
    inputs["input"] = args.input
    inputs["output"] = args.output
    inputs["size"] = size_by
    inputs["threads"] = args.threads
    inputs["overwrite"] = args.overwrite
    return inputs

def run_cmd_shell(command_str: str):
    try:
        output = subprocess.check_output(command_str, encoding='utf-8', shell=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Encountered an error executing the command: {command_str}")
        print(f"Error details:")
        print(f"Exit code={e.returncode}")
        print(f"Error message={e.output}")
        sys.exit(1)
    print(f"Command ran successfully")
    return output

def parse_size(size):
    """Calculate the size in bytes given human readable input"""
    units = {
        "B": 1, 
        "KB": 2**10,
        "MB": 2**20,
        "GB": 2**30,
        "TB": 2**40
        }
    size = size.upper()
    if not re.match(r' ', size):
        size = re.sub(r'([KMGT]?B)', r' \1', size)
    number, unit = [string.strip() for string in size.split()]
    return int(float(number)*units[unit])

def gen_partitions(inputs):
    """Run fpart to generate the partition files"""
    full_in_path = os.path.abspath(inputs["input"])
    full_out_path = os.path.abspath(inputs["output"])
    folder_name = os.path.basename(full_in_path)
    fpart_out = f"{os.path.join(full_out_path,folder_name)}_part"
    cmd = f"fpart -o {fpart_out} -s {inputs['size']} {full_in_path}"
    output = run_cmd_shell(cmd)

def run_tar(inputs):
    """
    ls | grep -E "*.[0-9]+" | xargs -P 8 tar -zcf {}.tar.gz -T {}
    """

    """
    p1 = subprocess.Popen(["find", inputs["output"], "-regex", "'.*/.*[0-9]+'"], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(
        ["xargs", 
        "-P", inputs['threads'],
        "-n", "1", 
        "-I", "file", 
        "tar", "-zcf", 
        "file.tar.gz", 
        "-T", "file"], 
        stdin=p1.stdout, stdout=subprocess.PIPE, shell=True)
    p1.stdout.close()
    output,err = p2.communicate()
    print(output, err)
    """
    cmd = f"find {inputs['output']} -regex '.*/.*[0-9]+' | xargs -P {inputs['threads']} -I file tar -zcf file.tar.gz -T file"
    output = subprocess.Popen(cmd, shell=True, encoding="utf-8", stderr=subprocess.PIPE)

def main():
    """ main code block """
    inputs = {
        'input' : None,
        'output' : "out",
        'size' : "20 Gb"
    }

    parser = get_args()
    args = parser.parse_args()
    inputs = check_inputs(args, inputs)
    gen_partitions(inputs)
    run_tar(inputs)

if __name__ == "__main__":
    main()

