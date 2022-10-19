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

def check_prog(prog_name: str) -> None:
    if shutil.which(prog_name) is None:
        print(f"{prog_name} is not in the $PATH and is required for this program")
        sys.exit(1)

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
    print(f"Command: {command_str}\nran successfully")
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
    cmd = f"find {inputs['output']} -regex '.*/.*[0-9]+' | xargs -P {inputs['threads']} -I file tar -zcf file.tar.gz -T file"
    proc = run_cmd_shell(cmd)
    #proc = subprocess.check_call(cmd, shell=True, encoding="utf-8", stderr=subprocess.PIPE)

def calc_tar_file_sizes(inputs:dict) -> list:
    file_sizes = [["file_name", "size_bytes", "size"]]
    abbrevs = (
        (1<<50, 'PB'),
        (1<<40, 'TB'),
        (1<<30, 'GB'),
        (1<<20, 'MB'),
        (1<<10, 'kB'),
        (1, 'bytes')
    )
    for file in os.listdir(inputs["output"]):
        if file[-2:] == "gz":
            bytes = os.path.getsize(os.path.join(inputs["output"],file))
            if bytes == 1:
                file_sizes.append([file, bytes, f"1 byte"])
            for factor, suffix in abbrevs:
                if bytes >= factor:
                    break
            file_sizes.append([file, bytes, f"{round(bytes / factor, 1)} {suffix}"])
            #file_sizes[file] = f"{bytes}\t{round(bytes / factor, 1)} {suffix}"
    return file_sizes

def print_align(array:list, summary_file:str) -> None:
    max_widths = []
    with open(summary_file, "w") as outfp:
        for column in zip(*array):
            max_widths.append(max(map(len, map(str, column))))
        for row in array:
            for width, cell in zip(max_widths, row):
                print(str(cell).ljust(width+5), end="", file=outfp)
            print(file=outfp)

def main():
    """ main code block """
    check_prog("fpart")
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
    file_sizes = calc_tar_file_sizes(inputs)
    print_align(file_sizes, os.path.join(inputs["output"], "sizes_of_tar_files.txt"))

if __name__ == "__main__":
    main()

