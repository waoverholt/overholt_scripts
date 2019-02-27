#!/usr/bin/env python

import sys
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC


in_file = sys.argv[1]
out_file = sys.argv[2]
start = int(sys.argv[3]) - 1
end = int(sys.argv[4])

fp = os.path.abspath(in_file)

seq = next(SeqIO.parse(fp, "fasta"))

#print end
#print seq[start:end].format("fasta")
print seq[start:end].reverse_complement().format("fasta")

#print seq.format("fasta")



