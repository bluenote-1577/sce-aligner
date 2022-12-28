import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import sys
import pysam

sam_files = sys.argv[1:]
passing_reads = dict()
for sam_file in sam_files:
    sam = pysam.AlignmentFile(sam_file, "r")
    for read in sam.fetch():
        if read.has_tag('de'):
            divergence = read.get_tag('de')
            if divergence > 0.020 and divergence < 0.080 and not read.is_secondary:
                passing_reads[read.query_name] = divergence
    sam.close()

pass_reads_file = open('passing_reads.txt', 'w')
for read,div in passing_reads.items():
    pass_reads_file.write(f"{read}\t{div}\n")
    
