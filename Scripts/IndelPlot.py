import pandas
import matplotlib.pyplot as plt
import matplotlib 
import pysam
import os 
import numpy
from collections import defaultdict
import glob2
import argparse
import seaborn


matplotlib.rcParams.update({'savefig.bbox':'tight'})


parser = argparse.ArgumentParser(description='Indel profile plot generation for a directory of BAM files', prog = 'BAM Indel Plot')
parser.add_argument('-i', '--inputDirectory')
parser.add_argument('-r', '--chrom')
parser.add_argument('-c', '--cutSite')
parser.add_argument('-b', '--buffer')
parser.add_argument('-t', '--tag')

args = parser.parse_args()

files = glob2.glob(args.inputDirectory)


cut_site = int(args.cutSite)
target_chrom = str(args.chrom)

# === Parameters ===
for b in range(len(files)):
    bam_file = files[b]
    print('ON: ' + str(bam_file))

    plotbuffer = int(args.buffer)
    start_pos = cut_site -plotbuffer             # Region start (0-based)
    end_pos = cut_site + plotbuffer                # Region end


    # === Initialize counters ===
    indel_counts = defaultdict(int)

    # === Open BAM file ===
    bam = pysam.AlignmentFile(os.path.join(bam_file), "rb")
    indels = pandas.DataFrame(numpy.arange(-100, 101, 1), columns = ['Size'])
    indels['Counts'] = 0
    indels.set_index(['Size'], inplace = True)
    
    # === Iterate through reads in region ===
    for read in bam.fetch(target_chrom, start_pos, end_pos):
        if read.is_unmapped:
            continue

        ref_pos = read.reference_start
        for cigar_op, cigar_len in read.cigartuples:
            if cigar_op == 1:  # Insertion
                # Count insertion at the base before insertion
                indels.loc[cigar_len, 'Counts'] = indels.loc[cigar_len, 'Counts'] + 1
            elif cigar_op == 2:  # Deletion
                indels.loc[-1*cigar_len, 'Counts'] = indels.loc[-1*cigar_len, 'Counts'] + 1
                ref_pos += cigar_len
            elif cigar_op in [0, 7, 8]:  # Match or mismatch
                ref_pos += cigar_len
            else:
                # Skip soft/hard clips and other ops
                pass

    bam.close()

    # === Prepare data for plotting ===
    positions = list(range(start_pos, end_pos))
    frequencies = [indel_counts.get(pos, 0) for pos in positions]

    indels.reset_index(inplace = True)
    
    x_smooth = numpy.linspace(indels['Size'].min(), indels['Size'].max(), 200)
    y_smooth = numpy.interp(x_smooth, indels['Size'], indels['Counts'])
    
    # === Plotting ===
    plt.figure(figsize=(12, 6))
    plt.figure()
    seaborn.displot(indels, x = 'Size', weights = 'Counts', kind = 'kde', bw_adjust = 0.05)
    expanded_data = indels.loc[indels.index.repeat(indels['Counts'])]
    plt.axvline(x=0, color='red', linestyle='--', label='Cut Site')
    plt.xlabel("Position", fontsize = 12)
    plt.ylabel("Cumulative Fraction of Reads", fontsize = 12)
    plt.title("Indel Frequency Relative to Cut Site", fontsize = 16)
    plt.legend()

    ax = plt.gca()
    ax.set_xlim([-25, 25])
    ax.set_ylim([0, 0.3])
    ax.tick_params(axis = 'x', labelsize = 12)
    ax.tick_params(axis = 'y', labelsize = 12)
    locs, labels = plt.xticks()
    plt.xticks(rotation = 60)
    plt.tight_layout()

    plt.savefig(str(args.tag)+"_indelPlot.png")
    
