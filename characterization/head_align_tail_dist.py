#!/usr/bin/env python
"""
Written by Saber HafezQorani on March, 2018
To get the length of head, aligned, and tail regions of an alignment.

Major changes from NanoSim to use minimap output sam file and also develop 2d length distribution.
"""

from __future__ import with_statement
import sys
import getopt
import numpy
import matplotlib.pyplot as plt
try:
    from six.moves import xrange
except ImportError:
    pass


def flex_bins(num_of_bins, ratio_dict, num_of_reads):
    count_reads = num_of_reads / num_of_bins
    k_of_bin = 0
    k_of_ratio = 0
    ratio_keys = sorted(ratio_dict.keys())
    num_of_keys = len(ratio_keys)

    ratio_bins = {}
    while k_of_bin < num_of_bins:
        if k_of_ratio >= num_of_keys:
            break

        start = k_of_ratio
        count = len(ratio_dict[ratio_keys[k_of_ratio]])
        k_of_ratio += 1

        while k_of_ratio < num_of_keys:
            tmp_count = count + len(ratio_dict[ratio_keys[k_of_ratio]])
            if abs(tmp_count - count_reads) >= abs(count - count_reads):
                break
            else:
                count = tmp_count
                k_of_ratio += 1

        k = (ratio_keys[start] if start else 0,
             ratio_keys[k_of_ratio] if k_of_ratio < num_of_keys else ratio_keys[k_of_ratio - 1] + 1)
        ratio_bins[k] = []
        for i in xrange(start, k_of_ratio):
            ratio_bins[k].extend(ratio_dict[ratio_keys[i]])

        k_of_bin += 1

    if k_of_ratio < num_of_keys - 1:
        k = (ratio_keys[k_of_ratio], ratio_keys[num_of_keys - 1] + 1)
        ratio_bins[k] = []
        for i in xrange(k_of_ratio, num_of_keys - 1):
            ratio_bins[k].extend(ratio_dict[ratio_keys[i]])

    return ratio_bins

def get_head_tail(cigar_string):

    head_info = cigar_string[0]
    tail_info = cigar_string[-1]

    if head_info.type == "S":
        head = head_info.size
    else:
        head = 0

    if tail_info.type == "S":
        tail = tail_info.size
    else:
        tail = 0

    return head, tail

def head_align_tail(outfile, num_of_bins, alnm_ftype):

    out1 = open(outfile + '_read_rellen_ecdf', 'w')
    out2 = open(outfile + '_read_totallen_ecdf', 'w')
    out3 = open(outfile + '_ht_ratio', 'w')
    out4 = open(outfile + "_align_ratio", 'w')
    out5 = open(outfile + "_unaligned_length_ecdf", "w")

    aligned = [] # we do not need this in TransNanoSim but I am just including for test.
    total = []
    list_ref_len = []
    unaligned_length = []
    count_aligned = 0
    count_unaligned = 0

    dict_ht_ratio = {}
    dict_align_ratio = {}
    dict_rellen = {}

    if alnm_ftype == "sam":
        sam_reader = HTSeq.SAM_Reader
        alnm_file_sam = outfile + "_primary.sam"
        t_alignments = sam_reader(alnm_file_sam)
        for alnm in t_alignments:
            if alnm.aligned:
                count_aligned += 1

                ref = alnm.iv.chrom
                ref_len = dict_ref_len[ref]
                ref_aligned = alnm.iv.length
                aligned.append(ref_aligned)

                read_len_total = len(alnm.read.seq)
                total.append(read_len_total)
                head, tail = get_head_tail(alnm.cigar)
                middle = read_len_total - head - tail

                #ratio aligned part over total length of the read
                alignment_ratio = float(middle) / read_len_total
                if middle not in dict_align_ratio:
                    dict_align_ratio[middle] = [alignment_ratio]
                else:
                    dict_align_ratio[middle].append(alignment_ratio)

                head_and_tail = head + tail
                if head != 0:
                    ht_ratio = float(head) / head_and_tail
                    if head_and_tail not in dict_ht_ratio:
                        dict_ht_ratio[head_and_tail] = [ht_ratio]
                    else:
                        dict_ht_ratio[head_and_tail].append(ht_ratio)

                #relative_length : total len of read over total len of reference transcriptome
                relative_length = float(read_len_total) / ref_len
                if ref_len not in dict_rellen:
                    dict_rellen[ref_len] = [relative_length]
                else:
                    dict_rellen[ref_len].append(relative_length)


    else:
        besthit_out = outfile + "_besthit.maf"
        with open(besthit_out, 'r') as f:
            for line in f:
                ref_line = line
                query_line = next(f)
                count_aligned += 1

                ref = ref_line.strip().split()
                aligned_ref = int(ref[3])
                aligned.append(aligned_ref)

                query = query_line.strip().split()
                head = int(query[2])
                middle = int(query[3])
                tail = int(query[5]) - int(query[2]) - int(query[3])
                total.append(int(query[5]))
                head_and_tail = int(query[5]) - int(query[3])
                alignment_ratio = float(query[3]) / float(query[5])
                if middle in dict_align_ratio:
                    dict_align_ratio[middle].append(alignment_ratio)
                else:
                    dict_align_ratio[middle] = [alignment_ratio]
                if head != 0:
                    ht_ratio = float(head) / head_and_tail
                    if head_and_tail in dict_ht_ratio:
                        dict_ht_ratio[head_and_tail].append(ht_ratio)
                    else:
                        dict_ht_ratio[head_and_tail] = [ht_ratio]


    # ecdf of length of aligned regions (2d length distribution) editted this part - Approach 2 relative length of ONT total over total length of the reference transcriptome it aligned to.
    rel_len_bins = flex_bins(num_of_bins, dict_rellen, count_aligned)
    rel_len_cum = dict.fromkeys(rel_len_bins.keys(), [])
    for key, value in rel_len_bins.items():
        hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ratio * 0.001)
        rel_len_cum[key] = cdf

    out1.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(rel_len_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out1.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t")
        for key in sorted(rel_len_cum.keys()):
            out1.write(str(rel_len_cum[key][i]) + "\t")
        out1.write("\n")
    out1.close()



    # ecdf of length of aligned reads
    max_length = max(total)
    hist_reads, bin_edges = numpy.histogram(total, bins=numpy.arange(0, max_length + 50, 50), density=True)
    cdf = numpy.cumsum(hist_reads * 50)
    out2.write("bin\t0-" + str(max_length) + '\n')
    for i in xrange(len(cdf)):
        out2.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t" + str(cdf[i]) + '\n')
    out2.close()



    # ecdf of head/total ratio
    ht_ratio_bins = flex_bins(num_of_bins, dict_ht_ratio, count_aligned)
    ht_cum = dict.fromkeys(ht_ratio_bins.keys(), [])
    for key, value in ht_ratio_bins.items():
        hist_ht, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ht * 0.001)
        ht_cum[key] = cdf

    out3.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(ht_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out3.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t")
        for key in sorted(ht_cum.keys()):
            out3.write(str(ht_cum[key][i]) + "\t")
        out3.write("\n")
    out3.close()



    # ecdf of align ratio
    align_ratio_bins = flex_bins(num_of_bins, dict_align_ratio, count_aligned)
    align_cum = dict.fromkeys(align_ratio_bins.keys(), [])
    for key, value in align_ratio_bins.items():
        hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ratio * 0.001)
        align_cum[key] = cdf

    out4.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(align_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out4.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t")
        for key in sorted(align_cum.keys()):
            out4.write(str(align_cum[key][i]) + "\t")
        out4.write("\n")
    out4.close()


    return count_aligned