#!/usr/bin/env python
"""
Written by Chen Yang on Mar 25th, 2015
To get the length of head, aligned, and tail regions of an alignment.

Major change in Apr 22nd

Updated in Nov 25th
"""

from __future__ import with_statement
import sys
import getopt
import numpy

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


def head_align_tail(outfile, num_of_bins):
    out5 = open(outfile + '_reftransc_totalONT_allbins_ecdf_v2', 'w')
    out1 = open(outfile + '_reftransc_totalONT_allbins_ecdf', 'w')
    out2 = open(outfile + '_aligned_reads_ecdf', 'w')
    out3 = open(outfile + '_ht_ratio', 'w')
    out4 = open(outfile + "_align_ratio", 'w')

    besthit_out = outfile + "_besthit.maf"

    aligned = []
    total = []
    ht_ratio = {}
    align_ratio = {}

    list_ref_len = []
    dict_reflen = {}

    with open(besthit_out, 'r') as f:
        for line in f:
            ref = line.strip().split()
            ref_len = int(ref[5])
            aligned_ref = int(ref[3])
            aligned.append(aligned_ref)

            query = next(f).strip().split()
            head = int(query[2])
            middle = int(query[3])
            tail = int(query[5])-int(query[2])-int(query[3])
            len_ONT = int(query[5])
            total.append(len_ONT)
            ht = int(query[5])-int(query[3])
            ratio = float(query[3])/float(query[5])
            rel_len = len_ONT / float(ref_len)
            #if rel_len > 1:
                #print("errorrrrrrrrr>>>>", ref[1], query[1],len_ONT, ref_len)
            if middle in align_ratio:
                align_ratio[middle].append(ratio)
            else:
                align_ratio[middle] = [ratio]
            if ht != 0:
                r = float(head) / ht
                if ht in ht_ratio:
                    ht_ratio[ht].append(r)
                else:
                    ht_ratio[ht] = [r]

            list_ref_len.append(ref_len)
            if ref_len not in dict_reflen:
                dict_reflen[ref_len] = [rel_len]
            else:
                dict_reflen[ref_len].append(rel_len)

            '''
            len_ONT = int(query[5])
            if ref_len not in dict_reflen:
                if len_ONT < ref_len:
                    dict_reflen[ref_len] = [len_ONT]
                    list_ref_len.append(ref_len)
            else:
                if len_ONT < ref_len:
                    dict_reflen[ref_len].append(len_ONT)
                    list_ref_len.append(ref_len)
            '''

    max_length = max(total)
    num_aligned = len(aligned)

    # ecdf of length of aligned regions (2d length distribution) editted this part - Approach 2 relative length of ONT total over total length of the reference transcriptome it aligned to.
    rel_len_bins = flex_bins(num_of_bins, dict_reflen, len(list_ref_len))

    rel_len_cum = dict.fromkeys(rel_len_bins.keys(), [])
    for key, value in rel_len_bins.items():
        hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ratio * 0.001)
        rel_len_cum[key] = cdf

    out5.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(rel_len_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out5.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t")
        for key in sorted(rel_len_cum.keys()):
            out5.write(str(rel_len_cum[key][i]) + "\t")
        out5.write("\n")

    out5.close()




    # ecdf of length of aligned regions (2d length distribution) editted this part - Approach 1 considering only length of ONT (total)
    hist_reflen_bins = flex_bins(num_of_bins, dict_reflen, len(list_ref_len))
    dict_value = dict.fromkeys(hist_reflen_bins.keys(), [])

    for key, value in hist_reflen_bins.items():
        dict_value[key] = value

    for item in sorted(dict_value.keys()):
        out1.write(str(item[0]) + "-" + str(item[1]) + "\t" + "\t".join([str(x) for x in dict_value[item]]) + "\n")

    # ecdf of length of aligned reads
    hist_reads, bin_edges = numpy.histogram(total, bins=numpy.arange(0, max_length + 50, 50), density=True)
    cdf = numpy.cumsum(hist_reads * 50)
    out2.write("bin\t0-" + str(max_length) + '\n')
    for i in xrange(len(cdf)):
        out2.write(str(bin_edges[i]) + '-' + str(bin_edges[i + 1]) + "\t" + str(cdf[i]) + '\n')






    # ecdf of head/total ratio
    # there needs to be at least one bin

    ht_ratio_bins = flex_bins(num_of_bins, ht_ratio, len(total))

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

    # ecdf of align ratio
    align_ratio_bins = flex_bins(num_of_bins, align_ratio, len(total))

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

    out1.close()
    out2.close()
    out3.close()
    out4.close()

    return num_aligned