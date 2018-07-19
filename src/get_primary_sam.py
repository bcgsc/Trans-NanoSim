#!/usr/bin/env python

from __future__ import with_statement
import HTSeq


def primary_and_unaligned(g_alnm, t_alnm, outfile):

    outsam_g = outfile + "_genome_alnm.sam"
    outsam_t = outfile + "_transcriptome_alnm.sam"
    if outsam_g == g_alnm:
        outsam_g = outfile + "_genome_alnm_processed.sam"
    if outsam_t == t_alnm:
        outsam_t = outfile + "_transcriptome_alnm_processed.sam"

    out1 = open (outsam_g, 'w')
    out2 = open(outsam_t, 'w')
    out3 = open(outfile + "_transcriptome_alnm_primary.sam", 'w')

    sam_reader = HTSeq.SAM_Reader
    g_alignments = sam_reader(g_alnm)
    for alnm in g_alignments:
        out1.write(alnm.original_sam_line)

    unaligned_list = []
    t_alignments = sam_reader(t_alnm)
    for alnm in t_alignments:
        out2.write(alnm.original_sam_line)
        if alnm.aligned:
            if not alnm.not_primary_alignment and not alnm.supplementary:
                out3.write(alnm.original_sam_line)
        else:
            unaligned_list.append(len(alnm.read.seq))

    out1.close()
    out2.close()
    out3.close()
    return unaligned_list
