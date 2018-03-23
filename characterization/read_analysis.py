#!/usr/bin/env python

"""
Created on March 2018

@author: Saber HafezQorani

This script generates read profiles from Oxford Nanopore 2D transcriptome reads.

"""


from __future__ import print_function
from __future__ import with_statement
from subprocess import call
from time import strftime
try:
    from six.moves import xrange
except ImportError:
    pass
import sys
import os
import getopt
import HTSeq
import pysam
import numpy
import head_align_tail_dist as align
import get_besthit
import besthit_to_histogram as error_model



# Usage information
def usage():
    usage_message = "./read_analysis.py <options>\n" \
                    "<options>: \n" \
                    "-h : print usage message\n" \
                    "-i : training ONT real reads, must be fasta files\n" \
                    "-r : reference genome of the training reads\n" \
                    "-t : reference transcriptome of the training reads\n" \
                    "-a : reference GTF/GFF3 annotation files\n" \
                    "-sg : User can provide their own genome alignment file, with sam extension\n" \
                    "-st : User can provide their own transcriptome alignment file, with sam extension\n" \
                    "-b : number of bins (for development), default = 20\n" \
                    "-o : The prefix of output file, default = 'training'\n"

    sys.stderr.write(usage_message)


def main(argv):

    # Parse input and output files
    infile = ''
    outfile = 'training'
    ref = ''
    maf_file = ''
    model_fit = True
    num_bins = 20
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:m:b:", ["infile=", "ref=", "outfile=", "no_model_fit"])
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-r", "--ref_g"):
            ref_g = arg
        elif opt in ("-t", "--ref_t"):
            ref_t = arg
        elif opt in ("-a", "--annotation"):
            annot = arg
        elif opt == "-sg":
            alignment_genome = arg
        elif opt == "-st":
            alignment_transcriptome = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt == "--no_model_fit":
            model_fit = False
        elif opt == "-b":
            num_bins = max(int(arg), 1)
        else:
            usage()
            sys.exit(1)

    if infile == '' or ref_g == '' or ref_t == '' or annot == '':
        print("Please specify the training reads and its reference genome and transcriptome along with the annotation GTF/GFF3 files!")
        usage()
        sys.exit(1)

    # READ PRE-PROCESS AND UNALIGNED READS ANALYSIS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")

    # Read the annotation GTF/GFF3 file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
    # If gtf provided, convert to GFF3 (gt gtf_to_gff3)
    annot_filename = annot.split(".")[0]
    annot_type = annot.split(".")[1]
    if annot_type.upper() == "GTF":
        call("gt gtf_to_gff3 -tidy -o " + annot_filename + ".gff3" + annot, shell=True)


    # Next, add intron info into gff3:
    call("gt gff3 -tidy -retainids -checkids -addintrons -o " + annot_filename + "_addedintron.gff3 " + annot_filename + ".gff3", shell=True)
    features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    dict_intron_info = {}
    annot_gff3_withintron = annot_filename + "_addedintron.gff3"
    gff_file = HTSeq.GFF_Reader(annot_gff3_withintron, end_included=True)
    for feature in gff_file:
        if "Parent" in feature.attr:
            info = feature.attr["Parent"].split(':')
            if info[0] == "transcript":
                feature_id = info[1]
                if feature_id not in dict_intron_info:
                    dict_intron_info[feature_id] = []
        if feature.type == "intron":
            features[feature.iv] += feature_id
            dict_intron_info[feature_id].append((feature.iv.start, feature.iv.end, feature.iv.length))

    #Read the length of reference transcripts from the reference transcriptome
    dict_ref_len = {}
    with open (ref_t) as f:
        for line in f:
            if line.startswith(">"):
                ref_id = line.split()[0][1:]
                dict_ref_len[ref_id] = 0
            else:
                dict_ref_len[ref_id] += len(line.strip())


    # If either of the alignment files are not provided (sam files):
    if alignment_genome == '' or alignment_transcriptome == '':
        alignment_genome = outfile + "_genome_alnm.sam"
        alignment_transcriptome = outfile + "_transcriptome_alnm.sam"

        # Alignment to reference genome
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference genome\n")
        call("minimap2 -ax splice " + ref_g + " " + infile + " > " + alignment_genome, shell=True)
        # Alignment to reference transcriptome
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
        call("minimap2 --cs -ax splice " + ref_t + " " + infile + " > " + alignment_transcriptome, shell=True)


    # Read the genome alignments to memory:
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading the genome alignments\n")
    dict_genome_alignment = {}
    SAM_or_BAM_Reader = HTSeq.SAM_Reader
    read_alignment_genomne = SAM_or_BAM_Reader(alignment_genome)
    for read in read_alignment_genomne:
        qname = read.read.name
        if qname not in dict_genome_alignment:
            dict_genome_alignment[qname] = read #Read the primary alignment only (ignores the secondary and supp alignments)

    # Read the transcriptome alignments to memory:
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading the transcriptome alignments\n")
    dict_trx_alignment = {}
    SAM_or_BAM_Reader = HTSeq.SAM_Reader
    read_alignment_trx = SAM_or_BAM_Reader(alignment_transcriptome)
    for read in read_alignment_trx:
        qname = read.read.name
        if qname not in dict_trx_alignment:
            dict_trx_alignment[qname] = read #Read the primary alignment only (ignores the secondary and supp alignments)

    # Aligned reads analysis
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
    num_aligned, error_dict = align.head_align_tail(outfile, num_bins, dict_trx_alignment, dict_ref_len)
    #num_aligned = head_align_tail(outfile, num_bins, dict_trx_alignment, dict_ref_len)

    # Un-aligned reads analysis
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Un-aligned reads analysis\n")

    unaligned_length = []
    num_unaligned = 0
    for qname in dict_trx_alignment:
        t_read = dict_trx_alignment[qname]
        if not t_read.aligned:
            num_unaligned += 1
            unaligned_length.append(len(t_read.read.seq))

    # Length distribution of unaligned reads
    out1 = open(outfile + "_unaligned_length_ecdf", 'w')
    if num_unaligned != 0:
        max_length = max(unaligned_length)
        hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length + 50, 50),
                                                          density=True)
        cdf = numpy.cumsum(hist_unaligned * 50)
        out1.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
        out1.write("bin\t0-" + str(max_length) + '\n')
        for i in xrange(len(cdf)):
            out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    else:
        out1.write("Aligned / Unaligned ratio:\t100%\n")

    out1.close()
    del unaligned_length


    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile, error_dict)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        #path = sys.argv[0].split("/")
        #r_path = '/'.join(path[:-1]) + '/' + "model_fitting.R"
        r_path = '/projects/btl/shafez/trans_nanosim/trans_nanosim-master/src/model_fitting.R'
        if os.path.isfile(r_path):
            call("R CMD BATCH '--args prefix=\"" + outfile + "\"' " + r_path, shell=True)
        else:
            sys.stderr.write("Could not find 'model_fitting.R' in ../src/\n" +
                  "Make sure you copied the whole source files from Github.")

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main(sys.argv[1:])

