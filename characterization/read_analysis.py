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
    ref_g = ''
    ref_t = ''
    annot = ''
    alignment_genome = ''
    alignment_transcriptome = ''
    aligner = ''
    alnm_file = ''
    model_fit = True
    num_bins = 20
    try:
        opts, args = getopt.getopt(argv, "hi:rg:rt:a:ga:ta:o:b:", ["infile=", "ref=", "outfile=", "no_model_fit"])
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-rg", "--ref_g"):
            ref_g = arg
        elif opt in ("-rt", "--ref_t"):
            ref_t = arg
        elif opt in ("-annot", "--annotation"):
            annot = arg
        elif opt in ("-a", "--aligner"):
            aligner = arg
        elif opt == "-ga":
            g_alnm = arg
        elif opt == "-ta":
            t_alnm = arg
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

    if aligner != '' and (alignment_genome != '' or alignment_transcriptome != ''):
        print("Please specify either an alignment file (-ga and -ta ) OR an aligner to use for alignment (-a )")
        usage()
        sys.exit(1)

    # READ PRE-PROCESS AND UNALIGNED READS ANALYSIS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")

    # Read pre-process
    in_fasta = outfile + ".fasta"
    if in_fasta == infile:
        in_fasta = outfile + "_processed.fasta"
    out_fasta = open(in_fasta, 'w')
    dic_reads = {}
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '>':
                name = '-'.join(line.strip()[1:].split())
                dic_reads[name] = ""
            else:
                dic_reads[name] += line.strip()
    for k, v in dic_reads.items():
        out_fasta.write('>' + k + '\n' + v + '\n')
    out_fasta.close()

    del dic_reads

    # Read the annotation GTF/GFF3 file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
    # If gtf provided, convert to GFF3 (gt gtf_to_gff3)
    annot_filename, annot_file_extension = os.path.splitext(annot)
    annot_file_extension = annot_file_extension[1:]
    if annot_file_extension.upper() == "GTF":
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

    #If both alignment files are provided:
    if g_alnm != "" and t_alnm != "":
        g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
        t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
        g_alnm_ext = g_alnm_ext [1:]
        t_alnm_ext = t_alnm_ext[1:]
        if g_alnm_ext != t_alnm_ext:
            print("Please provide both alignments in a same format: sam OR maf\n")
            usage()
            sys.exit(1)
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing the alignment files: " + g_alnm_ext + "\n")
        if t_alnm_ext == "maf":
            outmaf_g = outfile + "_genome_alnm.maf"
            outmaf_t = outfile + "_transcriptome_alnm.maf"
            if outmaf_g == g_alnm:
                outmaf_g = outfile + "_genome_alnm_processed.maf"
            if outmaf_t == t_alnm:
                outmaf_t = outfile + "_transcriptome_alnm_processed.maf"

            call("grep '^s ' " + g_alnm + " > " + outmaf_g, shell=True)
            call("grep '^s ' " + t_alnm + " > " + outmaf_t, shell=True)

            unaligned_length = list(get_besthit_maf.besthit_and_unaligned(in_fasta, outmaf_t, outfile))

        elif t_alnm_ext == "sam":

            unaligned_length = list(get_primary_sam.primary_and_unaligned(t_alnm, outfile))

    else:
        if aligner == "minimap2" or aligner == "":  # Align with minimap2 by default
            g_alnm_ext = "sam"
            t_alnm_ext = "sam"
            outsam_g = outfile + "_genome_alnm.sam"
            outsam_t = outfile + "_transcriptome_alnm.sam"
            # Alignment to reference genome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference genome\n")
            call("minimap2 --cs --MD -ax map-ont " + ref_g + " " + in_fasta + " > " + outsam_g, shell=True)
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
            call("minimap2 --cs -ax splice " + ref_t + " " + in_fasta + " > " + outsam_t, shell=True)

            unaligned_length = list(get_primary_sam.primary_and_unaligned(outsam_t, outfile))

        elif aligner == "LAST":
            g_alnm_ext = "maf"
            t_alnm_ext = "maf"
            outmaf_g = outfile + "_genome_alnm.maf"
            outmaf_t = outfile + "_transcriptome_alnm.maf"
            # Alignment to reference genome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference genome\n")
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 ref_genome " + in_fasta + " | grep '^s ' > " + outmaf_g, shell=True)
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference transcriptome\n")
            call("lastdb ref_transcriptome " + ref_t, shell=True)
            call("lastal -a 1 ref_transcriptome " + in_fasta + " | grep '^s ' > " + outmaf_t, shell=True)

            unaligned_length = list(get_besthit_maf.besthit_and_unaligned(in_fasta, outmaf_t, outfile))

        else:
            print("Please specify an acceptable aligner (minimap2 or LAST)\n")
            usage()
            sys.exit(1)


    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reads length distribution analysis\n")
    # Aligned reads length distribution analysis
    count_aligned = align.head_align_tail(outfile, num_bins, dict_ref_len, t_alnm_ext)

    # Unaligned reads length distribution analysis
    out1 = open(outfile + "_unaligned_length_ecdf", 'w')
    count_unaligned = len(unaligned_length)
    if count_unaligned != 0:
        max_length = max(unaligned_length)
        hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length + 50, 50),
                                                          density=True)
        cdf = numpy.cumsum(hist_unaligned * 50)
        out1.write("Aligned / Unaligned ratio:" + "\t" + str(count_aligned * 1.0 / count_unaligned) + '\n')
        out1.write("bin\t0-" + str(max_length) + '\n')
        for i in xrange(len(cdf)):
            out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    else:
        out1.write("Aligned / Unaligned ratio:\t100%\n")
    out1.close()

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile, t_alnm_ext)

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

