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
        elif opt == "-galnm":
            alignment_genome = arg
        elif opt == "-talnm":
            alignment_transcriptome = arg
        elif opt == "-alnmf":
            alnm_ftype = arg.upper()
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
        call("minimap2 --MD -ax map-ont " + ref_g + " " + infile + " > " + alignment_genome, shell=True)
        # Alignment to reference transcriptome
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
        call("minimap2 --cs -ax splice " + ref_t + " " + infile + " > " + alignment_transcriptome, shell=True)
        alnm_ftype = "SAM"

    else:

        if alnm_ftype == "MAF":
            out_maf_genome = outfile + "_genome_alnm.maf"
            out_maf_trx = outfile + "_transcriptome_alnm.maf"
            if out_maf_genome == alignment_genome:
                out_maf_genome = outfile + "_genome_alnm_processed.maf"
            if out_maf_trx == alignment_transcriptome:
                out_maf_trx = outfile + "_transcriptome_alnm_processed.maf"

            call("grep '^s ' " + alignment_genome + " > " + out_maf_genome, shell=True)
            call("grep '^s ' " + alignment_transcriptome + " > " + out_maf_trx, shell=True)
        elif alnm_ftype == "SAM":
            # do any preproccessing on SAM alnm files if necessary
            pass


    # Read the genome alignments to memory:
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading the genome alignments\n")

    if alnm_ftype == "SAM":
        dict_genome_alignment = {}
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
        read_alignment_genomne = SAM_or_BAM_Reader(alignment_genome)
        for alnm in read_alignment_genomne:
            qname = alnm.read.name
            if qname not in dict_genome_alignment: #primary alnms and unaligned reads (ignores the secondary and supp alignments)
            #test with flags too.
                dict_genome_alignment[qname] = alnm

        # Read the transcriptome alignments to memory:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading the transcriptome alignments\n")
        dict_trx_alignment = {}
        SAM_or_BAM_Reader = HTSeq.SAM_Reader
        read_alignment_trx = SAM_or_BAM_Reader(alignment_transcriptome)
        for alnm in read_alignment_trx:
            qname = alnm.read.name
            if qname not in dict_trx_alignment: #primary alnms and unaligned reads (ignores the secondary and supp alignments)
                # test with flags too.
                dict_trx_alignment[qname] = alnm

    elif alnm_ftype == "MAF":
        # read the MAF alnment files to the dictionaries
        unaligned_length, dict_trx_alignment = get_besthit.besthit_and_unaligned(in_fasta, out_maf_trx, outfile)
        num_unaligned = len(unaligned_length)

    else:
        print("Please specify an acceptable alignment files: MAF or SAM")
        usage()
        sys.exit(1)

    # Reads length distribution analysis
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reads length distribution analysis\n")
    num_aligned, num_unaligned, unaligned_length = align.head_align_tail(outfile, num_bins, dict_trx_alignment, dict_ref_len, alnm_ftype)
    #num_aligned = head_align_tail(outfile, num_bins, dict_trx_alignment, dict_ref_len)


    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile, dict_trx_alignment)

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

