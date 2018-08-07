#!/usr/bin/env python

from __future__ import with_statement

stranded = "no" # think about it. Should I input this info for cDNA ONT data or not?

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def intron_retention(gff_file, talnm_file, galnm_file, ref_t):
    gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)
    #t_alnm = "/projects/btl/shafez/analysis/testing_alignment_tools/minimap2/SRR5286960_trx_alm.sam"
    #g_alnm = "/projects/btl/shafez/analysis/testing_alignment_tools/minimap2/SRR5286960_genome_alm.sam"
    #ref_t = "/projects/btl_scratch/datasets/ref_transcriptome/mus_musculus/Mus_musculus.GRCm38.cdna.all.fa"

    #read the reference transcriptome to get their length.
    dict_ref_len = {}
    with open(ref_t) as f:
        for line in f:
            if line.startswith(">"):
                ref_id = line.split()[0][1:]
                dict_ref_len[ref_id] = 0
            else:
                dict_ref_len[ref_id] += len(line.strip())


    #read intron information from GFF file
    features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    c = 0
    dict_intron_info = {}
    for feature in gff_features:
        c += 1
        if c % 100000 == 0:
            print(str(c) + " features proccessed!")
        if "Parent" in feature.attr:
            info = feature.attr["Parent"].split(':')
            if info[0] == "transcript":
                feature_id = info[1]
                if feature_id not in dict_intron_info:
                    dict_intron_info[feature_id] = []
        if feature.type == "intron":
            # feature_id_2 = feature.name.split(':')[1] #feature_id_2 is same as feature_id above if feature is intron, I was just checking and testing it. then removed this line.
            features[feature.iv] += feature_id
            dict_intron_info[feature_id].append((feature.iv.start, feature.iv.end, feature.iv.length))


    dict_g_alnm = {}
    sam_reader = HTSeq.SAM_Reader
    g_alignments = sam_reader(galnm_file)
    for alnm in g_alignments:
        qname = alnm.read.name
        if alnm.aligned and not alnm.not_primary_alignment and not alnm.supplementary:
            dict_g_alnm[qname] = alnm
        if alnm.supplementary and qname in dict_g_alnm:
            del dict_g_alnm[qname]  # delete chimeric reads


    dict_t_alnm = {}
    sam_reader = HTSeq.SAM_Reader
    t_alignments = sam_reader(talnm_file)
    for alnm in t_alignments:
        qname = alnm.read.name
        if alnm.aligned and not alnm.not_primary_alignment and not alnm.supplementary:
            dict_t_alnm[qname] = alnm
        if alnm.supplementary and qname in dict_t_alnm:
            del dict_t_alnm[qname]  # delete chimeric reads


    #count the length of Intron retention events
    list_IR = []
    list_not_IR = []
    dict_ir_len = {}
    for qname in dict_g_alnm:
        galnm = dict_g_alnm[qname]
        if qname in dict_t_alnm:
            primary_trx = dict_t_alnm[qname].iv.chrom.split(".")[0]
            if stranded != "reverse":
                iv_seq = (co.ref_iv for co in galnm.cigar if (co.type in ('M', '=', 'X', 'D') and co.size > 0))
                #iv_seq = (co.ref_iv for co in galnm.cigar if co.type in ('M', 'D') and co.size > 0) #tested. test the above cases too to make sure about it.
            else:
                iv_seq = (invert_strand(co.ref_iv) for co in galnm.cigar if (co.type in ('M', '=', 'X', 'D') and co.size > 0))
            length_IR_total = 0
            length_IR_full = 0
            length_IR_intron = 0
            length_IR_read = []
            try:
                fs = set()
                fs2_temp = set()
                for iv in iv_seq:
                    for iv2, fs2 in features[iv].steps():
                        if fs2.intersection(set([primary_trx])):
                            length_IR_total += iv2.length
                            length_IR_intron += iv2.length
                            fs2_temp = fs2.intersection(set([primary_trx]))
                        else:
                            if length_IR_intron != 0:
                                length_IR_read.append(length_IR_intron)
                                for intron in dict_intron_info[primary_trx]:
                                    if length_IR_intron == intron[2]:
                                        length_IR_full += length_IR_intron
                                        fs = fs.union(fs2_temp)
                                length_IR_intron = 0
                                fs2_temp = set()
                if fs is None or len(fs) == 0:
                    list_not_IR.append(qname)
                elif len(fs) >= 1:
                    list_IR.append(qname)
            except UnknownChrom:
                list_not_IR.append(qname)

            dict_ir_len[qname] = [length_IR_total, length_IR_full]

    return dict_ir_len