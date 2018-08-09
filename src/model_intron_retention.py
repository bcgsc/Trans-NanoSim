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

def intron_retention(outfile, gff_file, talnm_file, galnm_file, ref_t):

    outfile = open(outfile + "_intron_retention_probs", 'w')
    gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)

    #read the reference transcriptome to get their length.
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading reference transcript lengths\n")
    dict_ref_len = {}
    with open(ref_t) as f:
        for line in f:
            if line.startswith(">"):
                ref_id = line.split()[0][1:]
                dict_ref_len[ref_id] = 0
            else:
                dict_ref_len[ref_id] += len(line.strip())


    #read intron information from GFF file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading intron coordinates from GFF file\n")
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

    #read primary genome alignment for each read
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read primary genome alignment for each read\n")
    dict_g_alnm = {}
    sam_reader = HTSeq.SAM_Reader
    g_alignments = sam_reader(galnm_file)
    for alnm in g_alignments:
        qname = alnm.read.name
        if alnm.aligned and not alnm.not_primary_alignment and not alnm.supplementary:
            dict_g_alnm[qname] = alnm
        if alnm.supplementary and qname in dict_g_alnm:
            del dict_g_alnm[qname]  # delete chimeric reads

    #read primary transcriptome alignment for each read
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read primary transcriptome alignment for each read\n")
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
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Calculating probabilites for each intron retention event\n")
    dict_first_intron_state = {False: 0, True: 0}
    dict_states = {(False, False): 0, (False, True): 0, (True, False): 0, (True, True): 0}
    for qname in dict_g_alnm:
        galnm = dict_g_alnm[qname]
        if qname in dict_t_alnm:
            talnm = dict_t_alnm[qname]
            primary_trx = talnm.iv.chrom.split(".")[0]
            if stranded != "reverse":
                iv_seq = (co.ref_iv for co in galnm.cigar if (co.type in ('M', '=', 'X', 'D') and co.size > 0))
                #iv_seq = (co.ref_iv for co in galnm.cigar if co.type in ('M', 'D') and co.size > 0) #tested. test the above cases too to make sure about it.
            else:
                iv_seq = (invert_strand(co.ref_iv) for co in galnm.cigar if (co.type in ('M', '=', 'X', 'D') and co.size > 0))

            list_IR_positions = []
            pos = []
            ir_info = False
            try:
                length_IR = 0
                for iv in iv_seq:
                    for iv2, fs2 in features[iv].steps():
                        if fs2.intersection(set([primary_trx])):
                            length_IR += iv2.length
                            pos.append(iv2.start)
                            pos.append(iv2.end)
                        else:
                            if length_IR != 0:
                                list_IR_positions.append(min(pos))
                                list_IR_positions.append(max(pos))
                                for intron in dict_intron_info[primary_trx]:
                                    if length_IR == intron[2]:
                                        ir_info = True
                                length_IR = 0
                                pos = []
            except UnknownChrom:
                ir_info = False
                pass

            if ir_info == False:
                if primary_trx in dict_intron_info:
                    if len(dict_intron_info[primary_trx]) >= 1: #if there is a intron
                        dict_first_intron_state[False] += 1
                    for i in range(1, len(dict_intron_info[primary_trx])):
                        dict_states[(False, False)] += 1
            else:
                # Now, go over all introns and check with the IR events
                # First we need to determine the state of first intron:
                first_intron = dict_intron_info[primary_trx][0]
                first_intron_spos = first_intron[0]
                first_intron_epos = first_intron[1]
                flag = False
                for IR_pos in list_IR_positions:
                    if first_intron_spos <= IR_pos <= first_intron_epos:
                        flag = True
                        break
                if flag == True:
                    dict_first_intron_state[True] += 1
                    previous_state = True
                else:
                    dict_first_intron_state[False] += 1
                    previous_state = False

                # Then we will go over other introns:
                current_state = False
                for i in range (1, len(dict_intron_info[primary_trx])):
                    intron = dict_intron_info[primary_trx][i]
                    current_state = False
                    intron_spos = intron[0]
                    intron_epos = intron[1]
                    for IR_pos in list_IR_positions:
                        if intron_spos <= IR_pos <= intron_epos:
                            current_state = True
                            break
                    #print(intron_spos, intron_epos, previous_state, current_state)
                    dict_states[(previous_state, current_state)] += 1
                    previous_state = current_state

    num_first_introns = dict_first_intron_state[True] + dict_first_intron_state[False]
    outfile.write("first intron IR event probability: " + str(dict_first_intron_state[True] / float(num_first_introns)) + "\n")
    outfile.write("previous_state\tcurrent_state\tprob\n")
    outfile.write("False\tFalse\t" + str(dict_states[(False, False)]) + "\n")
    outfile.write("False\tTrue\t" + str(dict_states[(False, True)]) + "\n")
    outfile.write("True\tFalse\t" + str(dict_states[(True, False)]) + "\n")
    outfile.write("True\tTrue\t" + str(dict_states[(True, True)]) + "\n")