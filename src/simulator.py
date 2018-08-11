# !/usr/bin/env python
"""
Created on Jan 11, 2018

@author: Saber HafezQorani

This script generates simulated Oxford Nanopore 2D transcriptome reads.

"""

from __future__ import print_function
from __future__ import with_statement
from subprocess import call
import sys
import glob
import getopt
import argparse
import random
import numpy
import time
from time import sleep
import operator
import re
from time import strftime
import mixed_model as mm

try:
    from six.moves import xrange
except ImportError:
    pass
try:
    import numpy as np
except ImportError:
    sys.exit("""You need numpy!
                install it from http://www.numpy.org/""")


PYTHON_VERSION = sys.version_info
VERSION = "1.0.0"
PRORAM = "Trans-NanoSim"
AUTHOR = "Saber HafezQorani (UBC & BCGSC)"
CONTACT = "shafezqorani@bcgsc.ca"

BASES = ['A', 'T', 'C', 'G']


# Usage information
def usage():
    usage_message = "./simulator.py [command] <options>\n" \
                    "<options>: \n" \
                    "-h : print usage message\n" \
                    "-r : reference genome in fasta file, specify path and file name, REQUIRED\n" \
                    "-c : The prefix of training set profiles, same as the output prefix in read_analysis.py, default = training\n" \
                    "-o : The prefix of output file, default = 'simulated'\n" \
                    "-n : Number of generated reads, default = 20,000 reads\n" \
                    "--max_len : Maximum read length, default = Inf\n" \
                    "--min_len : Minimum read length, default = 50\n" \
                    "--perfect: Output perfect reads, no mutations, default = False\n" \
                    "--KmerBias: prohibits homopolymers with length >= n bases in output reads, default = 6\n"

    sys.stderr.write(usage_message)


def list_to_range(input_list, min_l):
    l = [min_l]
    l.extend(input_list)
    output_list = []
    for i in range (0, len(l) - 1):
        r = (l[i], l[i+1])
        output_list.append(r)
    return output_list


def read_ecdf_test(list_cdf_range, list_value_range):
    lanes = 1
    ecdf_dict = {}
    l_cdf = [0.0] * lanes
    l_exp = [0.0] * lanes
    for j in range (len(list_cdf_range)):
        cdf = [float(x) for x in [list_cdf_range[j]]]
        value_range = [float(x) for x in list_value_range[j]]
        for i in xrange(lanes):
            if cdf[i] == l_cdf[i]:
                continue
            else:
                if l_cdf[i] != 0:
                    ecdf_dict[(l_cdf[i], cdf[i])] = (l_exp[i], value_range[1])
                else:
                    ecdf_dict[(l_cdf[i], cdf[i])] = (max(l_exp[i], value_range[1] - 10 * (value_range[1] - value_range[0])), value_range[1])
                l_exp[i] = value_range[1]
                l_cdf[i] = cdf[i]

        last_key = sorted(ecdf_dict.keys())[-1]
        last_value = ecdf_dict[last_key]
        ecdf_dict[last_key] = (last_value[0], value_range[1])

    return ecdf_dict


def read_ecdf(profile):
    # We need to count the number of zeros. If it's over 10 zeros, l_len/l_ratio need to be changed to higher.
    # Because it's almost impossible that the ratio is much lower than the lowest heuristic value.
    header = profile.readline()
    header_info = header.strip().split()
    ecdf_dict = {}
    lanes = len(header_info[1:])

    for i in header_info[1:]:
        boundaries = i.split('-')
        ecdf_dict[(int(boundaries[0])), int(boundaries[1])] = {}

    ecdf_key = sorted(ecdf_dict.keys())
    l_prob = [0.0] * lanes
    l_ratio = [0.0] * lanes

    for line in profile:
        new = line.strip().split('\t')
        ratio = [float(x) for x in new[0].split('-')]
        prob = [float(x) for x in new[1:]]
        for i in xrange(lanes):
            if prob[i] == l_prob[i]:
                continue
            else:
                if l_prob[i] != 0:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] = (l_ratio[i], ratio[1])
                else:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] \
                        = (max(l_ratio[i], ratio[1] - 10 * (ratio[1] - ratio[0])), ratio[1])
                l_ratio[i] = ratio[1]
                l_prob[i] = prob[i]

    for i in xrange(0, len(ecdf_key)):
        last_key = sorted(ecdf_dict[ecdf_key[i]].keys())[-1]
        last_value = ecdf_dict[ecdf_key[i]][last_key]
        ecdf_dict[ecdf_key[i]][last_key] = (last_value[0], ratio[1])

    return ecdf_dict


def make_cdf(dict_exp, dict_len):
    sum_exp = 0
    list_value = []
    dict_tpm = {}
    for item in dict_exp:
        if item in dict_len:
            sum_exp += dict_exp[item]
    for item in dict_exp:
        if item in dict_len:
            value = dict_exp[item] / float(sum_exp)
            list_value.append((item, value))


    sorted_value_list = sorted(list_value, key=lambda x: x[1])
    sorted_only_values = [x[1] for x in sorted_value_list]
    list_cdf = numpy.cumsum(sorted_only_values)
    ranged_cdf_list = list_to_range(list_cdf, 0)

    ecdf_dict = {}
    for i in xrange(len(ranged_cdf_list)):
        cdf_range = ranged_cdf_list[i]
        #ecdf_dict[cdf_range] = (sorted_value_list[i][0], dict_len[sorted_value_list[i][0]])
        ecdf_dict[cdf_range] = dict_len[sorted_value_list[i][0]]

    return ecdf_dict


def select_ref_transcript(input_dict):
    length = 0
    while True:
        p = random.random()
        for key, val in input_dict.items():
            if key[0] <= p < key[1]:
                length = val
                break
        if length != 0:
            break
    return length


def get_length_ratio(dict_ratio, len_ref):

    for item in dict_ratio.keys():
        if item[0] <= len_ref < item[1]:
            p = random.random()
            for k_r, v_r in dict_ratio[item].items():
                if k_r[0] <= p < k_r[1]:
                    r = (p - k_r[0]) / (k_r[1] - k_r[0]) * (v_r[1] - v_r[0]) + v_r[0]
                    ref_aligned = int(round(len_ref * r))
                    break
            return item, ref_aligned


def get_length_2d(len_dict, len_ref, max_l, min_l):

    for item in len_dict.keys():
        if item[0] < len_ref <= item[1]:
            key = item

    value = [int(x) for x in len_dict[key]]
    sorted_value = sorted(value)
    if max(value) > len_ref:
        value_range = sorted_value[:next(x[0] for x in enumerate(sorted_value) if x[1] > len_ref)]
    else:
        value_range = sorted_value

    if len(value_range) > 0:
        ratio_value = [float(1) / len(value_range) for x in value_range]
        list_cdf = numpy.cumsum(ratio_value)

        ranged_value_list = list_to_range(value_range, min_l)
        ranged_cdf_list = list_to_range(list_cdf, 0)

        ecdf_dict = {}
        for i in xrange(len(ranged_value_list)):
            ont_length_range = ranged_value_list[i]
            cdf_range = ranged_cdf_list[i]
            ecdf_dict[cdf_range] = ont_length_range

        middle_ref = 0
        while True:
            p = random.random()
            for cdf, length in ecdf_dict.items():
                if cdf[0] <= p < cdf[1]:
                    middle_ref = ((abs(cdf[1] - p) * length[0]) + (abs(cdf[0] - p) * length[1])) / abs(cdf[1] - cdf[0])
                    middle_ref = int(round(middle_ref))
                    #local_diff = abs(p - cdf[0]) / abs(cdf[1] - cdf[0])
                    #select_index = value_range.index(length[0]) * (1 - local_diff) + value_range.index(length[1]) * local_diff
                    #middle_ref = value_range[int(round(select_index))]
                    break
            if middle_ref != 0:
                break

        return key, middle_ref
    else:
        return key, 0


def get_length(len_dict, num, max_l, min_l):
    length_list = []
    for i in xrange(num):
        middle_ref = 0
        key = tuple(len_dict.keys())[0]
        while middle_ref <= min_l or middle_ref > max_l:
            p = random.random()
            for k_p, v_p in len_dict[key].items():
                if k_p[0] <= p < k_p[1]:
                    middle_ref = int(round((p - k_p[0]) / (k_p[1] - k_p[0]) * (v_p[1] - v_p[0]) + v_p[0]))
                    break
        length_list.append(middle_ref)

    return length_list


def read_profile(number, model_prefix, per, max_l, min_l):
    global unaligned_length, number_aligned, aligned_dict, reftotal_dict
    global first_match_hist, align_ratio, ht_dict, error_model_profile
    global error_markov_model, match_markov_model
    global dict_head, dict_tail
    global dict_ir_states

    # Read model profile for match, mismatch, insertion and deletions
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read error profile\n")
    sys.stdout.flush()
    error_model_profile = {}
    model_profile = model_prefix + "_model_profile"
    with open(model_profile, 'r') as mod_profile:
        mod_profile.readline()
        for line in mod_profile:
            new_line = line.strip().split("\t")
            if "mismatch" in line:
                error_model_profile["mis"] = [float(x) for x in new_line[1:]]
            elif "insertion" in line:
                error_model_profile["ins"] = [float(x) for x in new_line[1:]]
            else:
                error_model_profile["del"] = [float(x) for x in new_line[1:]]

    error_markov_model = {}
    with open(model_prefix + "_error_markov_model", "r") as error_markov:
        error_markov.readline()
        for line in error_markov:
            info = line.strip().split()
            k = info[0]
            error_markov_model[k] = {}
            error_markov_model[k][(0, float(info[1]))] = "mis"
            error_markov_model[k][(float(info[1]), float(info[1]) + float(info[2]))] = "ins"
            error_markov_model[k][(1 - float(info[3]), 1)] = "del"

    with open(model_prefix + "_first_match.hist", 'r') as fm_profile:
        first_match_hist = read_ecdf(fm_profile)

    with open(model_prefix + "_match_markov_model", 'r') as mm_profile:
        match_markov_model = read_ecdf(mm_profile)

    # Read length of unaligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read ECDF of unaligned reads\n")
    sys.stdout.flush()
    unaligned_length = [] #remove this line. no need for it.
    with open(model_prefix + "_unaligned_length_ecdf", 'r') as u_profile:
        new = u_profile.readline().strip()
        rate = new.split('\t')[1]
        # if parameter perfect is used, all reads should be aligned, number_aligned equals total number of reads.
        if per or rate == "100%":
            number_aligned = number
        else:
            number_aligned = int(round(number * float(rate) / (float(rate) + 1)))
        number_unaligned = number - number_aligned
        unaligned_dict = read_ecdf(u_profile)

    unaligned_length = get_length(unaligned_dict, number_unaligned, max_l, min_l)
    unaligned_dict.clear()

    # Read profile of aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read ECDF of aligned reads\n")
    sys.stdout.flush()

    # Read align ratio profile
    with open(model_prefix + "_align_ratio", 'r') as a_profile:
        align_ratio = read_ecdf(a_profile)

    # Read head/unaligned region ratio
    with open(model_prefix + "_ht_ratio", 'r') as ht_profile:
        ht_dict = read_ecdf(ht_profile)

    # Read length of aligned reads
    # If "perfect" is chosen, just use the total length ecdf profile, else use the relative length of ONT 2 reference
    if per:
        length_profile = model_prefix + "_aligned_reads_ecdf"
    else:
        length_profile = model_prefix + "_read_rellen_ecdf"

    with open(length_profile, 'r') as align_profile:
        aligned_dict = read_ecdf(align_profile)


    with open(model_prefix + "_reflen_total_ecdf", "r") as reftotal_profile:
        reftotal_dict = read_ecdf(reftotal_profile)

    dict_head = {}
    with open(model_prefix + "_head_sequences", "r") as f_head:
        for line in f_head:
            hseq = line.strip("\n")
            if len(hseq) not in dict_head:
                dict_head[len(hseq)] = [hseq]
            else:
                dict_head[len(hseq)].append(hseq)

    dict_tail = {}
    with open(model_prefix + "_tail_sequences", "r") as f_tail:
        for line in f_tail:
            tseq = line.strip("\n")
            if len(tseq) not in dict_tail:
                dict_tail[len(tseq)] = [tseq]
            else:
                dict_tail[len(tseq)].append(tseq)

    dict_ir_states = {}
    with open(model_prefix + "_intron_retention_probs", "r") as  ir_probs:
        fline = ir_probs.readline()
        f_intron_prob = float(fline.split("\t")[1])
        dict_ir_states["first"] = f_intron_prob
        sec_line = ir_probs.readline()
        for line in ir_probs:
            line_parts = line.split("\t")
            dict_ir_states[(line_parts[0], line_parts[1])] = float(line_parts[-1])


def get_ht_sequence(dict_ht, length):
    return random.choice(dict_ht[length])


def collapse_homo(seq, k):
    read = re.sub("A" * k + "+", "A" * (k - 1), seq)
    read = re.sub("C" * k + "+", "C" * (k - 1), read)
    read = re.sub("T" * k + "+", "T" * (k - 1), read)
    read = re.sub("G" * k + "+", "G" * (k - 1), read)

    return read


# Taken from https://github.com/lh3/readfq
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def simulation(ref, out, per, kmer_bias, max_l, min_l, exp):
    global unaligned_length, number_aligned, aligned_dict, reftotal_dict
    global genome_len, seq_dict, seq_len, dict_exp, ecdf_dict_ref
    global first_match_hist, align_ratio, ht_dict, match_markov_model
    global error_markov_model, error_model_profile
    global dict_head, dict_tail
    global dict_ir_states

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference transcriptome (length and expression)\n")
    sys.stdout.flush()

    seq_dict = {}
    seq_len = {}
    dict_exp = {}

    # Read in the reference transcriptome
    with open(ref, 'r') as infile:
        for seqN, seqS, seqQ in readfq(infile):
            info = re.split(r'[_\s]\s*', seqN)
            transcript_id = "-".join(info)
            seq_dict[transcript_id] = seqS
            #sys.stdout.write(".")
            #sys.stdout.flush()
    #sys.stdout.write("\n")
    #sys.stdout.flush()

    for key in seq_dict.keys():
        seq_len[key] = len(seq_dict[key])

    #Read the expression profile.
    #I need to improve it so that it reads any input and generates errors
    with open (exp, 'r') as exp_file:
        header = exp_file.readline()
        for line in exp_file:
            parts = line.split("\t")
            transcript_id = parts[0]
            tpm = float(parts[2])
            if transcript_id.startswith("ENS") and tpm > 0:
                dict_exp[transcript_id] = tpm

    #create the ecdf dict considering the expression profiles
    ecdf_dict_ref_exp = make_cdf(dict_exp, seq_len)

    # Start simulation
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of reads\n")
    sys.stdout.flush()
    out_reads = open(out + "_reads.fasta", 'w')
    out_error = open(out + "_error_profile", 'w')
    out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")

    # Simulate unaligned reads
    #sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of unaligned reads\n")
    #sys.stdout.flush()
    num_unaligned_length = len(unaligned_length)
    for i in xrange(num_unaligned_length):
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(i) + "\r")
        sys.stdout.flush()
        sleep(0.02)

        unaligned = unaligned_length[i]
        unaligned, error_dict = unaligned_error_list(unaligned, error_model_profile)
        new_read, new_read_name = extract_read(unaligned)
        new_read_name = new_read_name + "_unaligned_" + str(i)
        # Change lowercase to uppercase and replace N with any base
        new_read = case_convert(new_read)
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias, False)

        # Reverse complement half of the reads
        p = random.random()
        if p < 0.5:
            read_mutated = reverse_complement(read_mutated)
            new_read_name += "_R"
        else:
            new_read_name += "_F"
        out_reads.write(">" + new_read_name + "_0_" + str(unaligned) + "_0" + '\n')
        out_reads.write(read_mutated + "\n")

    del unaligned_length

    # Simulate aligned reads
    #sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    #sys.stdout.flush()

    if per:
        read_total_len = get_length(aligned_dict, number_aligned, max_l, min_l)
        del aligned_dict

        for i in xrange(number_aligned):
            new_read, new_read_name = extract_read(read_total_len[i])
            new_read_name = new_read_name + "_perfect_" + str(i)

            # Reverse complement half of the reads
            p = random.random()
            if p < 0.5:
                new_read = reverse_complement(new_read)
                new_read_name += "_R"
            else:
                new_read_name += "_F"

            out_reads.write(">" + new_read_name + "_0_" + str(read_total_len[i]) + "_0" + '\n')

            # Change lowercase to uppercase and replace N with any base
            new_read = case_convert(new_read)

            out_reads.write(new_read + "\n")
        out_reads.close()
        out_error.close()
        return

    i = 0
    while i < number_aligned:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(i + num_unaligned_length) + "\r")
        sys.stdout.flush()
        sleep(0.02)
        while True:
            #pick a reference to simulate a read out of it.
            ref_len_total = max(50, select_ref_transcript(ecdf_dict_ref_exp))
            #ref_len_total = get_length(reftotal_dict, 1, max_l, min_l)[0] # just for testing
            #get the align region ratio out of this ref len total
            key_range, ref_len_aligned = get_length_ratio(aligned_dict, ref_len_total)
            if ref_len_aligned != 0:
                middle_read, middle_ref, error_dict = error_list(ref_len_aligned, match_markov_model, first_match_hist, error_model_profile, error_markov_model)
                if middle_ref < key_range[1]:
                    break

        for k_align in sorted(align_ratio.keys()):
            if k_align[0] <= middle_read < k_align[1]:
                break

        p = random.random()
        for k_r, v_r in align_ratio[k_align].items():
            if k_r[0] <= p < k_r[1]:
                a_ratio = (p - k_r[0])/(k_r[1] - k_r[0]) * (v_r[1] - v_r[0]) + v_r[0]
                total = int(round(middle_read / a_ratio))
                remainder = total - int(round(middle_read))
                break

        if total > max_l:
            continue

        if remainder == 0:
            head = 0
            tail = 0
        else:
            for k_ht in sorted(ht_dict.keys()):
                if k_ht[0] <= remainder < k_ht[1]:
                    p = random.random()
                    for k_h, v_h in ht_dict[k_ht].items():
                        if k_h[0] <= p < k_h[1]:
                            ratio = (p - k_h[0])/(k_h[1] - k_h[0]) * (v_h[1] - v_h[0]) + v_h[0]
                            head = int(round(remainder * ratio))
                            tail = remainder - head
                            break
                    break

            # if remainder is larger than any empirical value, then randomly divide it into head and tail
            try:
                head
            except NameError:
                p = random.random()
                head = int(round(remainder * p))
                tail = remainder - head

        try:
            # Extract error introduced middle region from reference transcript
            new_read, new_read_name = extract_read_withrange(middle_ref, key_range)
            new_read_name = new_read_name + "_aligned_" + str(i + num_unaligned_length)
        except:
            print(i, ref_len, ont_total, middle, middle_ref, head, tail, key_range)
            #the problem is that middle_ref length is larger than anything in the key range. Actually we should take care of this thing at error_list function !
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": ERROR >> " + str(i))
            sys.exit(1)

        # Mutate read
        new_read = case_convert(new_read)
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)

        # Reverse complement half of the reads
        p = random.random()
        if p < 0.5:
            read_mutated = reverse_complement(read_mutated)
            new_read_name += "_R"
        else:
            new_read_name += "_F"

        #Add head and tail sequences
        if head != 0:
            nearest_head = min(dict_head, key=lambda x: abs(x - head))
            head_sequence = get_ht_sequence(dict_head, nearest_head)
        else:
            head_sequence = ""

        if tail != 0:
            nearest_tail = min(dict_tail, key=lambda x: abs(x - tail))
            tail_sequence = get_ht_sequence(dict_tail, nearest_tail)
        else:
            tail_sequence = ""

        read_mutated = head_sequence + read_mutated
        read_mutated += tail_sequence

        if kmer_bias:
            read_mutated = collapse_homo(read_mutated, kmer_bias)

        out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle_ref) + "_" +
                        str(tail) + '\n')
        out_reads.write(read_mutated + '\n')

        i += 1


    out_reads.close()
    out_error.close()

    align_ratio.clear()
    ht_dict.clear()


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def extract_read(length):
    global seq_dict, seq_len

    while True:
        new_read = ""
        key = random.choice(seq_len.keys())
        if length < seq_len[key]:
            ref_pos = random.randint(0, seq_len[key] - length)
            new_read = seq_dict[key][ref_pos: ref_pos + length]
            new_read_name = key + "_" + str(ref_pos)
            break
    return new_read, new_read_name


def extract_read_withrange(length, key_range):
    global seq_dict, seq_len

    seq_len_temp = {}
    for key, value in seq_len.items():
        if key_range[0] < value <= key_range[1]:
            seq_len_temp[key] = value

    while True:
        new_read = ""
        key = random.choice(seq_len_temp.keys())
        if length < seq_len_temp[key]:
            ref_pos = random.randint(0, seq_len_temp[key] - length)
            new_read = seq_dict[key][ref_pos: ref_pos + length]
            new_read_name = key + "_" + str(ref_pos)
            break
    return new_read, new_read_name


def unaligned_error_list(length, error_p):
    e_dict = {}
    error_rate = {(0, 0.4): "match", (0.4, 0.7): "mis", (0.7, 0.85): "ins", (0.85, 1): "del"}
    pos = 0
    last_is_ins = False
    while pos < length:
        p = random.random()
        for k_error in error_rate.keys():
            if k_error[0] <= p < k_error[1]:
                error_type = error_rate[k_error]
                break

        if error_type == "match":
            step = 1

        elif error_type == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
            e_dict[pos] = ["mis", step]

        elif error_type == "ins":
            step = mm.wei_geom(error_p["ins"][0], error_p["ins"][1], error_p["ins"][2], error_p["ins"][3])
            if last_is_ins:
                e_dict[pos + 0.1][1] += step
            else:
                e_dict[pos + 0.1] = ["ins", step]
                last_is_ins = True

        else:
            step = mm.wei_geom(error_p["del"][0], error_p["del"][1], error_p["del"][2], error_p["del"][3])
            e_dict[pos] = ["del", step]

        if error_type != "ins":
            pos += step
            last_is_ins = False

        if pos > length:
            length = pos

    return length, e_dict


def error_list(m_ref, m_model, m_ht_list, error_p, trans_p):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    middle_ref = m_ref
    prev_error = "start"

    # The first match come from m_ht_list
    p = random.random()
    k1 = list(m_ht_list.keys())[0]
    for k2, v2 in m_ht_list[k1].items():
        if k2[0] < p <= k2[1]:
            prev_match = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
            if prev_match < 2:
                prev_match = 2
    pos += prev_match

    # Select an error, then the step size, and then a match and so on so forth.
    while pos < middle_ref:
        # pick the error based on Markov chain
        p = random.random()
        for k in trans_p[prev_error].keys():
            if k[0] <= p < k[1]:
                error = trans_p[prev_error][k]
                break

        if error == "mis":
            step = mm.pois_geom(error_p[error][0], error_p[error][2], error_p[error][3])
        elif error == "ins":
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new += step
        else:
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new -= step

        if error != "ins":
            e_dict[pos] = [error, step]
            pos += step
            if pos >= middle_ref:
                l_new += pos - middle_ref
                middle_ref = pos
        else:
            e_dict[pos - 0.5] = [error, step]

        prev_error = error

        # Randomly select a match length
        for k1 in m_model.keys():
            if k1[0] <= prev_match < k1[1]:
                break
        p = random.random()
        for k2, v2 in m_model[k1].items():
            if k2[0] < p <= k2[1]:
                step = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
                break
        # there are no two 0 base matches together
        if prev_match == 0 and step == 0:
            step = 1

        prev_match = step
        if pos + prev_match > middle_ref:
            l_new += pos + prev_match - middle_ref
            middle_ref = pos + prev_match

        pos += prev_match
        if prev_match == 0:
            prev_error += "0"
    return l_new, middle_ref, e_dict


def mutate_read(read, read_name, error_log, e_dict, k, aligned=True):
    search_pattern = "A" * k + "+|" + "T" * k + "+|" + "C" * k + "+|" + "G" * k
    for key in sorted(e_dict.keys(), reverse=True):
        val = e_dict[key]
        key = int(round(key))

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    tmp_bases = list(BASES)
                    tmp_bases.remove(read[key + i])
                    new_base = random.choice(tmp_bases)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key + val[1]: key + val[1] + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key + val[1]:]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    new_base = random.choice(BASES)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key: key + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key:]

        read = new_read

        if aligned and val[0] != "match":
            error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            "\t" + ref_base + "\t" + new_bases + "\n")

    # If choose to have kmer bias, then need to compress homopolymers to 5-mer
    if k:
        read = collapse_homo(read, k)

    return read


def case_convert(seq):
    base_code = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
                 'N': ['A', 'T', 'C', 'G'], 'X': ['A', 'T', 'C', 'G']}

    up_string = seq.upper()
    up_list = list(up_string)
    for i in xrange(len(up_list)):
        if up_list[i] in base_code:
            up_list[i] = random.choice(base_code[up_list[i]])
    out_seq = ''.join(up_list)

    return out_seq


def main():

    parser = argparse.ArgumentParser(
        description='Given the read profiles from characterization step, ' \
                    'simulate transcriptome ONT reads and output error profiles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    reads = ''
    exp = ''
    model_prefix = ''
    out = ''
    number = ''
    max_readlength = ''
    min_readlength = ''
    kmer_bias = ''
    perfect = False

    parser.add_argument('-r', '--read', help='Input reads to quanity expression profiles', type = str)
    parser.add_argument('-f', '--ref', help='Input reference transcriptome', type = str, required= True)
    parser.add_argument('-e', '--expression', help='Expression profile in the specified format provided with documentation', type = str)
    parser.add_argument('-c', '--model_prefix', help='Address for profiles created in characterization step (model_prefix)', type = str, default= "training")
    parser.add_argument('-o', '--output', help='Output address for simulated reads', type = str, default= "simulated")
    parser.add_argument('-n', '--number', help='Number of reads to be simulated', type = int, default = 20000)
    parser.add_argument('-i', '--insertion_rate', help='Insertion rate (optional)', type = float, default= 1)
    parser.add_argument('-d', '--deletion_rate', help='Deletion rate (optional)', type = float, default= 1)
    parser.add_argument('-m', '--mismatch_rate', help='Mismatch rate (optional)', type = float, default= 1)
    parser.add_argument('-max', '--max_len', help='The maximum length for simulated reads', type=int, default= float("inf"))
    parser.add_argument('-min', '--min_len', help='The minimum length for simulated reads', type=int, default= 50)
    parser.add_argument('-k', '--KmerBias', help='Determine whether to considert Kmer Bias or not', type = int, default= 0)
    parser.add_argument('--perfect', help='Ignore profiles and simulated perfect reads', action='store_true')

    args = parser.parse_args()

    ref = args.ref
    model_prefix = args.model_prefix
    out = args.output
    number = args.number

    if args.read:
        reads = args.read
    if args.expression:
        exp = args.expression
    else:
        sys.stdout.write('The expression profile of the simulated reads will be determined by input reads used in characterization step. \n')
        exp = model_prefix
    if args.insertion_rate:
        ins_rate = float(args.insertion_rate)
    if args.deletion_rate:
        del_rate = float(args.deletion_rate)
    if args.mismatch_rate:
        mis_rate = float(args.mismatch_rate)
    if args.max_len:
        max_readlength = args.max_len
    if args.min_len:
        min_readlength = args.min_len
    if args.perfect:
        perfect = True
    if args.KmerBias:
        kmer_bias = args.KmerBias


    print ("Running the simulation step with following arguments: \n")
    print ("reads: ", reads)
    print ("ref: ", ref)
    print ("expression: ", exp)
    print ("model_prefix: ", model_prefix)
    print ("output: ", out)
    print ("number: ", number)
    print ("insertion_rate: ", ins_rate)
    print ("deletion: ", del_rate)
    print ("mismatch_rate: ", mis_rate)
    print ("max_readlength: ", max_readlength)
    print ("min_readlength: ", min_readlength)
    print("kmer_bias: ", kmer_bias)
    print ("perfect: ", perfect)


    # Generate log file
    sys.stdout.log = open(out + ".log", 'w')

    # Record the command typed to log file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
    sys.stdout.log.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
    sys.stdout.flush()

    if max_readlength < min_readlength:
        print("maximum read length must be longer than minimum read length!")
        sys.exit(1)

    # Read in reference transcriptome and generate simulated reads
    read_profile(number, model_prefix, perfect, max_readlength, min_readlength)

    simulation(ref, out, perfect, kmer_bias, max_readlength, min_readlength, exp)

    call("find . -name \*.pyc -delete", shell=True)
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
    sys.stdout.log.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
    sys.stdout.log.close()


if __name__ == "__main__":
    main()
