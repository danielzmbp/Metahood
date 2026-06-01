#!/usr/bin/env python3
import sys
import argparse
from collections import defaultdict

def parse_ko_cutoffs(cutoff_file):
    cut_ko = {}
    with open(cutoff_file) as f:
        f.readline()
        for line in f:
            sline = line.strip().split('\t')
            if sline[1] == '-':
                cut_ko[sline[0]] = [0.0, "NA"]
            else:
                cut_ko[sline[0]] = [float(sline[1]), sline[2]]
    return cut_ko

def from_hmmsearch_results(hmmsearch_output_path, cutoff_file,
                           evalue_cutoff=1e-05, bitscore_cutoff=0,
                           percent_aln_query_cutoff=0, percent_aln_reference_cutoff=0,
                           acc=True):
    specific_cutoffs = parse_ko_cutoffs(cutoff_file) if cutoff_file else {}

    for line in open(hmmsearch_output_path):
        if line.startswith('#'): continue

        sline = line.strip().split()
        if len(sline) < 22: continue
        
        seqname = sline[0]
        tlen = sline[2]
        ko_hmm = sline[3]
        accession = sline[4]
        qlen = sline[5]
        score = sline[7]
        i_evalue = sline[12]
        dom_score = sline[13]
        seq_from = sline[17]
        seq_to = sline[18]

        seq_list = [int(seq_from), int(seq_to)]
        perc_seq_aln = (max(seq_list)-min(seq_list))/float(tlen)
        perc_hmm_aln = (max(seq_list)-min(seq_list))/float(qlen)

        if acc:
            accession = ko_hmm

        if specific_cutoffs and ko_hmm in specific_cutoffs:
            if specific_cutoffs[ko_hmm][1] == 'full':
                if float(score) < specific_cutoffs[ko_hmm][0]:
                    continue
            elif specific_cutoffs[ko_hmm][1] == 'domain':
                if float(dom_score) < specific_cutoffs[ko_hmm][0]:
                    continue

        if(float(i_evalue) <= evalue_cutoff and
           float(score) >= bitscore_cutoff and
           perc_seq_aln >= percent_aln_query_cutoff and
           perc_hmm_aln >= percent_aln_reference_cutoff):
            yield seqname, [accession], float(i_evalue), (min(seq_list), max(seq_list))

def intervals_overlap(first, second):
    return max(first[0], second[0]) < min(first[1], second[1])

def generate_orf_annotation(file, cutoff_file):
    orf_to_annotations = defaultdict(list)
    for orf, [accession], evalue, r in from_hmmsearch_results(file, cutoff_file):
        orf_to_annotations[orf].append([accession, evalue, r])

    orf_to_ann = {orf:[val[0][0]] for orf, val in orf_to_annotations.items() if len(val)==1}
    issues = set(orf_to_annotations.keys()) - set(orf_to_ann.keys())
    
    for orf in issues:
        accessions = sorted(orf_to_annotations[orf], key=lambda x:x[1])
        real_acc = []
        to_del = set()
        for index, (acc, evalue, range_best) in enumerate(accessions):
            if acc in to_del:
                continue 
            real_acc.append(acc)
            for (acc_cand, eval_cand, range_cand) in accessions[index:]:
                if intervals_overlap(range_best, range_cand):
                    to_del.add(acc_cand)
        orf_to_ann[orf] = real_acc
    return orf_to_ann

def parse_results(hmmsearch_file, cutoff_file, out_file):
    orf_to_ann = generate_orf_annotation(hmmsearch_file, cutoff_file)
    with open(out_file, "w") as handle:
        for orf, ann in orf_to_ann.items():
            handle.write(f"{orf}\t{chr(9).join(ann)}\n")

def combine_results(input_files, cutoff_file, out_file):
    ko_to_def = {}
    with open(cutoff_file) as f:
        for line in f:
            if not line.startswith("KO"):
                parts = line.rstrip().split("\t")
                ko_to_def[parts[0]] = parts[-1]
                
    with open(out_file, "w") as handle:
        handle.write("orf\tKO\tKO definition\n")
        for file in input_files:
            for line in open(file):
                sline = line.rstrip().split("\t")
                orf = sline[0]
                for KO in sline[1:]:
                    handle.write(f"{orf}\t{KO}\t{ko_to_def.get(KO, '')}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    p_parse = subparsers.add_parser("parse")
    p_parse.add_argument("--hmmsearch", required=True)
    p_parse.add_argument("--cutoffs", required=True)
    p_parse.add_argument("--out", required=True)
    
    p_combine = subparsers.add_parser("combine")
    p_combine.add_argument("--inputs", nargs="+", required=True)
    p_combine.add_argument("--cutoffs", required=True)
    p_combine.add_argument("--out", required=True)
    
    args = parser.parse_args()
    
    if args.command == "parse":
        parse_results(args.hmmsearch, args.cutoffs, args.out)
    elif args.command == "combine":
        combine_results(args.inputs, args.cutoffs, args.out)
