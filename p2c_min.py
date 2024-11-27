#!/usr/bin/python
# -*- coding: utf-8 -*-

import errno
import os
import sys
import shutil
import subprocess
import datetime
import argparse
import numpy as np

class FileHandler :
    """ """
    def __init__(self, paf, outdir, bed, templates, log) :
        self.paf = paf
        self.bed = bed
        self.outdir = outdir
        self.log = log

        if os.path.isfile(paf) :
            self.paf = os.path.abspath(paf)
        else :
            raise Exception("ERROR: Alignment .paf file does not exist!")
            
        try:
            self.bed = os.path.abspath(bed)
        except:
            self.bed = None

        if templates == None :
            self.templates = os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates")
        else :
            self.templates = os.path.abspath(templates)

        if not os.path.isdir(self.templates) :
            raise Exception("ERROR: Invalid template file directory")

        if not os.path.isfile(os.path.join(self.templates, "circos.conf")) or not os.path.isfile(os.path.join(self.templates, "ideogram.conf")) or not os.path.isfile(os.path.join(self.templates, "ticks.conf")) :
            raise Exception("ERROR: Missing configuration files in template directory")

        self.make_outdir()

    def make_outdir(self) :
        try :
            os.makedirs(self.outdir)
            self.outdir = os.path.abspath(self.outdir)
        except OSError as e :
            if e.errno != errno.EEXIST :
                raise Exception("ERROR: Cannot create output directory!")
            else :
                self.outdir = os.path.abspath(self.outdir)
                print("WARNING: Directory already exists")

def parse_args() :
        parser = argparse.ArgumentParser(description='Circos plot all alignments against each target in the .paf file.')
        parser.add_argument('PAF',nargs=1,type=str,help="An alignment .paf formatted")
        parser.add_argument('Output',nargs=1,type=str,help="A directory path")
        parser.add_argument('--templates','-t',nargs=1,required=False,type=str,default=[None],help="A directory path containing the template files")
        parser.add_argument('--min-length','-m',nargs=1,type=int,default=[1000],required=False,help="Minimum alignment length. Default: 1000.")
        parser.add_argument('--min-query-length','-ql',nargs=1,type=int,default=[0],required=False,help="Minimum query length. Default: 0.")
        parser.add_argument('--min-target-length','-tl',nargs=1,type=int,default=[0],required=False,help="Minimum target length. Default: 0.")
        parser.add_argument('--bed','-b',nargs=1,type=str,default=[None],required=False,help="A .bed file with regions to highlight.")
        parser.add_argument('--log','-l',nargs=1,type=str,default=[None],required=False,help="A .txt file for the log output of circos.")
        parser.add_argument('-s','--self', action="store_true", help="Reference and targets are the same file. Default: %(default)s")
        return parser.parse_args()

def read_and_filter_alignments(paf, min_query_length, min_ref_length, is_self=False) :
    
    target_contigs = {}
    if not is_self:
        query_contigs = {}
    
    alignments = []
    circos_ids = {}
    contigs_with_valid_alignments = []
    
    f = open(paf)
    cid = 0
    for line in f :
        s=line.strip().split("\t")
        q = s[0]
        ql = int(s[1])
        t = s[5]
        tl = int(s[6])
        
        if ql < min_query_length:
            continue
        if tl < min_ref_length:
            continue
        
        if q not in circos_ids.keys():
            circos_ids[q] = cid
            cid += 1
        if t not in circos_ids.keys():
            circos_ids[t] = cid
            cid += 1
            
        target_in_target_contigs = t in target_contigs
        query_in_target_contigs = q in target_contigs
        if not target_in_target_contigs :
            target_contigs[circos_ids[t]] = tl
        
        if is_self and not query_in_target_contigs:
            target_contigs[circos_ids[q]] = ql
        
        if not is_self:
            if q not in query_contigs:
                query_contigs[circos_ids[q]] = ql
            
        qs = int(s[2])
        qe = int(s[3])
        ts = int(s[7])
        te = int(s[8])
        al = int(s[10]) # aligned length

        if tl <= min_ref_length :
            continue

        if ql <= min_query_length :
            continue

        alignments.append((circos_ids[q], circos_ids[t], qs, qe, ts, te))
        contigs_with_valid_alignments.append(circos_ids[q])
        contigs_with_valid_alignments.append(circos_ids[t])
        
    f.close()
    
    # In case is self
    if is_self:
        query_contigs = {}
        
    return alignments, contigs_with_valid_alignments, circos_ids, query_contigs, target_contigs

def make_karyotype(out, contigs_with_valid_alignments, reverse_circos_ids, contigs_lengths, query_contigs) :
    
    # Add each target line and its queries aligned to the karyotype
    o_t = open(os.path.join(out, "karyotype.txt"), "w")
    for cid in set(contigs_with_valid_alignments) :
        color = "chr2" if cid in query_contigs else "chr1"
        # Add target line to karyotype
        name = reverse_circos_ids[cid]
        length = contigs_lengths[cid]
        o_t.write(f"chr - c_{cid} {name} 0 {length} {color}\n")
    o_t.close()

def make_links(outdir, alignments) :
    o_t = open(os.path.join(outdir, "links.txt"), "w")
    for qid, tid, qs, qe, ts, te in alignments:
        o_t.write(f"c_{qid} {qs} {qe} c_{tid} {ts} {te}\n")
    o_t.close()

def make_highlights(o_t, bed, ):
    f = open(bed, "r")
    for line in f:
        s = line.strip().split("\t")
        
        try:
            cid = circos_ids[s[0]]
        except:
            print(f"WARNING: {s[0]} not found in corresponding .PAF (maybe too small and filtered out?)")
            continue
            
        start = int(s[1])
        end = int(s[2])

        if len(s) > 3:
            color = s[3]
        else:
            color = "orange"
        
        o_t.write(f"{cid} {start} {end} fill_color={color}\n")
        
    f.close()

def plot_paf(alignments, circos_ids, contigs_with_valid_alignments, query_contigs, target_contigs, outdir, templates, bed, log):
    
    # Get cid:name from name:cid
    reverse_circos_ids = {cid:name for name, cid in circos_ids.items()}
    
    # Get chromosome lengths
    contigs_lengths = target_contigs
    contigs_lengths.update(query_contigs)
    
    cwd = os.getcwd()
    # 1. Copy all template files into the subdir
    print("Copying files...")
    for name in ["circos.conf", "ticks.conf", "ideogram.conf"] :
        shutil.copy(os.path.join(templates, name), outdir)

    # 2. Make karyotype
    print("Building karyotype.txt...")
    make_karyotype(outdir, contigs_with_valid_alignments, reverse_circos_ids, contigs_lengths, list(query_contigs.keys()))
    ideograms = ";".join(f"c_{cid}" for cid in contigs_with_valid_alignments)
    replace_pattern(os.path.join(outdir, "circos.conf"), "<CHROMOSOMES_WILL_GO_HERE>", ideograms)

    # 3. Make links
    print("Building links.txt...")
    make_links(outdir, alignments)
    
    # 4. Highlights
    o_t = open(os.path.join(outdir, "highlight.txt"), "w")
    if bed is not None :
        make_highlights(o_t, bed, circos_ids)
    o_t.close()
    
    # For compatibility with v1
    o_t = open(os.path.join(outdir, "samtools_coverage.txt"), "w")
    o_t.close()
    o_t = open(os.path.join(outdir, "circos_coverage.txt"), "w")
    o_t.close()
    
    # 4. Run circos
    os.chdir(os.path.join(outdir))
    cmd = "circos"
    run(cmd, log=log)
    os.chdir(cwd)

def run(cmd, log=None) :
    print("Running circos in {}".format(os.getcwd()))
    if log is not None:
        flog = open(log, "w")
        proc = subprocess.Popen(cmd, stdout=flog, stderr=flog)
        proc.communicate()
        flog.close()
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        proc.communicate()
    print("Done!")

def replace_pattern(input_file, pattern, replacement) :
    # Read in the file
    with open(input_file, 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace(pattern, replacement)
    # Write the file out again
    with open(input_file, 'w') as file:
        file.write(filedata)

def main() :
    args = parse_args()
    templates = args.templates[0]
    paf = args.PAF[0]
    bed = args.bed[0]
    outdir = args.Output[0]
    min_len = args.min_length[0]
    min_query_length = args.min_query_length[0]
    min_ref_length = args.min_target_length[0]
    is_self = args.self
    log = args.log[0]
    
    # Fix file paths
    iFH = FileHandler(paf, outdir, templates, bed, log) # File handler
    
    # Get useful dictionaries
    alignments, valid_contigs, circ_ids, query_contigs, target_contigs = read_and_filter_alignments(iFH.paf, min_query_length, min_ref_length, is_self=is_self)
    
    # Prepares files and circos plot
    plot_paf(alignments, circ_ids, valid_contigs, query_contigs, target_contigs, iFH.outdir, iFH.templates, bed=iFH.bed, log=iFH.log)

if __name__ == '__main__':
        main()
