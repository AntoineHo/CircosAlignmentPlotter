#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import errno
import sys
import shutil
import subprocess
import datetime
import argparse

import pandas as pd
import numpy as np
from Bio import SeqIO # Need BIOPYTHON SEQ/IO

class Target :
    def __init__(self, name, tag, start, end, color, length) :
        """ """
        self.name = name
        self.tag = tag
        self.start = start
        self.end = end
        self.color = color
        self.length = length

    def __str__(self) :
        return "Target(name={}, tag={}, start={}, end={}, color={}, length={})".format(self.name,self.tag,self.start,self.end,self.color,self.length)

class Coords :
    def __init__(self, name_1, name_2, start_1, start_2, end_1, end_2) :
        """ """
        self.n1 = name_1
        self.s1 = start_1
        self.e1 = end_1
        self.n2 = name_2
        self.s2 = start_2
        self.e2 = end_2

    def __eq__(self, co) :
        return (self.n1 == co.n1 or self.n1 == co.n2) and (self.s1 == co.s1 or self.s1 == co.s2) and (self.e1 == co.e1 or self.e1 == co.e2)

class Alignment :
    """ """
    def __init__(self, query, Q_start, Q_end, T_start, T_end, length) :
        self.query = query
        self.Q_start = Q_start
        self.Q_end = Q_end
        self.T_start = T_start
        self.T_end = T_end
        self.length = length

class FH :
    """ """
    def __init__(self, query, targets, reference, paf, outdir, bed, scov, ccov, snps, templates) :
        self.paf = None
        self.query = None
        self.targets = None
        self.reference = None
        self.bed = None
        self.scov = None
        self.ccov = None
        self.snps = None
        self.outdir = outdir

        if os.path.isfile(paf) :
            self.paf = os.path.abspath(paf)
        else :
            raise Exception("ERROR: Alignment .paf file does not exist!")

        if os.path.isfile(query) :
            self.query = os.path.abspath(query)
        else :
            raise Exception("ERROR: Query .fa file does not exist!")

        if os.path.isfile(targets) :
            self.targets = os.path.abspath(targets)
        else :
            raise Exception("ERROR: Targets .bed file does not exist!")

        if os.path.isfile(reference) :
            self.reference = os.path.abspath(reference)
        else :
            raise Exception("ERROR: Reference .fasta file does not exist!")

        if bed != None :
            if os.path.isfile(bed) :
                self.bed = os.path.abspath(bed)
            else :
                raise Exception("ERROR: Highlights .bed file does not exist!")

        if scov != None :
            if os.path.isfile(scov) :
                self.scov = os.path.abspath(scov)
            else :
                raise Exception("ERROR: Coverage .cov file does not exist! {}".format(scov))

        if ccov != None :
            if os.path.isfile(ccov) :
                self.ccov = os.path.abspath(ccov)
            else :
                raise Exception("ERROR: Coverage .cov file does not exist! {}".format(ccov))

        if snps != None :
            if os.path.isfile(snps) :
                self.snps = os.path.abspath(snps)
            else :
                raise Exception("ERROR: SNPs .vcf file does not exist!")

        if templates == None :
            self.templates = os.path.join(os.path.dirname(os.path.realpath(__file__)), "templates")
        else :
            self.templates = os.path.abspath(templates)

        if not os.path.isdir(self.templates) :
            raise Exception("ERROR: Invalid template file directory")

        if not os.path.isfile(os.path.join(self.templates, "circos.conf")) or not os.path.isfile(os.path.join(self.templates, "ideogram.conf")) or not os.path.isfile(os.path.join(self.templates, "ticks.conf")) :
            raise Exception("ERROR: Missing configuration files in template directory")

        self.makeOutdir()

    def makeOutdir(self) :
        try :
            os.makedirs(self.outdir)
            self.outdir = os.path.abspath(self.outdir)
        except OSError as e :
            if e.errno != errno.EEXIST :
                raise Exception("ERROR: Cannot create output directory!")
            else :
                self.outdir = os.path.abspath(self.outdir)
                print("WARNING: Directory already exists")

    def __str__(self) :
        return "Query: {}\nReference: {}\nTargets: {}\nPAF: {}\nOutdir: {}\nTemplates: {}".format(self.query, self.reference, self.targets, self.paf, self.outdir, self.templates)

def parseArgs() :
        parser = argparse.ArgumentParser(description='Circos plot all alignments against each target in the .paf file.')
        parser.add_argument('PAF',nargs=1,type=str,help="An alignment .paf formatted")
        parser.add_argument('Query',nargs=1,type=str,help="A fasta file containing the target sequences")
        parser.add_argument('Reference',nargs=1,type=str,help="A fasta file containing the query sequences")
        parser.add_argument('BED',nargs=1,type=str,help="A .bed file containing the target contigs regions to plot")
        parser.add_argument('Output',nargs=1,type=str,help="A directory path")
        parser.add_argument('--templates','-t',nargs=1,required=False,type=str,default=[None],help="A directory path containing the template files")
        parser.add_argument('--min-length','-m',nargs=1,type=int,default=[1000],required=False,help="Minimum length. Default: 1000.")
        parser.add_argument('--bed','-b',nargs=1,type=str,default=[None],required=False,help="A .bed file with regions to highlight.")
        parser.add_argument('--scov','-sc',nargs=1,type=str,default=[None],required=False,help="A file formatted like the output of samtools depth.")
        parser.add_argument('--ccov','-cc',nargs=1,type=str,default=[None],required=False,help="A file formatted for circos.")
        parser.add_argument('--snps','-s',nargs=1,type=str,default=[None],required=False,help="A .vcf file of the query and/or target.")
        return parser.parse_args()

def readTargets(targets, correspondance, contig_lengths) :
    targets_to_plot = []
    for line in open(targets, "r") :
        s = line.strip().split("\t")
        #print("s: " + "_".join(x for x in s)) # DEBUG
        if len(s) == 0 or s[0] == "" :
            continue
        ctg = s[0]
        start = int(s[1])
        end = int(s[2])
        color = "orange"
        try :
            color = s[3]
        except :
            pass
        targets_to_plot.append(Target(ctg, correspondance[ctg], start, end, color, contig_lengths[ctg]))

    return targets_to_plot

def readPAF(paf, targets_to_plot) :
    dAlignments = {}
    names_of_targets_to_plot = [t.name for t in targets_to_plot]
    f = open(paf)
    for line in f :
        s=line.strip().split("\t")
        target = s[5]
        if target not in names_of_targets_to_plot :
            continue
        query = s[0]
        q_start = int(s[2])
        q_end = int(s[3])
        t_start = int(s[7])
        t_end = int(s[8])
        length = int(s[10])

        if target not in dAlignments.keys() :
            dAlignments[target] = [Alignment(query, q_start, q_end, t_start, t_end, length)]
        else :
            dAlignments[target].append(Alignment(query, q_start, q_end, t_start, t_end, length))
    f.close()
    return dAlignments

def buildCorrespondance(reference, query) :
    correspondance = {}
    contigs_lengths = {}

    n = 0
    for record in SeqIO.parse(reference, "fasta"):
        if record.id not in correspondance :
            correspondance[record.id] = "av{}".format(n)
        if record.id not in contigs_lengths :
            contigs_lengths[record.id] = len(record.seq)
        n += 1

    for record in SeqIO.parse(query, "fasta"):
        if record.id not in correspondance :
            correspondance[record.id] = "av{}".format(n)
        if record.id not in contigs_lengths :
            contigs_lengths[record.id] = len(record.seq)
        n += 1

    return correspondance, contigs_lengths

def makeKaryotype(out, dAln, targets_to_plot, target_fasta, query_fasta, correspondance, contigs_lengths, min_align_length) :
    show_up_ideograms = ""
    karyotype = []
    names = []
    o_t = open(os.path.join(out, "karyotype.txt"), "w")

    # Add each target line and its queries aligned to the karyotype
    for target in targets_to_plot :
        if target not in karyotype :
            karyotype.append(target)
            names.append(target.name)
            # Add target line to karyotype
            chrom = "chr - {} {} 0 {} chr2".format(target.tag, target.name, target.length)
            o_t.write(chrom+"\n")
        else :
            continue

        # Get the list of all long enough queries aligned to this target
        queries_aligned = []
        for aln in dAln[target.name] : # For each alignment in dict of alignments
            if aln.length < min_align_length :
                continue

            #
            #   =========    target region
            # -------------- query region
            #
            elif (aln.T_start <= target.start and aln.T_end >= target.end) :
                queries_aligned.append(aln.query)

            #
            #   =========    target region
            # ----- or ----- query region
            #
            elif (aln.T_start <= target.start and aln.T_end > target.start and aln.T_end <= target.end) or (aln.T_start >= target.start and aln.T_start <= target.end and aln.T_end >= target.end) :
                queries_aligned.append(aln.query)
                #print(aln.T_start, aln.T_end) # DEBUG

            #
            #   ========= target region
            #     ----    query region
            #
            elif (aln.T_start >= target.start and aln.T_end <= target.end) :
                queries_aligned.append(aln.query)
            else :
                #print(aln.T_start, aln.T_end) # DEBUG
                continue

        # Builds a show-up string to use in the template
        show_up_ideograms +=  ";" + target.tag + ":" + str(target.start) + "-" + str(target.end)

    # Add queries lines to karyotype
    for query in queries_aligned :
        if query not in names :
            names.append(query)
        else :
            continue # stop here
        tag = correspondance[query]
        query_length = contigs_lengths[query]
        chrom = "chr - {} {} 0 {} chr4".format(tag, query, query_length)
        o_t.write(chrom+"\n")

        # Builds a show-up string to use in the template
        show_up_ideograms += ";" + tag

    o_t.close()

    return show_up_ideograms

def makeLinks(outdir, dAln, targets_to_plot, min_align_length, correspondance) :
    o_t = open(os.path.join(outdir, "links.txt"), "w")
    # For self alignments: must remember coords and check if reverted
    pairs = []
	# All targets
    for p, target in enumerate(targets_to_plot) :
        for n, aln in enumerate(dAln[target.name]) :
            if aln.length < min_align_length :
                continue
            link1 = "link{} {} {} {}".format(str(p)+"_"+str(n), target.tag, aln.T_start, aln.T_end)
            o_t.write(link1+"\n")
            link2 = "link{} {} {} {}".format(str(p)+"_"+str(n), correspondance[aln.query], aln.Q_start, aln.Q_end)
            o_t.write(link2+"\n")
            pairs.append(Coords(target.name, aln.T_start, aln.T_end, aln.query, aln.Q_start, aln.Q_end))
    o_t.close()

def plotPAF(dAln, outdir, tpl, query_fasta, target_fasta, correspondance, min_align_length, contigs_lengths, targets_to_plot, bed=None, scov=None, ccov=None, snps=None) :
    cwd = os.getcwd()
    # 1. Copy all template files into the subdir
    print("Copying files...")
    for name in ["circos.conf", "ticks.conf", "ideogram.conf"] :
        shutil.copy(os.path.join(tpl, name), outdir)

    # 2. Make karyotype
    print("Building karyotype.txt...")
    show_up_ideograms = makeKaryotype(outdir, dAln, targets_to_plot, target_fasta, query_fasta, correspondance, contigs_lengths, min_align_length)
    replacePattern(os.path.join(outdir, "circos.conf"), "<CHROMOSOMES_WILL_GO_HERE>", show_up_ideograms[1:]) # [1:] used to remove the first ; char

    # 3. Make links
    print("Building links.txt...")
    makeLinks(outdir, dAln, targets_to_plot, min_align_length, correspondance)

    # OPTIONAL : make highlights
    print("Building highlights.txt...")
    if bed != None :
        o_t = open(os.path.join(outdir, "highlight.txt"), "w")

        for line in open(bed, "r") :
            s = line.strip().split("\t")
            ctg = s[0]
            start = int(s[1])
            end = int(s[2])
            color = "orange"
            try :
                color = s[3]
            except :
                pass
            if ctg in correspondance.keys() :
                region = "{} {} {} fill_color={}".format(correspondance[ctg], start, end, color)
                o_t.write(region + "\n")
            else :
                continue
        add_small_alignments(dAln, targets_to_plot, o_t, min_align_length, correspondance)
        o_t.close()
    else :
        o_t = open(os.path.join(outdir, "highlight.txt"), "w")
        add_small_alignments(dAln, targets_to_plot, o_t, min_align_length, correspondance)
        o_t.close()

    # OPTIONAL : add coverage information formatted by samtools
    print("Adding samtools formatted coverage...")
    if scov != None :
        df = pd.read_csv(cov, sep="\t", usecols=["REF","POS","COV"])
        dfs = []
        for target in targets_to_plot :
            dft = df.loc[df.loc[:,("REF")] == target.name]
            bins = list(np.linspace(0, contigs_lengths[target.name], 1001))
            dft["BIN"] = pd.cut(x=dft.loc[:,("POS")], bins=bins)
            gb = dft.groupby(["BIN"]).agg({"COV": 'mean', "POS":['min', 'max'], "REF":'first'})
            gb.replace({"REF":correspondance}, inplace=True)
            gb.columns = ['_'.join(col).strip() for col in gb.columns.values]
            gb.reset_index()
            #print(gb) # DEBUG
            #print("----") # DEBUG
            dfs.append(gb)
        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_csv(os.path.join(outdir, "coverage.txt"), sep=" ", columns=["REF_first", "POS_min", "POS_max", "COV_mean"], index=False, header=False)
    else :
        f = open(os.path.join(outdir, "samtools_coverage.txt"), "w")
        f.close()

    print("Adding circos formatted coverage...")
    if ccov != None :
        shutil.copy(ccov, os.path.join(outdir, "circos_coverage.txt"))
    else :
        f = open(os.path.join(outdir, "circos_coverage.txt"), "w")
        f.close()

    # OPTIONAL : add SNPs
    # TODO

    # 4. Run circos
    os.chdir(os.path.join(outdir))
    cmd = "circos"
    run(cmd)
    os.chdir(cwd)

def run(cmd) :
    print("Running circos in {}".format(os.getcwd()))
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    proc.communicate()
    print("Done!")

def add_small_alignments(dAln, targets_to_plot, outfile, min_align_length, correspondance) :
    for target in targets_to_plot :
        for alignment in dAln[target.name] :
            if alignment.length > min_align_length :
                continue
            start = alignment.T_start
            end = alignment.T_end
            color = "black"
            region = "{} {} {} fill_color={}".format(target.tag, start, end, color)
            outfile.write(region + "\n")

def replacePattern(input_file, pattern, replacement) :
    # Read in the file
    with open(input_file, 'r') as file :
        filedata = file.read()
    # Replace the target string
    filedata = filedata.replace(pattern, replacement)
    # Write the file out again
    with open(input_file, 'w') as file:
        file.write(filedata)

def main() :
    args = parseArgs()
    templates = args.templates[0]
    query = args.Query[0]
    targets = args.BED[0]
    reference = args.Reference[0]
    paf = args.PAF[0]
    bed = args.bed[0]
    scov = args.scov[0] # samtools formatted coverage
    ccov = args.ccov[0] # circos formatted coverage
    snps = args.snps[0]
    outdir = args.Output[0]
    min_len = args.min_length[0]
    iFH = FH(query, targets, reference, paf, outdir, bed, scov, ccov, snps, templates) # File handler
    print(iFH)

    # Get useful dictionaries
    correspondance, contigs_lengths = buildCorrespondance(iFH.reference, iFH.query)

    # Get targets object list
    targets_to_plot = readTargets(iFH.targets, correspondance, contigs_lengths)

    # Extract alignments in the .paf if these are in the targets bed file
    dAlns = readPAF(iFH.paf, targets_to_plot)

    # Prepares files and circos plot
    plotPAF(dAlns, iFH.outdir, iFH.templates, iFH.query, iFH.reference, correspondance, min_len, contigs_lengths, targets_to_plot, bed=iFH.bed, scov=iFH.scov, ccov=iFH.ccov, snps=iFH.snps)

if __name__ == '__main__':
        main()
