"""A script for optimisation"""

import rnafold, random
from genomestats import GenomeStats
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet

import argparse, os.path, os, numpy as np
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import optim, PCA, util

def get_arguments():
    parser = argparse.ArgumentParser(description="Perform codon optimisation")
    parser.add_argument("stats", help="Root name of preprocessed stats files")
    parser.add_argument("gene_sequence", nargs='*', help="Sequence file(s) to optimise")
    parser.add_argument("-s", "--scheme", 
            required=False, 
            default="simple",
            choices=["simple",
                "exact",
                "group",
                "second_rand",
                "second_maximum",
                "second_minimum",
                "demo",],
            help="Which optimisation scheme to use. Valid option are "+
            "\'simple\' (default): replace each codon with a "+
            "randomised codon with probability equal to background"+
            "(subject to --ignore-rare) "+
            "\'exact\': choose codons to match the average codon "+
            "bias as closely as possible, but order the codons "+
            "randomly "+
            "\'second_*\': Choose overall codon use using the "+
            "exact method, then choose the ordering of those "+
            "codons using second order statistics from the genome "+
            "\'second_rand\': choose second order codons randomly, "+
            "distributed according to the distribution found in the "+
            "genome "+
            "\'second_maximum\': choose the most likely possible "+
            "codon given the previous codon "+
            "\'second_minimum\': choose the least likely possible "+
            "codon given the previous codon - useful for testing "+
            "second order hypothesis. "+
            "\'demo\': generate good and bad second order sequences"
            )
    parser.add_argument("-p", "--prior_weight",
                        required=False,
                        type=float,
                        default=1.0,
                        help="Specify the prior weight for the PCA analysis")
    parser.add_argument("--group",
            required=False,
            type=lambda x: parse_group(x, parser),
            default='',
            help="Specify the group for optimisation based on a "+
            "subset of the genes: NAME,x,y,r")
    parser.add_argument("-r", "--ignore-rare", 
            required=False,
            type=float,
            default=0.,
            help="Do not include codons with frequency below given "+
            "percentage in output. "+
            "For example \'AUU\', \'AUC\', and \'AUA\' make up "+
            "49%%, 37%%, and 13%% of Isoleucine codons in Bacillus "+
            "Subtilis; setting --ignore-rare=15 will discard "+
            "\'AUA\' and emit \'AUU\' and \'AUC\' with "+
            "probabilities 0.57 and 0.43")
    parser.add_argument("-v", "--versions",
            required=False,
            type=int,
            default=1,
            help="Number of versions of the gene to produce, "+
            "useful for screening multiple codon optimised variants "+
            "for synthesis")
    parser.add_argument("-o", "--output-folder",
            required=False,
            default='',
            help="Folder to save output in, defaults to the same as "+
            "the input files")
    parser.add_argument("--amino",
            action='store_true',
            help="Input sequence to optimise is an amino acid "+
            "not a DNA sequence")
    parser.add_argument("--exclude",
            type=lambda x: parse_DNA_list(x, parser),
            help="Comma separated list of sequences to exclude from "+
                 "optimised sequences. Includes reverese translation "+
                 "and takes --upstream and --downstream into account. "+
                 "Useful for excluding sequences with particular cut "+
                 "sites.")
    parser.add_argument("--upstream",
            type=lambda x: parse_DNA(x, parser),
            default="",
            help="Sequence to prepend to the optimised sequences")
    parser.add_argument("--downstream",
            type=lambda x: parse_DNA(x, parser),
            default="",
            help="Sequence to append to the optimised sequences")
    parser.add_argument("--override",
            type=lambda x: parse_codons(x, parser),
            default=[],
            help="Comma separated <codon_number>:<codon> pairs used to "+
                 "override codon choice at particular locations, e.g. "+
                 "'--override 5:GGT,6:ATA' sets codons 5 and 6, regardless "+
                 "of the amino acids present in the input sequence")



    args = parser.parse_args()

    if args.ignore_rare:
        if args.ignore_rare > 100 or args.ignore_rare < 0:
            raise ValueError("--ignore-rare must be within [0,100] not {}".format(args.ignore_rare))

    if args.scheme == 'group':
        if not args.group:
            raise ValueError("argument --group is required for --scheme=group")

    return args

def parse_group(s, parser):
    if not s:
        return []
    ret = []

    v = s.split(',')
    if len(v) != 4:
        raise parser.error("pca-groups: group spec should include name,x,y,radius \"{}\"".format(g))
    try:
        ret.extend((v[0], float(v[1]), float(v[2]), float(v[3]),))
    except ValueError as e:
        raise parser.error("pca-groups: {}".format(e.args[0]))
    return ret

def parse_DNA(s, parser):
    invalid = [char for char in s.upper() 
               if char not in ['A','T','G','C']]
    if len(invalid) > 1:
        raise parser.error(
            "Unknown DNA bases '{}'".format("', '".join(invalid)))
    elif len(invalid) == 1:
        raise parser.error("Unknown DNA base '{}'".format(*invalid))
    return s

def parse_DNA_list(s, parser):
    r = [parse_DNA(x, parser) for x in s.split(',')]
    return r

def parse_codons(s, parser):
    r = []
    for l,cdn in [x.split(':') for x in s.split(',')]:
        v = (int(l), parse_DNA(cdn, parser))
        if len(v[1]) != 3:
            raise parser.error("codons must be of length 3")
        r.append(v)
    return r

def override(o, seq):
    if not o:
        return seq
    for i, cdn in o:
        seq = seq[:3*(i-1)] + cdn + seq[3*i:]
    return seq

def valid(exclude, seq):
    if not exclude:
        return True
    for e in exclude:
        if (seq.find(e) >= 0 or 
            seq.find(util.reverse_complement(e)) >= 0):
            return False
    return True

def save_score_figs(gs, head, title, scheme, seqs):
    ax = gs.plot_score([s[0] for s in seqs], [s[1] for s in seqs])
    ax.figure.savefig(os.path.join(head, "{}.s1.{}.png".format(title, scheme)))
    ax = gs.plot_score([s[0] for s in seqs], [s[1] for s in seqs], order=2)
    ax.figure.savefig(os.path.join(head, "{}.s2.{}.png".format(title, scheme)))

def main():

    args = get_arguments()

    gs = GenomeStats.from_file(args.stats)

    print(gs)

    if args.output_folder != '':
        if not os.path.isdir(args.output_folder):
            resp = ''
            while resp.upper() not in ['Y','N',]:
                resp = raw_input(("Output directory \'{}\' doesn't exist. "+
                    "Create it? [Y/N]: ").format(args.output_folder))
                if resp.upper() == 'Y':
                    os.mkdir(args.output_folder)
            else:
                return

    if not args.gene_sequence:
        #save a PCA plot
        pca = PCA.PrincipalComponentAnalysis.from_GenomeStats(gs, prior_weight=args.prior_weight)
        ax = pca.plot(colors="gray")
        ax.figure.savefig(os.path.join(args.output_folder, "stats.PCA.png"))
        print("No sequence to optimise")
        return

    optimisations = []

    for filename in args.gene_sequence:
        print("loading {}".format(filename))
        seq_records = util.load_sequence(filename)
        head,tail = os.path.split(filename)
        title,ext = os.path.splitext(tail)
        for sr in seq_records:
            title = sr.id if len(sr.id) > len(sr.name) else sr.name
            print("Optimising \'{}\'".format(title))
            if args.output_folder != '':
                head = args.output_folder


            #get a translation
            if not args.amino:
                str_seq = util.translate(str(sr.seq).upper())
            else:
                str_seq = str(sr.seq).upper()

            o = optim.Optimisation(gs, title, str_seq, args.scheme, args.group,
                                   args.ignore_rare, args.versions,
                                   args.exclude, args.upstream,
                                   args.downstream, args.override,
                                   args.prior_weight)

            o.run()
            optimisations.append(o)

    for o in optimisations:
        for name, out_str in zip(o.get_names(), o.get_results()):
            ofile = os.path.join(args.output_folder, name + ".fasta")
            out_sr = SeqRecord(Seq(out_str, Alphabet.generic_dna),
                               id = name,
                               name=name,
                               description="optimised")
            SeqIO.write(out_sr, ofile, "fasta")

if __name__ == '__main__':
    main()

