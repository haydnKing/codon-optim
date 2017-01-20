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
                "auto_PCA",
                "PCA",
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
    parser.add_argument("--pca-groups",
            required=False,
            type=lambda x: parse_pca_groups(x, parser),
            default='',
            help="Specify the groups for PCA optimisation in the "+
            "format: G1NAME,G1x,G1y[;G2NAME,G2x,G2y[;G3 ...]]")
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
    parser.add_argument("-t", "--save-table",
            required=False,
            default='',
            help="Save the codon table to the given file")
    parser.add_argument("-K", "--clusters",
            required=False,
            type=int,
            default=3,
            help="Cluseter the PCA scores of each available gene "+
            "into K clusters using GMM-EM and produce 1st order "+
            "codon tables based on those cluster weights")
    parser.add_argument("-p", "--prior-weight",
            required=False,
            type=float,
            default=5.0,
            help="Weight to give the priors for PCA analysis, "+
            "should be non-zero")
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

    if args.scheme == 'PCA':
        if not args.pca_groups:
            raise ValueError("argument --pca-groups is required for --scheme=PCA")

    return args

def parse_pca_groups(s, parser):
    if not s:
        return []
    ret = []
    groups = s.split(';')
    for g in groups:
        v = g.split(',')
        if len(v) != 4:
            raise parser.error("pca-groups: group spec should include name,x,y,radius \"{}\"".format(g))
        try:
            ret.append((v[0], float(v[1]), float(v[2]), float(v[3]),))
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
        v = (int(l), parse_DNA(cdn))
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
    if args.save_table != '':
        with open(args.save_table, 'w') as f:
            f.write(str(gs))

    if args.output_folder != '':
        if not os.path.isdir(args.output_folder):
            resp = ''
            while resp.upper() not in ['Y','N',]:
                resp = input(("Output directory \'{}\' doesn't exist. "+
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

    for filename in args.gene_sequence:
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

            sequences = []
            for v in range(args.versions):
                current_seqs = []
                if args.scheme == "simple":
                    current_seqs.append(optim.simple(gs, str_seq, args.ignore_rare/100.))
                elif args.scheme == "exact":
                    current_seqs.append(optim.exact(gs, str_seq, args.ignore_rare/100.))
                elif args.scheme == "second_rand":
                    current_seqs.append(optim.second(gs, str_seq, args.ignore_rare/100., mode='rand'))
                elif args.scheme == "auto_PCA":
                    current_seqs.extend(optim.auto_PCA(gs, 
                        str_seq, 
                        args.ignore_rare/100., 
                        args.clusters, 
                        args.prior_weight))
                elif args.scheme == "PCA":
                    current_seqs.extend(optim.by_PCA(gs, 
                        str_seq, 
                        args.pca_groups,
                        args.ignore_rare/100.,
                        args.prior_weight))
                elif args.scheme == "second_maximum":
                    current_seqs.append(optim.second(gs, str_seq, args.ignore_rare/100., mode='maximum'))
                elif args.scheme == "second_minimum":
                    current_seqs.append(optim.second(gs, str_seq, args.ignore_rare/100., mode='minimum'))
                elif args.scheme == "demo":
                    current_seqs.extend(optim.second_demo(gs, str_seq, args.ignore_rare/100.))

                if args.versions > 1:
                    current_seqs = [("{}.v{}".format(n, v), s) for n,s in current_seqs]

                current_seqs = [(n, args.upstream + override(args.override, s) + args.downstream) for
                                n,s in current_seqs]
                
                current_seqs = [(n,s) for n,s in current_seqs if
                                valid(args.exclude,s)]

                sequences.extend(current_seqs)

        if not args.amino:
            save_score_figs(gs, head, title, args.scheme, 
                    [('original',str(sr.seq)),] + sequences)
        else:
            save_score_figs(gs, head, title, args.scheme, 
                    sequences)

        pca = PCA.PrincipalComponentAnalysis.from_GenomeStats(gs, prior_weight=args.prior_weight)
        for name, seq in sequences:
            pca.add_sequence(name, seq)
        if not args.amino:
            pca.add_sequence('original', str(sr.seq))
            order = [gs.name(), 'original',]+[s[0] for s in sequences]
        else:
            order = [gs.name(), ]+[s[0] for s in sequences]

        #dirty hack
        cmap = plt.get_cmap()
        cols = [cmap(i/float(len(order)-1)) for i in range(len(order)-1)]
        ax = pca.plot(order=order,
                colors=['gray',] + cols)
        if args.scheme == 'PCA':
            #draw the circles
            for (n,x,y,r),c in zip(args.pca_groups, cols):
                ax.add_patch(mpatches.Circle((x,y), r, fill=False, color=c))
        ax.figure.savefig(os.path.join(head, title + ".PCA.png"))

        print("Generated {} sequences".format(len(sequences)))
        for name, out_str in sequences:
            ofile = os.path.join(head, title + "." + name + ".fasta")
            out_sr = SeqRecord(Seq(out_str, Alphabet.generic_dna),
                    id = sr.id,
                    name=sr.name,
                    description=sr.description)
            SeqIO.write(out_sr, ofile, "fasta")

if __name__ == '__main__':
    main()

