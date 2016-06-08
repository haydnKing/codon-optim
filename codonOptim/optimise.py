"""A script for optimisation"""

import rnafold, random
from genomestats import GenomeStats
import Bio.SeqIO as SeqIO
from util import load_sequence

import argparse, os.path, os, numpy as np
import matplotlib.pyplot as plt
import optim

def get_arguments():
	parser = argparse.ArgumentParser(description="Perform codon optimisation")
	parser.add_argument("stats", help="Root name of preprocessed stats files")
	parser.add_argument("gene_sequence", nargs='*', help="Sequence file(s) to optimise")
	parser.add_argument("-s", "--scheme", 
											required=False, 
											default="simple",
											choices=["simple",
															 "exact",
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
	parser.add_argument("-r", "--ignore-rare", 
											required=False,
											type=float,
											default=0.,
											help="Do not include codons with frequency below given "+
												"percentage in output. "+
												"For example \'AUU\', \'AUC\', and \'AUA\' make up "+
												"49%, 37%, and 13% of Isoleucine codons in Bacillus "+
												"Subtilis; setting --ignore-rare=15 will discard "+
												"\'AUA\' and emit \'AUU\' and \'AUC\' with "+
												"probabilities 0.57 and 0.43")
	parser.add_argument("-t", "--save-table",
											required=False,
											default='',
											help="Save the codon table to the given file")
	parser.add_argument("-K", "--cluster",
											required=False,
											type=int,
											default=0,
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
											
#	parser.add_argument("--UTR",
#											nargs=1,
#											required=False,
#											help="Sequence file providing the 5'UTR which will be "+
#											"used for the optimised gene. Secondary structure "+
#											"in the bases around the RBS (estimated as [-19,+12] "+
#											"from the start codon) will be minimised while "+
#											"respecting --ignore-rare if given")

	args = parser.parse_args()

	if args.ignore_rare:
		if args.ignore_rare > 100 or args.ignore_rare < 0:
			raise ValueError("--ignore-rare must be within [0,100] not {}".format(args.ignore_rare))

	return args

def save_score_figs(gs, head, title, scheme, names, seqs):
	ax = gs.plot_score(names, seqs)
	ax.figure.savefig(os.path.join(head, "{}.s1.{}.png".format(title, scheme)))
	ax = gs.plot_score(names, seqs, order=2)
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
				resp = raw_input(("Output directory \'{}\' doesn't exist. "+
													"Create it? [Y/N]: ").format(args.output_folder))
			if resp.upper() == 'Y':
				os.mkdir(args.output_folder)
			else:
				return

	if not args.gene_sequence:
		print("No sequence to optimise")
		ax = gs.plot_pca([], 
										 [],
										 prior_weight = args.prior_weight)
		ax.figure.savefig(os.path.join(args.output_folder, gs.name()+".PCA.png"))
		return

	for filename in args.gene_sequence:
		print("seq = load_sequence({})".format(filename))
		seq = load_sequence(filename)
		head,tail = os.path.split(filename)
		title,ext = os.path.splitext(tail)
		print("Optimising \'{}\'".format(title))
		if args.output_folder != '':
			head = args.output_folder
		
		names = []
		sequences = []
		name_fmt = "v{{:0{}d}}".format(int(np.ceil(np.log10(args.versions))))
		for v in range(args.versions):
			if args.versions > 1:
				print("\tversion {} / {}".format(v+1, args.versions))
			print("Using optimisation scheme \'{}\'".format(args.scheme))
			if args.scheme == "simple":
				sequences.append(optim.simple(gs, seq, args.ignore_rare/100.))
				names.append(name_fmt.format(v))
			elif args.scheme == "exact":
				sequences.append(optim.exact(gs, seq, args.ignore_rare/100.))
				names.append(name_fmt.format(v))
			elif args.scheme == "second_rand":
				sequences.append(optim.second(gs, seq, args.ignore_rare/100., mode='rand'))
				names.append(name_fmt.format(v))
			elif args.scheme == "second_maximum":
				sequences.append(optim.second(gs, seq, args.ignore_rare/100., mode='maximum'))
				names.append(name_fmt.format(v))
			elif args.scheme == "second_minimum":
				sequences.append(optim.second(gs, seq, args.ignore_rare/100., mode='minimum'))
				names.append(name_fmt.format(v))
			elif args.scheme == "demo":
				sequences.append(optim.second(gs, seq, args.ignore_rare/100., mode='rand'))
				names.append(name_fmt.format(v) + ".max")
				sequences.append(optim.second(gs, seq, args.ignore_rare/100., mode='irand'))
				names.append(name_fmt.format(v) + ".min")


			save_score_figs(gs, head, title, args.scheme, 
											['original',] + names, 
											[seq.seq,] + [s.seq for s in sequences])

			for o,n in zip(sequences, names):
				ofile = os.path.join(head, title + ".optim." + n + ".fasta")
				o.description = o.description + " v{}".format(v)
				SeqIO.write(o, ofile, "fasta")

		ax = gs.plot_pca(['original',] + names, 
										 [seq,] + sequences,
										 prior_weight = args.prior_weight)
		ax.figure.savefig(os.path.join(head, title + ".PCA.png"))


if __name__ == '__main__':
	main()

