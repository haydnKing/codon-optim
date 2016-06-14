"""A script for optimisation"""

import rnafold, random
from genomestats import GenomeStats
import Bio.SeqIO as SeqIO
from util import load_sequence

import argparse, os.path, os, numpy as np
import matplotlib.pyplot as plt
import optim, PCA

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
															 "second_PCA",
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
											


	args = parser.parse_args()

	if args.ignore_rare:
		if args.ignore_rare > 100 or args.ignore_rare < 0:
			raise ValueError("--ignore-rare must be within [0,100] not {}".format(args.ignore_rare))

	return args

def save_score_figs(gs, head, title, scheme, seqs):
	ax = gs.plot_score([s[0] for s in seqs], [str(s[1].seq) for s in seqs])
	ax.figure.savefig(os.path.join(head, "{}.s1.{}.png".format(title, scheme)))
	ax = gs.plot_score([s[0] for s in seqs], [str(s[1].seq) for s in seqs], order=2)
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
		print("No sequence to optimise")
		return

	for filename in args.gene_sequence:
		seq = load_sequence(filename)
		print("seq = load_sequence({}) name = {}".format(filename, seq.name))
		head,tail = os.path.split(filename)
		title,ext = os.path.splitext(tail)
		print("Optimising \'{}\'".format(title))
		if args.output_folder != '':
			head = args.output_folder
		
		sequences = []
		for v in range(args.versions):
			current_seqs = []
			if args.versions > 1:
				print("\tversion {} / {}".format(v+1, args.versions))
			print("Using optimisation scheme \'{}\'".format(args.scheme))
			if args.scheme == "simple":
				current_seqs.append(optim.simple(gs, seq, args.ignore_rare/100.))
			elif args.scheme == "exact":
				current_seqs.append(optim.exact(gs, seq, args.ignore_rare/100.))
			elif args.scheme == "second_rand":
				current_seqs.append(optim.second(gs, seq, args.ignore_rare/100., mode='rand'))
			elif args.scheme == "second_PCA":
				current_seqs.extend(optim.second_PCA(gs, 
																					seq, 
																					args.ignore_rare/100., 
																					args.clusters, 
																					args.prior_weight))
			elif args.scheme == "second_maximum":
				current_seqs.append(optim.second(gs, seq, args.ignore_rare/100., mode='maximum'))
			elif args.scheme == "second_minimum":
				current_seqs.append(optim.second(gs, seq, args.ignore_rare/100., mode='minimum'))
			elif args.scheme == "demo":
				current_seqs.append(optim.second(gs, seq, args.ignore_rare/100., mode='rand'))
				current_seqs.append(optim.second(gs, seq, args.ignore_rare/100., mode='irand'))

			if args.versions > 1:
				current_seqs = [("{}.v{}".format(n, v), s) for n,s in current_seqs]

			sequences.extend(current_seqs)

		save_score_figs(gs, head, title, args.scheme, 
										[('original',seq),] + sequences)

		pca = PCA.PrincipalComponentAnalysis.from_GenomeStats(gs, prior_weight=args.prior_weight)
		#for name, seq in sequences:
			#pca.add_sequence(name, seq)
		ax = pca.plot()
		ax.figure.savefig(os.path.join(head, title + ".PCA.png"))

		for name, seq in sequences:
			ofile = os.path.join(head, title + "." + name + ".fasta")
			SeqIO.write(seq, ofile, "fasta")

if __name__ == '__main__':
	main()

