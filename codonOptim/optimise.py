import rnafold, random
from genomestats import GenomeStats, inv_codon_table
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet
from util import load_sequence

import argparse, os.path, os, numpy as np
import matplotlib.pyplot as plt

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
															 "second_minimum",],
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
												"second order hypothesis"
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
		return

	for filename in args.gene_sequence:
		print("seq = load_sequence({})".format(filename))
		seq = load_sequence(filename)
		head,tail = os.path.split(filename)
		title,ext = os.path.splitext(tail)
		print("Optimising \'{}\'".format(title))
		if args.output_folder != '':
			head = args.output_folder
		
		plot_names = ['original',]
		plot_sequences = [seq,]
		for v in range(args.versions):
			if args.versions > 1:
				print("\tversion {} / {}".format(v+1, args.versions))
			print("Using optimisation scheme \'{}\'".format(args.scheme))
			if args.scheme == "simple":
				out = simple_optimise(gs, seq, args.ignore_rare/100.)
			elif args.scheme == "exact":
				out = exact_optimise(gs, seq, args.ignore_rare/100.)
			elif args.scheme == "second_rand":
				out = second_optimise(gs, seq, args.ignore_rare/100., mode='rand')
			elif args.scheme == "second_maximum":
				out = second_optimise(gs, seq, args.ignore_rare/100., mode='maximum')
			elif args.scheme == "second_minimum":
				out = second_optimise(gs, seq, args.ignore_rare/100., mode='minimum')

			ax = gs.plot_score(['original','optimised',], [seq.seq, out.seq])
			ax.figure.savefig(os.path.join(head, "{}.s1.{}.png".format(title, args.scheme)))
			ax = gs.plot_score(['original','optimised',], [seq.seq, out.seq], order=2)
			ax.figure.savefig(os.path.join(head, "{}.s2.{}.png".format(title, args.scheme)))

			plot_names.append('v{}'.format(v))
			plot_sequences.append(out)
			if args.versions > 1:
				cols = int(np.ceil(np.log10(args.versions)))
				ofile = os.path.join(head, 
														 title + 
														 (".optim.v{:0"+str(cols)+"d}.fasta").format(v))
				out.description = out.description + " v{}".format(v)
			else:
				ofile = os.path.join(head, title + ".optim.fasta")

			print("|1st order log-prob| = {:.2f} -> {:.2f}".format(gs.score(seq.seq),
																														 gs.score(out.seq)))
			print("|2nd order log-prob| = {:.2f} -> {:.2f}".format(gs.so_score(seq.seq),
																														 gs.so_score(out.seq)))

			SeqIO.write(out, ofile, "fasta")

		ax = gs.plot_pca(plot_names, plot_sequences, prior_weight=args.prior_weight)
		ax.figure.savefig(os.path.join(head, title + ".PCA.png"))



def translate(seq):
	return [inv_codon_table[str(seq[i:i+3]).upper()] for i in range(0, len(seq), 3)]

def exact_optimise(gs, sr, rare_codon_cutoff=0.):
	
	AA = translate(str(sr.seq))
	codons = gs.generate_codons(AA, cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AA:
		cdn = codons[aa].pop(np.random.randint(0, len(codons[aa])))
		oseq.append(cdn)
	oseq = ''.join(oseq)

	oAA = translate(oseq)
	if AA != oAA:
		print(AA)
		print(oAA)
		raise ValueError("Translations don't match")

	return SeqRecord(Seq(oseq, Alphabet.generic_dna),
									 id = sr.id,
									 name=sr.name,
									 description=sr.description + ". codon optimised")

def simple_optimise(gs, sr, rare_codon_cutoff=0.):

	AA = translate(str(sr.seq))
	oseq = gs.emit(AA, rare_codon_cutoff)

	oAA = translate(oseq) 
	if AA != oAA:
		print(AA)
		print(oAA)
		raise ValueError("Translations don't match")

	return SeqRecord(Seq(oseq, Alphabet.generic_dna),
									 id = sr.id,
									 name=sr.name,
									 description=sr.description + ". codon optimised")

def second_optimise(gs, sr, rare_codon_cutoff=0., mode='rand'):

	if mode not in ['rand','maximum','minimum']:
		raise ValueError("Unknown mode \'{}\'".format(mode))
	
	AA = translate(str(sr.seq))
	codons = gs.generate_codons(AA, cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AA:
		if not oseq:
			oseq.append(codons[aa].pop(np.random.randint(0, len(codons[aa]))))
		else:
			p = gs.so().loc[oseq[-1], codons[aa]]
			#if there are no instances, make the distribution uniform
			if p.sum() == 0:
				p[:] = np.ones(len(p))
			p = p / float(p.sum())
			if mode == 'rand':
				r = random.random()
				for i,s in enumerate(p.cumsum()):
					if r < s:
						oseq.append(codons[aa].pop(i))
						break
			elif mode == 'maximum':
				p.index = range(len(p))
				oseq.append(codons[aa].pop(np.argmax(p)))
			elif mode == 'minimum':
				p.index = range(len(p))
				oseq.append(codons[aa].pop(np.argmin(p)))

			
	oseq = ''.join(oseq)

	oAA = translate(oseq)
	if AA != oAA:
		print(AA)
		print(oAA)
		raise ValueError("Translations don't match")

	return SeqRecord(Seq(oseq, Alphabet.generic_dna),
									 id = sr.id,
									 name=sr.name,
									 description=sr.description + ". codon optimised")
	

def compare_priors(gs, priors=[0.1, 0.5, 1.0, 5.0, 10., 50., 100,]):
	cols = 3
	rows = np.ceil(len(priors)/float(cols))

	fig = plt.figure()
	for i,p in enumerate(priors):
		print("{}/{}".format(i+1, len(priors)))
		gs.plot_pca(prior_weight=p, ax=fig.add_subplot(rows, cols, i+1))

	fig.show()

if __name__ == '__main__':
	main()

