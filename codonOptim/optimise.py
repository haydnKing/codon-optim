import bias, rnafold
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet

import argparse, os.path

def get_arguments():
	parser = argparse.ArgumentParser(description="Perform codon optimisation")
	parser.add_argument("genome", help="Genome file to build codon table from")
	parser.add_argument("gene_sequence", nargs='+', help="Sequence file(s) to optimise")
	parser.add_argument("-s", "--scheme", 
											nargs=1, 
											required=False, 
											default="simple",
											choices=["simple",],
											help="Which optimisation scheme to use. Valid option are "+
											  "\'simple\' (default): replace each codon with a "+
												"randomised codon with probability equal to background"+
												"(subject to --ignore-rare)"
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

def load_sequence(filename):

	try:
		return SeqIO.read(filename, "fasta")
	except ValueError:
		pass

	try:
		return SeqIO.read(filename, "genbank")
	except ValueError:
		pass
	
	raise ValueError("Could not load file \'{}\'.".format(filename))

def main():

	args = get_arguments()

	genome = load_sequence(args.genome)

	b = bias.Bias(genome)

	print(b.stats())
	print(b)


	for filename in args.gene_sequence:
		seq = load_sequence(filename)
		if args.scheme == "simple":
			out = codon_optimise(b, seq, args.ignore_rare/100.)

		SeqIO.write(out, os.path.splitext(filename)[0] + ".optim.fasta", "fasta")

def translate(seq):
	return [bias.inv_codon_table[str(seq[i:i+3]).upper()] for i in range(0, len(seq), 3)]

def codon_optimise(b, sr, rare_codon_cutoff=0.1):

	AA = translate(str(sr.seq))
	oseq = []
	for aa in AA:
		oseq.append(b.emit(aa, rare_codon_cutoff))
	oseq = ''.join(oseq)

	oAA = translate(oseq) 
	if AA != oAA:
		raise ValueError("Translations don't match")


	return SeqRecord(Seq(oseq, Alphabet.generic_dna),
									 id = sr.id,
									 name=sr.name,
									 description=sr.description)
	


if __name__ == '__main__':
	main()

