import bias, rnafold
import Bio.SeqIO as SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.Alphabet as Alphabet


def main():

	marchantia = SeqIO.read("./data/m.polymorpha.gb", "genbank")

	b = bias.Bias(marchantia)

	print(b.stats())
	print(b)

	test = SeqIO.read("./data/efor-red.fasta", "fasta")
	UTR = SeqIO.read("./data/UTR.fasta", "fasta")

	out = codon_optimise(b, test, 0.1, UTR)

	SeqIO.write(out, "./data/efor-red-optim.fasta", "fasta")


def codon_optimise(b, sr, rare_codon_cutoff=0.1, UTR=None):

	oseq = []
	for icdn in (sr.seq[i:i+3].upper() for i in range(0, len(sr), 3)):
		oseq.append(b.emit(bias.inv_codon_table[icdn], rare_codon_cutoff))

	if UTR:
		start = str(UTR.seq)[-22:] + ''.join(oseq)[:16]
		print(start)
		print(rnafold.fold(start))

	return SeqRecord(Seq(''.join(oseq), Alphabet.generic_dna),
									 id = sr.id,
									 name=sr.name,
									 description=sr.description)
	


if __name__ == '__main__':
	main()

