from util import load_sequence
import Bio.SeqIO as SeqIO
from genomestats import GenomeStats
import argparse, os.path, os, pandas as pd


def get_arguments():
	parser = argparse.ArgumentParser(description="Preprocess a genome and extract all necessary stats")

	parser.add_argument("genome", help="Genome to build statistics from")
	parser.add_argument("-f", "--feature",
											required=False,
											default="CDS",
											help="Type of feature to use as coding sequence")
	parser.add_argument("-o", "--output-folder",
											required=False,
											default=os.getcwd(),
											help="Folder to store output files in")

	#parse the arguments
	args = parser.parse_args()


	return args

def main():

	args = get_arguments()

	name = os.path.splitext(os.path.split(args.genome)[1])[0]

	print("Loading genome file...")
	genome = load_sequence(args.genome)
	print("Calculating Statistics...")
	gs = GenomeStats.from_seqrecord(genome, featuretype=args.feature, name=name)
	print("Saving to {} as {}".format(args.output_folder, name))
	gs.save_stats(args.output_folder, name=name)


if __name__ == '__main__':
	main()
