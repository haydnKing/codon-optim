from genomestats import GenomeStats, inv_codon_table
from util import load_sequence

import argparse, os.path, os, numpy as np
import matplotlib.pyplot as plt

def ps_arguments(parser, args):
	if not args:
		return []
	a = args.split(',')
	if not len(a) == 3:
		raise parser.error("prior-scan: expected 3 arguments, got {}".format(len(a)))
	types = [float, float, int,]
	try:
		a = [t(x) for t,x in zip(types, a)]
	except ValueError:
		raise parser.error("prior-scan: could not parse {.__name__} from \'{}\'".format(t,x))
	return a

def get_arguments():
	parser = argparse.ArgumentParser(description="Perform Principal Component Analysis (PCA) on a genome, showing extra genes")
	parser.add_argument("-o", "--output-file",
											required=False,
											default='out.png',
											help="File to save output plot to")
	parser.add_argument("stats", help="Root name of preprocessed stats files")
	parser.add_argument("gene_sequence", nargs='*', help="Extra genes to plot alongside PCA")
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
												"could result in error if set to zero")
	parser.add_argument("--prior-scan",
											required=False,
											type=lambda x: ps_arguments(parser, x),
											default='',
											help="Perform a scan of prior weights. "+
											"Argument should be in the format START,STOP,NUMBER. "+
											"If this argument is given, prior-weight and cluster "+
											"are ignored")
	parser.add_argument("--columns",
											required=False,
											type=int,
											default=3,
											help="Number of columns to use when plotting the results "+
											"of prior-scan. A value <= 0 results in as many columns "+
											"as there are plots")

	args = parser.parse_args()

	return args

def main():

	args = get_arguments()

	gs = GenomeStats.from_file(args.stats)

	plot_names = []
	plot_seqs = []
	for filename in args.gene_sequence:
		seq = load_sequence(filename)
		head,tail = os.path.split(filename)
		title,ext = os.path.splitext(tail)

		plot_names.append(title)
		plot_seqs.append(str(seq.seq).upper())

	if not args.prior_scan:
		print("PCA 1 of 1")
		f = gs.plot_pca(plot_names, plot_seqs, prior_weight=args.prior_weight).figure
	else:
		f = plt.figure()
		rows = int(np.ceil(args.prior_scan[2] / float(args.columns)))
		for i,ps in enumerate(np.linspace(*args.prior_scan)):
			print("PCA {} of {}".format(i+1, args.prior_scan[2]))
			ax = f.add_subplot(rows, args.columns, i+1)
			gs.plot_pca(plot_names, plot_seqs, prior_weight=args.prior_weight, ax=ax)
			ax.set_title("prior weight = {:.1f}".format(ps))

	f.savefig(args.output_file)

if __name__ == '__main__':
	main()

