import Bio.SeqIO as SeqIO

AA = ['F', 'L', 'S', 'Y', '*', 'C', 'W', 'P', 'H', 'Q', 'R', 'I', 'M', 'T', 'N', 'K', 'V', 'A', 'D', 'E', 'G']
codon_table = {
	'A': ['GCT', 'GCC', 'GCA', 'GCG'], 
	'C': ['TGT', 'TGC'], 
	'E': ['GAA', 'GAG'], 
	'D': ['GAT', 'GAC'], 
	'G': ['GGT', 'GGC', 'GGA', 'GGG'], 
	'F': ['TTT', 'TTC'], 
	'I': ['ATT', 'ATC', 'ATA'], 
	'H': ['CAT', 'CAC'], 
	'K': ['AAA', 'AAG'], 
	'*': ['TAA', 'TAG', 'TGA'], 
	'M': ['ATG'], 
	'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 
	'N': ['AAT', 'AAC'], 
	'Q': ['CAA', 'CAG'], 
	'P': ['CCT', 'CCC', 'CCA', 'CCG'], 
	'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 
	'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 
	'T': ['ACT', 'ACC', 'ACA', 'ACG'], 
	'W': ['TGG'], 
	'V': ['GTT', 'GTC', 'GTA', 'GTG'], 
	'Y': ['TAT', 'TAC']
}
inv_codon_table = {}
for aa,cdn_list in codon_table.iteritems():
	for cdn in cdn_list:
		inv_codon_table[cdn] = aa

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

def check_folder_exists(name):
	if not os.path.isdir(name):
		resp = ''
		while resp.upper() not in ['Y','N',]:
			resp = raw_input(("Output directory \'{}\' doesn't exist. "+
												"Create it? [Y/N]: ").format(args.output_folder))
		if resp.upper() == 'Y':
			os.mkdir(args.output_folder)
		else:
			return False
	return True

