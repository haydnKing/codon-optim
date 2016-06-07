import Bio.SeqIO as SeqIO

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

