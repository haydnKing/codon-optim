"""Library for optimisations"""
import numpy as np, pandas as pd, copy
import util, PCA

def _verify(in_aa, out_str):
	out_aa= util.translate(out_str)

	if out_aa != in_aa:
		print("in : {}".format(in_aa))
		print("out: {}".format(out_aa))
		raise ValueError("Optimisation failed: translations don't match")

	return out_str

def exact(gs, AAseq, rare_codon_cutoff=0.):
	codons = _generate_codons(AAseq, gs.norm_bias(), cutoff=rare_codon_cutoff)

	oseq = []
	for aa in AAseq:
		cdn = codons[aa].pop(np.random.randint(0, len(codons[aa])))
		oseq.append(cdn)
	oseq = ''.join(oseq)

	return _verify(AAseq, ret)

def simple(gs, AAseq, rare_codon_cutoff=0.):

	oseq = gs.emit(AAseq, rare_codon_cutoff)

	return _verify(AAseq, oseq)

def second(gs, AAseq, rare_codon_cutoff=0., mode='rand'):

	if mode not in ['rand','maximum','minimum','irand']:
		raise ValueError("Unknown mode \'{}\'".format(mode))
	
	codons = _generate_codons(AAseq, gs.norm_bias(), cutoff=rare_codon_cutoff)
	oseq = _second(gs.so(), AAseq, codons, mode)
	return _verify(AAseq, oseq)


def _second(so, AAseq, codons, mode='rand'):

	oseq = []
	for aa in AAseq:
		if not oseq:
			oseq.append(codons[aa].pop(np.random.randint(0, len(codons[aa]))))
		else:
			p = so.loc[oseq[-1], codons[aa]]
			#if there are no instances, make the distribution uniform
			if p.sum() == 0:
				p[:] = np.ones(len(p))
			p = p / float(p.sum())
			if mode in ['rand', 'irand']:
				r = np.random.random()
				if mode == 'irand' and len(p) > 1:
					p = 1-p
					p = p / p.sum()
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
	return ''.join(oseq)


def auto_PCA(gs, 
						 AAseq, 
						 rare_codon_cutoff=0., 
						 GMM_components=3, 
						 prior_weight=1., 
						 PCA_components=3,
						 mode='rand'):

	ret = []
	#Perform PCA and GMM/EM clustering
	pca = PCA.PrincipalComponentAnalysis.from_GMM(gs.fo(),
																								K=GMM_components,
																								PCA_components=PCA_components,
																								prior_weight=prior_weight)

	#for each cluster
	for name, indexes in pca.labels().items():

		#calculate first order bias
		data = gs.fo().loc[indexes]

		bias = util.normalise(data.sum(0))

		#generate codon lists given fo bias	
		codons = _generate_codons(AAseq, bias, cutoff=rare_codon_cutoff)

		#order codons according to whole genome so preference
		oseq = _second(gs.so(), AAseq, codons, mode)
		
		seq = _verify(AAseq, oseq)

		ret.append(seq)
	
	return ret


def _generate_codons(AAseq, bias, cutoff=0.):
    """Generate unordered codons to use for each amino acid present in AAseq 
    such that the codon usage is as close as possible to bias."""
    bias = util.normalise(bias)

    out = {}
    for aa in util.AA:
        #list all codons which could be used for this aa
        cdn_list = [c for c in util.codon_table[aa] if bias[c] > cutoff]
        #how many codons do we need for this aa?
        count = len([1 for aas in AAseq if aas == aa])
        #what number of each codon should we have?
        counts = (bias[cdn_list] / np.sum(bias[cdn_list]))*count

        #sort by smallest residual
        counts = pd.DataFrame({'c':counts, 
                               'r':np.abs(counts-np.around(counts))
                              }).sort_index(by='r')['c']
        #assign integers
        overflow = 0.
        icounts = pd.Series(np.zeros(len(counts), dtype=int), index=counts.index)
        for i in range(len(counts)):
            icounts[i] = int(np.round(counts[i]+overflow))
            overflow = overflow + counts[i] - icounts[i]

        #list of codons
        out[aa] = []
        for cdn,count in icounts.iteritems():
            out[aa] = out[aa] + [cdn,]*count
            #shuffle the list (in some schemes, the codons are taken in list order
            #when the genome lacks information)
            np.random.shuffle(out[aa])

    return out

def by_PCA(gs, 
           AAseq, 
           pca_group,
           rare_codon_cutoff=0.,
           prior_weight=1.,
           mode='rand'):

    pca = PCA.PrincipalComponentAnalysis.from_GenomeStats(gs, prior_weight=prior_weight)

    name, x, y, r = pca_group
    indexes = pca.label_from_circle(name,x,y,r)
    codons = _generate_codons(AAseq, gs.get_bias(indexes), cutoff=rare_codon_cutoff)
    oseq = _second(gs.so(), AAseq, codons, mode)
    return _verify(AAseq, oseq)

def second_demo(gs, AAseq, rare_codon_cutoff, repeat=500):

	nb = gs.norm_bias()

	#pre-calculate codon probabilities
	codon_probs = {}
	for aa,ct in util.codon_table.items():
		codon_probs[aa] = (nb[ct] / nb[ct].sum()).cumsum()


	best_delta = 0.
	best_good = ''
	best_bad = ''
	for i in range(repeat):
		#randomly generate codons
		codons = {}
		for aa,ct in util.codon_table.items():
			c = []
			for i in range(len([1 for aa_ in AAseq if aa_ == aa])):
				r = np.random.random()
				for cdn,p in codon_probs[aa].items():
					if p > r:
						c.append(cdn)
						break
			codons[aa] = c

		bc = copy.deepcopy(codons)
		good = _second(gs.so(), AAseq, codons, 'rand')
		bad  = _second(gs.so(), AAseq, bc, 'irand')
		delta = gs.so_score(good) - gs.so_score(bad)
		print("delta = {}".format(delta))

		if delta > best_delta:
			best_delta = delta
			best_good = good
			best_bad = bad

	print("best_delta = {}".format(best_delta))
	return [('good', _verify(AAseq, best_good),), ('bad', _verify(AAseq, best_bad))]

class Optimisation:
    """Hold all the details of an optimisation run"""

    def __init__(self, genomestats, name, seq, amino, scheme, group=[], rare_cutoff=0., versions=1, 
                 exclude=[], upstream="", downstream="", override=[],
                 prior_weight=1.0):
        self.name = name
        if not amino:
            self.original = seq
            self.aseq = util.translate(seq)
        else:
            self.original = ""
            self.aseq = seq
        self.scheme = scheme
        self.group = group
        self.rare_cutoff = rare_cutoff
        self.versions = versions
        self.exclude = exclude
        self.upstream = upstream
        self.downstream = downstream
        self.override = override
        self.prior_weight = prior_weight
        self.gs = genomestats

    def run(self):
        self.output = []
        attempts = 0
        while len(self.output) < self.versions:
            seq = self._optimise()
            if self._validate(seq):
                self.output.append(seq)

            attempts = attempts + 1
            if attempts > 10 * self.versions:
                print(("Failed to generate {} versions of '{}' "+
                       "after {} attempts, got {} versions instead")
                      .format(self.versions, 
                              self.name, 
                              attempts,
                              len(self.output)))
                break

    def _optimise(self):
        seq = ""
        if self.scheme == "simple":
            seq = simple(self.gs, 
                         self.aseq, 
                         self.rare_cutoff/100.)
        elif self.scheme == "exact":
            seq = exact(self.gs, 
                        self.aseq, 
                        self.rare_cutoff/100.)
        elif self.scheme == "second_rand":
            seq = second(self.gs, 
                         self.aseq, 
                         self.rare_cutoff/100., 
                         mode='rand')
        elif self.scheme == "group":
            seq = by_PCA(self.gs, 
                         self.aseq, 
                         self.group,
                         self.rare_cutoff/100.,
                         self.prior_weight)
        elif self.scheme == "second_maximum":
            seq = second(self.gs, 
                         self.aseq, 
                         self.rare_cutoff/100., 
                         mode='maximum')
        elif self.scheme == "second_minimum":
            seq = second(self.gs, 
                         self.aseq, 
                         self.rare_cutoff/100., 
                         mode='minimum')

        return self._override(seq)

    def _override(self, seq):
        if not self.override:
            return seq
        for i, cdn in self.override:
            seq = seq[:3*(i-1)] + cdn + seq[3*i:]
        return seq

    def _validate(self, seq):
        if not self.exclude:
            return True
        for e in self.exclude:
            test_seq = (self.upstream[-len(e)+1:] + seq +
                        self.downstream[:len(e)-1])
            if (test_seq.find(e) >= 0 or 
                test_seq.find(util.reverse_complement(e)) >= 0):
                return False
        return True

    def get_results(self):
        return [self.upstream + o + self.downstream
                for o in self.output]

    def get_names(self):
        return ["{}.{}.v{}".format(self.name, self.scheme, i+1) 
                for i in range(len(self.output))]

