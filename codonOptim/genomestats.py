import numpy as np, pandas as pd
import Bio.SeqIO as SeqIO
import matplotlib.pyplot as plt
import os.path
import util

class GenomeStats:
    """Represent statistical information about a genome or collection of genes"""
    def __init__(self, _name, _data, _second_order, _scores, _ndata=None, _nso=None):
        """Calculate codon bias in CDS annotations.
        sr: genome seqrecord"""

        self._name = _name
        self._data = _data
        self._second_order = _second_order
        self._scores = _scores
        self._bias = self._data.sum(0)
        if _ndata is None:
            self._normed = util.normalise(self._bias)
        else:
            self._normed = _ndata
            if _nso is None:
                self._so_normed = util.so_normalise(self._second_order)
            else:
                self._so_normed = _nso

    @classmethod
    def from_seqrecord(cls, sr, featuretype='CDS', name=None):

        if not name:
            name = sr.name

        CDS = [f for f in sr.features if f.type == featuretype]

        _data = pd.DataFrame(np.zeros((len(CDS), 64), dtype=int), 
                             columns = util.list_codons())

        _second_order = pd.DataFrame(np.zeros((64,64), dtype=int),
                                     index = util.list_codons(),
                                     columns = util.list_codons())

        _scores = pd.DataFrame(np.zeros((len(CDS), 2)), 
                               columns = ['first', 'second',])

        _seqs = [util._extract(sr, cds) for cds in CDS]
        for i,seq in enumerate(_seqs):
            _data.loc[i,:] = util.get_bias(seq)
            util.add_second_order(_second_order, seq)

        #calculate scores
        _nd = util.normalise(_data.sum(0))
        _nso= util.so_normalise(_second_order)
        for i,seq in enumerate(_seqs):
            _scores.at[i,'first'] = util.score(_nd, seq)
            _scores.at[i,'second'] = util.so_score(_nso, seq)

        return cls(name, _data, _second_order, _scores, _nd, _nso)

    def save_stats(self, folder, name=None):
        if not name:
            name = self._name

        self._data.to_csv(os.path.join(folder, name + ".1.csv"))
        self._second_order.to_csv(os.path.join(folder, name + ".2.csv"))
        self._scores.to_csv(os.path.join(folder, name + ".scores.csv"))


    @classmethod
    def from_file(cls, folder, name = None):
        if not name:
            folder, name = os.path.split(folder)

        _data   = pd.read_csv(os.path.join(folder, name + ".1.csv"), index_col=0)
        _so     = pd.read_csv(os.path.join(folder, name + ".2.csv"), index_col=0)
        _scores = pd.read_csv(os.path.join(folder, name + ".scores.csv"), index_col=0)

        return cls(name, _data, _so, _scores)


    def raw_bias(self):
        return self._bias

    def norm_bias(self):
        return self._normed

    def fo(self):
        return self._data

    def fo_normed(self):
        return self._normed

    def so(self):
        return self._second_order

    def so_normed(self):
        return self._so_normed

    def name(self):
        return self._name

    def get_bias(self, indexes):
        return util.normalise(self._data.loc[indexes].sum(0))

    def emit(self, AAseq, cutoff=0.):
        """Emit a codon to code for the amino acid 'aa'. Ignore codons with 
        frequency below cutoff (e.g. 0.1 = 10%)
        """
        values = {}
        for aa in util.AA:
            cdn_list = [c for c in util.codon_table[aa] if self._normed[c] > cutoff]
            if not cdn_list:
                raise ValueError("No possible codons for aa=\'{}\'. Try a lower cutoff".format(aa))

            #emit a codon with probability equal to the remains
            values[aa] = self._bias[cdn_list] / float(self._bias[cdn_list].sum())

        oseq = []
        for aa in AAseq:
            r = np.random.random()
            for cdn,s in values[aa].cumsum().iterkv():
                if r < s:
                    oseq.append(cdn)
                    break
        return ''.join(oseq)

    def score(self, seq):
        return util.score(self._normed, seq)

    def so_score(self, seq):
        return util.so_score(self._so_normed, seq)

    def plot_score(self, names=[], seqs=[], colors=None, order=1):
        if order == 1:
            scores = self._scores['first']
        elif order == 2:
            scores = self._scores['second']
        else:
            raise ValueError("order must be 1 or 2, not {}".format(order))

        ax = plt.figure().gca()

        if colors is None:
            cmap = plt.get_cmap()
            colors = [cmap(i/float(len(names))) for i in range(len(names))]

        #histogram of the scores
        ax.hist(scores, bins=20, label = 'background', color='gray')

        for name,seq,color in zip(names,seqs,colors):
            if order == 1:
                score = self.score(seq)
            elif order == 2:
                score = self.so_score(seq)
                ax.axvline(score, color=color, label=name)

        ax.legend(loc=0)

        if order == 1:
            ax.set_xlabel("First Order Average log-probability")
        elif order == 2:
            ax.set_xlabel("Second Order Average log-probability")

        return ax

    def __str__(self):
        return util.fo_to_string(name=self._name, data=self._data)



