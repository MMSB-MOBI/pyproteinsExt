import sys, getopt, numpy
import json
import os
import subprocess
from cStringIO import StringIO
from Bio import SearchIO
import pyproteins.sequence.peptide

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"hi:o:",[])

    except getopt.GetoptError:
        print 'pfamCutter.py -i <inputfile> -o <tag for outputfiles>'
        sys.exit(22)
    for opt, arg in opts:
        if opt == '-h':
            print 'pfamCutter.py -i <inputfile> -o <tag>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt == "-o":
            outputTag = arg

    c = Cutter()
    c.cut(fileName=inputFile)
    c.write()


class Cutter():
    def __init__(self, hmmrPath="/Users/guillaumelaunay/work/bin/hmmer-3.1b2-macosx-intel/binaries", hmmProfilePath="/Users/guillaumelaunay/work/data/hmm/profiles"):
        self.path = {'profiles' : hmmProfilePath, 'binaries' : hmmrPath }
        self._index()
    def _index(self):
        self.profiles = [ f.replace(".h3i", "") for f in os.listdir(self.path["profiles"]) if f.endswith(".h3i") ]

    def cut(self, **kwargs):
        self.ruler="Q8NF91"
        if 'ruler' in kwargs:
            self.ruler = kwargs['ruler']

        if self.ruler not in self.profiles:
            raise ValueError(ruler + " is not a registred database of hmm profiles")

        if 'fileName' in kwargs:
            print self.path["profiles"] + '/' + self.ruler + ' ' + kwargs['fileName']
            p = subprocess.Popen([ self.path["binaries"] + '/hmmscan', self.path["profiles"] + '/' + self.ruler, kwargs['fileName'] ], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            out, err = p.communicate()
            self._parser(out, err)

    def write(self):
        for i, peptide in enumerate(self):
            with open (peptide.id.replace(" ", "_") + "_" + str(i) + ".fasta", "w") as f:
                f.write(peptide.fasta)

    def _parser(self, stdout, stderr):
            f = StringIO(stdout)
            self.results = SearchIO.read(f, "hmmer3-text")

    def __iter__(self):
        for hits in self.results:
            for i, h in enumerate(hits):
                r = h.query_range
                yield pyproteins.sequence.peptide.Entry(id=self.results.id + " " + hits.id + " " + str(r), seq=str(h.query.seq).replace("-", ""))

if __name__ == "__main__":
    main(sys.argv[1:])