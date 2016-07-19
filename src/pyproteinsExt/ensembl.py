import pyproteins.container.Core
import pyproteins.sequence.peptide
import json


#geneID="FBgn0261836" # MSP300 default

class Gene(object):
    def __init__(self, name=None,endPoint="http://rest.ensembl.org/"):
        self.endPoint = endPoint
        self.polypepSet = None
        if not name:
            raise ValueError("must specify a gene name")
        self.id = name
        # Check if it exists
    @property
    def ppSet(self):
        if not self.polypepSet:
            url = self.endPoint + "sequence/id/" + self.id + '?content-type=application/json;multiple_sequences=1;type=protein'
            self.polypepSet = polyPeptideSet(id=self.id, url=url)
        return self.polypepSet

class polyPeptideSet(pyproteins.container.Core.Container):
    def __init__(self, id, url):
        if not id:
            raise TypeError('identifier is empty')
        pyproteins.container.Core.Container.__init__(self, id, url=url) # smt

        self.data = json.loads(self.raw)

    def __iter__(self):
        for e in self.data:
            yield pyproteins.sequence.peptide.Entry(datum=e)
    def __getitem__(self, i):
        return pyproteins.sequence.peptide.Entry(datum=self.data[i])

