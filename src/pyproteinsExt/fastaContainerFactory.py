import pyproteinsExt.proteinContainer
import gzip


def parse(inputFile=None):
    mainContainer = Container()
    fnOpen = open
    t = 'r'
    if inputFile.endswith('gz'):
        fnOpen = gzip.open
        t = 'rt'
    try:
        f = fnOpen(inputFile, t)
    except IOError:
        print("Could not read file:", inputFile)
        return Container()

    instance = []
    with f:
        for l in f:
            if l.startswith(">"):
                if instance:
                    mainContainer = mainContainer.addParsing(Container(input=instance))
                instance = []
            instance.append(l.rstrip())
        mainContainer = mainContainer.addParsing(Container(input=instance))
    return mainContainer


def _parseBuffer(input):
    header = input[0].lstrip(">")
    prot = header.split(" ")[0]
    seq = "".join(input[1:])
    obj = Fasta(prot, header, seq)
    return {prot: obj}


class Container(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None):
        super().__init__(_parseBuffer, input)


class Fasta():
    def __init__(self, prot, header, seq):
        self.prot = prot
        self.seq = seq
        self.header = header

    def get_subsequence(self, start, end):
        return self.seq[start-1:end]
