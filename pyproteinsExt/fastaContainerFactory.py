import pyproteinsext.proteinContainer
import gzip


def parse(inputFile=None):
    """Will parse fasta file to create Container
    
    :param inputFile: path of fasta file to parse. 
    :return: fastaContainerFactory.Container
    """
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


class Container(pyproteinsext.proteinContainer.Container):
    """ Iterable container, that stores all Fasta objects """
    def __init__(self, input=None):
        super().__init__(_parseBuffer, input)


class Fasta():
    """
    Fasta object stores informations about one fasta entry. 

    :param prot: Protein id, first element of header, before first space
    :type prot: str
    :param seq: Amino acids or nucleotides sequence
    :type seq: str
    :param header: Fasta header
    :type header: str
    """

    def __init__(self, prot, header, seq):
        """ Constructor method """
        self.prot = prot
        self.seq = seq
        self.header = header

    def get_subsequence(self, start, end):
        """
        Get subsequence of entry, from start to end position
        
        :param start: Start position
        :type start: int
        :param end: End position
        :type end: int
        :return: subsequence
        :rtype: str
        """
        return self.seq[start-1:end]
