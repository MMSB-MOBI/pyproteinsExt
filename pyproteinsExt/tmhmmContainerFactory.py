import pyproteinsext.proteinContainer
import re
import gzip


def parse(inputFile=None):
    """
    Parse tmhmm output to create Container

    :param inputFile: path to tmhmm output
    :return: tmhmmContainerFactory.Container
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
            if l.startswith("#") and "Length" in l:
                if instance:
                    mainContainer = mainContainer.addParsing(Container(input=instance))
                instance = []
            instance.append(l.rstrip())
        mainContainer = mainContainer.addParsing(Container(input=instance))
    return mainContainer


def _parseBuffer(input):
    def parseFragments(list_fragments):
        fragmentsObj = []
        reLocation = re.compile("[\s]+([\d]+)[\s]+([\d]+)")
        for f in list_fragments:
            f_split = f.split("\t")
            cellular_location = f_split[2]
            start = f_split[3].split(" ")
            location = reLocation.findall(f_split[3])[0]
            start = int(location[0])
            end = int(location[1])
            fragment = Fragment(cellular_location, start, end)
            fragmentsObj.append(fragment)
        return fragmentsObj

    prot = input[0].split(" ")[1]
    nb_helix = [l for l in input if "Number of predicted TMHs" in l][0].split(":")[1].strip()
    prot_length = [l for l in input if "Length" in l][0].split(":")[1].strip()
    fragments = parseFragments([l for l in input if not l.startswith("#")])
    obj = TMHMM_Obj(prot, int(prot_length), int(nb_helix), fragments)
    return {prot: obj}


class Container(pyproteinsext.proteinContainer.Container):
    """ Iterable container, that stores TMHMM_Obj objects """
    def __init__(self, input=None):
        super().__init__(_parseBuffer, input)


class TMHMM_Obj():
    """
    TMHMM_Obj object stores informations about one TMHMM entry. 

    :param prot: Protein id
    :type prot: str
    :param prot_length: Protein length
    :type prot_length: int
    :param nb_helix: Number of predicted transmembrane helix by TMHMM
    :type nb_helix: int
    :param fragment: List of proteins fragments (inside, outside or helix)
    :type fragments: Fragment[]
    """

    def __init__(self, prot, prot_length, nb_helix, fragments):
        self.prot = prot
        self.prot_length = prot_length
        self.nb_helix = nb_helix
        self.fragments = fragments

    @property
    def topology_seq(self):
        """
        Protein sequence with 'o' for outside loop, 'i' for inside loop and helix number for helixes

        .. warning:: obsolete if we have more than 9 helixes
        """
        # /!\ WARNING : obsolete if we have more than 9 helixes !!! Change that.
        topology_seq = ["*"]*self.prot_length
        helix_number = 1
        for f in self.fragments:
            if f.cellular_location == 'inside':
                letter = "i"
            elif f.cellular_location == "outside":
                letter = "o"
            elif f.cellular_location == "TMhelix":
                letter = str(helix_number)
                helix_number += 1
            for i in range(f.start-1, f.end):
                topology_seq[i] = letter
        if "*" in topology_seq:
            print("WARNING : * in topology seq")
        return "".join(topology_seq)


class Fragment():
    """
    :param cellular_location: TMHMM cellular location of fragment (inside|outside|TMHelix)
    :type cellular_location: str
    :param start: fragment start position
    :type start: int
    :param end: fragment end position
    :type end: int
    """

    def __init__(self, cellular_location, start, end):
        self.cellular_location = cellular_location
        self.start = start
        self.end = end
