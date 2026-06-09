import pyproteins.container.customCollection
import pyproteins.container.Core
import re
from os.path import expanduser


refseqEntrySet = None


def isValidID(string):
    '''TO CHANGE'''
    return True


def strip(string):
    '''PROBABLY TO CHANGE'''
    subString = re.search("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string)
    if subString:
        return subString.group()


def getRefseqCollection():
    global refseqEntrySet
    if not refseqEntrySet:
        refseqEntrySet = EntrySet()

    return refseqEntrySet


class EntrySet(pyproteins.container.customCollection.EntrySet):
    def __init__(self, **kwargs):
        home = expanduser("~")
        cachePath = home
        super().__init__(collectionPath=cachePath, constructor=Entry,
                         typeCheck=isValidID, indexer=strip)


class Entry(pyproteins.container.Core.Container):
    def __init__(self, id, baseUrl="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", fileName=None):
        if not id:
            raise TypeError('identifier is empty')

        #url = baseUrl+"efetch.fcgi?db=nucleotide&id="+id+"&retmode=xml"
        url = baseUrl+"efetch.fcgi?db=protein&id="+id+"&retmode=xml"

        super().__init__(id, url=url, fileName=fileName)

        self.xmlHandler = self.getXmlHandler()
        if not self.xmlHandler:
            return None

        self.length = int(self.xmlHandler.find("GBSeq_length").text)
        self.parseFeatures()

    def __hash__(self):
        return hash(self.id)

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __eq__(self, other):
        return self.id == other.id

    def parseFeatures(self):
        self.Features = []
        for e in self.xmlHandler.find_all("GBFeature"):
            self.Features.append(Feature(e))

    def searchCDS(self, protein_id):
        features = [f for f in self.Features if f.type == "CDS" and f.qualifiers.get("protein_id", None) == protein_id]
        if len(features) == 1:
            return features[0]
        elif len(features) > 1:
            print("More than one feature. How to handle it ?")
            return None
        else:
            print("Protein not found")
            return None

    def getNeighborhood(self, feature, neighborhood_size):
        '''Get list of feature's neighbors, only CDS, from feature start - neighborhood size to feature end + neighborhood size. Only features completely
         contained in this limits are reported.'''
        def isNeighbor(location):
            if other_start_limit:
                if location.start >= other_start_limit:
                    return True

            if other_end_limit:
                if location.end <= other_end_limit:
                    return True

            if location.start >= start_limit and location.end <= end_limit:
                return True

            return False

        self.neighborhood = []
        if len(feature.locations) > 1:
            self.neighborhood.append("To check")
            return

        start = feature.locations[0].start
        end = feature.locations[0].end

        if start < neighborhood_size:
            start_limit = 0
            other_start_limit = self.length-(neighborhood_size-start)
        else:
            start_limit = start-neighborhood_size
            other_start_limit = None

        if end > self.length-neighborhood_size:
            end_limit = self.length
            other_end_limit = self.length-end
        else:
            end_limit = end+neighborhood_size
            other_end_limit = None

        for f in self.Features:
            if f == feature:
                continue
            if f.type != "CDS":
                continue
            neighbor = True
            for l in f.locations:
                if not isNeighbor(l):
                    neighbor = False

            if neighbor:
                self.neighborhood.append(f)


class Feature:
    def __init__(self, e):
        self.type = e.find("GBFeature_key").text
        self.location = e.find("GBFeature_location").text
        self.parseQualifiers(e)

    def parseQualifiers(self, e):
        self.qualifiers = {}
        for q in e.find_all("GBQualifier"):
            name = q.find("GBQualifier_name").text
            try:
                value = q.find("GBQualifier_value").text
            except AttributeError:
                value = ''
            self.qualifiers[name] = value


class Location:
    def __init__(self, start, end, strand, overlap0):
        self.start = start
        self.end = end
        self.strand = strand
        self.overlap0 = overlap0
