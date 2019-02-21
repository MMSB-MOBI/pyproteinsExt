#import urllib2
from urllib.request import urlopen
from bs4 import BeautifulSoup
import re
import pyproteins.container.customCollection
import pyproteins.container.Core
import pyproteinsExt.pfam as pfam
import json

from os.path import expanduser

PfamEntrySet = None
uniprotEntrySet = None

from json import JSONEncoder

import pyproteins.utils.make_json_serializable
#def _default(self, obj):
#    return getattr(obj.__class__, "toJSON", _default.default)(obj)
#
#_default.default = JSONEncoder.default
#JSONEncoder.default = _default


## Returns a list of all keyword in collection, along with the list of uniprotObj featuring them
def keyWordChart(uniprotObjIter, kwType='GO'):
    def kwMapper(obj, _type) : 
        if _type == 'GO':
            return obj.GO
        raise TypeError("implement other KW plz")
    
    kwChart = {}
    for uniprotObj in uniprotObjIter:
        for kwObj in  kwMapper(uniprotObj, kwType):
            if kwObj not in kwChart:
                kwChart[kwObj] = []
            kwChart[kwObj].append(uniprotObj)
    
    return sorted([ (k,v) for k,v in kwChart.items() ], key=lambda x : len(x[1]), reverse=True)

# Give link to uniprot Collection to allow proxy settings
# and cache setting for it and pfam
def getUniprotCollection ():
    global uniprotEntrySet
    if not uniprotEntrySet:
        uniprotEntrySet = EntrySet()

    return uniprotEntrySet

#def setCache(location):
#    print location
#    uniprotEntrySet.setCache(location=location)

#def proxySetting(**kwargs):
#    proxySetting(**kwargs)

def getPfamCollection ():
    global PfamEntrySet
    if not PfamEntrySet:
        PfamEntrySet = pfam.EntrySet()

    return PfamEntrySet

def proxySetting(**param):
    pyproteins.container.Core.proxySetting(param)

'''

    TODO Isoform, minimal -> affects the fasta sequence
                 Need to check isoform data xml structure, sequence variant specs of Uniprot

'''


def strip(string):
    subString = re.search("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string)
    if subString:
        return subString.group()

    return None

def capture(string):
    subString = re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string)
    if subString:
        return subString.group()

    return None

def isValidID(string):
    if re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", string):
        return True
    return False

class fetchEntries():
    def __init__(self, list):
        self.list = list
        self.entrySet = [Entry(id = id) for id in self.list]
    def __iter__(self):
        for entry in self.entrySet:
            yield entry



''' Collection of uniprot entries
customary cache location is "/Users/guillaumelaunay/work/data/uniprot"
'''
class EntrySet(pyproteins.container.customCollection.EntrySet):
    def __init__(self, **kwargs):
        home = expanduser("~")
        cachePath = home
        if 'collectionPath' in kwargs:
            cachePath = kwargs['collectionPath']

        super().__init__(collectionPath=cachePath, constructor=Entry, typeCheck=isValidID, indexer=strip)

    def serialize(self, **kwargs):
        print ("serializing uniprot collection")
        super().serialize(kwargs)
        print ("serializing pfam collection")
        getPfamCollection().serialize(kwargs)

class Entry(pyproteins.container.Core.Container):
    def __init__(self, id, baseUrl="http://www.uniprot.org/uniprot/", fileName=None):
        if not id:
            raise TypeError('identifier is empty')
        c_id = capture(id)
        if not c_id:
            raise ValueError('could not extract uniprot identifier from provided id parameter ' + id)
        super().__init__(c_id, url=baseUrl + str(c_id) + '.xml', fileName=fileName)
        #pyproteins.container.Core.Container.__init__(self, id, url=baseUrl + str(id) + '.xml', fileName=fileName)

        #print id + '-->' + str(fileName)

        self.xmlHandler = self.getXmlHandler()
        if not self.xmlHandler:
            return None
        
        self.name = self.xmlHandler.find('name').text
        self.fullName = self.xmlHandler.find('fullName').text
        self.geneName = self.xmlHandler.find('gene').find('name').text if self.xmlHandler.find('gene') else None
        self.parseKW()
        self.parseGO()
        self.parseLineage()
        self.parseAC()
        self.parseDomain()
        self.parseSse()
        self.parseSequence()
        self.parsePDB()
        self.parseMIM()
        self.parseDI()
        self.parseORPHA()

    def __hash__(self):
        return hash(self.id)

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __eq__ (self, other):
        return self.id == other.id

    def toJSON(self):
        #print ('toJSON')
        #asDict = {}
        #for k in self.__dict__.keys():
        #    if k == 'xmlHandler':
        ##        continue
        #    if k != 'GO':
        #        continue
        #    asDict[str(k)] = getattr(self, k)
        #    print(asDict)

        container = {}
        for k, v in self.__dict__.items():
            if k == 'name':
                container[k] = v
            if k == 'GO':
                container[k] = [ go.__dict__ for go in v ]
            if k == 'id':
                container[k] = v
            if k == 'geneName':
                container[k] = v
            if k == 'fullName':
                container[k] = v
        return container

    def parseAC(self):
        self.AC = [e.text for e in self.xmlHandler.find_all('accession')]

    def parseLineage(self):
        self.lineage = [ e.text for e in self.xmlHandler.find('lineage').find_all('taxon')]
        self.taxid = self.xmlHandler.find('organism').find_all('dbReference', type='NCBI Taxonomy')[0]['id']
    def parseMIM(self):
        self.MIM = []
        for e in self.xmlHandler.find_all("dbReference", type="MIM"):
            if str(e.parent.name) == 'entry':
                self.MIM.append(MimKW(e))

    def parseDI(self):
        self.DI = []
        for e in self.xmlHandler.find_all("disease"):
            self.DI.append(DI(e))

    def parseGO(self):
        self.GO = []
        for e in self.xmlHandler.find_all("dbReference", type="GO"):
            self.GO.append(GoKW(e))

    def parseORPHA(self):
        self.ORPHA = []
        for e in self.xmlHandler.find_all("dbReference", type="Orphanet"):
            if str(e.parent.name) == 'entry':
                self.ORPHA.append(OrphaKW(e))

    def parsePDB(self):
        self.pdbRef = []
        for e in self.xmlHandler.find_all("dbReference", type="PDB"):
            self.pdbRef.append(PDBref(e))

    def parseDomain(self):
        self.domains = []
        for e in self.xmlHandler.find_all("feature", type="domain"):
            buf = Domain(e, self.id)
            if buf.description:
                self.domains.append(Domain(e, self.id))
        for e in self.xmlHandler.find_all("feature", type="repeat"):
            buf = Domain(e, self.id)
            if buf.description:
                self.domains.append(Domain(e, self.id))
        if not self.domains:
        #    print "No domain data found for " + self.id + ", attempting pfam"
            try :
                self.domains = getPfamCollection().map(uniprotID=self.id)
            except ValueError as msg:
                print ("Could not bind uniprot to its pfam ressources reason\n" + str(msg))

    def parseSse(self):
        self.sse = []
        for e in self.xmlHandler.find_all("feature", type="strand"):
            self.sse.append(Sse(e))
        for e in self.xmlHandler.find_all("feature", type="turn"):
            self.sse.append(Sse(e))
        for e in self.xmlHandler.find_all("feature", type="helix"):
            self.sse.append(Sse(e))

    def parseKW(self):
        self.KW = []
        for e in self.xmlHandler.find_all("keyword"):
            self.KW.append(UniprotKW(e))
    def parseSequence(self):
        self.sequence = Sequence(self.xmlHandler.find("sequence", {"length" : True}))
    #    pass


    @property
    def fasta(self):
        return '>' + str(self.id) + ' ' + str(self.name) + '\n' + str(self.sequence)

    def peptideSeed(self):
        return {
            'id' : self.id,
            'desc' : self.name,
            'seq' : str(self.sequence)
        }



        #self.parseIsoform()

    def hasKW(self, keyword):
        if keyword.upper() in (kw.id.upper() for kw in self.KW):
            return True
        return False

    def hasGO(self, keyword):
        if keyword.upper() in (kw.id.upper() for kw in self.GO):
            return True
        return False

    def hasMIM(self, keyword):
        if keyword.upper() in (kw.id.upper() for kw in self.MIM):
            return True
        return False

    def hasORPHA(self, keyword):
        if keyword.upper() in (kw.id.upper() for kw in self.ORPHA):
            return True
        return False

    def hasDI(self, keyword):
        if keyword.upper() in (kw.id.upper() for kw in self.DI):
            return True
        return False

    # returns a position information container
    def pos(self, i, lookup=False):
        if  i < 1 or i >= len(self.sequence):
            raise IndexError("sequence index \""  + str(i) + "\" is out of range (" + str(len(self.sequence)) + ")")
        d = self._domainFlyCast(position=i)
        if len(d) > 1:
            raise ValueError("\"" + self.id + "\" lays more than one domain \"" + str(len(d)) +
                             "\" at position " + str(i)  )
        domain="SomeDefault"
        if d:
            domain = d[0]
        else:
            l = self._domainFlyCast(NterLookup=i)
            r = self._domainFlyCast(CterLookup=i)
            if l and r:
                domain = { 'outOfBounds' : 'hinge' }
            elif l:
                domain = { 'outOfBounds' : 'Cter' }
            elif r:
                domain = { 'outOfBounds' : 'Nter' }
            else:
                print ("Not a hinge a Cter or a Nter at " + str(i) + " in " + self.id)
        return Position(i, domain=domain, sse=self._getSse(i), letter=self.sequence[i])


    def _domainFlyCast(self, **kwargs):
        if not self.domains:
            raise ValueError("No domain data available for " + self.id)
            #return None
        if 'position' in kwargs:
            results = [ d._dict for d in self.domains if d.owns(kwargs['position']) ]
            return results
        if 'NterLookup' in kwargs:
            for i in range (int(kwargs['NterLookup']), 0, -1):
                d = self._domainFlyCast(position=i)
                if d:
                    return d#d[0]
        if 'CterLookup' in kwargs:
            for i in range (int(kwargs['CterLookup']), len(self.sequence) + 1):
                d = self._domainFlyCast(position=i)
                if d:
                    return d#d[0]


        return None

    def _getSse(self, position):
        i = int(position)
        if not self.sse:
            return None
        for d in self.sse:
            if d.begin <= i <= d.end:
                return d.type
        return "coil"

# Custom encoder for uniprot entity
class EntryEncoder(json.JSONEncoder):
    def default(self, entryObj):
        if isinstance(entryObj, pyproteins.container.Core.Container):
            container = {}
            for k, v in entryObj.__dict__.items():
                if k == 'name':
                    container[k] = v
                if k == 'GO':
                    container[k] = [ go.__dict__ for go in v ]
                if k == 'id':
                    container[k] = v
                if k == 'geneName':
                    container[k] = v
                if k == 'fullName':
                    container[k] = v
            return container
        # Error
        return json.JSONEncoder.default(self, entryObj)



class Position ():
    def __init__(self, number, **kwargs):
        self.number = number
        self.domain = kwargs['domain'] if 'domain' in kwargs else None

        self.sse = kwargs['sse'] if 'sse' in kwargs else None
        self.letter = kwargs['letter'] if 'letter' in kwargs else None

    def __repr__(self):
        return self.letter + str(self.number) + "(" + str(self.sse) + ") : " + str(self.domain)

class Sequence():
    def __init__(self, e):
        if not e:
            raise ValueError("Cant parse sequence")
        self.string = e.text.replace("\n", "").replace(" ", "")
        self.mass =e['mass']
    def __len__(self):
        return len(self.string)
    def __getitem__(self, i):
        if isinstance(i, slice):
            newSlice = slice(i.start - 1, i.stop - 1)
            return self.string[newSlice]
        else:
            if i < 1 or i > len(self):
                raise IndexError(str(i) + " out of range [1-" + self(len(self)) + "]")
            return self.string[i - 1]



    def __repr__(self):
        return self.string


class Domain(object):
    def __init__(self, xmlHandler, id):
        self.begin = [int(e['position']) for e in xmlHandler.find_all('begin')]
        self.end = [int(e['position']) for e in xmlHandler.find_all('end')]
        if len(self.begin) != len(self.end):
            raise ValueError("Number of end/stop differs " + str(self.begin) + "/" + str(self.end))
        if 'description' not in xmlHandler:
            #print "Warning Domain w/in \"" + id + "\" does not feature any description"
            self.description = None
        else:
            self.description = xmlHandler['description']
        self.carriedBy = id

    def __eq__(a, b):
        if a.id != b.id:
            return False
        if a.begin != b.begin:
            return False
        if a.end != b.end:
            return False
        if a.description != b.description:
            return False

        return True

    def __repr__(self):
        b = [ str (self.begin[i]) + "-" + str(self.end[i]) for i,e  in enumerate (self.begin) ]
        return "\"" + self.description + "\"\t" + ",".join(b) + "\n"

    def owns(self, position):
        try :
            j = int(position)
        except ValueError as err:
            print ("Improper amino acid position \"" + str(position) + "\"")
            return False

        for i, e  in enumerate (self.begin):
            if self.begin[i] <= j <= self.end[i]:
                return True

        return False



    @property
    def _dict(self):
        return self.__dict__

class Sse():
    def __init__(self, e):
        self.begin = int(e.find('begin')['position'])
        self.end =  int(e.find('end')['position'])
        self.type = e['type']
    def __repr__(self):
        return "SSE:" + self.type + " " + str(self.begin) + "-" + str(self.end)


class annotTerm:
    def __init__(self):
        pass
    def __hash__(self):
        return hash(str(self.id))
    def __eq__(self, other):
        return hash(self) == hash(other)

class GoKW(annotTerm):
    def __init__(self, e):
        self.id = e['id']
        self.term = e.find('property', type='term')['value']
        self.evidence = e.find('property',type='evidence')['value']
    def __repr__(self):
        return self.id + ":" + self.term + "{" + self.evidence + "}"

class MimKW(annotTerm):
    def __init__(self, e):
        self.id = e['id']
        self.value = e.find('property', type='type')['value']
    def __repr__(self):
        return self.id + ":" + self.value

class OrphaKW(annotTerm):
    def __init__(self, e):
        self.id = e['id']
        self.type = e.find("property")['type']
        self.value = e.find('property')['value']

    def __repr__(self):
        return self.id + ": (" + self.type + ")" + self.value

class DI(annotTerm):
    def __init__(self, e):
        self.id = e['id']
        self.name = e.find('name').string
        acronym = e.find('acronym')
        self.acronym = acronym.string if acronym else 'NA'
        self.description = e.find('description').string
    def __repr__(self):
        return self.id + ":" + self.name + " (" + self.acronym + ") {" + self.description + "}"



class PDBref():
    def __init__(self, e):
        self.id = e['id']
        self.method = e.find('property', type='method')['value'] if e.find('property', type='method') else None
        self.resolution = e.find('property',type='resolution')['value'] if e.find('property', type='resolution') else None
        self.chains = e.find('property',type='chains')['value'] if e.find('property', type='chains') else None

    def __repr__(self):
        string = "PDB_id: " + self.id
        elements = []
        if self.method:
            elements.append("method: " + self.method)
        if self.resolution :
            elements.append("resolution: " + self.resolution)
        if self.chains:
            elements.append("chains: " + self.chains)
        if elements:
            return string + ':{' + ','.join(elements)  + '}'

        return string

class UniprotKW():
    def __init__(self, e):
        self.id = e['id']
        self.term = e.text
    def __repr__(self):
        return self.id + ":" + self.term

