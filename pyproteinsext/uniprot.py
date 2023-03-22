from urllib.request import urlopen
import re
import pyproteins.container.customCollection
import pyproteins.container.Core
import pyproteinsext.pfam as pfam
import json

import gzip

from xml.etree.ElementTree import parse, dump, fromstring, register_namespace, ElementTree, tostring


from os.path import expanduser

PfamEntrySet = None
uniprotEntrySet = None

PfamCache = None

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
        home = expanduser("~")
        PfamCacheDefault = home
        PfamEntrySet = pfam.EntrySet(PfamCacheDefault)

    return PfamEntrySet

def proxySetting(**param):
    pyproteins.container.Core.proxySetting(param)

'''
    TODO Isoform, minimal -> affects the fasta sequence
                 Need to check isoform data xml structure, sequence variant specs of Uniprot
'''


def strip(string):
    subString = re.search(".xml$", string)
    if subString:
        return string.split(".")[0]

    return None

def capture(string):
    subString = re.search("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", string)
    if subString:
        return subString.group()

    return None

def isValidID(string):
    if re.match("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", string):
        return True
    if re.match("^[A-Z]{3}[0-9]{5}$", string):
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
        self.ns = '{http://uniprot.org/uniprot}'

        home = expanduser("~")
        cachePath = home     
        self.isXMLCollection = False
        self.isRedisCollection = False
        if 'streamXML' or 'collectionXML' in kwargs:
            self.isXMLCollection = True
            if 'streamXML' in kwargs:
                self.etree = parse(kwargs['streamXML'])
            elif 'collectionXML' in kwargs:
                if kwargs['collectionXML'].endswith('.gz'):
                    with gzip.open(kwargs['collectionXML'], 'r') as gz:
                        self.etree = parse(gz)
                else:
                    self.etree = parse(kwargs['collectionXML'])

            self.etree_root = self.etree.getroot()
            self.index = None
          #  print(f"==> {type(self.etree_root)} {type(self.etree)} <==")
        
        super().__init__(collectionPath=cachePath, constructor=Entry, typeCheck=isValidID, indexer=strip)
        if 'collectionXML' in kwargs:
            print (f"Acknowledged {len(self)} entries {kwargs['collectionXML']}")
      
    def __len__(self):
        if self.isXMLCollection:
            return len(list([ _ for _ in self ]))
        return super().__len__()

    def keys(self):
        """ Returns uniprot identifiers within collection as a generator """
        if self.isXMLCollection:
            for entry in self.etree_root.findall(f"{self.ns}entry"):
                for alt_acc in entry.findall(f"{self.ns}accession"):
                    yield(alt_acc.text)
        else :
            return super().keys()
    
    def has(self, id):
        if self.isXMLCollection:
            if self.index is None:
                self.index = {}
                for acc in self.keys():
                    self.index[acc] = True
            return id in self.index
        else :
            return super().has()

    def __iter__(self):
        if self.isXMLCollection:
            for entry in self.etree_root.findall(f"{self.ns}entry"):
                uniprotID = entry.find(f"{self.ns}accession").text
                yield Entry(uniprotID, xmlEtreeHandler = entry, xmlNS = self.ns)

    def get(self, uniprotID):
        if self.isXMLCollection:
            for entry in self.etree_root.findall(f"{self.ns}entry"):
                for acc in entry.findall(f"{self.ns}accession"):
                    if acc.text == uniprotID: # entry is the node matching provided UNIPROT accessor
                        return Entry(uniprotID, xmlEtreeHandler = entry, xmlNS = self.ns)
            return None
        else :
            print(f"Looking in XML files/dir collection for {uniprotID}")
            return super().get(uniprotID, xmlNS = self.ns)
        
        
    def serialize(self, ext=''):
        global PfamCache
        print ("serializing uniprot collection")
        super().serialize(ext=ext)
        if PfamCache:
            print ("serializing pfam collection")
            getPfamCollection().serialize(ext=ext)
    
    @property
    def taxids(self):
        taxids = set()
        for e in self:
            taxids.add(e.taxid)
        return list(taxids)

class Entry(pyproteins.container.Core.Container):
    ## Wrap and split kwargs
    # Re-encode etree xpath like search 
    # https://docs.python.org/3/library/xml.etree.elementtree.html#example
    # //div[@id='..' and @class='...]

    def __repr__(self):
        asStr = f"{self.id}:{self.AC}\n" 
        asStr += f"{self.name}:{self.fullName}({self.geneName})\n"
        
        if self.STRING_ID:
            asStr += f"STRING_ID:{self.STRING_ID}\n"
        
        asStr += f"taxid:{self.taxid}:{self.lineage}\n"
        asStr += f"KW:{self.KW}\n"
        asStr += f"GO:{self.GO}\n"
        
        return asStr

    def _buildXpath(self, xpath, **kwargs):
        
        _xpath_ = xpath.replace("/",f"/{self._ns}")
        _xpath_ += ' and '.join( [ f"[@{k}='{v}']" for k,v in kwargs.items() ] )

        if not ( _xpath_.startswith('/') or _xpath_.startswith('./') ) :
            _xpath_ = f"{self._ns}{_xpath_}"

        #print(f"building xpath from {xpath} to {_xpath_}")
        return _xpath_

    def _xmlAttrib(self, key, elem=None):
        if elem is None:
            elem = self.xmlHandler
        return elem.attrib[key]

    def _xmlFind(self, tag, elem=None, **kwargs):
        if elem is None:
            elem = self.xmlHandler
        if len(kwargs.keys()) > 1:
            raise(f"Asking for Too many attributes ({len(kwargs.keys())}) in tag xpath search")

        xpathStr = self._buildXpath(tag, **kwargs) 
        return elem.find(xpathStr)
        
    def _xmlFindAll(self, tag, elem=None, **kwargs):
        if elem is None:
            elem = self.xmlHandler
        
        xpath = self._buildXpath(tag, **kwargs) 
        return elem.findall(xpath)

    def __init__(self, id, baseUrl="http://www.uniprot.org/uniprot/", fetchable= True, fileName=None, xmlEtreeHandler=None, xmlNS=None):
        if not id:
            raise TypeError('identifier is empty')
        super().__init__(id, url=baseUrl + str(id) + '.xml', fileName=fileName)
        #pyproteins.container.Core.Container.__init__(self, id, url=baseUrl + str(id) + '.xml', fileName=fileName)
        
        if not xmlNS is None:
            self._ns = xmlNS

        if not xmlEtreeHandler is None:
            self.xmlHandler = xmlEtreeHandler
        else:
            print("Search for", f"{xmlNS}entry")
            self.xmlHandler = self.getXmlHandler(fetchable=fetchable).find(f"{xmlNS}entry")
        
        if self.xmlHandler is None:
            return None
        
        self.name = self._xmlFind("./name").text
        
        if not self._xmlFind("./protein/recommendedName/fullName") is None:
            self.fullName = self._xmlFind("./protein/recommendedName/fullName").text
        else:
            self.fullName = self.name
        
        self.geneName =  None
        e = self._xmlFind("gene")
        if not e is None:
           self.geneName = self._xmlFind("./name", elem=e).text

        self.STRING_ID = None
        e = self._xmlFind("./dbReference", type="STRING")
        if not e is None:
            self.STRING_ID = self._xmlAttrib('id', elem=e)
        
        self.parseAC()
        self.parseLineage()
        self.parseKW()
        self.parseGO()
        self.parseSequence()
        self.parse_subcellular_location()


# Following oarsing stages are disabled since we got rid of bs4
# Need to port them  to lxml, making use of xpath syntax
    def PARSER_TO_PORT_TO_ETREE(self):

       
        self.parseSse()
        self.Ensembl = self.parseEnsembl()
        self.GeneID = self.parseGeneID()
        self.parseSequence()
        self.parsePDB()
        self.parseMIM()
        self.parseDI()
        self.parseORPHA()
        self.xref = self.get_xref()
        self.parseInterpro()

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
            if k == 'taxid':
                container[k] = v
        return container

    def parseAC(self):
        self.AC = [ e.text for e in  self._xmlFindAll("./accession") ]
          
    def parseLineage(self):
        self.lineage = [ e.text for e in  self._xmlFindAll("./organism/lineage/taxon") ]
        e = self._xmlFind('organism/dbReference', type='NCBI Taxonomy')
        self.taxid = self._xmlAttrib('id', elem=e)
    
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
        for e in self._xmlFindAll("dbReference", type="GO"):
            self.GO.append(GoKW(
                self._xmlAttrib('id', elem=e),
                self._xmlAttrib('value', elem = self._xmlFind("property", elem=e, type="term") ),
                self._xmlAttrib('value', elem = self._xmlFind("property", elem=e, type="evidence") )
            ))
        
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
        try:
            self.domains=getPfamCollection().map(uniprotID=self.id)
        except:
            self.domains=[]

        #self.domains = []
        #for e in self.xmlHandler.find_all("feature", type="domain"):
        #    buf = Domain(e, self.id)
        #    if buf.description:
        #        self.domains.append(Domain(e, self.id))
        #for e in self.xmlHandler.find_all("feature", type="repeat"):
        #    buf = Domain(e, self.id)
        #    if buf.description:
        #        self.domains.append(Domain(e, self.id))
        #if not self.domains:
        #    print "No domain data found for " + self.id + ", attempting pfam"
        #    try :
        #        self.domains = getPfamCollection().map(uniprotID=self.id)
        #    except ValueError as msg:
        #        print ("Could not bind uniprot to its pfam ressources reason\n" + str(msg))

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
        
        for e in self._xmlFindAll("./keyword"):
            self.KW.append( UniprotKW(self._xmlAttrib('id', elem=e), e.text) ) 

    def parseSequence(self):
        self.sequence = self._xmlFind("./sequence").text.replace("\n", "").replace(" ", "")
        #self.sequence = Sequence(self.xmlHandler.find("sequence", {"length" : True}))
    #    pass
    
    def parseEnsembl(self):
        #Search for Ensembl id : ENSGXXXXXXXXXXX
        Ensembl_id = []
        for e in self.xmlHandler.find_all("dbReference", type="Ensembl"):
            for e_gene_id in e.find_all('property',type='gene ID'):
                if e_gene_id["value"] not in Ensembl_id:
                    Ensembl_id.append(e_gene_id["value"])
        return Ensembl_id
        
    def parseGeneID(self):
        GeneID= []
        for e_gene_id in self.xmlHandler.find_all("dbReference", type="GeneID"):
            if e_gene_id["id"] not in GeneID:
                GeneID.append(e_gene_id["id"])
        return GeneID

    def parse_subcellular_location(self):
        self.subcellular_location = []
        for citation in self._xmlFindAll('comment', type="subcellular location"):
            for sl in self._xmlFindAll('subcellularLocation', elem = citation):
                for location in self._xmlFindAll('location', elem = sl):
                    self.subcellular_location.append(location.text)
        # except : 
        #     print("WARNING", self.AC, "subcellular location can't be parsed")
        #     print(citation, sl, location)

        # self.subcellular_location = []
        # citations = self._xmlFindAll("comment", type="subcellular location")
        # for citation in citations:
        #     annotation = self._xmlFind('subcellularLocation/location', elem=citation)
           

        #     # if not annotation:
        #     #     print("no annotation", self.AC, annotation)
        #     # else:
        #     try:
        #         self.subcellular_location.append(annotation.text)
        #     except:
        #         print("WARNING", self.AC, "subcellular location can't be parsed")
        #         self.subcellular_location = ['WARNING']
        
    def get_xref(self):
        # print("GET XREF")
        dic_xref = {'EMBL': {}, 'RefSeq': {}}
        # Search EMBL
        for e in self.xmlHandler.find_all("dbReference", type="EMBL"):
            if str(e.parent.name) == 'entry':
                for e_prot_id in e.find_all('property',type='protein sequence ID'):
                    dic_xref["EMBL"][e["id"]] = e_prot_id["value"]
        # Search RefSeq
        for e in self.xmlHandler.find_all("dbReference", type="RefSeq"):
            if str(e.parent.name) == 'entry':
                for e_prot_id in e.find_all('property',type='nucleotide sequence ID'):
                    dic_xref["RefSeq"][e_prot_id["value"]] = e["id"]
        return dic_xref

    def parseInterpro(self):
        self.Interpro = []
        for e in self.xmlHandler.find_all("dbReference", type = "InterPro"):
            for name in e.find_all("property", type = "entry name"):
                self.Interpro.append(Interpro(e))

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
    
    @property
    def isGOannot(self):
        return not len(self.GO) == 0
        
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
        self.mass =e._xmlAttrib('mass')
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
    def __init__(self, id,term, evidence):
        self.id = id
        self.term = term
        self.evidence = evidence
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
    def __init__(self, id, text):
        self.id = id
        self.term = text
    def __repr__(self):
        return self.id + ":" + self.term     

class Interpro():
    def __init__(self,e):
        self.id = e["id"]
        self.getName(e)

    def __repr__(self):
        return self.id + ":" + self.name

    def getName(self,e): 
        for name in e.find_all("property", type = "entry name"):
           self.name = name["value"]


class Genome():
    def __init__(self,xmlHandler):
        self.searchEMBL(xmlHandler)
        self.searchRefSeq(xmlHandler)

    def searchEMBL(self,xmlHandler):
        self.EMBLRef=[]
        self.EMBLProteinRef=[]
        for e in xmlHandler.find_all("dbReference", type="EMBL"):
            if str(e.parent.name) == 'entry':
                self.EMBLRef.append(e['id'])
                for e_prot_id in e.find_all('property',type='protein sequence ID'):
                    self.EMBLProteinRef.append(e_prot_id['value'])

    def searchRefSeq(self,xmlHandler):
        self.RefSeqRef=[]
        self.RefSeqProteinRef=[]
        for e in xmlHandler.find_all("dbReference", type="RefSeq"):
            if str(e.parent.name) == 'entry':
                self.RefSeqProteinRef.append(e['id'])
                for e_prot_id in e.find_all('property',type='nucleotide sequence ID'):
                    self.RefSeqRef.append(e_prot_id['value'])
