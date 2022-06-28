#import urllib2
from bs4 import BeautifulSoup
import re
import xml.etree.ElementTree as ET
import json


BIOGRID_KEY="8877654116411665fc885b0fc3014cd2"
BIOGRID_UNIPROT_MAPPER_URL="http://thebiogrid.org/downloads/archives/External%20Database%20Builds/UNIPROT.tab.txt"
'''
BIOGRID_ORDERED_JSON_KEYS=["BIOGRID_INTERACTION_ID", "ENTREZ_GENE_A", "ENTREZ_GENE_B",
"BIOGRID_ID_A", "BIOGRID_ID_B", "SYSTEMATIC_NAME_A", "SYSTEMATIC_NAME_B", "OFFICIAL_SYMBOL_A",
"OFFICIAL_SYMBOL_B", "SYNONYMS_A", "SYNONYMS_B", "EXPERIMENTAL_SYSTEM", "EXPERIMENTAL_SYSTEM_TYPE",
"PUBMED_AUTHOR", "PUBMED_ID", "ORGANISM_A", "ORGANISM_B", "THROUGHPUT", "QUANTITATION",
"MODIFICATION", "PHENOTYPES", "QUALIFICATIONS", "TAGS", "SOURCEDB"]
# Fields description at "http://wiki.thebiogrid.org/doku.php/downloads:biogrid_json"
'''
BIOGRID_ORDERED_JSON_KEYS=["BIOGRID_ID_A", "BIOGRID_ID_B", "ENTREZ_GENE_A", "ENTREZ_GENE_B",
 "SYSTEMATIC_NAME_A", "SYSTEMATIC_NAME_B", "OFFICIAL_SYMBOL_A",
"OFFICIAL_SYMBOL_B", "SYNONYMS_A", "SYNONYMS_B", "EXPERIMENTAL_SYSTEM", "EXPERIMENTAL_SYSTEM_TYPE",
"PUBMED_AUTHOR", "PUBMED_ID", "ORGANISM_A", "ORGANISM_B", "THROUGHPUT", "QUANTITATION",
"MODIFICATION", "PHENOTYPES", "QUALIFICATIONS", "TAGS", "SOURCEDB"]


class BIOGRID_DATUM(object):
    def __init__(self, id, data, mapper):
        self.biogridMapper=mapper
        self.data = data
        self.biogridID = id

    def __repr__(self):
        string = ''
        for field in BIOGRID_ORDERED_JSON_KEYS:
            v = self.data[field] if field in self.data else '-'
            if field.startswith("BIOGRID_ID_"):
                w = self.biogridMapper(biogridId=v)
                v = v if not w else "uniprotkb:" + w
            string = string + str(v) + "\t"
        string += str(self.biogridID)
        return string

    @property
    def species(self):
        return (self.data['ORGANISM_A'], self.data['ORGANISM_B'])

    @property
    def uniprotPair(self):
        pp = self.interactors
        if pp[0][0] == "uniprokb:" and pp[1][0] == "uniprokb:":
            a = pp[0][1]
            b = pp[1][1]
            (a,b) = (b,a) if b < a else (a,b)
            return (a, b)
        return None

    @property
    def interactors(self):

        x = self.data['BIOGRID_ID_A']
        y = self.data['BIOGRID_ID_B']
        u = self.biogridMapper(biogridId=x)
        v = self.biogridMapper(biogridId=y)

        if u:
            x = ("uniprokb:", u)
        else:
            x = ("biogrid:", x)

        if v:
            y = ("uniprokb:", v)
        else:
            y = ("biogrid:", y)
        d = ([x],[y]) ## Overkill but required for ducktyping w/ psq.interactors
        #print d
        #return (x, y)
        return d

class BIOGRID(object):
    def __init__(self, **kwargs):  #mapperUrl=BIOGRID_UNIPROT_MAPPER_URL, uniprotMapFile=None
        self.biogridMapper = BIOGRIDMAPPER()
        self.loadBiogridMapper(**kwargs)
        self.data = {}

    def __iter__(self):
        for k in self.data:
            yield BIOGRID_DATUM(k, self.data[k], self.biogridMapper)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return "\n".join( [str(biogridDatum) for biogridDatum in self ] )


    def getBiomolecules(self, _type="uniprot"):

        moleculeList = []
        if _type == "uniprot":
            #for datum in self:
                #print datum
            #    b = datum.interactors
                #print b

            return list(set([ mol[0][1] for datum in self for mol in datum.interactors if mol[0][0] == "uniprokb:" ]))

    def readFile(self, fileName):
        data = ''

        with open(fileName, 'r') as f:
            data = f.read()

        self.load(data)

    # biogridordered keys is the output sequence of dump method
    def load(self, stream, type="biogridOrderedKeys"):
        if not stream:
            print("You must provide a mitab input")
            return
        bufferStr = []

        for line in stream.split("\n"):
            bufferStr.append(line)
        print(len(bufferStr))
        self.tsvBiogridParser(bufferStr)

    def tsvBiogridParser(self, iBuffer):
        for rec in iBuffer:
            #print "----->"  + rec
            record = rec.split("\t")
            key = record.pop()
            self.data[key] = {} # recover the biogrid interaction identifier as primary key
            print(key)
            for l, y in zip(BIOGRID_ORDERED_JSON_KEYS, record):
             #   print '-->' + y

                if y.startswith('uniprotkb:'):
                    m = re.search('^uniprotkb:(.+)$', y)
                    if m:
                        y = self.biogridMapper(uniprotId = m.group(1))
                self.data[key][str(l)] = y #if y != '-' else None

    def clear(self):
        self.data = {}


    def zQuery(self, **kwargs):
        if 'uniprotA' in kwargs and 'uniprotB' in kwargs:
            hArray =(kwargs['uniprotA'], kwargs['uniprotB'])

            cloneOne = BIOGRID()
            chunks = [hArray[0][x:x+50] for x in xrange(0, len(hArray[0]), 50)]
            for p1 in chunks:
                cloneOne.uniprotQuery(p1)

            cloneTwo = BIOGRID()
            chunks = [hArray[1][x:x+50] for x in xrange(0, len(hArray[1]), 50)]
            for p2 in chunks:
                cloneTwo.uniprotQuery(p2)

            self.data = { x:cloneOne.data[x] for x in cloneOne.data if x in cloneTwo.data }

    def dump(self, file=None):
        if file:
            with open(file, 'w') as f:
                f.write("#" + "\t".join(BIOGRID_ORDERED_JSON_KEYS) + "\n")
                f.write(self.__repr__())
        else:
            return self.__repr__()
    def query(self, uniprotId=None, geneIdA=None, geneIdB=None, species=None):
        if uniprotId:
            self.uniprotQuery(uniprotId)
        elif geneIdA and geneIdB:
            self.genePairQuery(geneIdA, geneIdB)
        elif species:
           self.specieQuery(species)

    def specieQuery(self, taxid):
        # 1st check if organism is valid endopoit
        url = "http://webservice.thebiogrid.org/organisms/?format=json&accesskey=" + BIOGRID_KEY
        taxonDict = self._urlToDict(url)
        taxQueryList = []
        taxid = str(taxid)
        # Treat provide taxid 1st as a ncbi taxon identifier, then as a substring of a litteral name

        if taxid in taxonDict:
            taxQueryList.append(taxid)
        else:
            taxQueryList = [ tx for tx in taxonDict if taxid in taxonDict[tx] ]

        if not taxQueryList:
            print (taxid + ' is not a supported Biogrid specie')
            return
        else :
            print (taxQueryList)
        # get the stuff

        for tx in taxQueryList :
            url = "http://webservice.thebiogrid.org/interactions?taxid=" + str(tx) + "&format=json&accesskey=" + BIOGRID_KEY
            data = self._urlToDict(url)
            #print url
            #print data
            self.data.update(data)

    def genePairQuery(self, geneIdA, geneIdB):
        url = "http://webservice.thebiogrid.org/interactions/?geneList="
        url = url + geneIdA + "|" + geneIdB + "&format=json&accessKey=" + BIOGRID_KEY
        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError as error:
            print (url + "\nHTTP ERROR " + str(error.code))
            return None
        except urllib2.URLError as error:
            print (url + "\n" + str(error.reason))
            return None
        raw = response.read()
        response.close()
        self._parse(raw)
        self._filter(geneIdA, geneIdB)
        #http://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=cdc27|apc1|apc2&evidenceList=Affinity Capture-MS|Two-hybrid&accesskey=[ACCESSKEY]


    # Returns another instance of BG with filter data elements
    def filter(self, **kwargs):
        target = BIOGRID()
        buf = kwargs['uniprot']
        # For now we dont look into alias uniprot identifiers
        if 'uniprot' in kwargs:
            if isinstance(kwargs['uniprot'], list):
                buf = set(kwargs['uniprot'])
            elif isinstance(kwargs['uniprot'], basestring):
                buf = set([kwargs['uniprot']])
#            elif isinstance(kwargs['uniprot'], set):

            for bgData in self:
                up = bgData.uniprotPair
                if not up:
                    continue
                if set(up) & buf :
                    #print str(up) + ' <===> ' + str(buf)
                    target.data.update({bgData.biogridID : self.data[bgData.biogridID]})

        return target



    def _filter(self, idA, idB):
        dataBuffer = {}
        for interKey in self.data:
            datum = self.data[interKey]
            #print self.data[interKey]
            if idA not in [datum["SYNONYMS_A"], datum["SYNONYMS_B"], datum["SYNONYMS_A"], datum["SYSTEMATIC_NAME_B"], datum["SYSTEMATIC_NAME_A"], datum["OFFICIAL_SYMBOL_A"], datum["OFFICIAL_SYMBOL_B"]]:
                continue
            if idB not in [datum["SYNONYMS_A"], datum["SYNONYMS_B"], datum["SYNONYMS_A"], datum["SYSTEMATIC_NAME_B"], datum["SYSTEMATIC_NAME_A"], datum["OFFICIAL_SYMBOL_A"], datum["OFFICIAL_SYMBOL_B"]]:
                continue
            dataBuffer[interKey] = datum
        self.data = dataBuffer

    def uniprotQuery(self, uniprotId,includeInteractors=True):
        url = "http://webservice.thebiogrid.org/interactions/?geneList="
        qString = None
        if isinstance(uniprotId, list):
            buf = []
            for x in uniprotId:
                biogridId = self.biogridMapper(uniprotId=x)
                if not biogridId:
                    continue
                buf.append(biogridId)
            qString = '|'.join(buf)
        else:
            qString = self.biogridMapper(uniprotId=uniprotId)

        if not qString:
                return

        incBool = 'true' if includeInteractors else 'false'
        url = url + qString + "&searchbiogridids=true&includeInteractors=" + incBool + "&format=json&accessKey=" + BIOGRID_KEY

        data = self._urlToDict(url)
        self.data.update(data)

    def _urlToDict(self,url):

        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError as error:
            print (url + "\nHTTP ERROR " + str(error.code))
            return None
        except urllib2.URLError as error:
            print (url + "\n" + str(error.reason))
            return None
        raw = response.read()
        response.close()

        buf = json.loads(raw)
        if isinstance(buf, list):
            if len(buf) == 0:
                return {}
            raise TypeError('unexpected list on json parsing' + str(buf))
        # We expected dict w/ unique key per interactions
        return buf

    def loadBiogridMapper(self, uniprotMapFile=None):

        if uniprotMapFile:
            with open(uniprotMapFile,'r') as f:
                print ('reading uniprot mapping from ' + uniprotMapFile)
                data = f.read()
                self._loadMapper(data)
            return

        print ('reading uniprot mapping from ' + BIOGRID_UNIPROT_MAPPER_URL)
        try:
            response = urllib2.urlopen(BIOGRID_UNIPROT_MAPPER_URL)

        except urllib2.HTTPError as error:
            print ("loadBioGridMapper url resolution failed")
            print ("HTTP ERROR " + str(error.code))
            return None
        except urllib2.URLError as error:
            print ("loadBioGridMapper url resolution failed")
            print (error.reason)
            return None
        raw = response.read()
        response.close()
        self._loadMapper(raw)


    def _loadMapper(self, stream):
        self.biogridMapper.load(stream)
    def analyse(self):
        if len(self) == 0: return None
        container = { "pmids" : [], "experiments" : [], "experimentTypes" : [],
                     "highThroughput" : False, "lowThroughput" : False }
        for k, record in self.data.iteritems():
            if record["PUBMED_ID"] not in container["pmids"]:
                container["pmids"].append(record["PUBMED_ID"])
                container["experiments"].append(record["EXPERIMENTAL_SYSTEM"])
                container["experimentTypes"].append(record["EXPERIMENTAL_SYSTEM_TYPE"])
                if "Low Throughput" in record["THROUGHPUT"]: container["lowThroughput"] = True
                if "High Throughput" in record["THROUGHPUT"]: container["highThroughput"] = True
        return container
    def getExperimentalSystems(self):

        registredSystem  = []
        for key in self.data:
            if self.data[key]["EXPERIMENTAL_SYSTEM"] not in registredSystem:
                registredSystem.append(self.data[key]["EXPERIMENTAL_SYSTEM"])

        return registredSystem

class BIOGRIDMAPPER():
    def __init__(self):
        self.uniprotToBiogrid = {}
        self.biogridToUniprot = {}
    def __call__ (self, uniprotId=None, biogridId=None):
        if biogridId and uniprotId:
            self.add(str(biogridId), str(uniprotId))
            return
        if (biogridId):
            value = self.toUniprot(str(biogridId))
            return value
        if (uniprotId):
            value = self.toBiogrid(str(uniprotId))
            return value

    def read(self, fileName):
        s=''
        with open(fileName, 'r') as f:
            s = f.read()
        self.load(s)

    def load(self, stream):
        for line in stream.split("\n"):
            if line.startswith("#") or not line:continue
            array = line.split()
            self(uniprotId=array[0], biogridId=array[1])

    def toUniprot (self, id):
        if id in self.biogridToUniprot:
            return self.biogridToUniprot[id]
        #print "Biogrid mapper failed to convert \"" + str(id) + "\" to a uniprot ID"
        return None

    def toBiogrid (self, id):
        if id in self.uniprotToBiogrid:
            return self.uniprotToBiogrid[id]
        #print "Biogrid mapper failed to convert \"" + str(id) + "\" to a biogrid ID"
        return None

    def add(self, biogridId=None, uniprotId=None):
        self.biogridToUniprot[str(biogridId)] = str(uniprotId)
        self.uniprotToBiogrid[str(uniprotId)] = str(biogridId)
