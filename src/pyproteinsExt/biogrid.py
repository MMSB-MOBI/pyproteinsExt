import urllib2
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
class BIOGRID():
    def __init__(self, mapperUrl=BIOGRID_UNIPROT_MAPPER_URL):
        self.biogridMapper = BIOGRIDMAPPER()
        self.loadBiogridMapper()
        self.data = {}

    def __iter__(self):
        for k in self.data:
            yield self.data[k]

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        string = ""
        for key in self.data:
            data = self.data[key]
            for field in BIOGRID_ORDERED_JSON_KEYS:
                v = data[field]
                if field.startswith("BIOGRID_ID_"):
                    v = self.biogridMapper(biogridId=v)
                    if not v: v = "-"
                string = string + str(v) + "\t"
            string = string.rstrip() + "\n"
        return string

    def clear(self):
        self.data = {}

    def dump(self, file=None):
        if file:
            with open(file, 'w') as f:
                f.write("#" + "\t".join(BIOGRID_ORDERED_JSON_KEYS) + "\n")
                f.write(self.__repr__())
        else:
            print self.__repr__()
    def query(self, uniprotId=None, geneIdA=None, geneIdB=None):
        if uniprotId:
            self.uniprotQuery(uniprotId)
        elif geneIdA and geneIdB:
            self.genePairQuery(geneIdA, geneIdB)

    def genePairQuery(self, geneIdA, geneIdB):
        url = "http://webservice.thebiogrid.org/interactions/?geneList="
        url = url + geneIdA + "|" + geneIdB + "&format=json&accessKey=" + BIOGRID_KEY
        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError as error:
            print url + "\nHTTP ERROR " + str(error.code)
            return None
        except urllib2.URLError as error:
            print url + "\n" + str(error.reason)
            return None
        raw = response.read()
        response.close()
        self._parse(raw)
        self._filter(geneIdA, geneIdB)
        #http://webservice.thebiogrid.org/interactions/?searchNames=true&geneList=cdc27|apc1|apc2&evidenceList=Affinity Capture-MS|Two-hybrid&accesskey=[ACCESSKEY]

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

    def uniprotQuery(self, uniprotId):
        url = "http://webservice.thebiogrid.org/interactions/?geneList="
        biogridId = self.biogridMapper(uniprotId=uniprotId)
        if not biogridId: return
        url = url + biogridId + "&searchbiogridids=true&includeInteractors=true&format=json&accessKey=" + BIOGRID_KEY
        try:
            response = urllib2.urlopen(url)
        except urllib2.HTTPError as error:
            print url + "\nHTTP ERROR " + str(error.code)
            return None
        except urllib2.URLError as error:
            print url + "\n" + str(error.reason)
            return None
        raw = response.read()
        response.close()
        self._parse(raw)

    def _parse(self, raw):
        self.data = json.loads(raw)

    def loadBiogridMapper(self):
        try:
            response = urllib2.urlopen(BIOGRID_UNIPROT_MAPPER_URL)

        except urllib2.HTTPError as error:
            print "loadBioGridMapper url resolution failed"
            print "HTTP ERROR " + str(error.code)
            return None
        except urllib2.URLError as error:
            print "loadBioGridMapper url resolution failed"
            print error.reason
            return None
        raw = response.read()
        response.close()
        for line in raw.split("\n"):
            if line.startswith("#") or not line:continue
            array = line.split()
            self.biogridMapper(uniprotId=array[0], biogridId=array[1])

    def analyse(self):
        if len(self) == 0: return None
        container = { "pmids" : [], "experiments" : [], "experimentTypes" : [],
                     "highThroughput" : False, "lowThroughput" : False }
        for record in self:
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

    def toUniprot (self, id):
        if id in self.biogridToUniprot:
            return self.biogridToUniprot[id]
        print "Biogrid mapper failed to convert \"" + str(id) + "\" to a uniprot ID"
        return None

    def toBiogrid (self, id):
        if id in self.uniprotToBiogrid:
            return self.uniprotToBiogrid[id]
        print "Biogrid mapper failed to convert \"" + str(id) + "\" to a biogrid ID"
        return None

    def add(self, biogridId=None, uniprotId=None):
        self.biogridToUniprot[str(biogridId)] = str(uniprotId)
        self.uniprotToBiogrid[str(uniprotId)] = str(biogridId)
