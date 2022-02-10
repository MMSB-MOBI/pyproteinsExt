import pyproteins.container.customCollection
import pyproteins.container.Core
import re

''' Collection of Pfam entries '''

def strip(fileName):
    string = re.sub("\.[^\.]+$", "", fileName)
    return string

class EntrySet(pyproteins.container.customCollection.EntrySet, object):
    def __init__(self, collectionPath):
        if not collectionPath:
            raise ValueError("Please specify a path to pfam cache argument required")
        super().__init__(collectionPath=collectionPath, constructor=Entry, typeCheck=None, indexer=strip)
        
    def map(self, uniprotID=None):
        if uniprotID:
            e = self.get(uniprotID)
            return e.matches

class Entry(pyproteins.container.Core.Container):
    def __init__(self, id, baseUrl="http://pfam.xfam.org/protein/", fileName=None):
        if not id:
            raise TypeError('identifier is empty')
        #super(Entry, self).__init__(id, url=baseUrl + str(id) + '?output=xml', fileName=fileName)
        super().__init__(id, url=baseUrl + str(id) + '?output=xml', fileName=fileName)
        #pyproteins.container.Core.Container.__init__(self, id, url=baseUrl + str(id) + '?output=xml', fileName=fileName)
        self._parser()
    #@property
    #def xmlH(self):
    #    return self.getXmlHandler()


    def _parser(self):
        e = self.getXmlHandler()
        if not e:
            raise ValueError("not a valid Pfam xml Handler")
        if not e.find('description'):
            raise ValueError("not a valid Pfam ressource")

        self.description = re.sub("^[\n\s]*([\S]+.*[\S]+)[\n\s]*$",r"\1", e.find("description").text)

        self.matches = [ Match(x) for x in e.find_all('match') ]



#<match accession="PF00143" id="Interferon" type="Pfam-A">
#<location start="30" end="185" ali_start="31" ali_end="185" hmm_start="2" hmm_end="156" evalue="2.5e-62" bitscore="221.60"/>
#</match>
class Match():
    def __init__(self, xmlHandler):
        if xmlHandler:

            self.id = xmlHandler['id']
            self.accession = xmlHandler['accession']
            self.type = xmlHandler['type']
            self.locations = [ { k : v  for e in xmlHandler.find_all('location') for k, v in e.attrs.items()} ]

    def __repr__(self):
        string = self.id + " " + self.accession + " " + self.type + ":"
        for loc in self.locations:
            string += "{start:" + loc['start'] + " end:" + loc['end'] + "}"
        return string


# Following methods ensure duck-typing w/ Uniprot.Entry.Domain class (ie common interface)
    def owns(self, position):
        i = int (position)
        for loc in self.locations:
            if i >= int(loc['start']) and i <= int(loc['end']):
                return True
        return False

    @property
    def _dict(self):
        begin = [ int(l['start']) for l in self.locations ]
        end = [ int(l['end']) for l in self.locations ]
        return { 'description' : self.id + ":" + self.accession + ":" + self.type, 'begin' : begin, 'end' : end }
