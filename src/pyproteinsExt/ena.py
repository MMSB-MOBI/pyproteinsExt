import pyproteins.container.customCollection
import pyproteins.container.Core
from os.path import expanduser
import re
import os 
import io 
from skbio import DNA

enaEntrySet=None

def getENACollection():
    print("Get ENA Collection")
    global enaEntrySet
    enaEntrySet=EntrySet()
    return enaEntrySet    


class EntrySet(pyproteins.container.customCollection.EntrySet):
    def __init__(self, **kwargs):
        home = expanduser("~")
        cachePath = home
      #  if 'collectionPath' in kwargs:
      #      cachePath = kwargs['collectionPath']
      #  if 'pfamCollectionPath' in kwargs:
      #      PfamCache = kwargs['pfamCollectionPath']

        super().__init__(collectionPath=cachePath, constructor=Entry, indexer=strip)

    def serialize(self, **kwargs):
        print ("serializing uniprot collection")
        super().serialize(kwargs)

class Entry(pyproteins.container.Core.Container):
    def __init__(self, id, baseUrl="https://www.ebi.ac.uk/ena/data/view/", fileName=None):
        if not id:
            raise TypeError('identifier is empty')
        super().__init__(id, url=baseUrl + str(id) + '&display=xml', fileName=fileName)

        self.xmlHandler=self.getXmlHandler()
        self.CDS=[e for e in self.xmlHandler.find_all("feature") if e["name"]=="CDS"]

        
def strip(string):
    subString = re.search(".xml$", string)
    if subString:
        return string.split(".")[0]
    return None
