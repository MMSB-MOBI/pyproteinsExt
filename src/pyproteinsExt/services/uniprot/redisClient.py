from pyrediscore import *
from .entryProxy import EntrySchema
import progressbar


entryProxyLoader = EntrySchema().load

# Basic storage no account for alternativ ID
# receive uniprot ID
# hash mainID
@connect
@store
def storeEntry(uniprotObj, *args, **kwargs):
    return f"uniprot:{uniprotObj.id}", uniprotObj.toJSON()

@connect
@storeMany
def storeManyEntries(entryList, *args, **kwargs):
    return { f"uniprot:{e.id}" : e.toJSON() for e in entryList}
    
@connect
@delete
def removeEntries(uniprotIDs, *args, **kwargs):
    return [f"uniprot:{_}" for _ in uniprotIDs]
#batch=50
#    for key in listUniprotKey():

@connect
@get
def getUniProtEntry(uniprotID, *args, raw=False, **kwargs):
    return (f"uniprot:{uniprotID}", entryProxyLoader)

@connect
@mget
def mgetUniProtEntry(uniprotIDs, *args, raw=False, **kwargs):
    _ = [ f"uniprot:{uniprotID}" for uniprotID in uniprotIDs ]
    return (_, entryProxyLoader)

@connect
@listKey
def listUniprotKey(*args, prefix=False, **kwargs):
    return ('uniprot:*', 'uniprot:')