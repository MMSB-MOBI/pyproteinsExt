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
@delete
def removeEntries(uniprotIDs, *args, **kwargs):
    return [f"uniprot:{_}" for _ in uniprotIDs]
#batch=50
#    for key in listUniprotKey():
    pass

#@connect
#@delete
#def uniprotClean(*args, **kwargs):
#batch=50
#    for key in listUniprotKey():
#    pass

@connect
@get
def getUniProtEntry(uniprotID, *args, raw=False, **kwargs):
    return (f"uniprot:{uniprotID}", entryProxyLoader)

@connect
@listKey
def listUniprotKey(*args, prefix=False, **kwargs):
    return ('uniprot:*', 'uniprot:')