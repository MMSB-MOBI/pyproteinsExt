from pyrediscore import *
from ..models.entryProxy import EntrySchema

entryProxyLoader = EntrySchema().load

# Basic storage no account for alternativ ID
# receive uniprot ID
# hash mainID
@connect
@store
def storeEntry(e, as_dict=False, *args, **kwargs):
    return f"uniprot_id:{e.id}", e.toJSON() if not as_dict else e

@connect
@storeMany
def storeManyEntries(entryList, as_dict=False,  *args, **kwargs):
    return { f"uniprot_id:{e.id}" : e.toJSON() if not as_dict else e for e in entryList}
    
@connect
@delete
def removeEntries(uniprotIDs, *args, **kwargs):
    return [f"uniprot_id:{_}" for _ in uniprotIDs]
#batch=50
#    for key in listUniprotKey():

@connect
@get
def getUniProtEntry(uniprotID, *args, raw=False, **kwargs):
    return (f"uniprot_id:{uniprotID}", entryProxyLoader)

@connect
@mget
def mgetUniProtEntry(uniprotIDs, *args, raw=False, **kwargs):
    _ = [ f"uniprot_id:{uniprotID}" for uniprotID in uniprotIDs ]
    return (_, entryProxyLoader)

@connect
@listKey
def listUniprotKey(*args, prefix=False, **kwargs):
    return ('uniprot_id:*', 'uniprot_id:')