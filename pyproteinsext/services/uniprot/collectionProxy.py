from .redisClient import listUniprotKey, getUniProtEntry,\
    removeEntries, setDatabaseParameters, storeEntry,\
    storeManyEntries, mgetUniProtEntry
from progressbar import Percentage, Bar, ETA, AdaptiveETA, ProgressBar, UnknownLength, Counter, Timer
#from progressbar import ProgressBar, Bar, Counter, Timer, ETA, Percentage, RotatingMarker
def bootstrap(**kwargs):
    print("Bootstraping uniprot Entry redis Collection")
    setDatabaseParameters(**kwargs)
    print(f"{__len__()} uniprot entries were found in store")

def cleanup(displayCnt=True): 
    cnt = 0
    _ = []

    widgets = ['Deleting ', Counter(), ' ', Percentage(),\
               ' ', Bar(),\
               ' ', Timer(),\
               ' ', AdaptiveETA()]
    with ProgressBar(widgets=widgets, maxval= __len__()\
        if displayCnt else UnknownLength\
        , redirect_stdout=True) as bar:
        for key in listUniprotKey():
            cnt += 1
            _.append(key)
            if cnt % 50 == 0:
                try:
                    remove(_)
                except KeyError:
                    print("Missing keys to delete")
                _ = [] 
                bar.update(cnt)
        if _:
            try:
                remove(_)
            except KeyError:
                print("Missing keys to delete")                
            bar.update(cnt)

def convert(entryIterator, bulkSize=250):
    _ = []
    cnt=0
    #print(">>", entryIterator, "<<")
    widgets = ['Loading ',  Counter(), ' ', Percentage(),\
               ' ', Bar(),\
               ' ', Timer(),\
               ' ', AdaptiveETA()]
    with ProgressBar(widgets=widgets,maxval=len(entryIterator)\
        , redirect_stdout=True) as bar:
        for e in entryIterator:
            _.append(e)
            cnt+=1
            if cnt % bulkSize == 0:
                storeManyEntries(_)                
                bar.update(cnt)
                _ = []
        if _:
            storeManyEntries(_)                
            bar.update(cnt)

def getSliceIDs(cstart=0, cstop=None):
    cnt = 0
    results = [ ]
    for e in __iter__():
        if not cstop is None:
            if cnt == cstop:
                return results
        
        if cnt >= cstart:           
            results.append(e.id)
        cnt += 1
    
    return results

def add(e):
    """Add a single uniprot entry to the store"""
    try :
        storeEntry(e)
    except KeyError as err :
        print(f"No need to add {e.id} already in store")

# Proof of concept for deletion
# Improvment stage 1 : send slice to removeEtnries. current delete decorator is one by one
# Imorvement stage 2 : use pipeline or scriptiong in redisCLient
def remove(uniprotIDs=None):
    print(f"Removing {uniprotIDs}")
    
    if uniprotIDs: # slow delete all 
        if type(uniprotIDs) == list:        
            removeEntries(uniprotIDs)
        elif type(uniprotIDs) == str:
            removeEntries([uniprotIDs])
    
# Get Object
# deserialize
def get(uniprotID, raw=False):
    _ =  getUniProtEntry(uniprotID)
    #print(f"PWEEPWEE {_}")
    return _

def mget(uniprotIDs, raw=False):
    _ =  mgetUniProtEntry(uniprotIDs, raw=raw)
    return _

# Not sure its possible at module level
def __iter__():
    chunckSize = 250
    _ = []
    cnt = 0
    for _id in listUniprotKey():
        cnt += 1 
        _.append(_id)
        if cnt % chunckSize == 0:            
            for e in mget(_):
                yield e
            _ = []
    if _:
        for e in mget(_):
            yield e 
    print(f"Completed total {cnt} iterations")

# Not possible at module level
# TypeError: object of type 'module' has no len()
def __len__():
    cnt = 0
    for _id in listUniprotKey(): 
        cnt += 1
    return cnt

def length():
    cnt = 0
    for _id in listUniprotKey(): 
        cnt += 1
    return cnt


