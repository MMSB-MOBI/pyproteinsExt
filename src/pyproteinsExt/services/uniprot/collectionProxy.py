from .redisClient import listUniprotKey, getUniProtEntry,\
    removeEntries, setDatabaseParameters, storeEntry
import progressbar
#from progressbar import ProgressBar, Bar, Counter, Timer, ETA, Percentage, RotatingMarker
def bootstrap(**kwargs):
    print("Bootstraping uniprot Entry redis Collection")
    setDatabaseParameters(**kwargs)

def cleanup(): 
    print("Cleaning up")
   
    cnt = 0
    _ = []
    with progressbar.ProgressBar(max_value=progressbar.UnknownLength\
        , redirect_stdout=True) as bar:
        for key in listUniprotKey():
            cnt += 1
            _.append(key)
            if cnt % 50 == 0:
                remove(_)
                _ = [] 
                bar.update(cnt)
        if _:
            remove(_)
            bar.update(cnt)
#print(f"{cnt} uniprot entries deleteted")

def convert(entryIterator):
    cnt=0
    print(">>", entryIterator, "<<")
    with progressbar.ProgressBar(max_value=len(entryIterator)\
        , redirect_stdout=True) as bar:
        for e in entryIterator:
        #print(f"converting {e.id}")
            try :
                storeEntry(e)
            except KeyError as err :
                print(f"No need to add {e.id} already in store")
            cnt+=1
            bar.update(cnt)

# Proof of concept for deletion
# Improvment stage 1 : send slice to removeEtnries. current delete decorator is one by one
# Imorvement stage 2 : use pipeline or scriptiong in redisCLient
def remove(uniprotIDs=None):
    #print(f"Removing {uniprotIDs}")
    
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

def __iter__():
    for _id in listUniprotKey(): 
        yield get(_id)

def __len__():
    cnt = 0
    for _id in listUniprotKey(): 
        cnt += 1
    return cnt

