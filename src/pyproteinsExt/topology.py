import pyproteinsExt.hmmrContainerFactory as hmmr 
import pyproteinsExt.proteinContainer

def parse(hmmrOut=None):
def parse(hmmrOut,tmhmmOut):
    container=Container()
    hmmrContainer=hmmr.parse(hmmrOut)
    tmhmmContainer=tmhmm.parse(tmhmmOut)    

    dic_container={'hmmr':hmmrContainer,'tmhmm':tmhmmContainer}
    container.addParsing(Container(input=dic_container))   
    return container

class Container(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None):
        super().__init__(_parseBuffer,input)

class Topology(): 
    def __init__(self,prot,hmmr,tmhmm):
        self.prot=prot
        self.hmmr=hmmr
        self.tmhmm=tmhmm

def _parseBuffer(dic_container):
    hmmrContainer=dic_container['hmmr']
    tmhmmContainer=dic_container['tmhmm']
    dic_obj={}
    for p in hmmrContainer.pIndex: 
        hmmr=hmmrContainer.pIndex[p]
        tmhmm=tmhmmContainer.entries[p]
        obj=Topology(p,hmmrContainer.pIndex[p],tmhmm)
        dic_obj[p]=obj
    return dic_obj                   

