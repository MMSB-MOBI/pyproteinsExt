import pyproteinsExt.hmmrContainerFactory as hmmr 
import pyproteinsExt.proteinContainer

def parse(hmmrOut=None):
    container=Container()
    hmmrContainer=None
    if hmmrOut: 
        hmmrContainer=hmmr.parse(hmmrOut)

    dic_container={'hmmr':hmmrContainer}
    container.addParsing(Container(input=dic_container))   
    return container

class Container(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None):
        super().__init__(_parseBuffer,input)

class Topology(): 
    def __init__(self,prot,hmmr):
        self.prot=prot
        self.hmmr=hmmr

def _parseBuffer(dic_container):
    hmmrContainer=dic_container['hmmr']
    dic_obj={}
    for p in hmmrContainer.pIndex: 
        obj=Topology(p,hmmrContainer.pIndex[p])
        dic_obj[p]=obj
    return dic_obj                   

