import pyproteinsExt.hmmrContainerFactory as hmmr 
import pyproteinsExt.tmhmmContainerFactory as tmhmm
import pyproteinsExt.fastaContainerFactory as fasta
import pyproteinsExt.proteinContainer

def check_if_same_proteins(dic_container):
    list_proteins=[]
    if dic_container['hmmr']:
        hmmr_proteins=set([p for p in dic_container['hmmr'].pIndex])
        list_proteins.append(hmmr_proteins)
    if dic_container['tmhmm']:    
        tmhmm_proteins=set([e.prot for e in dic_container['tmhmm']])
        list_proteins.append(tmhmm_proteins)
    if dic_container['fasta']:
        fasta_proteins=set([e.prot for e in dic_container['fasta']])
        list_proteins.append(fasta_proteins)  
    if len(list_proteins)<=1: 
        return True
    for i in range(len(list_proteins)):
        for j in range(i+1,len(list_proteins)): 
            prot1=list_proteins[i]
            prot2=list_proteins[j]
            if len(prot1)!=len(prot2):
                return False
            if prot1.difference(prot2):
                return False
    return True            


def parse(hmmrOut,tmhmmOut,fastaOut):
    container=TopologyContainer()
    hmmrContainer=hmmr.parse(hmmrOut)
    tmhmmContainer=tmhmm.parse(tmhmmOut)
    fastaContainer=fasta.parse(fastaOut)    

    dic_container={'hmmr':hmmrContainer,'tmhmm':tmhmmContainer,'fasta':fastaContainer}
    if not check_if_same_proteins(dic_container):
        raise Exception("not same proteins ",hmmrOut,tmhmmOut,"Check !")
    container.addParsing(TopologyContainer(input=dic_container))   
    return container

class TopologyContainer(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None):
        super().__init__(_parseBuffer,input)  

class Topology(): 
    def __init__(self,prot,hmmr,tmhmm,fasta):
        self.prot=prot
        self.hmmr=hmmr
        self.tmhmm=tmhmm
        self.fasta=fasta

def _parseBuffer(dic_container):
    hmmrContainer=dic_container['hmmr']
    tmhmmContainer=dic_container['tmhmm']
    fastaContainer=dic_container['fasta']
    dic_obj={}
    for p in hmmrContainer.pIndex: 
        hmmr=hmmrContainer.pIndex[p]
        tmhmm=tmhmmContainer.entries[p]
        fasta=fastaContainer.entries[p]
        obj=Topology(p,hmmrContainer.pIndex[p],tmhmm,fasta)
        dic_obj[p]=obj
    return dic_obj                   

