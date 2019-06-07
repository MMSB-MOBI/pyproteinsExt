import pyproteinsExt.hmmrContainerFactory as hmmr 
import pyproteinsExt.tmhmmContainerFactory as tmhmm
import pyproteinsExt.fastaContainerFactory as fasta
import pyproteinsExt.proteinContainer
from collections import OrderedDict
import re
from ete3 import NCBITaxa
from ete3 import Tree
from statistics import mean
from igraph import Graph

def check_if_same_proteins(dic_container):
    hmmr_proteins=set()
    tmhmm_proteins=set()
    fasta_proteins=set()
    if dic_container['hmmr']:
        hmmr_proteins=set([h.prot for h in dic_container['hmmr'].hmmrEntries])
    if dic_container['tmhmm']:    
        tmhmm_proteins=set([e.prot for e in dic_container['tmhmm']])
    if dic_container['fasta']:
        fasta_proteins=set([e.prot for e in dic_container['fasta']])

    list_proteins=[hmmr_proteins,tmhmm_proteins,fasta_proteins]
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

    #if len(fastaContainer)==0: 
    #    return container 

    dic_container={'hmmr':hmmrContainer,'tmhmm':tmhmmContainer,'fasta':fastaContainer}
    if not check_if_same_proteins(dic_container):
        raise Exception("not same proteins ",hmmrOut,tmhmmOut,"Check !")
    container.addParsing(TopologyContainer(input=dic_container))   
    return container

class TopologyContainer(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None,ete3_tree=None,domain_entries=None):
        super().__init__(_parseBuffer,input)
        self.ete3_tree=ete3_tree
        self.domain_entries=domain_entries
          
    def filter(self,fPredicat,**kwargs): 
        new_container=TopologyContainer()
        for e in self : 
            if fPredicat(e,**kwargs):
                new_container.addEntry(e)
        return new_container        

    def get_domain_mfasta(self,domain):
        mfasta=''
        for e in self: 
            for hit in e.hmmr: 
                if hit.domain == domain:
                    if float(hit.hit.iEvalue)>1e-3:
                        print("OOOO")
                    header=">"+hit.prot+" "+hit.domain
                    #print(header)
                    seq=self.get_seq(hit.hit)
                    mfasta+=header+"\n"+seq+"\n"
        return mfasta                

    def get_seq(self,hit):
        seq=hit.aliStringLetters   
        seq=seq.replace("-","")
        seq=seq.upper()
        return seq    

    def proteins_mfasta(self):
        mfasta=''
        for e in self: 
            mfasta+=">"+e.fasta.header+"\n"+e.fasta.seq+"\n"
        return mfasta 

    '''Redefine addParsing parent method to include domain_entries'''

    def complete_hmmr(self,hmmscan_out):  
        container=hmmr.parse(hmmscan_out)
        new_proteins=set([h.prot for h in container.hmmrEntries])
        self_proteins=set(self.entries.keys())
        if len(new_proteins)!=len(self_proteins):
            raise Exception("full Pfam proteins and original proteins are not the same")
        else: 
            if new_proteins.difference(self_proteins):
                raise Exception("full Pfam proteins and original proteins are not the same")        
        for e in self: 
            hits=[h for h in container.hmmrEntries if h.prot==e.prot]
            e.hmmr+=hits

    def filter_hit(self,fPredicat,**kwargs):
        new_container=TopologyContainer()
        for e in self: 
            add=False
            hit_to_add=[]
            for hit in e.hmmr: 
                if fPredicat(hit,**kwargs):
                    add=True
                    hit_to_add.append(hit)
            if add :     
                new_e=Topology(e.prot,hit_to_add,e.tmhmm,e.fasta,e.taxo)    
                new_container.addEntry(new_e)
        return new_container 
    def compute_overlapped_domains(self,overlap_accept_size):
        self.reinitialize_overlapped_domains()
        for e in self: 
            for i in range(len(e.hmmr)):
                for j in range(i+1,len(e.hmmr)):
                    hit1=e.hmmr[i]
                    hit2=e.hmmr[j]
                    if hit1.is_overlapping(hit2,overlap_accept_size):
                        hit1.overlapped_hits.append(hit2)
                        hit2.overlapped_hits.append(hit1)

    def reinitialize_overlapped_domains(self):
        for h in [h for e in self for h in e.hmmr ]: 
            h.reinitialize_overlapped_hits()      
    def create_domain_entries(self):
        def initialize_domain_entries():
            domain_entries={}
            domains=set([h.domain for e in self for h in e.hmmr])
            for d in domains: 
                domainObj=Domain(d,set(),set(),set())
                domain_entries[d]=domainObj 
            return domain_entries    
        domain_entries=initialize_domain_entries()
        for e in self:
            for h in e.hmmr: 
                domain_entries[h.domain].hits.add(h)
                domain_entries[h.domain].proteins.add(e.prot)
                domain_entries[h.domain].taxo.add(e.taxo)
        domain_entries=OrderedDict(sorted(domain_entries.items(),key=lambda kv: len(kv[1].proteins),reverse=True))        
        self.domain_entries=domain_entries    
class Topology(): 
    def __init__(self,prot,hmmr,tmhmm,fasta,taxo=None):
        self.prot=prot
        self.hmmr=hmmr
        self.tmhmm=tmhmm
        self.fasta=fasta
        self.taxo=taxo

    def get_taxo(self,function_get_taxid):
        ncbi=NCBITaxa()
        taxname=None
        taxrank=None
        taxid=function_get_taxid(self)
        taxname_dic=ncbi.get_taxid_translator([taxid])
        if taxname_dic:
            taxname=taxname_dic[int(taxid)]
            taxrank_dic=ncbi.get_rank([taxid])
            if taxrank_dic : 
                taxrank=taxrank_dic[int(taxid)]
        self.taxo=Taxo(taxid,taxname,taxrank)  

class Domain(): 
    def __init__(self,name,hits,proteins,taxo):
        self.name=name
        self.hits=hits
        self.proteins=proteins
        self.taxo=taxo
        self.upper_node=None
        self.mean_distance=None

def _parseBuffer(dic_container):
    hmmrContainer=dic_container['hmmr']
    tmhmmContainer=dic_container['tmhmm']
    fastaContainer=dic_container['fasta']
    dic_obj={}
    for f in fastaContainer: 
        p=f.prot
        hmmr=[h for h in hmmrContainer.hmmrEntries if h.prot==p]
        tmhmm=tmhmmContainer.entries[p]
        fasta=fastaContainer.entries[p]
        obj=Topology(p,hmmr,tmhmm,fasta)
        dic_obj[p]=obj
    return dic_obj   

class Taxo():
    def __init__(self,taxid,taxname,taxrank):
        self.taxid=taxid
        self.taxname=taxname
        self.taxrank=taxrank
