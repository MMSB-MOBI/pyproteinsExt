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
from colour import Color
import pyproteinsExt.refseq as refseq
from os.path import expanduser
import pyproteinsExt.uniprot as uniprot
import copy 

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

    def __getitem__(self, index):
        return list(self.entries.values())[index]

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
                new_e=Topology(e.prot,hit_to_add,e.tmhmm,e.fasta,e.taxo,e.uniprot_entry,e.annotated_domains_fragments,e.Nter_UR_fragment,e.Cter_UR_fragment,e.helix_fragments,e.loop_fragments,e.unknown1_fragment,e.unknown2_fragment)    
                new_container.addEntry(new_e)
        return new_container 

    def filter_last_helix(self,distance=90):
        new_container=TopologyContainer()
        c=0
        for e in self:
            new_helix_fragments=copy.deepcopy(e.helix_fragments)
            new_loop_fragments=copy.deepcopy(e.loop_fragments)
            if len(e.helix_fragments)==7:
                c+=1
                new_helix_fragments=new_helix_fragments[:-1]
                new_loop_fragments=new_loop_fragments[:-1]
            elif e.helix_fragments[-1]["start"]-e.helix_fragments[-2]["end"]>distance:
                c+=1
                new_helix_fragments=new_helix_fragments[:-1]
                new_loop_fragments=new_loop_fragments[:-1]
            new_e=Topology(e.prot,e.hmmr,e.tmhmm,e.fasta,e.taxo,e.uniprot_entry,e.annotated_domains_fragments,e.Nter_UR_fragment,e.Cter_UR_fragment,new_helix_fragments,new_loop_fragments,e.unknown1_fragment,e.unknown2_fragment)        
            new_container.addEntry(new_e)
        print(c,"proteins have been filtered")    
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
    def create_ete3_tree(self):
        ncbi=NCBITaxa()
        taxids=set([e.taxo.taxid for e in self])
        if None in taxids: 
            raise Exception("Entries doesn't have taxids")
        tree=ncbi.get_topology(list(taxids))
        
        #Complete Tree object with list of domains and proteins for each node 
        node_list=[]
        for n in tree.traverse('postorder'): #Browse tree in postorder, starts from leaf and ascend to root 
            n.sameDomainNode=set()
            node_list.append(n)
            n.domains=set([h.domain for e in self for h in e.hmmr if e.taxo.taxid==n.name])
            n.proteins=set([e.prot  for e in self if e.taxo.taxid==n.name])
            if n.get_descendants():
                for child in n.children: 
                    n.domains.update(child.domains)
                    n.proteins.update(child.proteins)

        #Complete Tree object with list of nodes with same domains for each node            
        c=0
        for i in range(len(node_list)):
            c+=1
            for j in range(i+1,len(node_list)):
                n1=node_list[i]
                n2=node_list[j]
                if len(n1.domains)==len(n2.domains):
                    if not n1.domains.difference(n2.domains):
                        n1.sameDomainNode.add(n2)
                        n2.sameDomainNode.add(n1)
        self.ete3_tree=tree                   
            
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
    

    def compute_upper_node_and_distance(self,core_domains=[]):
        ncbi=NCBITaxa()
        if not self.domain_entries:
            raise Exception("Compute domain_entries first.")
        if not self.ete3_tree:
            raise Exception("Compute ete3_tree first.")

        for d in self.domain_entries.values(): 
            if d.name not in core_domains:  
                distances=[]
                if len(d.taxo)==1:
                    taxo=list(d.taxo)[0]
                    d.upper_node=self.ete3_tree.search_nodes(name=taxo.taxid)[0]
                    d.mean_distance=0
                else:
                    list_taxids=list(set([t.taxid for t in d.taxo]))
                    domain_tree=ncbi.get_topology(list_taxids)
                    traverse_generator=domain_tree.traverse()
                    d.upper_node=next(traverse_generator)   
                    for i in range(len(list_taxids)):
                        for j in range(i+1,len(list_taxids)):
                            dist=self.ete3_tree.get_distance(list_taxids[i],list_taxids[j])
                            distances.append(dist)
                    d.mean_distance=mean(distances)    

    def create_domain_graph(self,core_domains):

        def get_vertex_size(nb_domain,interval,min_size):
            if nb_domain==1: 
                return min_size 
            else: 
                size=min_size+interval*(nb_domain-1)
                return size

        def get_edge_size(nb_occ,interval): 
            return nb_occ*interval     

        g=Graph()
        dic_edges={}
        dic_nb_domain={}
        all_domains=set()
        for p in self: 
            domains=set([h.domain for h in p.hmmr])
            for cd in core_domains: 
                domains.discard(cd)
            related_domains=domains.copy()
            for d in domains:
                if d not in core_domains: 
                    if d not in dic_nb_domain :
                        dic_nb_domain[d]=0 
                    dic_nb_domain[d]+=1    
                related_domains.remove(d)
                all_domains.add(d)
                for d2 in related_domains : 
                    edge=tuple(sorted((d,d2)))
                    if edge not in dic_edges: 
                        dic_edges[edge]=0
                    dic_edges[edge]+=1  
                    
        list_edges=[]
        list_weight=[]
        for e in dic_edges: 
            list_edges.append(e)
            list_weight.append(dic_edges[e])  
        
        g=Graph()
        g.add_vertices(len(all_domains))
        g.vs["name"]=list(all_domains)
        g.vs["label"]=g.vs["name"]
        g.add_edges(list_edges)
        for vertex in g.vs : 
            vertex["weight"]=dic_nb_domain[vertex["name"]]
            vertex["size"]=get_vertex_size(vertex["weight"],2,5)
        for e in g.es : 
            source=[v for v in g.vs if v.index==e.source][0] 
            target=[v for v in g.vs if v.index==e.target][0]
            edge_tuple=tuple(sorted((source["name"],target["name"])))
            nb_occ=dic_edges[edge_tuple]
            min_domains=min(source["weight"],target["weight"])
            e["weight"]=nb_occ/min_domains
            e["label"]=round(nb_occ/min_domains,2)
            e['width']=e["weight"]*5

        return g
                        
    def separate_seq_into_fragments(self):        
        for e in self: 
            e.get_annotated_domains_fragments()
            e.get_Nter_UR_fragment()
            e.get_Cter_UR_fragment()
            e.get_helix_fragments()
            e.get_loop_fragments()
    def get_unknown_fragments(self):
        for e in self : 
            e.get_unknown1_fragment()
            e.get_unknown2_fragment()        
                        
class Topology(): 
    def __init__(self,prot,hmmr,tmhmm,fasta,taxo=None,uniprot_entry=None,annotated_domains_fragments=None,Nter_UR_fragment=None,Cter_UR_fragment=None,helix_fragments=None,loop_fragments=None,unknown1_fragment=None,unknown2_fragment=None):
        self.prot=prot
        self.hmmr=hmmr
        self.tmhmm=tmhmm
        self.fasta=fasta
        self.taxo=taxo
        self.uniprot_entry=uniprot_entry
        self.annotated_domains_fragments=annotated_domains_fragments
        self.Nter_UR_fragment=Nter_UR_fragment
        self.Cter_UR_fragment=Cter_UR_fragment
        self.helix_fragments=helix_fragments
        self.loop_fragments=loop_fragments
        self.unknown1_fragment=unknown1_fragment
        self.unknown2_fragment=unknown2_fragment

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

    def get_annotated_domains_fragments(self):
        annotated_domains_fragments=[]
        for hit in self.hmmr: 
            dic_core_domains={}
            if hit.domain in dic_core_domains: 
                raise Exception("Several hits for domain. Handle this part")
            dic_core_domains={'name':hit.domain,'seq':hit.get_sequence(),'start':hit.get_start(),'end':hit.get_end()}
            annotated_domains_fragments.append(dic_core_domains)       
        annotated_domains_fragments.sort(key=lambda r:r["start"])
        self.annotated_domains_fragments=annotated_domains_fragments
        return annotated_domains_fragments

    def get_Nter_UR_fragment(self):
        dic={'name':'N-ter_UR','start':1}
        Nter=self.tmhmm.fragments[0]
        #print(self.annotated_domains_fragments[0]["start"])
        if Nter.end > self.annotated_domains_fragments[0]["start"]: 
            dic["end"]=self.annotated_domains_fragments[0]["start"]
        else: 
            dic["end"]=Nter.end
        dic["seq"]=self.fasta.get_subsequence(dic["start"],dic["end"])  
        self.Nter_UR_fragment=dic
        #return dic       

    def get_Cter_UR_fragment(self):
        dic={'name':"C-ter_UR"} 
        Cter=self.tmhmm.fragments[-1]
        dic["end"]=Cter.end 
        if Cter.start < self.annotated_domains_fragments[-1]["end"] :
            dic["start"]=self.annotated_domains_fragments[-1]["end"]
        else: 
            dic["start"]=Cter.start
        dic["seq"]=self.fasta.get_subsequence(dic["start"],dic["end"])  
        self.Cter_UR_fragment=dic
        #return dic

    def get_unknown1_fragment(self):
        dic={'name':"Unknown1"}
        start=self.helix_fragments[-1]["end"]+1
        end=[d for d in self.annotated_domains_fragments if d["name"]=="fad_binding_prokaryotes"][0]["start"]-1
        if end < start:
            self.unknown1_fragment=[]
        else:     
            seq=self.fasta.get_subsequence(start,end)
            dic["start"]=start
            dic["end"]=end
            dic["seq"]=seq
            self.unknown1_fragment=[dic]

    def get_unknown2_fragment(self):
        dic={"name":'Unknown2'}
        start=[d for d in self.annotated_domains_fragments if d["name"]=="fad_binding_prokaryotes"][0]["end"]+1
        end=[d for d in self.annotated_domains_fragments if d["name"]=="nad_binding_prokaryotes"][0]["start"]-1
        seq=self.fasta.get_subsequence(start,end)
        dic["start"]=start
        dic["end"]=end
        dic["seq"]=seq
        self.unknown2_fragment=[dic]

    def get_helix_fragments(self):
        list_fragments=[]
        helixes=[f for f in self.tmhmm.fragments if f.cellular_location=="TMhelix"]
        helix_number=1
        for h in helixes: 
            dic={'name':"TMhelix_"+str(helix_number),'start':h.start,'end':h.end}
            helix_number+=1
            dic["seq"]=self.fasta.get_subsequence(dic["start"],dic["end"])
            list_fragments.append(dic)
        self.helix_fragments=list_fragments
        #return list_fragments

    def get_loop_fragments(self):
        list_fragments=[]
        loops=[f for f in self.tmhmm.fragments[1:-1] if f.cellular_location!="TMhelix"]
        count_inside=0
        count_outside=0
        for l in loops: 
            dic={}
            if l.cellular_location=="inside":
                count_inside+=1
                dic["name"]="inside_loop_"+str(count_inside)
                dic["start"]=l.start
                dic["end"]=l.end
            elif l.cellular_location=="outside":
                count_outside+=1
                dic["name"]="outside_loop_"+str(count_outside)
                dic["start"]=l.start
                dic["end"]=l.end    
            dic["seq"]=self.fasta.get_subsequence(dic["start"],dic["end"])    
            list_fragments.append(dic)
        self.loop_fragments=list_fragments    
        #return list_fragments            

    def set_uniprot_xref(self, uColl):
        p_id = self.prot.split("|")[1]
        try:
            uniprot_entry = uColl.get(p_id)
            self.uniprot_xref = uniprot_entry.xref
        except ValueError as ve:
            ve = str(ve)
            if ve != "Error, empty xmlHandler":
                raise Exception()
            self.uniprot_xref = None


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
