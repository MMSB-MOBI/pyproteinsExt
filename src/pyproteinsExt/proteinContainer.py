class Container(object):
    def __init__(self,fParsing=None,input=None):
       # self.queryHmmFile=None             #../../../data/fad_binding.hmm
       # self.targetSequenceDatabase=None   #./Trembl_50.fasta
       # self.query=None                    #PF08022_full  [M=104]
        #self.upType = upType
        self._transpose = None
        if not input:
           self.entries={}
        else:
            self.entries=self.parsing(input,fParsing)
    
    def addParsing(self, other):
        for k in other.entries: 
            if k not in self.entries: 
                self.entries[k]=other.entries[k]
        return self

    def addEntry(self,new_entry): 
        if not new_entry.prot in self.entries: 
            self.entries[new_entry.prot]=new_entry      

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        for protein in self.entries: 
            yield self.entries[protein]  

    def filter(self,fPredicat,**kwargs): 
        new_container=Container()
        for e in self : 
            if fPredicat(e,**kwargs):
                new_container.addEntry(e)
        return new_container

    def parsing(self,input,fParsing):
        return fParsing(input)        