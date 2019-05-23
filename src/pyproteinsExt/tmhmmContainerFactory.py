import pyproteinsExt.proteinContainer
import re

def parse(inputFile=None):
    bigBuffer = ''
    instances_list=[]
    mainContainer = Container()
    fnOpen = open
    t = 'r'
    if inputFile.endswith('gz'):  
        fnOpen = gzip.open
        t = 'rt'
    try:
        f = fnOpen(inputFile, t)
    except IOError:
        print ("Could not read file:", inputFile)
        return Container()

    instance=[]
    with f:
        for l in f:
            if l.startswith("#") and "Length" in l : 
                if instance : 
                    mainContainer=mainContainer.addParsing(Container(input=instance))        
                instance=[]
            instance.append(l.rstrip())    
        mainContainer=mainContainer.addParsing(Container(input=instance))  
    return mainContainer
    #print(instances_list[-1])
   

def filter_nb_helix(e,**kwargs):
    nb_helix_min=kwargs.get("nb_helix_min",None)
    nb_helix_max=kwargs.get("nb_helix_max",None)
    if not nb_helix_min:
        raise TypeError("filter_nb_helix needs nb_helix_min argument")         
    if not nb_helix_max:
        raise TypeError("filter_nb_helix needs nb_helix_max argument")   

    if e.nb_helix >= nb_helix_min and e.nb_helix <= nb_helix_max:
        return True 

    return False           

class Container(pyproteinsExt.proteinContainer.Container):
    def __init__(self, input=None):
        super().__init__(_parseBuffer,input)
       

def _parseBuffer(input):
    def parseFragments(list_fragments):
        fragmentsObj=[]
        reLocation=re.compile("[\s]+([\d]+)[\s]+([\d]+)")
        for f in list_fragments:
            f_split=f.split("\t")
            cellular_location=f_split[2]
            start=f_split[3].split(" ")
            location=reLocation.findall(f_split[3])[0]
            start=int(location[0])
            end=int(location[1])
            fragment=Fragment(cellular_location,start,end)
            fragmentsObj.append(fragment)
        return fragmentsObj    

    prot=input[0].split(" ")[1]
    nb_helix=[l for l in input if "Number of predicted TMHs" in l][0].split(":")[1].strip()
    prot_length=[l for l in input if "Length" in l][0].split(":")[1].strip()
    fragments=parseFragments([l for l in input if not l.startswith("#")])
    obj=TMHMM_Obj(prot,int(prot_length),int(nb_helix),fragments)
    return {prot:obj}


class TMHMM_Obj(): 
    def __init__(self,prot,prot_length,nb_helix,fragments):
        self.prot=prot
        self.prot_length=prot_length
        self.nb_helix=nb_helix
        self.fragments=fragments
        self.topology_seq=self.get_topology_seq()

    def get_topology_seq(self):
        topology_seq=["*"]*self.prot_length
        helix_number=1
        for f in self.fragments: 
            if f.cellular_location=='inside':
                letter="o"
            elif f.cellular_location=="outside":
                letter="i"
            elif f.cellular_location=="TMhelix":
                letter=str(helix_number) 
                helix_number+=1
            for i in range(f.start-1,f.end): 
                topology_seq[i]=letter
        if "*" in topology_seq : 
            print("WARNING : * in topology seq")       
        return "".join(topology_seq)
                             


class Fragment():
    def __init__(self,cellular_location,start,end):
        self.cellular_location=cellular_location
        self.start=start
        self.end=end  


