# pyproteinsext

[Documentation](https://cecilpert.github.io/my_pyproteinsext/)

## notebook test&examples

TO DO

## Generic container modules

TO DO

## Anotation modules

TO DO

## Specific container modules

### Multiple Sequence Alignment

#### Reading a file

```python
import pyproteins.sequence.msa as msaLib
oMsa = msaLib.Msa(fileName="/Users/guillaumelaunay/work/projects/MSA/clustalw.aln")
```

#### Accessing sequences

```python
print(oMsa[0])
```

will display

```python
{'header': 'sp|Q5SJH5|RIMM_THET8',
'sequence':'------MRLVEIGRFGAPYALKGGLRF--RGEP---VVLHLER----'}
```

##### Accessing matching sequences

Retrieve a sequence in the msa, by specifying a predicate function or a regular expression.
The predicate function will be applied to each record in turn, it will be passed a dictionary with the **index** and **record** keys.
Records matching the lookup will be returned in a list with sequence gap striped.

###### indexed-based search

```python
def f(d):
    return d["index"] < 10
recordList = oMsa.recordLookup(predicate=f)
print(len(recordList))
print(recordList[0])
```

will display

```python
10
{'header': 'sp|Q5SJH5|RIMM_THET8', 'sequence': 'MRLVEIGRFGAPYALKGGLRFRGEPVVLHLERVYVEGHGWRAIEDLYRVGEELVVHLAGVTDRTLAEALVGLRVYAEVADLPPLEEGRYYYFALIGLPVYVEGRQVGEVVDILDAGAQDVLIIRGVGERLRDRAERLVPLQAPYVRVEEGSIHVDPIPGLFD'}
```

###### header content based search

```python
def g(d):
    return re.search("THET", d["record"]['header'])
recordList = oMsa.recordLookup(predicate=g)
print(len(recordList))
print(recordList[0])
```

will display

```python
3
{'header': 'sp|Q5SJH5|RIMM_THET8', 'sequence': 'MRLVEIGRFGAPYALKGGLRFRGEPVVLHLERVYVEGHGWRAIEDLYRVGEELVVHLAGVTDRTLAEALVGLRVYAEVADLPPLEEGRYYYFALIGLPVYVEGRQVGEVVDILDAGAQDVLIIRGVGERLRDRAERLVPLQAPYVRVEEGSIHVDPIPGLFD'}
```

#### Transformations

The following methods all return a new **Alignment** object.

##### Column deletions

###### Purging gap

Specify a treshold of gap frequencies to filter out columns, default value is 0.5

```python
oMsa.gapPurge(gapRatio=0.5)
```

###### sequence based masking

Use any sequence of the alignment to delete all columns where this sequence features a gap. Default *master* sequence is number 0.

```python
oMsa.maskMaster(self, masterIndex=0)
```

###### sequence based filtering

Sequence within a MSA can be filtered according to their relationships with a master sequence. The predicate function will be applied to all sequence in turn.
Predicate will be passed 3 arguments :

* a **3-tuple statistic** wrapping (_sequence identity_, _sequence similarity_, _sequence coverage_) of the master with respect to current sequence.
* the **master sequence** object
* the **current sequence** object

Returned object is a MSA of at least one sequence (the master).

An optional named *masterIndex* can be pass to use an alternative sequence as reference.


```python
# Defining predicate, here minimal coverave of 85%
def f(stat, iSeq, jSeq):
    return stat[2] > 0.85
bMsa = oMsa.masterFilter(predicate=f)
print(bMsa.fastaDump())
```

will print,

```python
>sp|Q5SJH5|RIMM_THET8
---------MRLVEIGRFGAPYALKGGLRFRGEPVVLHLERVYVEGHGWRAIEDLYRVGE
ELVVHLAGVTDRTLAEALVGLRVYAEVADLPPLEEGRYYYFALIGLPVYVEGRQVGEVVD
ILDAGAQDVLIIRGVGERLRDRAERLVPLQAPYVRVEEGSIHVDPIPGLFD
>tr|H9ZRG5|H9ZRG5_THETH
---------MRLVEIGRFGAPYALKGGLRFRGEPVVLHLERVYVEGHGWRAIEDLYRVGE
ELVVHLAGVTDRTLAEALVGLRVYAEVADLPPLEEGRYYYFALIGLPVYVEGRQVGEVVD
ILDAGAQDVLIIRGVGERLRDRAERLVPLQAPYVRVEEGGIHVDPIPGLFD
>tr|E8PJQ1|E8PJQ1_THESS
MGLWHNGLGMRLVEIGRFGAPYALRGGLKFRGEPVVAHLERVYVEGHGWRAVEDLYQVGD
DLVVHLAGVSSRELAEPLVGLRVYAEVEELPPLEEGRYYYFALIGLPVYVGGLKMGEVVD
ILDAGAQDVLVIRGVGERLRDQTERLVPLQAPYVRVEEEGIHVEPIPGLFD
>tr|B7A7I3|B7A7I3_THEAQ
-------MAGRLVEIGRFGAPYALAGGLKFRGEPVVAHLTRIYVEGHGWRAVEDLYQVGE
ELVVHLAGVSTRELAEALVGLRVYAEVADLPPLEEGQYYYFALIGLPVYVEGQKVGEVAD
ILDAGAQDVLVIRGVGERLRDRAERLVPLQAPYVRVEAEGIHVEPIPGLFD
>tr|K7QWL8|K7QWL8_THEOS
---------MRLVEIGRFGAPYALKGGLRFRGEPVVLHLERVYVEGHGFRAVEDLYRVGE
VLILHLAGVSTRELAEALVGLRVYAEVEDLPPLEEGQYYYFALVGLPVYVGEEQVGEVAD
ILDAGAQDVLVIRGIGERLRDQRERLVPLQAPYVTVEEGRILVEPIPGLFD
```

###### free slicing

TO DO


#### Meta-data and statistics

* oMsa.shape: [sequenceNumber, columnNumber]


### PDB container


### Protein data container

#### Load a PDB file

```python
import pyproteinsext.structure.coordinates as PDB
parser = PDB.Parser()
pdbObj = parser.load(file="./1syq.pdb")
```

##### Display SEQRES
```python
pdbObj.SEQRES["A"]
```
`A` is chain name.

#### Aligning SEQRES and vald ATOM RECORD

##### Create wrapper peptide object

```python
import pyproteins.sequence.peptide as pep
p1 = {'id' : "SEQRES",
    'desc' : 'pdb file fasta translation',
    'seq' : pdbObj.SEQRES["A"]
}
pepSeqRes = pep.Entry(p1)
pepCoor = pep.Entry(pdbObj.chain("A").peptideSeed())
```

##### Align "peptide" sequences

```python
import pyproteins.alignment.nw_custom as N
import pyproteins.alignment.scoringFunctions as scoringFunctions
blosum = scoringFunctions.Needle().fScore
nw = N.nw(gapOpen=-10, gapExtend=-0.5, matchScorer=blosum)
aliResObj = nw.align(pepPDB, pepUniProt)
print(aliResObj)
```
Example with a sequence peptide from to a PDB file.
In this illustration, PDB file is *2vkn*.

```python
#parsing PDB
import pyproteinsext.structure.coordinates as PDB
parser = PDB.Parser()
pdbObj = parser.load(file="path/2ns7.pdb")
pdbObj.SEQRES
```
It’s possible to create a tuple with AA name and his position ; with command :
```python
import pyproteins.sequence.peptide as pep
AApdb = [(pep.threeToOne(aa.name), int(aa.num)) for aa in pdbObj.byres()]
```
`threeToOne` is a translater of AA name code at 3 letters to AA name code at 2 letters.

#### Uniprot container

You can access the content of any uniprot element. Corresponding XML file we ll be download locally if needed in a user defined cache directory.

```python
import pyproteinsext.uniprot as uniprot
uniprot.proxySetting(https="https://yourproxy:port", http="http://yourproxy:port")

uColl = uniprot.getUniprotCollection()
uColl.setCache(location='/Users/guillaumelaunay/work/data/uniprot')
uniprot.getPfamCollection().setCache(location='/Users/guillaumelaunay/work/data/pfam')

obj=uColl.get("P98160")


print(obj.GO)
print("\n")
print(obj.DI)
print("\n")
print(obj.peptideSeed())
print("\n")
print(obj.fasta)
```

will print

```bash
[GO:0005605:C:basal lamina{ECO:0000501}, ...]

[DI-02288:Schwartz-Jampel syndrome (SJS1) {Rare autosomal recessive disorder
characterized by permanent myotonia (prolonged failure of muscle relaxation)
and skeletal dysplasia, resulting in reduced stature, kyphoscoliosis,
bowing of the diaphyses and irregular epiphyses.},
...]

{'id': 'P98160', 'desc': 'PGBM_HUMAN', 'seq': 'MGWRAAGALLLALLLHGRLLAVTHGLRAYDGLSLPEDIETVTA...}

>P98160 PGBM_HUMAN
MGWRAAGALLLALLLHGRLLAVTHGLRAYDGLSLPED...
```

#### HMMR results container

You can give one or several file to parse method. For each protein, an entry `hmmrObj` is created

```python
import pyproteinsext.hmmrContainer as hm
hmmrContainer = hm.parse('hmmsearch_A.out', 'hmmsearch_B.out')
```

**hmmrObj** attributes and properties:  
* prot : protein name
* domain : domain name
* hit : dictionnary that contains hmmr hit informations, like score, evalue, alignment positions...
* sequence : give protein sequence that corresponds to domain
* start : give domain start in protein
* end : give domain end in protein

```python
for e in hmmrContainer:
    print(e.prot)
    print(e.domain)
    print(e.hit)
    print(e.sequence)
    print(e.start)
    print(e.end)
```
will print 
```text
tr|A0A1Q6ZBW6|A0A1Q6ZBW6_9ARCH
PF08022_full
{'hmmID': 'PF08022_full', 'aliID': 'tr|A0A1Q6ZBW6|A0A1Q6ZBW6_9ARCH', 'header': '1  score: 16.0 bits;  conditional E-value: 0.00014', 'score': '16.0', 'bias': '0.0', 'cEvalue': '0.00014', 'iEvalue': '0.24', 'hmmFrom': '26', 'hmmTo': '102', 'aliFrom': '37', 'aliTo': '99', 'envFrom': '27', 'envTo': '102', 'acc': '0.89', 'hmmStringLetters': 'fkwkpGqhvylsvpsisklllesHPFtiasapekddelslvirarggwtkrlaelaekseaesksklkvlieGPYGa', 'matchString': ' +++pGq+++++vp + ++     P ++  ++  ++e  +vi++ g ++++l e+              ++ GPYG+', 'aliStringLetters': 'HQANPGQFAMVWVPGVDEV-----PMSVLAIHG-KSEAGVVIKKGGPVSTALWEKKVGD--------IFFVRGPYGH', 'hmmSymbolStuff': {'CS': 'XXXXXXXXXXXXXXXXXXX--XXXXXXXXXXXX.XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX--TTSTTSH'}, 'aliSymbolStuff': {'PP': '5789*************88.....******999.******************9976555........69*******6'}}
HQANPGQFAMVWVPGVDEVPMSVLAIHGKSEAGVVIKKGGPVSTALWEKKVGDIFFVRGPYGH
37
99
...
```    

#### TMHMM results container

Container that store TMHMM results. You can only give one tmhmm result file to parse (for now).

```python
import pyproteinsext.tmhmmContainerFactory as tmhmm
tmhmmContainer = tmhmm.parse('tmhmm.out')
```
Container is a collection of TMHMM_Obj objects.

**TMHMM_Obj** attributes and properties:
* prot : protein name
* prot_length : protein length
* nb_helix : number of predicted helixes
* fragments : List of Fragment object. Each fragment have attributes cellular_location, start and end. 
* topology_seq : protein sequence with 'o' for outside loop, 'i' for inside loop and helix number for helixes. **WARNING : doesn't work if there are more than 9 helixes (use 2 characters instead of 1)**

```python
for e in tmhmmContainer:
    print(e.prot)
    print(e.prot_length)
    print(e.nb_helix)
    print(len(e.fragments),"fragments")
    for f in e.fragments:
        print(f.cellular_location, f.start, f.end)
```
will print 
```text
tr|A0A2E0DFV0|A0A2E0DFV0_9EURY
155
4
9 fragments
outside 1 5
TMhelix 6 25
inside 26 53
TMhelix 54 76
outside 77 90
TMhelix 91 109
inside 110 129
TMhelix 130 152
outside 153 155
...
```


#### Pfam container

### Protein-Protein interaction containers

#### psicquicData container API

Descriptions of the Psicquic service and related MITAB format can be found [here](https://psicquic.github.io/)

```python
import pyproteinsext.psicquic as psq
psqObj = psq.PSICQUIC(offLine=True)
psqObj.read(mitabFile)
psqAllInR6 = psqObj.filter(predicate=f)
```

#### mitabObject API

A __mitabObject__ stores a set of _psicquic.PSQDATA_ objects. Each _psicquic.PSQDATA_ object handles one mitab record. Managment of the _psicquic.PSQDATA_ set is done through the _pyproteins.container.Core.dnTree_ container model.

In practice, a __mitabObject__ can be used to

###### List of all interactors

```python
mitabTopologyObject.keys()
```

###### Obtain a one-level dictionary of all the partnairs of a query protein along with their list of _psicquic.PSQDATA_

```python
mitabTopologyObject["P38801"]
```

###### Obtain all the _psicquic.PSQDATA_ of a specific pair of proteins

```python
mitabTopologyObject["P38801"]["P24783"]
```

```text
[uniprotkb:P24783	uniprotkb:P38801	intact:EBI-5602|uniprotkb:Q05456|uniprotkb:D6W169	intact:EBI-1909|uniprotkb:D3DL32	psi-mi:dbp2_yeast(display_long)|uniprotkb:DBP2(gene name)|psi-mi:DBP2(display_short)|uniprotkb:YNL112W(locus name)|uniprotkb:N1945(orf name)|uniprotkb:DEAD box protein 2(gene name synonym)|uniprotkb:p68-like protein(gene name synonym)	psi-mi:lrp1_yeast(display_long)|uniprotkb:LRP1(gene name)|psi-mi:LRP1(display_short)|uniprotkb:Like an rRNA processing protein 1(gene name synonym)|uniprotkb:rRNA processing protein 47(gene name synonym)|uniprotkb:RRP47(gene name synonym)|uniprotkb:YC1D(gene name synonym)|uniprotkb:Yeast C1D domain-containing protein(gene name synonym)|uniprotkb:YHR081W(locus name)	psi-mi:"MI:0111"(dihydrofolate reductase reconstruction)	Tarassov et al. (2008)	pubmed:18467557|mint:MINT-6673767|imex:IM-14275	taxid:559292(yeast)|taxid:559292(Saccharomyces cerevisiae)	taxid:559292(yeast)|taxid:559292(Saccharomyces cerevisiae)	psi-mi:"MI:0915"(physical association)	psi-mi:"MI:0471"(MINT)	intact:EBI-6319091|imex:IM-14275-472	author score:99.00|intact-miscore:0.37]
```

The key order does not matter.

## Database modules

### DB FS

A [library](https://github.com/glaunay/pyproteinsExt) to build and query large database using the file system.

#### Building a database

##### Fetch [desired multifasta gziped file database](https://www.uniprot.org/downloads)

##### Split this file into smaller zip file

```sh
python uniprotFastaFS.py ~/tmp/vUniprot --cluster uniprot_sprot_2018_11.fasta.gz
```

Several `node_*.gz` multifasta files  were created.

##### Index the node files

Create a node folder foreach node file, uses the file system architecture to index their content. The `--nodes` argument is a regexp to specify the node number(s) to index.

```sh
python uniprotFastaFS.py ~/tmp/vUniprot --nodes '*'
```

All indexation processes are performed in parrallel and the resulting
subfolder organisations are independant. By independant folder architecture, we mean, as illustrated in the exemple below, that node subfolders (eg:`node_1` and `node_2`) may present similar _prefix_ have subfolders (eg:`A,B`).

```sh
vUniprot
    |____node_1
    | |____A
    | | |____index.txt
    | | |____data.gz
    | |____B
    | | |____index.txt
    | | |____data.gz
    |____node_2
    | |____A
    | | |____index.txt
    | | |____data.gz
    | |____B
    | | |____index.txt
    | | |____data.gz
```

#### Querying a database

You can fecth a fasta entry the following way, optionally dumping it to file.

```sh
python uniprotFastaFS.py ~/tmp/vUniprot --get P98160
```

#### Please consult CLI help for additional options
