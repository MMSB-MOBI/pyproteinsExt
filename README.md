# pyproteinsExt

## notebook test&examples

TO DO

## Generic container modules

TO DO

## Anotation modules

TO DO

## Specific container modules

### PDB container


### Protein data container

#### Load a PDB file

```python
import pyproteinsExt.structure.coordinates as PDB
parser = PDB.Parser()
pdbObj = parser.load(file="./1syq.pdb")
```

##### Display SEQRES
```python
pdbObj.SEQRES["A"]
```

#### Aligning SEQRES and vald ATOM RECORD

##### Create wrapper peptide object

```python
import pyproteins
p1 = {'id' : "SEQRES",
    'desc' : 'pdb file fasta translation',
    'seq' : pdbObj.SEQRES["A"]
}
pepSeqRes = pyproteins.sequence.peptide.Entry(p1)
pepCoor = pyproteins.sequence.peptide.Entry(pdbObj.chain("A").peptideSeed())
```

##### Align "peptide" sequences

```python
nw = pyproteins.alignment.nw_custom.nw(gapOpen=-10, gapExtend=-0.5)
aliResObj = nw.align(pepSeqRes, pepCoor)
print(aliResObjnw)
```

TO DO

#### Uniprot container

TO DO

```python
import json
uObj = uColl.get('Q8DR57')
json.dumps(uObj, cls=EntryEncoder)
```

#### Pfam container

### Protein-Protein interaction containers

#### psicquicData container API

TO DO

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