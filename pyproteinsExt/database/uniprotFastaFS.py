"""
    Build a database of fasta entries using file system directory architecture
    Usage:
        uniprotFastaFS.py <location> [--cluster <bigFile> --radius <radius>]
        uniprotFastaFS.py <location> [--nodes <volNumExp> --threads <nWorker> --size <size>]
        uniprotFastaFS.py <location> [--input <file> --size <size> --auto]
        uniprotFastaFS.py <location> [--get <proteinID> --out <fileOut>]
        uniprotFastaFS.py <location> (-i | --info)
        uniprotFastaFS.py --view <file>

    Options:
        -h --help                             Show this screen
        -d <file>, --input <file>             MultiFasta input file
        -s <size>, --size <size>              Volume size [default: 500]
        -a, --auto                            Build without asking confirmation
        -i, --info                            Display Database statistics
        -v <file>, --view <file>              Display Input file statistics
        -g <proteinID>, --get <proteinID>     Identifier to retrieve data
        -o <fileOut>, --out <fileOut>         output file [default: stdout]
        -c <file>, --cluster <file>           archive file to dispatch to nodes
        -r <radius>, --radius <radius>        maximum number or entries per nodes
        -n <volNumExp>, --nodes <volNumExp>   Building the volumes of marching regExp  number nodes
        -t <nWorker>, --threads <nWorker>     Set the number of workers for parrallel building of node volumes

"""
#
# Wrapper to a file system base database of fasta sequences
# NB: insert shoukld 1st stat nodes to modify smallest one
import re
import os
import glob
import gzip
import sys
import shutil
from multiprocessing import Pool, TimeoutError
#import zipfile


from docopt import docopt

# iterate over entries found in file
# forEach insertID 
# then
# forEach load
def batchBuild(rootFolder, fileName, Nsize=500, alreadyIndexed=False):

    if not alreadyIndexed:
        scratchInsert(rootFolder, fileName, Nsize)

    scratchLoad(rootFolder, fileName)

    print("Scratch indexing and data dumping complete, archiving (Vsize=" + str(Nsize) + ")")
    a=0
    for iFile in glob.iglob(rootFolder + '/**/data.raw', recursive=True):
        with open(iFile, 'r') as fInput:
            a += 1
            jFile = iFile.replace('.raw','.gz')
            with gzip.open(jFile, 'wb') as fOutput:
                for l in fInput:    
                    fOutput.write(l.encode('utf-8'))
    
        os.remove(iFile)
    print("Archiving of " + str(a) + " volumes successfull")

def preview(fileName):
    c = 0
    for elem in fileCrawl(fileName):
        c += 1
    print("Input File features " + str(c) + " entries ")


def scratchInsert(rootFolder, fileName, Nsize):
    cnt = 0
    
    for elem in fileCrawl(fileName):
        _insertID(elem['id'], rootFolder, Nsize, BuildStep=True)
        cnt += 1

    print(str(cnt) + " elements inserted from scratch")

def scratchLoad(rootFolder, fileName):
    print("Adding content from scratch")
    for elem in fileCrawl(fileName):
        _load(rootFolder, elem['id'], elem['content'], BuildStep=True)

def stat(rootFolder):
    stat = []
    print("Computing stats w/ file pattern " + rootFolder + '/**/index.txt')
    t=0
    for iFile in glob.iglob(rootFolder + '/**/index.txt', recursive=True):
        #print (iFile)
        c = 0
        with open(iFile, 'r') as f:
            for l in f:
                c += 1
        stat.append((iFile[len(rootFolder)+1:].rstrip('/index.txt').replace('/',''), c))
        t += c
    print("Total indexed entries " + str(t))
    return sorted(stat, key=lambda x:x[1], reverse=True) 

## Open flat file eventually zip compress and pass it to fastaStrea Stream yielder
def fileCrawl(fileName):
    fOpen = open
    mode = 'r'
    _args = {}
    if fileName.endswith('gz'):
        fOpen = gzip.open
        mode = 'rt'
        _args = { "encoding" : 'utf-8' }
       #print("Reading from zipped source")
        
    with fOpen(fileName, mode, **_args) as f:
        fs = fastaStream(f)
        for d in fs:
            yield d

# Stream Read and iterate other each fasta record
class fastaStream(object):
    def __init__(self, iStream):
        self.iStream = iStream

    def __iter__(self):
        datum = {'id' : None, 'content' : None}
        for l in self.iStream:
            if l.startswith(">"):
                if datum['id']:
                    yield datum
            
                datum['id'] = l.split('|')[1]
                datum['content'] = l
            else:
                datum['content'] += l
        yield datum



# Given a collection of fasta and a splitLvl
# Create the required directories 
# Dump the corresponding gziped fasta file and index txt file
# WE LL SEE

def _load(workDir, _id, content, BuildStep=False):
    vDir, lvl, vName, gotcha = _getElemDir(_id, workDir)
    if not vDir:
        raise ValueError("No suitable volume for " + _id + " content")
    
    vInfo = _dirIndex(vDir, vName)
    elemInfo = [ (e[0], e[2]) for e in vInfo if e[0] == _id ] # Filter out vInfo of elemnt to load
    if len(elemInfo) > 1:
        raise TypeError("Volume key duplicate")

    if not elemInfo:
        raise ValueError("Matching volume does not hold key  " + _id)
    elemInfo = elemInfo[0]
    
    if elemInfo[1] == '1' and not BuildStep:
        print (elemInfo[0] + ' data already loaded' )
        return True

    updateVolumeDatum(vDir, _id, content, BuildStep=BuildStep)
    return True


def updateVolumeDatum(vDir, _id, content, BuildStep=False):
    dataFile = vDir + '/data.gz' if not BuildStep else vDir + '/data.raw'
# Build call -> aappend to flat file    
    if BuildStep:
        with open(dataFile, 'a') as f:
            f.write(content)
        return True

# Updating archived data
    toDump = content
    if os.path.isfile(dataFile):
        with gzip.open(dataFile, 'rt', encoding='utf-8') as fRead :
            toDump = ('').join([ l for l in fRead ])
        toDump += content    

    with gzip.open(dataFile, 'wb') as f:
            f.write(toDump.encode('utf-8'))
    
    # Updating index file
    rawInfo = []
    with open(vDir + '/index.txt', 'r') as f:
        for l in f:
            _buf = l.split()
            if _buf[0] ==  _id:
                _buf[1] = '1'
            rawInfo.append(_buf)
   
    with open(vDir + '/index.txt', 'w') as f:
       f.write( '\n'.join([ '\t'.join(e) for e in rawInfo ])  + '\n' )
    
    return True

# Create desired new directory w/ index file of elem provided list
def createVolume(volumeID, rootDir, elem=[], data=[]):
    newVolumePath = rootDir + '/' + ("/").join(volumeID)
    
    os.mkdir(newVolumePath)
    with open(newVolumePath + '/index.txt' , 'w') as f:
        for e in elem:
            f.write('\t'.join(e) + '\n')

    if data:
        with gzip.open(newVolumePath + '/data.gz', 'wb') as f:
            f.write('\n'.join(data).encode('utf-8'))
    
    return newVolumePath

## If queryKeys is provided returns single key dictionary
def _dataDict(dirPath, queryKey=None):
    if not os.path.isfile(dirPath + '/data.gz'):
        return None
    _d = {}
    
    for d in fileCrawl(dirPath + '/data.gz'):
        if queryKey:
            if d['id'] == queryKey:
                return { queryKey : d['content'] }
            continue
        _d[d['id']] = d['content']
    
    return _d

def _insertID(elem, rootDir, N=2, BuildStep=False):
   # if exists(elem, rootDir) and not Force:
   #     print (elem + ' already in database')
   #     return True
    
    
    ## Look for deepest possible
    dirPath, dirLvl, vName, gotcha = _getElemDir(elem, rootDir)
    ## Create it if none suitable
    if not dirPath:
        dirPath = createVolume(elem[0], rootDir) # Special case back to root
    ## Check the size of the volume were about to insert into
    #vInfo = _dirIndex(dirPath, vName)
    #if elem in [ e[0] for e in vInfo ]:
    if gotcha:
        print(elem + ' already in DB no ID to insert')
        return True

    vInfo = _dirIndex(dirPath, vName)
    #print ("About to insert id in a " + str(len(dirElem)) + ' elements index')
    ## If size limits distribute the current volume content into subvolumes
    if len(vInfo) >= N:
        #print("Volume "+ dirPath + ' reached size limit (' + str(N) + ') spliting it' )
        # LOADING Fasta Record
        heavyData = _dataDict(dirPath)
        if shotgun(rootDir, vInfo, dirLvl, N, vName, heavyData):
    ## Redo the insertion attempt
            #print("Shotgun completed attempting inserting again")
            _insertID(elem, rootDir, N)
            return True
        else :# Is elem to add in its max suffix volume
            if len(elem) > len(vName) + 1: # Too long to be indexed here
            #    print(str(vName) + " No spread out of shotgun => creating volume " + elem[:len(vName) + 1] + " to index " + elem)
                createVolume(elem[:len(vName) + 1], rootDir, elem=[elem])
                return True
            #print("No spread of " + vName +  " possible appending terminal suffix element " + elem)
#else: # We were not able to shotgun current index, bc all indexed elements have max len suffix
# We inject           
#        return True
    ## Actual insertion -> index update
    with open(dirPath + '/index.txt', 'a') as f:
        status = 0 if not BuildStep else 1
        f.write(elem + '\t' + str(status) + '\n')
    #print(elem + ' appened to ' + dirPath + '/index.txt')
    return True
# Separate element within current directory 
# which are already stored in their longest possible suffix => setA
# From the other => setB
# setB will be split in subdirectories
# setA will overwrite index file of current directory 
def shotgun(rootDir, elemList, dirLvl, N, cVolName, hRecord):
    #print( str(dirLvl) + ' ' + cVolName + " , LOADING " + str(elemList))
    small = []
    long = []
    for e in elemList:
        t = small if e[1] else long
        t.append( (e[0], e[2]) ) ## Storin name and load status

    #print('small ' + str(small) + ' // long' + str(long))
    
    if len(long) == 0: # Cnat do anything
        return False
    
    #print('Spreading ' + str(len(long)) + ' non max suffix element' )
    #print('XX ' + str(dirLvl))
    
    
    newLots = {}
    for e in long:
        k = cVolName + e[0][dirLvl]
        if not k in newLots:
            newLots[k] = []
        newLots[k].append(e)


    for k in newLots:
        newData = []
        if hRecord:
            newData = [ hRecord[ e[0] ]for e in newLots[k] if e[0] in hRecord ]
        createVolume(k, rootDir, elem=newLots[k], data=newData)
    # Erase previous index with max suffix name only
    upIndex = rootDir + '/' + ('/').join(cVolName) + '/index.txt'
    upData = rootDir + '/' + ('/').join(cVolName) + '/data.gz'
    if small:
        with open(upIndex, 'w') as f:
            for e in small:
                f.write('\t'.join(e) + '\n')
        
        if hRecord:
            newDataSmall = [ hRecord[ e[0] ] for e in small if e[0] in hRecord ]
            with gzip.open(upData, 'wb') as f:
                f.write('\n'.join(newDataSmall).encode('utf-8'))

    else:
        os.remove(upIndex)
        if os.path.isfile(upData):
            os.remove(upData)
    return True
        
    #print("Should add " + str(newLetters) + ' to ' + rootDir)
    

def get(elem, rootDir):
    
    dirPath, dirLvl, vName, gotcha = _getElemDir(elem, rootDir)
    if not dirPath:
        raise ValueError("No suitable volume for " + elem + " content")
    
    if not gotcha:
        return None

    d = _dataDict(dirPath)
    return d[elem]


def exists(elem, rootDir):
    dirPath, dirLvl, vName, gotcha = _getElemDir(elem, rootDir)
    return gotcha
 #  if not dirPath:
 #       return False
 #   return elem in [ e[0] for e in _dirIndex(dirPath, vName) ]

# "(elem, isLeave)"
# Eg : A/B/index.txt
# ABC is a leave, no deeper directory suffixin possible
# ABFE is not

def _dirIndex(cDir, vName):
    iFile = cDir + '/index.txt'
    
    if not os.path.isfile(iFile):
        #print('No index file at ' + cDir)
        return []
        #raise ValueError('No index file at ' + cDir)
    data = []
    with open(iFile, 'r') as f:
        for l in f:
            [eName, eLoad] = l.split()
            elemSuf = eName[:-1] # Getting rid of carraige return and last char --> longest possible suffix
            isLeave = str(elemSuf) == str(vName)
            data.append( (eName, isLeave, eLoad) )
    
    return data
    
# Provided a fastaElem,
# Pool version, 
# returns
# - the directory where it should be located (deepest possible)
# - the subtree level
# - the list of entries in this volume
def _getElemDir(_id, rootDir):
    #print("Looking for a volume for request " + _id)
    defVal = [None, 0, None, False]

    nodePath = glob.iglob(rootDir + '/node_*/') if glob.glob(rootDir + '/node_*/') else [ rootDir + '/' ]    
    hitPaths = [ _getElemDirNode(iNode, _id) for iNode in nodePath  ]
   
    
    for iDir in sorted( hitPaths, key=lambda x : x[1] ) : 
        if not iDir[0]: # Current node do not lay a single lvl suitable directory
            continue
        if iDir[1] > defVal[1]: ### No check on max size of current Node
            defVal = list(iDir) + [False]

        #print (iDir + " :: content")
        if _id in [ e[0] for e in _dirIndex(iDir[0], iDir[2]) ]:
            ans = list(iDir) + [True]
            return ans
   # [ for iDir in hitPaths ]
    #print ( sorted( hitPaths, key=lambda x : x[1] ) )
    return defVal

## Given one node folder tree get deepest matching directory path
def _getElemDirNode(rootDir, _id):
    offset = -1
    while len(_id) + offset > 0:
        cPath = rootDir + ("/").join(_id[:offset])
        #print(cPath)
        if os.path.isdir(cPath):
            lvl = len(_id) + offset
            #print("Got following volume location : " + cPath + ' at lvl ' + str(lvl))
            return (cPath, lvl, str(_id[:offset]))
        offset-=1
    
    return (None, 0, None)

# Will split an input file into subfile of maxSize content, and zip them
def setNodes(rootFolder, inputFile, maxSize):
    nodeCounter = 0
    fCurrent = None
    fCurrName = None

    def _wrap(c, fName, fOut):
        fOut.close()
        with open(fName, 'rb') as f_in:
            with gzip.open(rootFolder + "/node_" + str(c) + ".gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        #gzip.GzipFile(rootFolder + "/node_" + str(c) + ".gz", mode='w').write(fName)
        os.remove(fName)

    for eNum, elem in enumerate( fileCrawl(inputFile) ):
        if eNum % int(maxSize) == 0:
            
            if fCurrent: # close and zip the current one, then remove the flat version
               _wrap(nodeCounter, fCurrName, fCurrent)
            # Open the new current one
            nodeCounter += 1
            fCurrName = rootFolder + "/node_" + str(nodeCounter) + ".raw"
            fCurrent = open(fCurrName, "w" )
        fCurrent.write(str( elem['content']) )

    if fCurrent:
        _wrap(nodeCounter, fCurrName, fCurrent)
       
    print ("needs " + str(nodeCounter) + " volumes for maxSize of " + str(maxSize))
    

def pClusterBuild(inputPack):
    fileSource = inputPack[0]
    kwargs = inputPack[1]

    rootFolder = os.path.dirname(fileSource)
    nodeID = re.findall(r'[\d]+',os.path.basename(fileSource))[0]
    nodeRootFolder = rootFolder + '/node_' + str(nodeID)
    os.makedirs(nodeRootFolder)
    batchBuild(nodeRootFolder, fileSource, **kwargs)


if __name__ == "__main__":
    arguments = docopt(__doc__, version='uniprotFastaFSDB 2.0')
    #print(arguments)
    if arguments['--input']:
        if arguments['--auto']:
            print ("Automatic Building confirmation")
        else:
            answer = None
            while answer not in ("y", "n"):
                answer = input("Building database at " + arguments['<location>'] + " ? [y/n]")
                if answer == "y":
         # Do this.
                    pass
                elif answer == "n":
         # Do that.
                    pass
                else:
    	            print("Please enter y/n")

            if answer == "n":
                sys.exit(1)

        print( "Building database at " + arguments['<location>'] )
        batchBuild(arguments['<location>'], arguments['--input'], Nsize=int(arguments['--size']))

    if arguments['--info']:
        data = stat(arguments['<location>'])
        ff = None if arguments['--out'] == 'stdout' else open(arguments['--out'], "a")
        print(data, file=ff)

    if arguments['--get']:
        data = get(arguments['--get'], arguments['<location>'])
        ff = None if arguments['--out'] == 'stdout' else open(arguments['--out'], "a")
        print(data, file=ff)
    
    if arguments['--view']:
        preview(arguments['--view'])
    #arguments['location'])

    if arguments['--cluster']:
        print("Cluster mode")
        maxClusterIndex = arguments['--radius'] if arguments['--radius'] else 100000
        setNodes(arguments['<location>'], arguments['--cluster'], maxClusterIndex)
    
    if arguments['--nodes']:

        print (arguments)
        nWorker = arguments['--threads'] if arguments['--threads'] else 8
       
        
        opt = {}
        if arguments['--size']:
            opt['Nsize'] = int(arguments['--size'])
        
        reNumNode = arguments['--nodes'].replace('?', '.').replace('*', '.*')
        reNumNode = '_' + reNumNode + '\.gz$'
        #print(reNumNode)
        inputClusterIter =  [ (iFile, opt) for iFile in glob.iglob(arguments['<location>'] + '/node_*.gz' ) if re.search(reNumNode, iFile) ]
        print(inputClusterIter)
        with Pool(processes=nWorker) as pool:
            pool.map(pClusterBuild, inputClusterIter)
        
        print ("Done")
        #pClusterBuild(arguments['<location>'], number)

        '''
python ~/work/DVL/python3/pyproteinsExt/src/pyproteinsExt/database/uniprotFastaFS.py ~/tmp/vTrembl --node '??'

        '''