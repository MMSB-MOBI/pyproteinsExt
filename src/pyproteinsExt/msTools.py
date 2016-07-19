import re
import Uniprot
import copy
import math
from Matrisome import Matrisome
from Queue import Queue
from threading import Thread
from collections import OrderedDict


Matrisome = Matrisome()
uniprotEntrySet = Uniprot.EntrySet()
'''
Entry object attributes are defined as xls columns headers found in interact files

'''


def testEntry(uniProtObj):
    goWords = [
        { 'id' : 'GO:0016020', 'txt' : 'membrane'},
        { 'id' : 'GO:0005578', 'txt' : 'proteinaceous extracellular matrix'},
        { 'id' : 'GO:0005201', 'txt' : 'extracellular matrix structural constituent'},
        { 'id' : 'GO:0005604', 'txt' : 'basement membrane'},
        { 'id' : 'GO:0031012', 'txt' : 'extracellular matrix'},
        { 'id' : 'GO:0005605', 'txt' : 'basal lamina'},
        { 'id' : 'GO:0005576', 'txt' : 'extracellular region'},
        { 'id' : 'GO:0005615', 'txt' : 'extracellular space'}
    ]
    keyWords = [{ 'id' : 'KW-0472', 'txt' : 'Membrane'},
                { 'id' : 'KW-0272', 'txt' : 'Extracellular matrix'},
                { 'id' : 'KW-0084', 'txt' : 'Basement membrane'},
                { 'id' : 'KW-0964', 'txt' : 'Secreted'},
                { 'id' : 'KW-0339', 'txt' : 'Growth factor'},
                { 'id' : 'KW-0202', 'txt' : 'Cytokine'}
    ]
    badGoWords = [
        { 'id' : 'GO:0031965', 'txt' : 'nuclear membrane' },
        { 'id' : 'GO:0005634', 'txt' : 'nucleus' }
    ]

    badGoFound = []
    for goWord in badGoWords:
        if uniProtObj.hasGO(goWord['id']):
            #raise WordException(word = goWord['txt'])
            badGoFound.append(goWord['txt'])

    goFound = []
    for goWord in goWords:
        if uniProtObj.hasGO(goWord['id']):
            goFound.append(goWord['txt'])

    kwFound = []
    for kw in keyWords:
        if uniProtObj.hasKW(kw['txt']):
            kwFound.append(kw['txt'])

    if not kwFound and not goFound and not badGoFound : return None
    return {'goMatches' : goFound, 'kwMatches' : kwFound, 'badGoMatches' : badGoFound }


keymapRef = ["index", "probability", "spectrum", "precursor_neutral_mass", "massdiff", "specno", "mod", "peptide", "libra1", "libra2",
            "libra3", "libra4", "Symbol", "Accession", "Name", "Start_pre", "Start_mat", "Cleavage_site", "Nterm_annot", "Isoforms",
            "IAS", "nterm_rat_1", "nterm_rat_1_norm", "substr_p_1", "qcf_1"]

class EntrySet:

    def __iter__(self):
        for k, d in self.data.iteritems():
            yield d

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):

        l = [ d for k, d in self.data.iteritems()]
        if isinstance(index, int):
            #...    # process index as an integer
            return l[index]
        elif isinstance(index, slice):
            start, stop, step = index.indices(len(l))    # index is a slice
            return l[index]
            #...    # process slice
            #print start, stop, step
        else:
            raise TypeError("index must be int or slice")

    def __init__(self, fileName=None, entrySet=None, string=None, data=None, keymap=None, sigma=None):
        self.data = {}
        if entrySet:
            #self.keymap = copy.deepcopy(entrySet.keymap)
            #self.data = copy.deepcopy(entrySet.data)
            self.keymap = copy.deepcopy(entrySet.keymap)
            for d in entrySet:
                self.add( copy.deepcopy(d) )
            return
        if fileName:
            if not sigma:
                print "Warning you load a file w/out specifying normalisation parameters"
            buffer = None
            with open(fileName, 'r') as f:
                for line in f:
                    if re.match("^[\s]+$", line):
                        continue
                    if not buffer:
                        buffer = line
                    else:
                        buffer = buffer + line
            self._loadString(buffer)
        if data and keymap:
            self.keymap = keymap
            self.data = data
            return
        elif keymap:
            self.keymap = keymap
        elif string:
            self._loadString(buffer)

        if sigma:
            for e in self:
                e.normalize(sigma)

    @property
    def deepCount(self):
        cnt = len(self)
        for d in self:
            cnt += len(d.siblings)
        return cnt

    @property
    def clipCount(self):
        d = {}
        for e in self:
            clp = e.clipData
            if not clp:
                continue
            if clp['symbol'] in d:
                d[clp['symbol']].append(e)
            else:
                d[clp['symbol']] = [e]
        for k in d:
            print k + ' ' + str(len(d[k]))

        return d

    @property
    def peptides(self):
        def cmpSort (a, b):
            if not a.sigmaFold and not b.sigmaFold:
                return 0
            if not a.sigmaFold:
                return -1
            if not b.sigmaFold:
                return 1
            if max(a.sigmaFold) > max(b.sigmaFold):
                return 1
            if max(a.sigmaFold) == max(b.sigmaFold):
                return 0
            return -1

        l = [ v for k,v in self.data.iteritems() ]
        l.sort(cmp=cmpSort)
        lr = l[::-1]
        return [ [ d.getRefSeq(), d.sigmaFold ] for d in lr ]

    def __div__(a, b): # return a list of element peptideSeq and sigma ratio
        if not a.keymap == b.keymap:
            raise ValueError('Warning, Adding entrysets of different types\n' + str(a.keymap) + '\n' +  str(b.keymap) + '\n')

        results = []
        for e in a:
            eList = b.findall(e)
            if not eList:
                continue
            if len(eList) > 1:
                raise ValueError('unexpected number of peptide hit (' + str(len(eList)) + ') for sigma comparaison, should be 0 or 1')

            x = e['signal']
            y = eList[0]['signal']
            r = None
            if x and y:
                r = math.ceil((x/y)*1000)/1000 if y != 0 else 0
            results.append({'e' : [e, eList[0]], 'peptide' : e.getRefSeq(),
                'ratio' : r, 'max' : x if x > y else y })

        results.sort(cmp=lambda x,y:cmp(x['ratio'], y['ratio']))
        return results

    def __or__(a, b): # non redundant union
        if not a.keymap == b.keymap:
            raise ValueError('Warning, Adding entrysets of different types\n' + str(a.keymap) + '\n' +  str(b.keymap) + '\n')

        s = EntrySet(entrySet=a)
        cnt = 0
        for e in b:
            eMatch = s.findall(entry=e)
            if len(eMatch) > 1 :
                raise ValueError('Warning, Found multiple instance of identical sequences in stack')
            if len(eMatch) == 1 :
                cnt += 1
                eMatch[0].addSiblings(e)
            else :
                s.data[hash(e)] = e

        print str(cnt) + ' peptides discarded, ' + str(len(s))+ ' retained'
        return s

    def __and__(a, b): # intersection
        if not a.keymap == b.keymap:
            raise ValueError('Comparing entrysets of different types')

        eSet = EntrySet(data = {}, keymap=copy.deepcopy(a.keymap))
        for k in a.data :
            if k in b.data:
                eA = copy.deepcopy(a.data[k])
                eSet.add(eA)
                eB = b.data[k]
                if eA != eB:
                    eA.addSiblings(b.data[k])

        return eSet

    def __xor__(a, b): # immplementation is obviously non optimal
        if not a.keymap == b.keymap:
            raise ValueError('Comparing entrysets of different types')

        eSet = EntrySet(data = {}, keymap=copy.deepcopy(a.keymap))
        for k in a.data :
            if not k in b.data:
                eA = copy.deepcopy(a.data[k])
                eSet.add(eA)

        for k in b.data :
            if not k in a.data:
                eB = copy.deepcopy(b.data[k])
                eSet.add(eB)

        return eSet


    def __eq__(a, b):
        status =  True if len (set(a.data.keys()) & set(b.data.keys())) == 0 else False
        return status

    def clone (self, fullPeptideHash=False): # return a copy of the set without siblings references
        eSet = EntrySet(data={}, keymap = copy.deepcopy(self.keymap))
        for e in self:
            eNew = copy.deepcopy(e)
            eNew.siblings = []
            if fullPeptideHash:
                eNew.mark()
            eSet.add(eNew)
        return eSet

    def findall(self, entry=None, peptide=None):
        if not entry and not peptide and siblingsAbove is None:
            raise ValueError('undefined entry object to look for ?')
        eHits = []
        if entry:
                e = self.has(entry)
                if e:
                    eHits.append(e)
        elif peptide :
            for e in self:
                if e.getRefSeq() == peptide:
                    eHits.append(e)
        return eHits

    def has(self, entry):
        if not entry: return None
        if hash(entry) in self.data:
            return self.data[hash(entry)]
        return None

    def _threadPool(self):
        num_worker_threads = 10
        def worker():
            while True:
                e = q.get()
                e._boundUniprot()
                q.task_done()

        q = Queue()
        for i in range(num_worker_threads):
            t = Thread(target=worker)
            t.daemon = True
            t.start()

        for e in self:
            q.put(e)


        print '*** Main thread waiting ' + str(q.qsize())
        q.join()
        print '*** Done'
    # block until all tasks are done

    # Here we modify the content of each entry, typically b4 tsv dump
    # Carefull, siblings are not modified
    def enrich(self, termTypes):
        if not termTypes:
            raise ValueError("single term or list of enrichment type required")

        if "ECM" in termTypes:
                self.keymap.extend(['Matrisome', 'ECMkwMatches', 'ECMgoMatches','ECMbadGoMatches'])
        elif "matrisome" in termTypes:
              self.keymap.append("Matrisome")

        for e in self:
            if "ECM" in termTypes:
                e.isECM(annot=True)
            elif "matrisome" in termTypes:
                e.isMatrisome(annot=True)



    def filter(self, verbose=True, cleavageSite=None, annotatedAs=None, signalOver=None, signalBelow=None, sigmaFoldOver=None, sigmaFoldBelow=None, evalType="average", siblingsAbove=None):
        eSet = EntrySet(data = {}, keymap=copy.deepcopy(self.keymap))

        if siblingsAbove is not None:
            for e in self:
                if len(e.siblings) > siblingsAbove:
                    eSet.add(copy.deepcopy(e))

        if annotatedAs is "matrisome":
            for datum in self:
                if datum.isMatrisome():
                    eSet.add(copy.deepcopy(datum))
        elif annotatedAs is "ECM":
            self._threadPool()
            for datum in self:
                if datum.isECM():
                    if verbose:
                        print datum.getRefSeq()
                    eSet.add(copy.deepcopy(datum))
        elif signalOver or signalBelow or sigmaFoldOver or sigmaFoldBelow:
            if not 'nterm_rat_1' in self.keymap:
                raise ValueError('no terminal ratio to assess signal')
            if signalOver:
                lo_tresh = float(signalOver)
            if signalBelow:
                up_tresh = float(signalBelow)
            if sigmaFoldOver:
                lo_tresh = float(sigmaFoldOver)
            if sigmaFoldBelow:
                up_tresh = float(sigmaFoldBelow)

            for datum in self:
                v = datum.signal if signalBelow or signalOver else datum.sigmaFold
                if not v:
                    continue
                ok = True
                x = sum(v) / len(v)
                if evalType == "max":
                    x = max(v)
                elif evalType == "min":
                    x = min(v)

                if signalOver or sigmaFoldOver:
                    if x < lo_tresh:
                        ok = False
                if signalBelow or sigmaFoldBelow:
                    if x > up_tresh:
                        ok = False

                if ok:
                    eSet.add(copy.deepcopy(datum))

        elif cleavageSite:  # Should have 4 values, inDomain//Nter//Cter//hinge
            '''
            cleavageSite must be an array of strings ()

            outDomain on Nter
            outDomain on Cter
            outDomain in hinges
            1st letter status of Cleavage site left position
            2nd letter status of Cleavage site right position
            N: not in domain and in Nter protein region
            C: not in domain and in Cter protein region
            H: not in domain and in an hinge (between two domains) protein region
            d: in domain

            NN, CC, Nd, dC, Hd, dH

            '''
            for e in self:
                clipData = e.clipData
                if not clipData:
                    continue
                if clipData['symbol'] in cleavageSite:
                    eSet.add(copy.deepcopy(e))

        return eSet

    def listUniprotID(self):
        for datum in self:
            yield datum['uniprotID']

    def __repr__(self):
        string = "N terminome entry set\nkeymap:\n" + str(self.keymap) + "\ndataArray\n"
        for e in self:
            string = string +  str(e._asList())
        return string + "\n"

    def tsvDump(self, enrichWith=None, deep=False, fileName=None): ## union case
        if enrichWith:
            term = []
            if isinstance(enrichWith, basestring):
                term.append(enrichWith)
            elif isinstance(enrichWith, list):
                term.extend(enrichWith)
            self.enrich(term)

        asString = ("\t").join(self.keymap) + "\n"
        for e in self:
            asString += e.join("\t") + "\n"
            if deep:
                for eS in e.siblings:
                    asString += eS.join("\t") + "\n"
        if fileName:
            with open(fileName, "w") as text_file:
                text_file.write(asString)
            print "tsv content written to \"" + fileName + "\""
            return

        return asString

    def bingoDump(self, fileName='bingo.cy'): ## union case
        asString = ("\n").join(self.listUniprotID()) + "\n"
        with open(fileName, "w") as text_file:
            text_file.write(asString)

    def _loadString(self, rawData):
        buffer = rawData.splitlines()
        self.keymap = buffer.pop(0).split("\t")
        for line in buffer:
            eNew = Entry(self.keymap, line)
            self.add(eNew)

    def add(self, eNew): # check hash if present add to element sibling
            e = self.has(eNew)
            if e:
                e.addSiblings(eNew)
            else :
                self.data[hash(eNew)] = eNew

    def amineRatio(self):
        d = {'nAminoSub' : 0, 'nAminoGrp' : 0, 'nEpsilonSub' : 0}
        for e in self:
            aminoSubList = e.getPtmList(ptmFilter=['145.11', '272.20'])
            epsilonAminoSubList = e.getPtmList(ptmFilter=['272.20'])
            nAminoGrp = e.count(type="amineGroup")
            d['nAminoSub'] = d['nAminoSub'] + len(aminoSubList)
            d['nEpsilonSub'] = d['nEpsilonSub'] + len(epsilonAminoSubList)
            d['nAminoGrp'] = d['nAminoGrp'] + nAminoGrp
        return d

'''

'''
class Entry:
    def __init__(self, keymap, rawData, fullPeptideHash=False):
        self.data = rawData.split("\t")
        self.keymap = list(keymap)

        if not fullPeptideHash :
            string = re.sub('\[[^\]]+\]', '', self['peptide'])
            s1=re.sub("^[^\.]{0,1}\.n{0,1}","", string)
            self._refSeq = s1
        else :
            self._refSeq = self['peptide']

        self._uniprotBound = None
        self._sigmaFold = None
        self.siblings = [] # Stores Entry reference of merged objects

    def __hash__(self):
        return hash(self._refSeq)
#    def __cmp__(self):
#        return object.__cmp__(self)


    @property
    def uniprot(self):
        if self._boundUniprot():
            return self._uniprotBound
        return None

    @property
    def clipData(self):
        if not self['Cleavage_site']:
            return None
        if not self.uniprot:
            return None

        m = re.search("\(([\d]+)\)\.\(([\d]+)\)", self['Cleavage_site'])
        #print m.group()
        if not m:
            return None
        clipData = { 'symbol' : '' }

        for i,e in enumerate (["Nter", "Cter"]):
            try :
                clipData[e] = self.uniprot.pos( int(m.group(i + 1)) )
            except ValueError as msg:
                clipData[e] = None
            except IndexError as msg:
                clipData[e] = None

            if not clipData[e] : # append left born of cleavage site being 0 ..
                clipData["symbol"] += 'n'
            elif 'description' in clipData[e].domain   :
                clipData["symbol"] += 'd'
            elif 'outOfBounds' in clipData[e].domain:
                if clipData[e].domain['outOfBounds'] == 'Cter':
                    clipData["symbol"] += 'C'
                elif clipData[e].domain['outOfBounds'] == 'Nter':
                    clipData["symbol"] += 'N'
                elif clipData[e].domain['outOfBounds'] == 'hinge':
                    clipData["symbol"] += 'h'
                else:
                    raise ValueError("Unexpected domain status symbol in " + str(clipData[e]))
            else:
                raise ValueError("Unexpected domain status symbol in " + str(clipData[e]))

        if not clipData["Cter"] and not clipData["Nter"]:
            return None

        return clipData

    @property
    def sigmaFold(self):
        if not self.siblings:
            if self._sigmaFold:
                return [self._sigmaFold]
            else:
                return None
        return [e._sigmaFold for e in self.flattenSiblingsOut() if e._sigmaFold]

    @property
    def signal(self):
        l = []

        try :
            s = self['nterm_rat_1_norm'].replace(",", ".")
            v = float(s)
            l.append(v)
        except ValueError as err:
            pass
        for e in self.siblings:
            try :
                s = e['nterm_rat_1_norm'].replace(",", ".")
                v = float(s)
                l.append(v)
            except ValueError as err:
                pass
        return l

    def addSiblings(self, value):
        if self == value:
            return

        eSet = self.siblings

        if isinstance (value, list):
            eSet.extend(value)
        elif value:
            eSet.append(value)
        self.siblings = list(OrderedDict.fromkeys(eSet))

    def flattenSiblingsOut(self) :
        l = [self]
        for e in self.siblings:
            l.extend(e.flattenSiblingsOut())
        m = list(OrderedDict.fromkeys(l))
        return m

    def __eq__(a, b):
        if a.keymap != b.keymap:
            return False
        for k in a.keymap:
            if a[k] != b[k]:
                return False
        return True

    def normalize(self, sigma):
        value = pow(2,sigma)
        self._sigmaFold = self.signal[0] / value if self.signal else None

    def _boundUniprot(self):
        if not self._uniprotBound:
            id = self['uniprotID']
            if not id:
                return False
            try :
                self._uniprotBound = uniprotEntrySet.get(id)#Uniprot.Entry(id)
            except TypeError as msg:
                print "Following Ms element could not be bound to uniprot entity, reason " + str(msg)
                print self
                return False
        return True

    def isMatrisome(self, annot=False):
        matrisomeTerm = Matrisome.get(uniprotID=self['uniprotID'])
        if annot:
            self['matrisomeTerm'] = matrisomeTerm[0]['Category'] + ":" + matrisomeTerm[0]['Division'] if matrisomeTerm else 'NA'
        if matrisomeTerm:
            return True
        return False

    def isECM(self, annot=False):
        if not self._boundUniprot():
            return False
        matches = testEntry(self._uniprotBound)
        isMatrisome = self.isMatrisome(annot=annot)

        if annot:
            if not matches:
                self['ECMkwMatches'] = "NA"
                self['ECMgoMatches'] = "NA"
                self['ECMbadGoMatches'] = "NA"
            else :
                self['ECMkwMatches'] = (',').join(matches['kwMatches']) if matches['kwMatches'] else 'NA'
                self['ECMgoMatches'] = (',').join(matches['goMatches']) if matches['goMatches'] else 'NA'
                self['ECMbadGoMatches'] = (',').join(matches['badGoMatches']) if matches['badGoMatches'] else 'NA'
        if isMatrisome:
            return True
        if not matches:
            return False
        if matches['kwMatches'] or matches['goMatches']:
            return True
        return False
# TO DO SMTG WITH NuCLEUS related keywords
        #e['badGoMatches'] = (',').join(matches['badGoMatches'])

    def __repr__(self):
        asDict = self._asDict()
        return str(asDict)

    def __len__(self):
        s = self.unmark(copy=True)
        return len(s)

    def __getitem__(self, key):


        if key == "uniprotID":
            return self._guessUniprotID()

        if key not in self.keymap: # Raise Error here
            raise ValueError("Improper key \"" + key + "\"\n" + str(self.keymap))
        return self.data[self.keymap.index(key)]

    def __setitem__(self, key, value):
        #self.vector[key] = item
        if not key in self.keymap:
            self.keymap.append(key)
            self.data.append(value)
        else :
            i = self.keymap.index(key)
            self.data[i] = self.data[i] + " " + value

    def join (self, FS):
        if not FS:
            FS = ","
        return FS.join(self.data)

    def getRefSeq(self):
        return self._refSeq

    def unmark(self, copy=False):
        #/\[[^\]+]\]//g
        string = re.sub('\[[^\]]+\]', '', self['peptide'])
        s1=re.sub("^[^\.]{0,1}\.n{0,1}","", string)
        #s2=re.sub(".[^\.]{0,1}$","", s1)
        if copy:
            return s1
        self._refSeq = s1

    def getPtmList(self, ptmFilter=None):
        ptmList = self._listPtm()
        if not ptmFilter:
            return ptmList

        subPtmList = [ ptm for ptm in ptmList if ptm in ptmFilter ]
        return subPtmList

    def _listPtm(self):
        m = re.findall("\[([^\]]+)\]", self['peptide'])
        return m

    def count(self,type=None):
        if not type:
            seq = self.unmark(copy=True)
            return len(seq) - 1
        if type == "amineGroup":
            cnt = self._countAmine()
            return cnt
        raise ValueError("Unexpected type to look/count for in sequences" + "\"" + str(type) + "\"")

    def _countAmine(self):
        seq = self.unmark(copy=True)
        seq = re.sub("\.[^\.]$","", seq)
        subn = re.findall("K", seq)

        print seq + " => " + str(len(subn))

        return len(subn)

    def mark(self): # changes the hash code
        self._refSeq = self['peptide']

    def _guessUniprotID(self):
        if 'Accession' in self.keymap:
            return self['Accession']
        if 'protein' in self.keymap:
            m = re.search("^(sp|tr)\|([^\|]+)|", self['protein'])
            if m:
                return m.group(2)
        return None
    def _asDict (self):
        asDict = { self.keymap[i] : d for i, d in enumerate(self.data) }
        return asDict

    def _asList (self):
        return self.data

class Sequence:
    def __init__(self, rawString):
        pass
    def __eq__(a,b):
        pass
