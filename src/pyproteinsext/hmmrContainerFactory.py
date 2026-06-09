import re, io
import gzip

reInstance = re.compile('(\# hmm([\S]+) :: search ([\S]+) against a ([\S]+) database.+?(?=\[ok\])\[ok\])', re.DOTALL)
reType = re.compile('^\# hmm([\S]+) :: search ([\S]+) against a ([\S]+) database[\s]*$')
# hmmsearch :: search profile(s) against a sequence database
reQuery = re.compile('Query:[\s]+([\S].*)[\s]+\[[LM]=[\d]+\]')
reStart = re.compile('^Scores for complete sequences \(score includes all domains\):')
reSum = re.compile('(([\d\.e-]+[\s]+){7})([0-9]+)[\s]+([\S]+)[\s]+([\S].*)$')
reInclude = re.compile('------ inclusion threshold ------')
reDetail = re.compile('Domain annotation for each ([\S]+) \(and alignments\):')
reStop = re.compile('Internal pipeline statistics summary')
reAliLongName = re.compile('^>>[\s]*([\S]+.*)$')
reDomainDetail = re.compile('  == domain [\d]+(?!.* == domain [\d])')
reAlignDef = re.compile('^[\s]*[\d]+[\s]*[\S][\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]]{2}[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]]{2}[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]{2}[\s]+([\.\d]+)')
reBothStringLetter = re.compile('^([\s]+([\S]+)[\s]+([\S]+)[\s]+)([\S]+)[\s]+([\S]+)')
reBothStringSymbol = re.compile('^([\s]+)([\S]+)')
reEmpty = re.compile('^[\s]*$')
reStrandsLine = re.compile("^([\s]*([\S]+)[\s]+([\d]+)[\s])([\S]+)[\s]([\d]+)[\s]*$")

#reInstance = re.compile('(# hmmsearch.+?(?=\[ok\])\[ok\])', re.DOTALL)

def parse(*inputFiles):
    """
    Parse one or several hmmr results file and store them into container. 

    :param *inputFiles: One or several path to input files, comma-separated
    :return: hmmrContainerFactory.Container
    """

    upType = None
    bigBuffer = ''
    mainContainer = Container()
    for inputFile in inputFiles:
        if inputFile:
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

            with f:
                for l in f:
                    bigBuffer += l

                    #print(l)

                runReports = reInstance.findall(bigBuffer)
                #print(runReports)
                if not runReports:
                    raise ValueError("Could not parse upper level in provided hmm(scan/search) file")

                for d in runReports:
                    mainContainer.addParsing(Container(input=io.StringIO(d[0]), upType=d[1]))
    return mainContainer


class Container(object):
    """
    Container that stores hmmr results under different data structure

    :param hmmrEntries: 
    :param dIndex:
    :param pIndex:
    """

    def __init__(self, input=None, upType=None):
       # self.queryHmmFile=None             #../../../data/fad_binding.hmm
       # self.targetSequenceDatabase=None   #./Trembl_50.fasta
       # self.query=None                    #PF08022_full  [M=104]
        self.upType = upType
        if not input:
            self.dIndex={}
            self.pIndex={}
            self.hmmrEntries=set()
        else:
            self.dIndex,self.pIndex,self.hmmrEntries = _parseBuffer(input)


    def addParsing(self, other):
        if self.upType != other.upType:
            if self.upType == None:
                self.upType = other.upType
            elif other.upType == None:
                other.upType = self.upType
        if self.upType != other.upType:
            raise TypeError('Adding', self.upType,' and ', other.upType)

        #self.summary.extend(other.summary)
        #self.details.extend(other.details)
        self.updateParsing(other)
        return self

    def updateParsing(self,other):
        for d in other.dIndex:
            if d not in self.dIndex:
                self.dIndex[d]=set()
            self.dIndex[d].update(other.dIndex[d])

        for p in other.pIndex:
            if p not in self.pIndex:
                self.pIndex[p]=set()
            self.pIndex[p].update(other.pIndex[p])

        self.hmmrEntries.update(other.hmmrEntries)

    def addIndex(self,dic,index,match):
        if index not in dic:
            dic[index]=set()
        dic[index].add(match)

    def __iter__(self):
        '''Iter through matches'''
        for d in self.dIndex:
            for m in self.dIndex[d]:
                yield m

    def add(self,*matches):
        for match in matches:
            prot_id=match.aliShortID
            dom_id=set([hit.hmmID for hit in match.data])
            if len(dom_id)>1:
                raise Exception("dom_id can't have length > 1. Check your code.")
            dom_id=dom_id.pop()
            self.addIndex(self.pIndex,prot_id,match)
            self.addIndex(self.dIndex,dom_id,match)

    def filterProteins(self,fPredicat,**kwargs):
        new_container=Container()
        for protein in self.pIndex:
            matches=self.pIndex[protein]
            if fPredicat(matches,**kwargs):
                new_container.add(*matches)
        return new_container

def _parseBuffer(input):
    dIndex={}
    pIndex={}
    hmmrEntries=set()
    inclusionBool = True
    readBool = False
    detailBool = False
    detailBuffer = []
    summaryBool = False

    queryID = None
    emptyLinesInRow = 0

    for l in input:
       # print(l)
        if not queryID:
            m = reQuery.search(l)
            if m:
                queryID = m.groups()[0]

        if reStart.search(l):
            readBool = True
            continue

        if re.search('\-\-\- full sequence \-\-\-   \-\-\- best 1 domain \-\-\-    \-#dom\-', l):
            #print ("UP")
            inclusionBool = True
            continue

        if reDetail.search(l):
            readBool = False
            detailBool  = True
            continue
        #if reStop.search(l):
        #    break

        if reEmpty.match(l):
            emptyLinesInRow += 1
            if emptyLinesInRow > 1:
               # print('tweet')
                readBool = False
                detailBool = False
        else:
            emptyLinesInRow = 0

        if detailBool:
            if l.startswith(">>"):
                detailBuffer.append(l)
            else :
                if len(detailBuffer) == 0:
                    if re.search('No targets detected that satisfy reporting thresholds', l):
                        print ('No hit found')
                        return ([],[])
                    continue
                detailBuffer[-1] += l

    if not queryID:
        raise ValueError('Could not find query identifier in header')
    #print(detailBuffer)
    dIndex[queryID]=set()
    for rawData in detailBuffer:
        #print(rawData)
        match=Match(rawData,queryID)
        subjctID=match.aliShortID
        if not subjctID:
            continue
        if subjctID not in pIndex:
            pIndex[subjctID]=set()
        pIndex[subjctID].add(match)
        dIndex[queryID].add(match)
        #print(match.data)
        if match.data:
            for hit in match.data:
                domain=hit.hmmID
                prot=hit.aliID
                hmmrEntries.add(HMMObj(prot,domain,hit))
        else:
            hmmrEntries.add(HMMObj(subjctID,None,None))
    #details = [ Match(rawData, queryID) for rawData in detailBuffer ]
    return (dIndex,pIndex,hmmrEntries)

class Match (object):
    def __init__(self, bufferString, queryID):
        self.aliLongID = None
        self.aliShortID = None
        self.hmmID = None
        self.queryID = queryID
        self.parseDetailEntry(bufferString)

    #def __hash__(self):
    #    return hash(self.aliLongID + self.hmmID)

    def parseDetailEntry(self, bufferString):
        buff = bufferString.split('  == domain ')
        self.data = []
        if len(buff) == 1:
            print ("Warning:: " + buff[0])
            self.aliLongID=buff[0].split("\n")[0].lstrip(">> ")
            self.aliShortID=self.aliLongID.split(" ")[0]
            self.data =[]
            return

        for u in buff[0].split('\n'):
            m = reAliLongName.search(u)
            if m:
                self.aliLongID = m.groups()[0]

            m = reAlignDef.search(u)
            if m :
                self.data.append( Alignment(m.groups()) )

        if len( self.data ) != len(buff[1:]):
            raise ValueError('uneven alignment definitions (' + str(len(self.data)) + ') and details (' + str(len(buff[1:])) + ')\n\n' + str(buff))

        for alignmentObject, rawAli in zip( self.data, buff[1:]):
            alignmentObject.parse(rawAli)

        self.aliShortID =  self.data[0].aliID
        self.hmmID = self.data[0].hmmID


    def __repr__(self):
        return str( self.__dict__ )

    def _repr_html_(self):
        htmlContent = ''
        for d in self.data:
            htmlContent += '<table style="border:none">'
            for name, x, y, field in [ (self.hmmID, d.hmmFrom, d.hmmTo, d.hmmStringLetters), ("", "", "", d.matchString), (self.aliShortID, d.aliFrom, d.aliTo, d.aliStringLetters) ]:
                htmlContent += '<tr style="border:none"><td style="border:none">' + name + '</td>' \
                                + '<td style="border:none; color:red; font-size:0.8em">'+ x + '</td><td style="border:none; color:red; font-size:0.8em">' + y + '</td>' \
                                + '<td align="center" style="min-width:20px; text-align:center;border:1px solid rgba(0, 0, 0, .05);">' \
                                + '</td><td align="center" style="min-width:20px; text-align:center;border: 1px solid rgba(0, 0, 0, .05);">'.join(field) +  '</td></tr>'
            htmlContent += '</table>'

        return htmlContent

class Alignment(object):

    def __init__(self, args):

        self.hmmID = None
        self.aliID = None
        self.header = None
        self.score   = args[0]
        self.bias    = args[1]
        self.cEvalue = args[2]
        self.iEvalue = args[3]
        self.hmmFrom = args[4]
        self.hmmTo   = args[5]
        self.aliFrom = args[6]
        self.aliTo   = args[7]
        self.envFrom = args[8]
        self.envTo   = args[9]
        self.acc     = args[10]

        #self.hmmStringSymbols = ''
        self.hmmStringLetters = ''
        self.matchString      = ''
        self.aliStringLetters = ''
        self.hmmSymbolStuff = {}
        self.aliSymbolStuff = {}


    def parse(self, bufferString):

        buf = [ l for l in bufferString.split("\n") ] # if l != ''
        self.header = buf.pop(0)

        #print ('Parsing this:', buf)


        hitLines = [ i for i,l in enumerate(buf) if reStrandsLine.match(l) ]
       # print (hitLines)
       # print(len(hitLines))



        for i in range(0, len(hitLines), 2):
            # Update offset values
            m = reStrandsLine.match(buf[ hitLines[i] ])
            n = reStrandsLine.match(buf[ hitLines[i + 1] ])
            self.aliID = n.groups()[1]
            self.hmmID = m.groups()[1]

#            print ('strand parsing', m.groups())
            offset = len(m.groups()[0])
            length = len(m.groups()[3])
            #print('Line offset is', offset, 'char chunk length is', length)
            reMatchString = r".{" + str(offset) + r"}(.{" + str(length) + r"})"
            mMatch = re.match(reMatchString, buf[ hitLines[i] + 1 ])
            if not mMatch:
                raise ValueError("Could not parse alignment match record", buf[ hitLines[i] + 1 ], "with", reMatchString)
            self.matchString += mMatch.groups()[0]
            reStuffString = r".{" + str(offset) + r"}(.{" + str(length) + r"}) ([A-Z]{2})"
        # i, j a row number, for i we go up and assign to profile data
        #                    for j we go down and assigne to ali data
            self.hmmStringLetters, self.hmmSymbolStuff = plumbParse(buf, hitLines[i], self.hmmStringLetters, self.hmmSymbolStuff, reStuffString,goingUp=True)
            self.aliStringLetters, self.aliSymbolStuff = plumbParse(buf,  hitLines[i + 1], self.aliStringLetters, self.aliSymbolStuff, reStuffString, goingUp=False)
            #print('---->x', self.hmmStringLetters,'<-----x')

    def _parse(self, bufferString):

        buf = [ l for l in bufferString.split("\n") if l != '' ]
        self.header = buf.pop(0)

        iOffset=5
        if len(buf)%5 != 0:
            if len(buf)%4 != 0:
                raise ValueError('unexpected detail alignment data (n=' +  str(len(buf)) +'):\n\n'+str(buf))
            else:
                iOffset=4

        i = 0
        while i < len(buf):
            m = reBothStringSymbol.match(buf[i])
            if iOffset == 5:
                if m:
                    self.hmmStringSymbols +=  m.groups()[1]
            j = i + 1 if iOffset == 5 else i
            m = reBothStringLetter.match(buf[j])
            stringChop = ''
            if m:
                offset = m.groups()[0]
                self.hmmID = m.groups()[1]
                stringChop = m.groups()[3]
                self.hmmStringLetters += stringChop


            j = i + 2 if iOffset == 5 else i + 1
            blank = ' ' * len(offset)
            my_regex = r"^" + re.escape(blank) + r"(.{"+ str(len(stringChop)) + r"})"
            m = re.match(my_regex, buf[j], re.IGNORECASE)

            if m:
                self.matchString += m.groups()[0]

            j = i + 3 if iOffset == 5 else i + 2
            m = reBothStringLetter.match(buf[j])
            stringChop = ''
            if m:
                self.aliID = m.groups()[1]
                stringChop = m.groups()[3]
                self.aliStringLetters += stringChop
            j = i + 4 if iOffset == 5 else i + 3
            m = reBothStringSymbol.match(buf[j])
            if m:
                offset = m.groups()[0]
                self.aliStringSymbols +=  m.groups()[1]

            i += iOffset

        #if not (len(self.aliStringSymbols) == len(self.hmmStringLetters) == len(self.aliStringSymbols) == len(self.aliStringLetters) == len(self.matchString)):
        #    raise ValueError('uneven alignemnt strings')

    def __repr__(self):
        return str( self.__dict__ )

## Low -level alignment buffer parser
def plumbParse(buffer, index, mainStringAccumulator, stuffContainer, reSymbolStuff, goingUp=True):
    while(buffer[index]):
        line = buffer[index]
        #print (index, '::',line)

        if line == '':
            break;
        m = reStrandsLine.match(line)
        mAlt = re.match(reSymbolStuff, line)
        if m:
        #    print ('ACC', m.groups())
        #    print(mainStringAccumulator)
            mainStringAccumulator += m.groups()[3] # Check index
        #    print('Strand View internal here::', mainStringAccumulator)
        elif mAlt:
            stuffKey = mAlt.groups()[1]
            if not stuffKey in stuffContainer:
                stuffContainer[stuffKey] = ''
            stuffContainer[stuffKey] += mAlt.groups()[0]
        else:
            raise ValueError("Unknown alignment record", line)

        index += -1 if goingUp else 1
       # print ("-->", index)
    return mainStringAccumulator, stuffContainer

class HMMObj():
    """
    Object that stores information about one protein domain hit

    :param prot: 
    :param domain: 
    :param hit:
    """

    def __init__(self,prot,domain,hit):
        self.prot=prot
        self.domain=domain
        self.hit=hit
        self.overlapped_hits=[]

    @property    
    def sequence(self):
        """
        Domain sequence in protein, without indels.
        """
        return self.hit.aliStringLetters.replace("-","").upper()

    @property
    def start(self):
        """
        Domain start position in protein
        """
        return int(self.hit.aliFrom)

    @property
    def end(self):
        """
        Domain end position in protein
        """
        return int(self.hit.aliTo)    

    def is_overlapping(self,other_hit,accept_overlap_size):
        """
        Check if two HMMObj overlaps.

        :param other_hit: 
        :param accept_overlap_size: 

        :return: True if overlapping size is bigger than accept_overlap_size. False elsewhere. 
        """
        start1=self.start
        end1=self.end
        start2=other_hit.start
        end2=other_hit.end
        residues1=set([i for i in range(start1,end1+1)])
        residues2=set([i for i in range(start2,end2+1)])
        if len(residues1.intersection(residues2))>accept_overlap_size:
            return True
        return False

    def reinitialize_overlapped_hits(self):
        self.overlapped_hits=[]