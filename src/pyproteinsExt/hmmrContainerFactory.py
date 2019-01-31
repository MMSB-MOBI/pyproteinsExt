import re,sys, io
import gzip

#
# Factory to parse hmmscan // hmmsearch outputs
#
#


#
# Query - Atomic
# 

reInstance=re.compile('(\# hmm([\S]+) :: search ([\S]+) against a ([\S]+) database.+?(?=\[ok\])\[ok\])', re.DOTALL)
reType=re.compile('^\# hmm([\S]+) :: search ([\S]+) against a ([\S]+) database[\s]*$')
# hmmsearch :: search profile(s) against a sequence database
reQuery   = re.compile('Query:[\s]+([\S].*)[\s]+\[[LM]=[\d]+\]')
reStart   = re.compile('^Scores for complete sequences \(score includes all domains\):')
reSum     = re.compile('(([\d\.e-]+[\s]+){7})([0-9]+)[\s]+([\S]+)[\s]+([\S].*)$')
reInclude = re.compile('------ inclusion threshold ------')
reDetail  = re.compile('Domain annotation for each ([\S]+) \(and alignments\):')
reStop    = re.compile('Internal pipeline statistics summary')
reAliLongName  = re.compile('^>>[\s]*([\S]+.*)$')
reDomainDetail = re.compile('  == domain [\d]+(?!.* == domain [\d])')
reAlignDef = re.compile('^[\s]*[\d]+[\s]*[\S][\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]]{2}[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]]{2}[\s]+([\S]+)[\s]+([\S]+)[\s]+[\.\[\]{2}[\s]+([\.\d]+)')
reBothStringLetter = re.compile('^([\s]+([\S]+)[\s]+([\S]+)[\s]+)([\S]+)[\s]+([\S]+)')
reBothStringSymbol = re.compile('^([\s]+)([\S]+)')
reEmpty = re.compile('^[\s]*$')
reStrandsLine = re.compile("^([\s]*([\S]+)[\s]+([\d]+)[\s])([\S]+)[\s]([\d]+)[\s]*$")

#reInstance = re.compile('(# hmmsearch.+?(?=\[ok\])\[ok\])', re.DOTALL)

def parse(inputFile=None):
    upType = None
    bigBuffer = ''
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
            if not runReports:
                raise ValueError("Could not parse upper level in provided hmm(scan/search) file")
            

            mainContainer = Container()
            
            for d in runReports:           
                mainContainer += Container(input=io.StringIO(d[0]), upType=d[1])
            return mainContainer

class Container(object):
    def __init__(self, input=None, upType=None):
       # self.queryHmmFile=None             #../../../data/fad_binding.hmm
       # self.targetSequenceDatabase=None   #./Trembl_50.fasta
       # self.query=None                    #PF08022_full  [M=104]
        self.upType = upType
        self._transpose = None
        if not input:
            self.summary = []
            self.details = []
        else:
            self.summary, self.details = _parseBuffer(input)

    def __add__(self, other):
        if self.upType != other.upType:
            if self.upType == None:
                self.upType = other.upType
            elif other.upType == None:
                other.upType = self.upType
        if self.upType != other.upType:
            raise TypeError('Adding', self.upType,' and ', other.upType)

        self.summary.extend(other.summary)
        self.details.extend(other.details)
        return self

    def __len__(self):
        return len(self.details)

# returns a list of proteins referencing HMM hits
    def T(self):
        if not self._transpose:
            t = {}
            for d in self.details:
                if not d.aliShortID:
                    continue
                if d.aliShortID not in t:
                    t[d.aliShortID] = {}
                if d.hmmID in t[d.aliShortID]:
                    print('Warning(' + d.aliShortID + '), multiple hit w/ identical HMM query on a single protein')
                else :
                    t[d.aliShortID][d.hmmID] = []
        
                t[d.aliShortID][d.hmmID].append(d)    
            self._transpose = t
        
        return self._transpose
        


def _parseBuffer(input):
    summary = []
    details = []
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
            summaryBool = True
            inclusionBool = True
            continue

        if summaryBool and reEmpty.match(l):
            #print ("LOW")
            summaryBool = False
            continue

        #    print("?")
        #    if l.starstwith('Domain annotation for each model'):         
        #        print("NN")  
        #        summaryBool = False
        #        continue

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


        #if readBool:
            #buf = l.split("    ")
            #if len(buf) == 10:
            #print(str(len(buf))+ ':: ' + str(buf))

        if summaryBool:
            if reInclude.search(l):
                inclusionBool = False
            #print("Summary parsing attempt")
            m = reSum.search(l)        
            if m:
                scores = m.group(1).split()
                summary.append({
                    'fullSequence' : {
                        'E-value'    : scores[0],
                        'score'      : scores[1],
                        'bias'       : scores[2]
                    },
                    'bestOneDomain' : {
                        'E-value'    : scores[3],
                        'score'      : scores[4],
                        'bias'       : scores[5]
                    },
                    'exp'        : scores[6],
                    'N'          : m.group(3),
                    'Model'   : m.group(4),
                    'Description': m.group(5),
                    'included'   : inclusionBool
                })
                    #print (summary[-1])
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
    details = [ Match(rawData, queryID) for rawData in detailBuffer ]
    return (summary, details)

class Match (object):
    def __init__(self, bufferString, queryID):
        self.aliLongID = None
        self.aliShortID = None
        self.hmmID = None
        self.queryID = queryID
        self.parseDetailEntry(bufferString)
        
    def __hash__(self):
        return hash(self.aliLongID + self.hmmID)

    def parseDetailEntry(self, bufferString):
        buff = bufferString.split('  == domain ')
        self.data = []
        if len(buff) == 1:
            print ("Warning:: " + buff[0])
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