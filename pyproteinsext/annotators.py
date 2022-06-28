import pyproteins.services.utils

#
# Returns
# ("Fields", .... )
#
#
#


# Canonical annotator data structure

#ECMannotator = {
#    'tag' : "ECM",  #Name of the annotator
#        'positiveAnnotationList' : [ #set of properties that selects an entry
#        { 'name' : 'goWords', # Name of the column in csv
#          'target' : 'GO', # properties/Attribute to get from the tested object
#           'content' : [  # Values to test against got attributes
#            { 'id' : 'GO:0016020', 'txt' : 'membrane'}, # a test value, only id is effectively tested
#        ]
#      },
#        { 'name' : 'keyWords',
#          'target' : 'KW',
#          'content' : [
#                { 'id' : 'KW-0472', 'txt' : 'Membrane'},
#            ],
#        }
#    ],
#    'negativeAnnotationList' : [
#        { 'name' : 'badGoWords',
#          'target' : 'GO',
#          'content' : [
#            { 'id' : 'GO:0031965', 'txt' : 'nuclear membrane' }
#            ]
#        }
#    ]
#}

# Annotator type/sub properties checks
def _checkConstraints(datum):
    if 'tag' not in datum:
        raise ValueError("No \'tag\' values in annotator dictionary input")
    if not 'positiveAnnotationList' in datum or not 'negativeAnnotationList' in datum:
        raise ValueError("No \'positiveAnnotationList\' or \'negativeAnnotationList\' values in annotator dictionary input")

def _checkConstraint(data):
    for d in data:
        if not 'name' in data:
            raise ValueError("No \'name\' value in annotation element")
        if not 'target' in data:
            raise ValueError("No \'target\' value in annotation element")
        if not 'content' in data:
            raise ValueError("No \'content\' value in annotation element")
    for d in data['content']:
        if not 'id' in d:
            raise ValueError("No \'id\' value in annotation element content")



#
#
def stringifyContraintList(data):
    _str = ''
    for terms in data:
        #print '-->' + str(terms) + '<--'
        _str += '\tname : ' + terms['name'] + '\n\ttarget : ' + terms['target'] + '\n'
        _str += '\tcontent:\n' + '\n'.join(['\t\t' + x.id + ', ' + x.txt for x in terms['content'] ])
        _str += '\n'
    return _str

class AnnotationTerm(object):
    def __init__(self, datum):
        if 'id' not in datum:
            raise ValueError("missing id in annotation term specs")
        self.id = datum['id']
        self.txt = datum['txt'] if 'txt' in datum else 'N/A'

    def __str__(self):
        #return term['id'] + '(' + term['txt'].replace('(', '\(').replace(')', '\)') + ')'
        return self.id + '(' + self.txt + ')'
    def __repr__(self):
        return str(self)

class Uniprot(object):
    def __repr__(self):
        _str = 'tag : ' + self.tag + '\n'
        if 'positiveAnnotationList' in self.constraints:
            if self.constraints['positiveAnnotationList']:
                _str += 'positive annotation terms:\n'
                _str += stringifyContraintList(self.constraints['positiveAnnotationList'])
        #print _str
        if 'negativeAnnotationList' in self.constraints:
            if self.constraints['negativeAnnotationList']:
                _str += 'negative annotation terms:\n'
                _str += stringifyContraintList(self.constraints['negativeAnnotationList'])

        return _str

    def __init__(self, tag="default"):
        self.tag = tag
        self.constraints = {}

    def loadConstraints(self, datum):
        _checkConstraints(datum)
        if not self.constraints:
            self.constraints = datum

    def addPositive(self, datum):
        self._add(datum, 'positiveAnnotationList')
    def addNegative(self, datum):
        self._add(datum, 'negativeAnnotationList')

    def _add(self, datum, aType):
        if not isinstance(datum, list):
            datum = [datum]
        toAddList = []
        for d in datum:
            _checkConstraint(d)
            toAdd = {'name' : d['name'],
                     'target' : d['target'],
                     'content' :[]}
            toAdd['content'] = [ AnnotationTerm(x) for x in d['content']]
            toAddList.append(toAdd)

        if not aType in self.constraints:
            self.constraints[aType] = []
        self.constraints[aType] += toAddList

    #returns :
    #  termsNameX  : [[value, ..], ...]    #N element annotated
    #  termsNameY  : [[value, ..], ...]    #N element annotated
    #  status : [True, False ... ]         #N element status
    def annotateAll(self, elementList):
        tType = ['positiveAnnotationList', 'negativeAnnotationList']
        data = {}
        status = []
      
        for e in elementList:
           
            mBool, datum = self.annotate(e)
          
            status.append(mBool)
            for aType in tType:
                if not datum[aType]:
                    break
                for d in datum[aType]:
                    if d['name'] not in data:
                        data[d['name']] = []
                    data[d['name']].append(d['matches'])
     
        return (data, status)
            #annotationHits[annotationType]
            #d = { 'name' : annotation['name'], 'matches' : [] }

#   return the annotations found in the entry
    def annotate(self, elem):
        e = None
        if hasattr(elem, 'GO'):
            e = elem
        elif hasattr(elem, '_uniprotBound'):
            e = elem._uniprotBound
        else:
            raise ValueError('Cant find valid uniprot Object, cant annotate')

        if not self.constraints:
            raise ValueError( 'constraints required to perform annotation test')
        (matchBool, matchContent) = testEntry(e, self.constraints)
        return (matchBool, matchContent)
# Returns a boolean if the uniprot entry matches  specified annotation rules
    def isValid(self, e):
        # unroll dictWords, goWords keyWords badGoWords
        if not self.constraints:
            raise ValueError( 'constraints required to perform annotation test')

        (matchBool, matchContent) = testEntry(e._uniprotBound, self.constraints)
        return matchBool

    def pandify(self):
        pass
        #for annotationType in ['positiveAnnotationList', 'negativeAnnotationList']:
        #    for termsList in matchContent[annotationType]:
        #        tag = annotator['tag'] + '_' + termsList['name']
        #        self[tag] = 'NA' if not termsList['matches'] else (';').join(['(ID:' + d['id'] + ')' + d['txt'] for d in termsList['matches']])


def testEntry(uniProtObj, annotator):
    annotationHits = {
        "positiveAnnotationList" : [],
        "negativeAnnotationList" : []
        }
    blackListedTermFound = False
    wishedListedTermFound = False
    for annotationType in ['positiveAnnotationList', 'negativeAnnotationList']:
        if annotationType in annotator:
            for annotation in annotator[annotationType]:
                _testBoundMethod = None
                if annotation['target'] == 'GO':
                    _testBoundMethod = uniProtObj.hasGO
                elif annotation['target'] == 'KW':
                    _testBoundMethod = uniProtObj.hasKW
                elif annotation['target'] == 'MIM':
                    _testBoundMethod = uniProtObj.hasMIM
                elif annotation['target'] == 'DI':
                    _testBoundMethod = uniProtObj.hasDI
                elif annotation['target'] == 'ORPHA':
                    _testBoundMethod = uniProtObj.hasORPHA
                else:
                    raise ValueError("Unknown annotation to test on uniprot Entry " + str(annotation['target']))

                d = { 'name' : annotation['name'], 'matches' : [] }
                for term in annotation['content']:
                    if _testBoundMethod(term.id):
                        d['matches'].append(term)
                        if annotationType == 'positiveAnnotationList':
                            wishedListedTermFound = True
                        else:
                            blackListedTermFound = True
                annotationHits[annotationType].append(d)

            entryBool = False if blackListedTermFound else wishedListedTermFound
    return (entryBool, annotationHits)


class Matrisome:
    def __init__(self, masterFile=None):
        if not masterFile:
            raise ValueError ("You must provide a masterFile argument to Matrisome constructor")

        tsvBuffer = pyproteins.services.utils.tsvToDictList(fileName=masterFile)
        self.data = tsvBuffer['data']
        self.keymap = tsvBuffer['keymap']
        self.accessors = self._index()

    # index entry using relevant key found in their record
    def _index(self):
        accessors = {'uniprot' : {}, 'gene' : {}} # any other meaningfull accessor
        for d in self.data:
            #print "-->" + d['UniProt_IDs'] + '<--'
            if d['UniProt_IDs']:
                for k in d['UniProt_IDs'].split(':'):
                    if k in accessors['uniprot']:
                        print (k + ' found multiple time')
                    else :
                        accessors['uniprot'][k] = []
                    accessors['uniprot'][k].append(d)
            if d['Gene Symbol']:
                for k in d['Gene Symbol'].split(':'):
                    if k in accessors['gene']:
                        print (k + ' found multiple time')
                    else :
                        accessors['gene'][k] = []
                    accessors['gene'][k].append(d)
        return accessors

    def get(self, uniprotID=None):
        if uniprotID:
            if not uniprotID in self.accessors['uniprot']:
                return []
            return [{'Category' : d['Category'], 'Division' : d['Division']} for d in self.accessors['uniprot'][uniprotID]]
        return []