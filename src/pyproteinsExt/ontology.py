# This package provides utility function to manipulate the GO/MI ontology terms
# check performed on passed string is a bit clunky
# "obo:MI_0090" or "MI:0090" are OK
import re
from owlready2 import *

def isOboRegular(term):
   # print "Coucou " + term
  #  print re.match(r'^obo:[A-Z]{2}:[0-9]{4}$', term)
    return re.match(r'^obo:[A-Z]{2}_[0-9]{4}$', term)#

def isOboNamespaced(term):
    return re.match(r'[^:]+:[0-9]+', term)

# Manipulating the MolecularInteraciton ontology
#print model.classes
#aClass = model.getClass("MI_0090")
#aClass[0].children()
#model.printClassTree()


class Ontology():
    def __init__(self, file=None):
        #super(Ontology, self).__init__(ressource)
        if file:
            self.onto = get_ontology("file://" + file).load()

# find terms based on namespaced Id or text description in label
    def find(self, stringTerm):
        if isOboNamespaced(stringTerm):
            return self.onto.search(id=stringTerm)
        return self.onto.search(label=stringTerm)

    def findOne(self, stringTerm):
        if isOboNamespaced(stringTerm):
            return self.onto.search_one(id=stringTerm)
        return self.onto.search_one(label=stringTerm)

 # if term is a/list of string coherce it into matching ontology classes
    def _coherceIntoMany(self, data):
        if not isinstance(data, list):
            data = [data]
        res = []
        for e in data:
            res += self.find(e)
        return list(set(res))

# Given a list of term of termlike, (including regexp)
# Get the non redundant list of terms that match the supplied argument and all their sons

    def harvest(self, termLikeList):

        seedNodes = self._coherceIntoMany(termLikeList)
        print ('Total number of nodes matching request : ' + str(len(seedNodes) ) )

        secondaryNodes = []
        for seedNode in seedNodes:
            secondaryNodes += self.onto.search(subclass_of=seedNode)
        secondaryNodes = list(set(secondaryNodes))
        print ('Total number of potential nodes found under seed nodes : ' + str(len(secondaryNodes) ) )

        return list(set( seedNodes + secondaryNodes ) )

    def isSonOf(self, x, y):
        x, y = self._coherce([x, y])
        return y in self._getLineage(x)


    def _getLineage(self, term):
        term = self._coherce(term)[0]

        lineage = []

        tmpNode = self._rollupNode(term)
        if not tmpNode:
            print(str(term) + ' has no parent')

        while tmpNode:
            lineage.append(tmpNode)
            tmpNode = self._rollupNode(tmpNode)

        return lineage

    def _rollupNode(self, e):
        if not e.is_a:
            return None
        return e.is_a[0]

    # Given a list of terms regroup them in clusters based on the member list 
    # that are the highest in the hierarchy
    # INCOMPLETE
    def cluster(self, termsList, lvl=1):

        def _hasAnyDad(e, l):
            for _e in l:
                #print(str(_e) + ' ' + str(e))
                if self.isSonOf(e, _e):
                    return True
            return False
        # Create a list of term that have no parents

        rootNodes = [ {e : []} for e in termsList if not _hasAnyDad(e, termsList) ]
    
        return rootNodes



# Find and Count foreach provided term in
# the correesponding parents in range domain the number
# Not tested w/ owlReady !!!
    def project(self, domainTerms, rangeTerms, flatDic=False):

        domainTerms = self._coherce(domainTerms)
        res = { x : [] for x in self._coherce(rangeTerms) }
        others = []
        for query in domainTerms:
            known = False
            for rangeTerm in res:
                if self.isSonOf(query, rangeTerm):
                    res[rangeTerm].append(query)
                    known = True
            if not known:
                others.append(query)

        if not flatDic:
            res["others"] = others
            return res

        data = { "term" : [], "parent" : [] }
        for rangeTerm in res:
            for son in res[rangeTerm]:
                data["term"].append( str(son.bestLabel()) )
                data["parent"].append( str(rangeTerm.bestLabel()) )
        for query in others:
            data["term"].append( str(query.bestLabel()) )
            data["parent"].append( "others" )

        return data

    # if term is a string coherce it in a ontology class object
    def _coherce(self, data):
        if not isinstance(data, list):
            data = [data]
        return [ self._termCoherce(e) for e in data ]

    def _termCoherce(self, termLike):
        if isinstance(termLike, str):
            term = self.findOne(termLike)
            if not term:
                raise ValueError("No such term found : \"" + str(termLike) + "\"")
            termLike = term

        return termLike

