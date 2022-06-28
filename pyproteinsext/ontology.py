# This package provides utility function to manipulate the GO/MI ontology terms
# check performed on passed string is a bit clunky
# "obo:MI_0090" or "MI:0090" are OK
import re
from owlready2 import *


def isOboRegular(term):
    # print "Coucou " + term
    # print re.match(r'^obo:[A-Z]{2}:[0-9]{4}$', term)
    return re.match(r'^obo:[A-Z]{2}_[0-9]{4}$', term)


def isOboNamespaced(term):
    return re.match(r'[^:]+:[0-9]+', term)


# Manipulating the MolecularInteraciton ontology
# print model.classes
# aClass = model.getClass("MI_0090")
# aClass[0].children()
# model.printClassTree()


class Ontology():
    def __init__(self, owlFile=None, url=None, oType='GO'):
        self.type =oType
        # super(Ontology, self).__init__(ressource)
        if not owlFile is None :
            self.onto = get_ontology("file://" + owlFile).load()
        elif not url is None:
            self.onto = get_ontology(url).load()
        else:
            raise IOError("No ontology source specifed!")
    
    # find terms based on namespaced Id or text description in label
    def find(self, stringTerm):
        if isOboNamespaced(stringTerm):
            elem = self.onto.search(id=stringTerm)            
            if not elem:               
                elem = self.onto.search(hasAlternativeId=stringTerm)
            if not elem:
                raise KeyError(f"No element \"{stringTerm}\"")
            return elem
        return self.onto.search(label=stringTerm)

    def findOne(self, stringTerm):
    #    if isOboNamespaced(stringTerm):
    #        return self.onto.search_one(id=stringTerm)
    #    return self.onto.search_one(label=stringTerm)
        if isOboNamespaced(stringTerm):
            elem = self.onto.search_one(id=stringTerm)
            if not elem:
                #print("T2")
                elem = self.onto.search_one(hasAlternativeId=stringTerm)
            if not elem:
                #print("T3")
                raise KeyError(f"No element \"{stringTerm}\"")
            return  elem
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
        print('Total number of nodes matching request : ' + str(len(seedNodes)))

        secondaryNodes = []
        for seedNode in seedNodes:
            secondaryNodes += self.onto.search(subclass_of=seedNode)
        secondaryNodes = list(set(secondaryNodes))
        print('Total number of potential nodes found under seed nodes : ' + str(len(secondaryNodes)))

        return list(set(seedNodes + secondaryNodes))

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

    @staticmethod
    def _rollupNode(e):
        if not e.is_a:
            return None
        return e.is_a[0]

    # Given a list of terms regroup them in clusters based on the member list
    # that are the highest in the hierarchy
    # INCOMPLETE
    def cluster(self, termsList):

        def _hasAnyDad(e, l):
            for _e in l:
                # print(str(_e) + ' ' + str(e))
                if self.isSonOf(e, _e):
                    return True
            return False

        # Create a list of term that have no parents

        rootNodes = [{e: []} for e in termsList if not _hasAnyDad(e, termsList)]

        return rootNodes

    # Find and Count foreach provided term in
    # the correesponding parents in range domain the number
    # Not tested w/ owlReady !!!
    def project(self, domainTerms, rangeTerms, flatDic=False):

        domainTerms = self._coherce(domainTerms)
        res = {x: [] for x in self._coherce(rangeTerms)}
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

        data = {"term": [], "parent": []}
        for rangeTerm in res:
            for son in res[rangeTerm]:
                data["term"].append(str(son.bestLabel()))
                data["parent"].append(str(rangeTerm.bestLabel()))
        for query in others:
            data["term"].append(str(query.bestLabel()))
            data["parent"].append("others")

        return data

    # if term is a string coherce it in a ontology class object
    def _coherce(self, data):
        if not isinstance(data, list):
            data = [data]
        return [self._termCoherce(e) for e in data]

    def _termCoherce(self, termLike):
        if isinstance(termLike, str):
            term = self.findOne(termLike)
            if not term:
                raise ValueError("No such term found : \"" + str(termLike) + "\"")
            termLike = term

        return termLike

    def constructTreeFromRoot(self, root_id: str, things_to_add: dict = None):
        """
        Construct a Tree instance from root_id node (id=root_id).
        Returns a Tree object.

        things_to_add may be an iterable of string specified which attribute of the owl node should be imported with.
        Default imported attributes are id and label, stored in node.id and node.label.
        Every other key you need to add should be specified in things_to_add. An inexistant property in a node will be ignored.
        Data will be available in node.misc, in a dictionnary.

        things_to_add may be a function taking a OwlNode in parameter (also known as ThingClass) and returning None OR a list of string.
        List of string will be treated for this node like specified in the previous paragraph.

        Examples:
            # Classic

            my_tree = ontology.constructTreeFromRoot("MI:0001")

            # Let's adding some additionnals informations to the nodes

            # With iterables

            my_tree = ontology.constructTreeFromRoot("MI:0001", ['hasExactSynonym'])

            print(my_tree.root.misc['hasExactSynonym'])  # prints "interaction detect"

            # With functions

            my_tree = ontology.constructTreeFromRoot("MI:0001", lambda x: ['hasExactSynonym'] if len(x.label) > 1 else None)
        """

        def find(predicate, iterable):
            for e in iterable:
                if predicate(e):
                    return e
            return None

        def generateUntilSeed(predicate, iterable):
            for e in iterable:
                yield e
                if predicate(e):
                    break

        def getThingsToAdd(node):
            things = things_to_add

            if callable(things_to_add):
                things = things_to_add(node)

            if things:
                d = {}

                for k in things:
                    try:
                        d[k] = getattr(node, k)
                    except AttributeError:
                        # Does not exists
                        pass

                if len(d):
                    return d

            return None

        nodes = self.harvest(root_id)

        tree = Tree()

        rootnode = self.onto.search_one(id=root_id)

        if not rootnode:
            raise IndexError("Seed does not exists")

        tree.root = Node(rootnode.id[0], rootnode.label[0], misc=getThingsToAdd(rootnode))

        for c in nodes:
            # print(c, c.id, c.label)
            try:
                c_id, c_label, c_misc = c.id[0], c.label[0], getThingsToAdd(c)

                lineage = map(
                    lambda x: (x.id[0], x.label[0], getThingsToAdd(x)) if x.name != 'Thing' else ("__root__", None),
                    self._getLineage(c_id))

                # Filtre lineage jusqu'à trouver le root. Arrête ensuite
                lineage = list(generateUntilSeed(lambda x: x[0] == root_id, lineage))

                if not find(lambda x: x[0] == root_id, lineage):
                    continue

                tree.append(lineage, c_id, c_label, c_misc)
            except IndexError:
                pass

        return tree


# Needed for create tree for Ontology.constructTreeFromRoot()
class Tree:
    """Tree with unique nodes, using storage containing only IDs and labels."""

    def __init__(self):
        self.root = None
        self.node_dict = {}

    def toDict(self):
        """
        Get dictionnary version of a tree.

        Supposing an interface like:

        interface NodeChild { name: str /* Node label */, children: { [childId: str]: NodeChild }, misc?: any }


        Dictionnary will be organised like:


        {
            root_node_id: NodeChild
        }
        """
        return {
            self.root.id: self.root.toDict()
        }

    def append(self, parents, current_id, current_label, misc=None):
        """Append a node to the tree.

        Parents must be a list of tuples containing in [0]: id, [1]: label. Ordered by child proximity, the first one is the nearest.
        It can be a "misc" dict in [2] position.
        """
        # Si le noeud existe déjà
        if self.findInTree(current_id):
            return

        i = 0
        for p in parents:
            parent_id = p[0]
            parent_l = p[1]
            parent_misc = p[2] if len(p) > 2 else None

            if self.findInTree(parent_id):
                break

            self.append(parents[i + 1:], parent_id, parent_l, parent_misc)

            i += 1

        new_node = Node(current_id, current_label, misc=misc)

        if len(parents) == 0:
            self.root.parent = new_node
            new_node.addChild(self.root)
            self.root = new_node
        else:
            parent = self.findInTree(parents[0][0])
            # print("parent", parents, parents[0][0], parent, self.root.children)
            new_node.parent = parent
            parent.addChild(new_node)

        self.node_dict[new_node.id] = new_node

    def findInTree(self, id: str):
        if self.root is None:
            return None

        if id in self.node_dict:
            return self.node_dict[id]

        return self.root if self.root.id == id else self.root.findInNode(id)

    def clone(self, deep=True):
        new_tree = Tree()
        new_tree.root = self.root.clone(deep=deep)

        return new_tree

    def prune(self, seeds: list):
        new_tree = self.clone()

        # Prune new tree with specific seeds
        okay = new_tree.root.prune(set(seeds))

        if not okay:
            raise IndexError("Cant find any node matching given seeds")

        return new_tree

    def toEte3(self):
        return self.root.toEte3()


class Node:
    def __init__(self, id, label, parent=None, misc=None):
        self.parent = parent
        self.id = id
        self.label = label
        self.children = {}
        self.misc = misc

    def addChild(self, child):
        if not isinstance(child, Node):
            raise TypeError('Must be Node type')

        child.parent = self
        self.children[child.id] = child

    def removeChild(self, id):
        del self.children[id]

    def childExists(self, id):
        return id in self.children

    def childrens(self):
        return self.children.values()

    def isLeaf(self):
        return self.children.__len__() == 0

    def toDict(self):
        d = {
            'name': self.label,
            'children': {x.id: x.toDict() for x in self.childrens()}
        }

        if self.misc:
            d['misc'] = self.misc

        return d

    def findInNode(self, id: str):
        if id in self.children:
            return self.children[id]

        for c in self.childrens():
            n = c.findInNode(id)

            if n:
                return n

        return None

    # Need to import ete3 and give TreeNode to this func
    def toEte3(self, TreeNode):
        node = TreeNode(name=self.id)
        node.add_feature("label", self.label)

        for c in self.childrens():
            node.add_child(c.toEte3())

        return node

    def clone(self, deep=True):
        # misc will be swallow copied
        cur = Node(self.id, self.label, misc=self.misc)

        if deep:
            for c in self.childrens():
                child_clone = c.clone(deep=deep)
                cur.addChild(child_clone)

        return cur

    # seeds is an iterable, ideally a set
    def prune(self, seeds):
        one_node_is_okay = False

        children_to_remove = []

        for e in self.children:
            okay = self.children[e].prune(seeds)

            if not okay:
                children_to_remove.append(e)
            else:
                one_node_is_okay = True

        for c in children_to_remove:
            self.removeChild(c)

        if not one_node_is_okay:
            one_node_is_okay = self.id in seeds

        return one_node_is_okay
