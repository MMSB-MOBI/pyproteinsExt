import Uniprot
import os


class noIndexError(Exception):
    def __init__(self, id):
        self.id = id

class noIdError(Exception):
    pass
#  Guess using uniprot identifier

class TMHMM_ele():
    def __init__(self, name=None):
        self.name = name
        self.data = []
    def push(self, datum):
        self.data.append(datum)
    def __repr__(self):
        states = self.getStatusList()
        return self.name + "\n" + ('').join(states)
    def getStatusList(self):
        states = [d['status'] for d in self.data]
        return states


class Parser():
    def __init__(self, dirLocation=None):
        ## indexing file elements
        self.index = {}
        self.dirLocation = []
        self.dirLocation.append(dirLocation)
        for dirPath in self.dirLocation:
            self._index(dirPath)

    def _index(self, dirPath):
        for file in os.listdir(dirPath):
            if not file.endswith(".plp"):
                continue

            id = Uniprot.strip(file)
            if id :
                self.index[id] = dirPath + "/" + file

    def get(self, string = None, header=None, aminoAcidSeq=None):

        id = Uniprot.strip(string)
        try :
            if not id:
                raise noIdError()
            if not self.index[id]:
                raise noIndexError(id)

        except noIdError:
            print "No protein identifier found to search tmhmm index from string " + string
            return None
        except noIndexError as e:
            print "No file found in index for protein identifier " + e.id + " striped from string " + string
            return None

        stateList = self.parse(self.index[id], id)
        return stateList

    def parse(self, file=None, name=None):

        # open file
        with open (file) as inp:
            element = TMHMM_ele(name)
            for line in inp:
                if not line.startswith("#"):
                    buffer = line.split()
                    scores = [float(i) for i in buffer[2:5]]
                    m = max(scores)
                    i = scores.index(m)
                    status = 'I'
                    if i == 1 :
                        status = 'M'
                    elif i == 2:
                        status = 'E'

                    element.push({'aa' : buffer[1], 'position' : buffer[0], 'inside' : buffer[2], 'membrane' : buffer[3], 'outside' : buffer[4], 'status' : status})
        return element
