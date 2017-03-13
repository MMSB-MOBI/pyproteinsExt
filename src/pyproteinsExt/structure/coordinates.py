import re


'''
GL 09/03/2017
 This is minimu implementation of PDB parsing and coordinates mangament

'''

''' PDB format reminder

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

'''

class Parser(object):
    def __init__(self):
        pass

    def load(self, **kwargs):
        structureObj = Structure();

        if 'file' in kwargs:
            with open(kwargs['file'], 'r') as f:
                currentAtomList = structureObj.pop();
                for l in f:
                    if l.startswith("ATOM "):
                        currentAtomList.append(Atom(string=l))
                    elif l.startswith("ENDMDL"):
                        currentAtomList = structureObj.pop();

        return structureObj



class Model(object):
    def init(self):
        pass

class Structure(object):
    def __init__(self):
        self.model = []
        self.currModel = 1
        self.bufferSelection = []

    @property
    def residueNumber(self):
        return len([ res for res in self.byres() ])

    def byres(self):

        i = 0
        data = self.model[self.currModel - 1]
        w_resSeq = data[0].resSeq
        w_chainID = data[0].chainID
        i_start = 0
        i_stop = 0
        while i < len(data):
            #print str(i) + " -- " + str(len(data))
            cur_atom = data[i]
            #print str(w_resSeq) + ' vs ' + str(cur_atom.resSeq) + ' and ' + w_chainID + ' vs ' + cur_atom.chainID
            if w_resSeq != cur_atom.resSeq or w_chainID != cur_atom.chainID:
                #i_stop = i + 1
                #print "pulling at " + str(i_start) + ":" + str(i)
                yield Residue(data[i_start:i])

                w_resSeq  = cur_atom.resSeq
                w_chainID = cur_atom.chainID

                i_start = i
            i += 1
        yield Residue(data[i_start:])


    def pop(self): # Append an empty array to model list and return its referenc

        self.model.append([])
        return self.model[-1]

    def select(exp): # TO DO
        pass
        #
        # Should implment binary tree construction from the parsing of pymol-like expression
        # All tree relationship could then be mapped to set operation on atoms lists
        #
        #        http://pyparsing.wikispaces.com/
'''
        buf = []
        if 'chainID' in kwargs:
            for atom in self.bufferSelection:#self.model[currModel - 1]:
                if atom.chainID == kwargs['chainID']:
                    buf.append(atom)
        elif 'resName' in kwargs:
            p = re.compile('my pattern')
            for atom in self.bufferSelection:#self.model[currModel - 1]:
                if p.match(atom.resName):
                    buf.append(atom)


        self.bufferSelection = buf
        return self
'''

    #def push(self, **kwargs):
    #    if 'atomRecord' in kwargs:


class Residue(object):
    def __init__(self, atomArray):
        self.data = atomArray

    def __len__(self):
        return len(self.data)
    def __str__(self):
        return ''.join([ str(a) for a in self.data ])

    def __getitem__(self, key):
        return self.data[key]

    def __iter__(self):
        for d in self.data:
            yield d
        #        return self.data

    @property
    def id(self):
        return self.data[0].resName + str(self.data[0].resSeq) + ':' + self.data[0].chainID

class Atom(object):
    def __init__(self, **kwargs):
        if 'string' in kwargs:
            buf=kwargs['string']
            bufLen = len(buf)
            if bufLen < 55:
                raise TypeError("Unexpected atom record length \"" + str(buLen) + "\" on input at \n" + buf)

            self.recordName=buf[0:6].replace(" ", "")
            self.serial=int(buf[6:11])
            self.name=buf[12:16].replace(" ", "")
            self.altLoc= buf[16]
            self.resName=buf[17:20].replace(" ", "")
            self.chainID=buf[21]
            self.resSeq=int(buf[22:26]) # Not sure
            self.iCode=buf[26]
            self.x=float(buf[30:38])
            self.y=float(buf[38:46])
            self.z=float(buf[46:54])
            self.occupancy  = float(buf[54:60])           if bufLen > 54 else None
            self.tempFactor = float(buf[60:66])           if bufLen > 60 else None
            self.element    = buf[76:78].replace(" ", "") if bufLen > 76 else None
            self.charge     = buf[78:80].replace(" ", "") if bufLen > 78 else None


    def __hash__(self):
        tup = (self.recordName, self.serial, self.altLoc, self.resName, self.chainID, self.resSeq, self.iCode, self.x, self.y, self.z, self.occupancy, self.tempFactor, self.element, self.charge)
        return hash(tup)

    def __str__(self):
        if self.charge:
            return "%-5s%5d %4s%s%4s %s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s" % (self.recordName, self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq,self.iCode, self.x, self.y, self.z, self.occupancy, self.tempFactor, self.element, self.charge)
        if self.element:
            return "%-5s%5d %4s%s%4s %s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  " % (self.recordName, self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq,self.iCode, self.x, self.y, self.z, self.occupancy, self.tempFactor, self.element)
        if self.tempFactor:
            return "%-5s%5d %4s%s%4s %s%4s%s   %8.3f%8.3f%8.3f%6.2f%6.2f              " % (self.recordName, self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq,self.iCode, self.x, self.y, self.z, self.occupancy, self.tempFactor)
        if self.occupancy:
            return "%-5s%5d %4s%s%4s %s%4s%s   %8.3f%8.3f%8.3f%6.2f                 " % (self.recordName, self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq,self.iCode, self.x, self.y, self.z, self.occupancy)
        return "%-5s%5d %4s%s%4s %s%4s%s   %8.3f%8.3f%8.3f                          " % (self.recordName, self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq,self.iCode, self.x, self.y, self.z)



