import re
import copy
import numpy
import math
import os

aaCode = {
	'ALA' : 'A',
	'ARG' : 'R',
	'ASN' : 'N',
	'ASP' : 'D',
	'GLU' : 'E',
	'GLN' : 'Q',
	'GLY' : 'G',
	'HIS' : 'H',
	'CYS' : 'C',
	'LYS' : 'K',
	'MET' : 'M',
	'LEU' : 'L',
	'ILE' : 'I',
	'MET' : 'M',
	'PRO' : 'P',
	'SER' : 'S',
	'THR' : 'T',
	'TRP' : 'W',
	'TRY' : 'W',
	'TYR' : 'Y',
	'PHE' : 'F',
	'VAL' : 'V',
	'MSE' : 'M'
}

def translate(x):
	global aaCode

	for k,v in aaCode.items():
		if x == v:
			return k
		elif  x == k:
			return v

	if len(x) == 1:
		return 'UNK'
	return 'X'

'''
GL 09/03/2017
 This is minimu implementation of PDB parsing and coordinates mangament

'''

''' PDB format reminder

COLUMNS		DATA  TYPE	FIELD		DEFINITION
-------------------------------------------------------------------------------------
 1 -  6		Record name   "ATOM  "
 7 - 11		Integer	   serial	   Atom  serial number.
13 - 16		Atom		  name		 Atom name.
17			 Character	 altLoc	   Alternate location indicator.
18 - 20		Residue name  resName	  Residue name.
22			 Character	 chainID	  Chain identifier.
23 - 26		Integer	   resSeq	   Residue sequence number.
27			 AChar		 iCode		Code for insertion of residues.
31 - 38		Real(8.3)	 x			Orthogonal coordinates for X in Angstroms.
39 - 46		Real(8.3)	 y			Orthogonal coordinates for Y in Angstroms.
47 - 54		Real(8.3)	 z			Orthogonal coordinates for Z in Angstroms.
55 - 60		Real(6.2)	 occupancy	Occupancy.
61 - 66		Real(6.2)	 tempFactor   Temperature  factor.
77 - 78		LString(2)	element	  Element symbol, right-justified.
79 - 80		LString(2)	charge	   Charge  on the atom.

'''

class Parser(object):
	def __init__(self):
		pass

	def load(self, **kwargs):
		structureObj = Structure();

		def _read(stream):
			currentAtomList = structureObj.pop();
			newModel = False
			for l in stream:
				#print ">>" + l + "<<"
				if l.startswith("ATOM "):
					if newModel:
						currentAtomList = structureObj.pop();
						newModel = False
					currentAtomList.append(Atom(string=l))
				elif l.startswith("ENDMDL"):
					newModel = True

				if l.startswith("SEQRES "):
					structureObj.loadSEQRES(l)

		if 'file' in kwargs:
			with open(kwargs['file'], 'r') as f:
				_read(f)

		if 'stream' in kwargs:
			#print kwargs['stream'] + "<--"
			_read(kwargs['stream'].splitlines(True))

		if len(structureObj.model[0]) == 0:
			raise ValueError("Empty coordinates")

		# Assign a name
		if not structureObj.name:
			if 'file' in kwargs:
				base = os.path.basename(kwargs['file'])
				structureObj.name = os.path.splitext(base)[0]

		if not structureObj.name:
			structureObj.name = "a protein"

		return structureObj



class Model(object):
	def init(self):
		pass



class Structure(object):
	def __init__(self):
		self.name = None
		self.model = []
		self._trace = None
		self._fasta = None
		self.currModel = 1
		self.bufferSelection = []
		self._chainList = None
		self._residues = None
		self._resID = None
		self.SEQRES = {}

	def __len__(self):
		return len(self.model[self.currModel - 1])

	def loadSEQRES(self, line):

		aaSeq = line.split()[2:]
		segID = aaSeq.pop(0)
		aaSeq.pop(0)
		if segID not in self.SEQRES:
			self.SEQRES[segID] = ''
		
		self.SEQRES[segID] += ''.join([ translate(aa) for aa in aaSeq ])
		


	def setCoordinateFromDictorize(self, dictorizedSelf):
		atomList = self.model[self.currModel - 1]
		#print atomList
		#print len(atomList)
		#print len(dictorizedSelf['x'])
		for i in range(len(atomList)):
	#        print i
			self.model[self.currModel - 1][i].x = dictorizedSelf['x'][i]
			self.model[self.currModel - 1][i].y = dictorizedSelf['y'][i]
			self.model[self.currModel - 1][i].z = dictorizedSelf['z'][i]

	@property
	def atomDictorize(self):
		atomList = self.model[self.currModel - 1]
		(x, y, z, seqRes, chainID, resName, name) = self.atomVectorize;
		return { "x" : x,
				 "y" : y,
				 "z" : z,
				 "seqRes" : seqRes,
				 "chainID" : chainID,
				 "resName" : resName,
				 "name" : name
				}

	@property
	def atomVectorize(self):
		atomList = self.model[self.currModel - 1]
		seqRes = []
		chainID = []
		x = []
		y = []
		z = []
		resName = []
		name = []
		for a in atomList:
			seqRes.append(a.seqRes)
			chainID.append(a.chainID)
			x.append(a.x)
			y.append(a.y)
			z.append(a.z)
			name.append(a.name)
			resName.append(a.resName)
		return (x, y, z, seqRes, chainID, resName, name)

	@property
	def getResID(self):
		if not self._resID:
			self._resID = [residue.data[0].getResID for residue in self.byres()]
		return self._resID
	# Browse model and create a new pdbObject with a sinlge model that satifise the provided conditions
	# Eg: ask for the 1st model w/ oligomeric composition chain=["A","B","C"]
	def modelReduce(self, **kwargs):
		modelMatch = -1

		if 'chain' in kwargs:
			chainList = kwargs['chain']
			for i, model  in enumerate(self.model):
				for segID in chainList:
					found = False
					for p_atom in model:
						if p_atom.chainID == segID:
							print(p_atom)
							found = True
							break
					if not found:
						break
				if found:
					modelMatch = i
					break
		if modelMatch < 0:
			return None

		print(" Model ", str(modelMatch), "satisfies conditions")

		cloneObj = Structure()
		clonedModel =  self.model[modelMatch]
		currentAtomList = cloneObj.pop()
		for p_atom in clonedModel:
			currentAtomList.append(copy.deepcopy(p_atom))
		return cloneObj

	def peptideSeed(self): # A function to operate w/ pyprotein utilities like nw_custom

		return {
			'id' : self.name,
			'desc' : 'pdb file fasta translation',
			'seq' : self.fasta
		}

# move a structure along a vector
	def nudge(self, C):
		V = [ numpy.array( d.coordinates ) for d in self.atomRecord ]
		V += C
		for i, a in enumerate (V):
			self.model[self.currModel - 1][i].x = V[i,0]
			self.model[self.currModel - 1][i].y = V[i,1]
			self.model[self.currModel - 1][i].z = V[i,2]


	def centerOrigin(self):
		V = [ numpy.array( d.coordinates ) for d in self.atomRecord ]
	# move to centroid
		C = sum(V)/len(V)
		print(C)
		V -= C
		for i, a in enumerate (V):
			#print
			self.model[self.currModel - 1][i].x = V[i,0]
			self.model[self.currModel - 1][i].y = V[i,1]
			self.model[self.currModel - 1][i].z = V[i,2]

		return C



	def rotate(self, **kwargs):
# 1st we move centroid to origin
# 2nd we rotate
# 3rd we move centroid back to original
# U=RotationMatrix  OR (alpha=0, beta=0, gamma=0)
		centeringBool = True

		if ('nocenter' in kwargs):
			if kwargs['nocenter']:
				centeringBool = False

		if 'U' not in kwargs and ('alpha' not in kwargs or 'beta' not in kwargs or 'gamma' not in kwargs):
			raise ValueError("Specify a 3x3 rotation matrix or alpha,beta,angle for rotation around X,Y and Z axis");

		V = [ numpy.array( d.coordinates ) for d in self.atomRecord ]

		if centeringBool:
	# move to centroid
			C = sum(V)/len(V)
			V -= C

		if 'U' not in kwargs:
			alpha=kwargs['alpha']
			beta=kwargs['beta']
			gamma=kwargs['gamma']

	# Rotate around origin
			Rx = numpy.matrix([
						[1,          0,           0],
						[0, math.cos(alpha), -math.sin(alpha)],
						[0, math.sin(alpha),  math.cos(alpha)]
					   ])
			Ry = numpy.matrix([
						[math.cos(beta),  0,  math.sin(beta)],
						[0,          1,          0],
						[-math.sin(beta), 0,  math.cos(beta)]
						])
			Rz = numpy.matrix([
						[math.cos(gamma),  -math.sin(gamma), 0],
						[math.sin(gamma),  math.cos(gamma), 0],
						[0,                    0, 1]
						])

		#print V[0]


		for i, a in enumerate (V):

			v = numpy.matrix([[a[0]],
						   [a[1]],
						   [a[2]]])
			if 'U' not in kwargs:
				v_ = Rz*Ry*Rx*v
			else :
				v_ = kwargs['U']*v

			self.model[self.currModel - 1][i].x = v_[0,0]
			self.model[self.currModel - 1][i].y = v_[1,0]
			self.model[self.currModel - 1][i].z = v_[2,0]


	## Move back
	 #   t = numpy.matrix([[C[0]],
					   #[C[1]],
					   #[C[2]]])
	  #  tBack = t
		if centeringBool:
			for a in self.model[self.currModel - 1]:
				a.x += C[0]
				a.y += C[1]
				a.z += C[2]
		#print self.model[self.currModel - 1][0]


	def __str__(self):
		asString = '' # some header here
		for n, model in enumerate(self.model):
			self.currModel = n + 1
			asString += ''.join([ str(a) if str(a).endswith("\n") else str(a) + '\n' for a in self.atomRecord ])
			#asString += "ENDMDL\n"
		return asString

	def clone(self):
		cloneObj = Structure()

		for struct in self.model:
			currentAtomList = cloneObj.pop();
			for p_atom in struct:
				currentAtomList.append(copy.deepcopy(p_atom))

		return cloneObj

	@property
	def trace(self):
		if not self._trace:
			self._trace = []
			for r in self.byres():
				for a in r:
					if a.name == 'CA':
						self._trace.append(a)

		return self._trace

	@property
	def residueNumber(self):
		return len([ res for res in self.byres() ])

	# Return consecutive 1letter amino-acid for residues with defined Calpha coordinates
	@property
	def fasta(self):
		if not self._fasta:
			self._fasta = ''.join([ r.fasta for r in self.byres(strict=True) ])
		return self._fasta

	@property
	def atomRecord(self):
		for atom in self.model[self.currModel - 1]:
			yield atom

	@property
	def chainList(self):
		if not self._chainList:
			buf = [ r.chain for r in self.byres() ]
			self._chainList = list(set(buf))
		return self._chainList

	def byres(self, strict=False):

		i = 0
		data = self.model[self.currModel - 1]
		w_resSeq = data[0].resSeq
		w_chainID = data[0].chainID
		w_iCode = data[0].iCode
		i_start = 0
		i_stop = 0
		while i < len(data):
			#print str(i) + " -- " + str(len(data))
			cur_atom = data[i]
			#print str(w_resSeq) + ' vs ' + str(cur_atom.resSeq) + ' and ' + w_chainID + ' vs ' + cur_atom.chainID
			if w_resSeq != cur_atom.resSeq or w_chainID != cur_atom.chainID or w_iCode != cur_atom.iCode:
				#i_stop = i + 1
				#print "pulling at " + str(i_start) + ":" + str(i)
				x = Residue(data[i_start:i])
				if (strict and x.hasCalpha) or not strict:
					yield x
				w_resSeq  = cur_atom.resSeq
				w_chainID = cur_atom.chainID
				w_iCode = cur_atom.iCode

				i_start = i
			i += 1

		x = Residue(data[i_start:i])
		if (strict and x.hasCalpha) or not strict:
			yield x

	def pop(self): # Append an empty array to model list and return its referenc

		self.model.append([])
		return self.model[-1]

	# Return a pdb object limited to specified chain ID
	def chain(self, chain):
		chainList = []
		try:
			chainList = [e for e in chain]
		except TypeError:
			chainList.append(chain)

		cloneObj = Structure()

		for struct in self.model:
			currentAtomList = cloneObj.pop();
			for p_atom in struct:
				if p_atom.chainID in chainList:
					currentAtomList.append(copy.deepcopy(p_atom))

		cloneObj.name = self.name + ':' + ','.join(chainList)
		if len(cloneObj.model[0]) == 0:
			print("Empty selection for chain(s) \"", str(chain), "\"")
			return None

		return cloneObj


	def select(exp): # TO DO
		pass
		#
		# Should implment binary tree construction from the parsing of pymol-like expression
		# All tree relationship could then be mapped to set operation on atoms lists
		#
		#        http://pyparsing.wikispaces.com/

class Residue(object):
	def __init__(self, atomArray):
		self.data = atomArray

	def __len__(self):
		return len(self.data)

	def __repr__(self):
	   return '\n'.join([ str(a) for a in self.data ])

	def __getitem__(self, key):
		return self.data[key]

	def __iter__(self):
		for d in self.data:
			yield d

	def __str__(self):
		return "%3s %4d%c %s" %(self.name, self.num, self.iCode, self.chain)

	def __eq__(self, other):
		return (self.chain, self.seqRes) == (other.chain, other.seqRes)

	def __ne__(self, other):
		return not(self == other)

	@property
	def fasta(self):
		return translate(self.data[0].resName)

	@property
	def id(self):
		return self.data[0].resName + str(self.data[0].resSeq) + ':' + self.data[0].chainID

	@property
	def name(self):
		return self.data[0].resName

	@property
	def seqRes(self):
		return self.data[0].seqRes

	@property
	def num(self):
		return self.data[0].resSeq

	@property
	def chain(self):
		return self.data[0].chainID

	@property
	def iCode(self):
		return self.data[0].iCode

	@property
	def seqRes(self):
		return self.data[0].seqRes

	@property
	def hasCalpha(self):
		for a in self:
			if a.name == "CA":
				return True
		return False

	def asPdbRecord(self):
		return self.__repr__()

class Atom(object):
	def __init__(self, **kwargs):
		if 'string' in kwargs:
			buf=kwargs['string'].rstrip()
			bufLen = len(buf)
			if bufLen < 54:
				raise TypeError("Unexpected atom record length \"" + str(bufLen) + "\" on input at \n" + buf)
			self.recordName=buf[0:6].replace(" ", "")
			self.serial=int(buf[6:11])
			self.name=buf[12:16].replace(" ", "")
			self.altLoc= buf[16]
			self.resName=buf[17:20].replace(" ", "")
			self.chainID=buf[21]
			self.resSeq=buf[22:26] # Not sure
			self.iCode=buf[26]
			self.x=float(buf[30:38])
			self.y=float(buf[38:46])
			self.z=float(buf[46:54])
			self.occupancy  = float(buf[54:60])           if bufLen > 54 else None
			self.tempFactor = float(buf[60:66])           if bufLen > 60 else None
			self.element    = buf[76:78].replace(" ", "") if bufLen > 76 else None
			self.charge     = buf[78:80].replace(" ", "") if bufLen > 78 else None

	@property
	def seqRes(self):
		return str(self.resSeq) + str(self.iCode)

	@property
	def coordinates(self):
		return [self.x, self.y, self.z]

	@property
	def toVector(self):
		return [self.x, self.y, self.z, self.resSeq + self.iCode, self.chainID]

	@property
	def getResID(self):
	#   return self.resName + str(self.resSeq) + ':' + self.chainID + ':' + self.iCode
		return self.resName + self.seqRes + ':' + self.chainID

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
