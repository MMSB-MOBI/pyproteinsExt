#import pyproteinsExt.structure.coordinates
from multiprocessing import Pool
import pyproteins.sequence.msa
import pyproteins.alignment.nw_custom
import pyproteins.sequence.peptide

import math,numpy,subprocess

import xml.etree.ElementTree

def needle(pdbObjectOne, pdbObjectTwo):

    pepA = pyproteins.sequence.peptide.Entry(pdbObjectOne.peptideSeed())
    pepB = pyproteins.sequence.peptide.Entry(pdbObjectTwo.peptideSeed())

    #print pepA
    #print pepB
    nw = pyproteins.alignment.nw_custom.nw(gapOpen=-10, gapExtend=-0.5)
    ali = nw.align(pepA, pepB)
    return ali


def minDist(atomList1, atomList2):
    dMin = 9999999.99999
    for v1 in [ numpy.array((a.x ,a.y, a.z)) for a in atomList1 ]:
        for v2 in [ numpy.array((b.x ,b.y, b.z)) for b in atomList2 ]:
            d = numpy.linalg.norm(v1 - v2)
            if d < dMin:
                dMin = d
    return dMin


def euclidianDist(atom1, atom2):
    mtx_atom1 = numpy.array(atom1.coordinates)
    mtx_atom2 = numpy.array(atom2.coordinates)
    return numpy.linalg.norm(mtx_atom1 - mtx_atom2)


class Cell(object):
    def __init__(self, iLabel, jLabel, value):
        self.iLabel = iLabel
        self.jLabel = jLabel
        self.value = value

    def __str__(self):
        return '(' + self.iLabel + ', ' + self.jLabel + ') : ' + str(self.value)



class MeshMap(object):
    def __init__(self, s1, s2):
        pass


class ContactMap(object):
    def __init__(self, s1, s2):
        self.nb_raws = s1.residueNumber
        self.nb_columns = s2.residueNumber
        self.mtx = numpy.zeros((self.nb_raws, self.nb_columns))

        self._resArrayOne = [ x for x in s1.byres() ]
        self._resArrayTwo = [ x for x in s2.byres() ]

        self.rl = [r.id for r in self._resArrayOne]
        self.cl = [r.id for r in self._resArrayTwo]

        for i1, r1 in enumerate(self._resArrayOne):
            for i2, r2 in enumerate(self._resArrayTwo):
                self.mtx[i1, i2] = minDist(r1, r2)

    # matrix element accessor along with row/column label, mostly for inspection
    def __getitem__(self, tup):
        x, y = tup
        return Cell(self.rl[x], self.cl[y], self.mtx[x, y])

    def __str__(self):
        asString = 'Contact map ' + str(len(self.rl)) + 'x' + str(len(self.cl)) + '\n'
        asString += '\t\t' + ''.join(["%9s" % r for r in self.cl]) + '\n'
        for i, row in enumerate(self.mtx):
            asString += "%9s" % self.rl[i]
            asString += ' '.join(["%9.1f" % d for d in row])
            asString += '\n'
        return asString

    def Q(self, d=4.5):
        f = lambda x: 1 if x < d else 0
        DistToQ = numpy.vectorize(f)
        return DistToQ(self.mtx)

    def residuesInterfacialBool(self, d=4.5):
        q = self.Q(d)
        def f(row):
            if numpy.sum(row) > 0:
                return 1
            return 0

        return interfaceBoolList( [({ 'num' : r.num, 'name' : r.name, 'chain' : r.chain }, f(q[i])) for i,r in enumerate(self._resArrayOne)],
            [({ 'num' : r.num, 'name' : r.name, 'chain' : r.chain}, f(q.T[i])) for i,r in enumerate(self._resArrayTwo)] )

    def weighted_contact_number(self):
        mtx = numpy.zeros((self.nb_raws, self.nb_columns))

        for i in range(len(self._resArrayOne)):
            for j in range(len(self._resArrayTwo)):
                if i != j:
                    mtx[i, j] = 1/(self.mtx[i, j]**2)

        return mtx.sum(axis=1,dtype=float)


class ContactMap_intra(object):
    def __init__(self, struc, cutoff=5.0):
        self.nb = struc.residueNumber
        self.mtx = numpy.zeros((self.nb,self.nb))

        self._resArray = [ x for x in struc.byres() ]

        self.l = [r.id for r in self._resArray]

        self.counter_infcutoff = 0

        for i in range(len(self._resArray)):
            for j in range(i+1,len(self._resArray)):
                self.mtx[i, j] = minDist(self._resArray[i], self._resArray[j])

                if self.mtx[i, j] < cutoff:
                    self.counter_infcutoff += 1

    # matrix element accessor along with row/column label, mostly for inspection
    def __getitem__(self, tup):
        x, y = tup
        return Cell(self.l[x], self.l[y], self.mtx[x, y])

    def __str__(self):
        asString = 'Contact map intra ' + str(len(self.l)) + 'x' + str(len(self.l)) + '\n'
        asString += '\t\t' + ''.join(["%9s" % r for r in self.l]) + '\n'
        for i, row in enumerate(self.mtx):
            asString += "%9s" % self.l[i]
            asString += ' '.join(["%9.1f" % d for d in row])
            asString += '\n'
        return asString

    def weighted_contact_number(self):
        mtx = numpy.zeros((self.nb,self.nb))

        for i in range(len(self._resArray)):
            for j in range(i+1,len(self._resArray)):
                mtx[i, j] = 1/(self.mtx[i, j]**2)

        symetric_mtx = mtx + mtx.T

        return symetric_mtx.sum(axis=1,dtype=float)


class ContactMap_intra_grid(object):
    def __init__(self, struc, cutoff=5.0):
        self.struc = struc
        self._resArray = self.struc.getResID
        self._residuePairRegistry = {}
        self.cutoff = cutoff
        self._parsing()
        self._build_grid()
        self._list_next_unique_neighbors()
        self._calculate_distances()
        self._build_ContactMap()

    # matrix element accessor along with row/column label, mostly for inspection
    def __getitem__(self, tup):
        x, y = tup
        return Cell(self._resArray[x], self._resArray[y], self.mtx[x, y])

    def __str__(self):
        asString = 'Contact map intra grid ' + str(len(self._resArray)) + 'x' + str(len(self._resArray)) + '\n'
        asString += '\t\t' + ''.join(["%9s" % r for r in self._resArray]) + '\n'
        for i, row in enumerate(self.mtx):
            asString += "%9s" % self._resArray[i]
            asString += ' '.join(["%9.1f" % d for d in row])
            asString += '\n'
        return asString

    def _parsing(self):
        # ---
        # PARSING
        # ---

        self.counter_atoms = 0

        self.min_x,self.min_y,self.min_z = 9999999.99999,9999999.99999,9999999.99999
        self.max_x,self.max_y,self.max_z = -9999999.99999,-9999999.99999,-9999999.99999

        for atom in self.struc.atomRecord:
            # ---
            # GET MINIMUMS AND MAXIMUMS OF EACH COORDINATES X,Y,Z
            # ---

            self.counter_atoms += 1

            self.min_x = atom.x if atom.x < self.min_x else self.min_x
            self.min_y = atom.y if atom.y < self.min_y else self.min_y
            self.min_z = atom.z if atom.z < self.min_z else self.min_z

            self.max_x = atom.x if atom.x > self.max_x else self.max_x
            self.max_y = atom.y if atom.y > self.max_y else self.max_y
            self.max_z = atom.z if atom.z > self.max_z else self.max_z

        # ---
        # GET DIMENSIONS OF THE 3D GRID
        # ---

        self.dim_x = int(math.floor((self.max_x-self.min_x)/self.cutoff))+1
        self.dim_y = int(math.floor((self.max_y-self.min_y)/self.cutoff))+1
        self.dim_z = int(math.floor((self.max_z-self.min_z)/self.cutoff))+1

    def _coor_Cartesian2Grid(self,atom):
        # ---
        # GET AN ATOM AND RETURN ITS NORMALIZED COORDINATES
        # ---
        return [int(math.floor((atom.x-self.min_x)/self.cutoff)),int(math.floor((atom.y-self.min_y)/self.cutoff)),int(math.floor((atom.z-self.min_z)/self.cutoff))]

    def _build_grid(self):
        self.grid_3D = numpy.empty((self.dim_x,self.dim_y,self.dim_z),dtype=numpy.object_)
        self.grid_3D.fill([])
        self.grid_3D = numpy.frompyfunc(list,1,1)(self.grid_3D)

        for atom in self.struc.atomRecord:
            [i,j,k] = self._coor_Cartesian2Grid(atom)
            self.grid_3D[i,j,k].append(atom)

    def _iter_next_unique_neighbors(self, x, y, z):
        value_list = []

        for i in range(x, x+2):
            if i == self.dim_x:
                break

            for j in range(y-1, y+2):
                if j < 0:
                    continue
                if j == self.dim_y:
                    break

                for k in range(z-1, z+2):
                    if k < 0:
                        continue
                    if k == self.dim_z:
                        break

                    if ( i > x or j > y or (j > y-1 and k > z-1) ) and self.grid_3D[i,j,k] != []:
                        value_list.append(self.grid_3D[i,j,k])

        return value_list

    def _list_next_unique_neighbors(self):
        self.value_LIST = []

        for i in range(self.dim_x):
            for j in range(self.dim_y):
                for k in range(self.dim_z):

                    if self.grid_3D[i,j,k] != []:
                        self.value_LIST.append(self._iter_next_unique_neighbors(i,j,k))

    def _calculate_distances(self):
        self.dist_LIST = []
        self.counter_dist = 0

        for i in range(len(self.value_LIST)):
            current = self.value_LIST[i][0]
            neighbors_flat_LIST = [j for i2 in self.value_LIST[i][1:] for j in i2]

            for j in range(len(current)):

                for k in range(j+1,len(current)):
                    self.dist_LIST.append([[current[j],current[k]],euclidianDist(current[j],current[k])])
                    self.counter_dist += 1

                for l in range(len(neighbors_flat_LIST)):
                    self.dist_LIST.append([[current[j],neighbors_flat_LIST[l]],euclidianDist(current[j],neighbors_flat_LIST[l])])
                    self.counter_dist += 1

    def _build_ContactMap(self):
        self.mtx = numpy.empty((len(self._resArray),len(self._resArray)),dtype=float)
        self.mtx.fill(9999999.99999)

        for idx in self.dist_LIST:
            i = self._resArray.index(idx[0][0].getResID)
            j = self._resArray.index(idx[0][1].getResID)
            if idx[1] > self.cutoff:
                continue
            if idx[1] < self.mtx[i, j] and i != j:
                self.mtx[i, j] = idx[1]
                self.mtx[j, i] = idx[1]
                self._registerResiduePair(i, j)

    def _registerResiduePair(self, i, j):
        (resID_A, resID_B) = (self._resArray[i], self._resArray[j]) if i < j else (self._resArray[j], self._resArray[i])
        if resID_A not in self._residuePairRegistry:
            self._residuePairRegistry[resID_A] = [resID_B]
        elif resID_B not in self._residuePairRegistry[resID_A]:
            self._residuePairRegistry[resID_A].append(resID_B)


class ContactOrder(object):
    def __init__(self, struc_name, struc, cutoff=5.0):
        self._contiguous = None
        self.list_CO = []

        for curr_ch in struc.chainList:

            self.sub_struct = struc.chain(curr_ch)
            cm_intra = ContactMap_intra_grid(self.sub_struct,cutoff)
            counter,SUM = 0,0

            for i in range(len(cm_intra._resArray)):
                for j in range(i+1,len(cm_intra._resArray)):
                    if cm_intra.mtx[i, j] < cutoff:
                        SUM += abs(i-j)
                        counter += 1

            self.list_CO.append([struc_name,curr_ch,'{:.2f}'.format(float(SUM)/counter),self.contiguous])

    def __iter__(self):
        for item in self.list_CO:
            yield item

    @property
    def contiguous(self):
        if not self._contiguous:
            self._contiguous = True
            trace_list = [atom for atom in self.sub_struct.trace]
            for i in range(1,len(trace_list)):
                curr_atom = trace_list[i]
                prev_atom = trace_list[i-1]
                dist = minDist([curr_atom],[prev_atom])
                if float(dist) > 4.0:
                    self._contiguous = False
                    break
        return self._contiguous


class interfaceBoolList(object):
    def __init__(self, l1, l2):
        self.l1 = l1
        self.l2 = l2

    def __getitem__(self, k):
        if(k == 0):
            l = self.l1
        elif(k == 1):
            l = self.l2
        else :
            raise ValueError("index \"" + str(k) + "\"out of bounds ")
        s = self.asString(l)
        return s

    def toList(self, k):
        data = []

        if(k == 0):
            l = self.l1
        elif(k == 1):
            l = self.l2
        else :
            raise ValueError("index \"" + str(k) + "\"out of bounds ")
        for d in l:
            data.append({ 'res' : d[0], 'cc' : d[1] })
        return data

    def asString(self, l):
        return '\n'.join([ "%4s %s %s %d" % (d[0]['num'], d[0]['name'], d[0]['chain'], d[1]) for d in l ])


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = numpy.dot(P, U)
    return P


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = numpy.dot(numpy.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = numpy.linalg.svd(C)
    d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = numpy.dot(V, W)

    return U


def quaternion_rmsd(P, Q):
    """
    based on doi:10.1016/1049-9660(91)90036-O
    Rotate matrix P unto Q and calculate the RMSD
    """
    rot = quaternion_rotate(P, Q)
    P = numpy.dot(P,rot)
    return rmsd(P, Q)


def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3,:3]
    return rot


def makeW(r1,r2,r3,r4=0):
    """
    matrix involved in quaternion rotation
    """
    W = numpy.asarray([
             [r4, r3, -r2, r1],
             [-r3, r4, r1, r2],
             [r2, -r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return W


def makeQ(r1,r2,r3,r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = numpy.asarray([
             [r4, -r3, r2, r1],
             [r3, r4, -r1, r2],
             [-r2, r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return Q


def quaternion_rotate(X, Y):
    """
    Calculate the rotation
    """
    N = X.shape[0]
    W = numpy.asarray([makeW(*Y[k]) for k in range(N)])
    Q = numpy.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = numpy.asarray([numpy.dot(Q[k].T,W[k]) for k in range(N)])
    W_minus_Q = numpy.asarray([W[k] - Q[k] for k in range(N)])
    C1 = -numpy.sum(Qt_dot_W,axis=0)
    C2 = 0.5*N
    C3 = numpy.sum(W_minus_Q,axis=0)
    A = numpy.dot(C3.T,C3)*C2-C1
    eigen = numpy.linalg.eigh(A)
    r = eigen[1][:,eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """

    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return numpy.sqrt(rmsd/N)

def blastThem(pdbA, pdbB, mode):

    def _f(pdb, name, fname):
        with open (fname, 'w') as f:
            s = '>' + name + '\n' + ''.join([ x + '\n' if i%80 == 0 and i > 0 else x for i,x in enumerate (pdb.fasta) ])
            if len(pdb.fasta) % 81 != 0:
                s+='\n'
            f.write(s)

    _f(pdbA, 'proteinA', 'proteinA.fasta')
    _f(pdbB, 'proteinB', 'proteinB.fasta')

    def _blastpgp(name):
        j = 1
        print( 'basltpg stuff', name, str(j) )
        try:
            subprocess.call(['blastpgp', '-j', str(j), '-i', name + '.fasta', '-d', 'nr70i', '-m', '7', '-o', name + '.blast'])
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                print('blastpgp not found')
        # handle file not found error.
            else:
                print('Something wrong w/ blastpgp')
        return name + '.blast'


    if mode == "debug":
        print("no blast ran, debug mode")
        return 'blastDebugMode'

    print("running psiblast")
    with Pool(2) as p:
        x = p.map(_blastpgp, ['proteinA', 'proteinB'])
    
def clustThem(mode):

    if mode == "debug":
        print("no clustal ran, debug mode")
        pMsa = pyproteins.sequence.msa.Msa(fileName='proteins.aln')
        return pMsa

    root = xml.etree.ElementTree.parse('proteinA.blast')
    hits = [ elem.text for elem in root.findall('BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_id') ]
    mfastaRaw = ''
    with open('proteinA.fasta', 'r') as myfile:
        mfastaRaw=myfile.read()
    with open('proteinB.fasta', 'r') as myfile:
        mfastaRaw+=myfile.read()
    for name in list(set(hits)):
        output = subprocess.check_output(['fastacmd', '-d', 'nr70i', '-s', str(name)])
        mfastaRaw += output

    with open('proteins.mfasta', 'w') as f:
        f.write(mfastaRaw)

    subprocess.call(['clustalw2', '-INFILE=proteins.mfasta'])

    pMsa = pyproteins.sequence.msa.Msa(fileName='proteins.aln')

    return pMsa

def aliFit(structA, structB, aliArrayOne, aliArrayTwo):
    equiResNum = []
    i = -1
    j = -1

    for x,y in zip(aliArrayOne, aliArrayTwo):
        xBool = False
        yBool = False
        if x != '-' :
            xBool = True
            i += 1
        if y != '-' :
            yBool = True
            j += 1

        if not xBool or not yBool:
            continue
        equiResNum.append( { 'iRes' : { 'num' : i , 'name' : x }, 'jRes' : { 'num' : j , 'name' : y } } )

    print("Superimposing on", str(len(equiResNum)), 'CA coordinates')

    P = [ numpy.array( structA.trace[ d['iRes']['num']].coordinates) for d in equiResNum ]
    Q = [ numpy.array( structB.trace[ d['jRes']['num']].coordinates) for d in equiResNum ]

    normal_rmsd = rmsd(P, Q)

    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    U = kabsch(P, Q)
    #p_all -= Pc
    #p_all = numpy.dot(p_all, U)
    #write_coordinates(p_atoms, p_all, title="{} translated".format(args.structure_a))
    k_rmsd = kabsch_rmsd(P, Q)
    q_rmsd = quaternion_rmsd(P, Q)

    print ("Normal RMSD:", normal_rmsd)
    print ("Kabsch RMSD:", k_rmsd)
    print ("Quater RMSD:", q_rmsd)
    return (U, normal_rmsd, k_rmsd, q_rmsd)

def fit(structObj1, structObj2, **kwargs):

    mode = "blast"
    aliSeq1 = None
    aliSeq2 = None
    if 'mode' in kwargs:
        mode = kwargs['mode']

    if mode == "needle":
        ali = needle(structObj1, structObj2)
        (aliSeq1, aliSeq2) = ali.aaWords

    elif mode == "blast":
        blastThem(structObj1, structObj2, mode)
        msaObj = clustThem(mode)

        aliSeq1 = msaObj[0]['sequence']
        aliSeq2 = msaObj[1]['sequence']

    print (aliSeq1, "\n", aliSeq2)

    _tuple = aliFit(structObj1, structObj2, aliSeq1, aliSeq2)
    return _tuple


