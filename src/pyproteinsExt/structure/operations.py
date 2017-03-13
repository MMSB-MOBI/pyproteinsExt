import numpy as np


def minDist(atomList1, atomList2):
    dMin = 9999999.99999
    for v1 in [ np.array((a.x ,a.y, a.z)) for a in atomList1 ]:
        for v2 in [ np.array((b.x ,b.y, b.z)) for b in atomList2 ]:
            d = np.linalg.norm(v1 - v2)
            if d < dMin:
                dMin = d
    return d

class ContactMap(object):
    def __init__(self, s1, s2):
        self.mtx = np.zeros((s1.residueNumber, s2.residueNumber))

        self._resArrayOne = [ x for x in s1.byres() ]
        self._resArrayTwo = [ x for x in s2.byres() ]

        self.rl = [r.id for r in self._resArrayOne]
        self.cl = [r.id for r in self._resArrayTwo]
        for i1, r1 in enumerate(self._resArrayOne):
            for i2, r2 in enumerate(self._resArrayTwo):
                self.mtx[i1, i2] = minDist(r1, r2)

    def __str__(self):
        print 'Contact map ' + str(len(self.rl)) + 'x' + str(len(self.cl)) + '\n'
        asString = '         ' + ''.join(["%9s" % r for r in self.cl]) + '\n'
        for i, row in enumerate(self.mtx):
            asString += "%9s" % self.rl[i]
            asString += ' '.join(["%9.1f" % d for d in row])
            asString += '\n'
        return asString

    def Q(self, d=4.5):
        f = lambda x: 1 if x < d else 0
        DistToQ = np.vectorize(f)
        return DistToQ(self.mtx)

    def residuesInterfacialBool(self, d=4.5):
        q = self.Q(d)
        def f(row):
            if np.sum(row) > 0:
                return 1
            return 0

        return ( [(r.id, f(q[i])) for i,r in enumerate(self._resArrayOne)],
                 [(r.id, f(q.T[i])) for i,r in enumerate(self._resArrayTwo)])

