from __future__ import division
import numpy, scipy.sparse, scipy.sparse.linalg, time
import MBIlib


class MBI(object):

    def __init__(self, P, xs, ms0=None, ks0=None):
        nx = len(xs)
        nf = 1 if len(P.shape) is nx else P.shape[nx]

        ns = numpy.array(P.shape[:nx],order='F')
        ms = numpy.array([ms0[i] if ms0 is not None else int(ns[i]/3) for i in range(nx)],int)
        ks = numpy.array([min(ks0[i] if ks0 is not None else 4, ms[i]) for i in range(nx)],int)
        nT = numpy.prod(ns)

        self.xs = xs
        self.ns, self.ms, self.ks = ns, ms, ks
        self.nx, self.nf, self.nT = nx, nf, nT

        ts = []
        for i in range(nx):
            k, m, n = ks[i], ms[i], ns[i]
            d = MBIlib.knotopen(k, k+m)
            ts.append(MBIlib.paramuni(k+m, m, n, d))

        t = numpy.zeros(P.shape[:nx]+(nx,),order='F')
        for ind,x in numpy.ndenumerate(t):
            t[ind] = ts[ind[-1]][ind[ind[-1]]]
        t = t.reshape((nT,nx),order='F')
        self.t = t
        P = numpy.reshape(P,(nT,nf),order='F')

        B = self.getJacobian(0, 0)
        BT = B.transpose()
        BTB, BTP = BT.dot(B), BT.dot(P)

        nC = numpy.prod(ms)
        C = numpy.zeros((nC,nf),order='F')
        for i in range(nf):
            C[:,i] = scipy.sparse.linalg.cg(BTB,BTP[:,i])[0]

        Cx = []
        for i in range(nx):
            k, m, n = ks[i], ms[i], ns[i]
            B = self.assembleJacobian(0, 0, 1, n, k*n, k, m, ts[i])
            BT = B.transpose()
            BTB, BTP = BT.dot(B), BT.dot(xs[i])
            Cx.append(scipy.sparse.linalg.cg(BTB,BTP)[0])
            Cx[-1][0] = xs[i][0]
            Cx[-1][-1] = xs[i][-1]

        self.C, self.Cx = C, Cx

    def assembleJacobian(self, d1, d2, nx, nP, nB, ks, ms, t):
        Ba, Bi, Bj = MBIlib.computejacobian(d1, d2, nx, nP, nB, ks, ms, t)
        return scipy.sparse.csc_matrix((Ba,(Bi,Bj)))

    def getJacobian(self, d1, d2):
        return self.assembleJacobian(d1, d2, self.nx, self.nT, self.nT*numpy.prod(self.ks), self.ks, self.ms, self.t)

    def evaluate(self, x, d1=0, d2=0):
        nx, nf, nT = self.nx, self.nf, self.nT
        ns, ms, ks = self.ns, self.ms, self.ks
        nP = x.shape[0]

        minV = numpy.min
        maxV = numpy.max
        for i in range(nx):
            if minV(x[:,i]) < minV(self.Cx[i]):
                print 'MBI error: min value out of bounds', i, minV(x[:,i]), minV(self.Cx[i])
                #raise Exception('MBI evaluate error: min value out of bounds')
            if maxV(x[:,i]) > maxV(self.Cx[i]):
                print 'MBI error: max value out of bounds', i, maxV(x[:,i]), maxV(self.Cx[i])
                #raise Exception('MBI evaluate error: max value out of bounds')

        t = numpy.zeros((nP,nx),order='F')
        for i in range(nx):
            t[:,i] = MBIlib.inversemap(ks[i], ms[i], nP, x[:,i], self.Cx[i])

        i1, i2 = max(0, d1-1), max(0, d2-1)
        nC, nCx1, nCx2 = self.C.shape[0], self.Cx[i1].shape[0], self.Cx[i2].shape[0]
        return MBIlib.evaluate(d1, d2, nx, nf, nC, nCx1, nCx2, nP, ks, ms, t, self.C, self.Cx[i1], self.Cx[i2])
