"""Independent arbitrary-precision radial quadrature for frozen affine media."""
from __future__ import annotations
import mpmath as mp
import numpy as np


def exact(value):
    numerator, denominator = float(value).as_integer_ratio()
    return mp.mpf(numerator) / denominator


class AffineReference:
    def __init__(self, field, xyz, p, coefficient, frequency, dps=80, order=16):
        self.dps, self.order = dps, order
        with mp.workdps(dps):
            c = np.asarray(field['coefficients'])
            h = np.asarray(field['h'])
            if not (np.all(c[:2] == 0) and np.all(c[2] == c[2, 0])):
                raise ValueError('Reference requires exactly affine coefficients')
            g, n0, h0 = exact(c[2, 0]), exact(c[3, 0]), exact(h[0])
            if any(exact(c[3, i]) != n0 + g * (exact(h[i])-h0) for i in range(c.shape[1])):
                raise ValueError('Affine reduction is not algebraically exact')
            self.proof = {'exact_affine_reduction': True, 'segments': int(c.shape[1]), 'coefficient_conversion': 'binary64 exact integer ratio'}
            self.g, self.n0, self.h0 = g, n0, h0
            self.alpha = exact(coefficient)/(exact(frequency)*1000000)**2
            self.x0, self.p0 = mp.matrix([exact(x) for x in xyz]), mp.matrix([exact(x) for x in p])
            self.r0 = mp.norm(self.x0)
            self.er = self.x0/self.r0
            radial = mp.fdot(self.p0, self.er)
            if radial <= 0: raise ValueError('Local reference requires upward propagation')
            transverse = self.p0-radial*self.er
            self.et = transverse/mp.norm(transverse)
            self.L = self.r0*mp.norm(transverse)
            self.E = mp.fdot(self.p0, self.p0)+self.alpha*self.density(self.r0)
            self.nodes, self.weights = mp.gauss_quadrature(order, 'legendre')

    def density(self, radius):
        return self.n0+self.g*(radius-6371-self.h0)

    def local(self, radius):
        n2 = self.E-self.alpha*self.density(radius)
        q2 = n2-(self.L/radius)**2
        if n2 <= 0 or q2 <= 0: raise ValueError('Turning/cutoff outside local reference scope')
        return mp.sqrt(n2), mp.sqrt(q2)

    def integrals(self, lower, upper):
        half, middle = (upper-lower)/2, (upper+lower)/2
        sums = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        for node, weight in zip(self.nodes, self.weights):
            r = middle+half*node
            n,q = self.local(r)
            for i,v in enumerate((n/q, self.L/(r*r*q), 1/q)):
                sums[i] += weight*v
        return [half*x for x in sums]

    def query(self, values):
        with mp.workdps(self.dps):
            previous_s, previous_r, theta, group = mp.mpf(0), self.r0, mp.mpf(0), mp.mpf(0)
            rows = []
            for value in values:
                s = exact(value)
                if s < previous_s: raise ValueError('Queries must be nondecreasing')
                ds = s-previous_s
                n,q = self.local(previous_r)
                r = previous_r+ds*q/n
                for _ in range(12):
                    integral = self.integrals(previous_r, r)
                    n,q = self.local(r)
                    delta = (integral[0]-ds)*q/n
                    r -= delta
                    if abs(delta) < mp.power(10, -self.dps+15): break
                else: raise ArithmeticError('High precision radius inversion did not converge')
                integral = self.integrals(previous_r,r)
                theta += integral[1]; group += integral[2]
                er = mp.cos(theta)*self.er+mp.sin(theta)*self.et
                et = -mp.sin(theta)*self.er+mp.cos(theta)*self.et
                n,q = self.local(r)
                position = r*er; momentum = q*er+self.L/r*et
                rows.append(([float(v) for v in position],[float(v) for v in momentum],float(group),[mp.nstr(v, self.dps) for v in position]+[mp.nstr(group,self.dps)]))
                previous_s,previous_r = s,r
            return {'xyz':np.array([x[0] for x in rows]),'p':np.array([x[1] for x in rows]),'group':np.array([x[2] for x in rows]),'decimal': [x[3] for x in rows], 'g100': rows[-1][2]}


def self_check():
    field={'h':np.array([0.,200.]),'coefficients':np.zeros((4,1))}
    reference=AffineReference(field,[6471.,0.,0.],[.5,np.sqrt(.75),0.],80.61638494312905,10,80,16)
    result=reference.query([0.,10.,100.])
    expected=np.array([6471.,0.,0.])+np.array([0.,10.,100.])[:,None]*np.array([.5,np.sqrt(.75),0.])
    assert np.max(np.abs(result['xyz']-expected))<1e-11
    assert np.max(np.abs(result['group']-np.array([0.,10.,100.])))<1e-12
    with mp.workdps(80):
        assert exact(.1) != mp.mpf('0.1')
    return {'vacuum_straight':True,'exact_conversion':True}
