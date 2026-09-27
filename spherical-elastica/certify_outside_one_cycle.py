"""Exact fixed-point interval attempt for the remaining exterior-drop branch.

All reported bounds use integer arithmetic, outward rounding, integer square
roots, and a rational Machin enclosure of pi. Floating point is used only to
choose a trial spherical cap; any such nonnegative rational choice is valid.
An incomplete run is NOT a proof. A complete report retains every terminal
box, which must cover [0,2] x [0,1] without gaps.
"""
from pathlib import Path
from fractions import Fraction
from math import isqrt
import json
import time

ROOT = Path(__file__).parent
S = 1 << 44


def ceildiv(a, b):
    return -((-a)//b)


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo, self.hi = lo, lo if hi is None else hi
        if self.lo > self.hi:
            raise ValueError('Inverted interval')

    @classmethod
    def rat(cls, n, d=1):
        return cls(n*S//d, ceildiv(n*S, d))

    def __add__(self, y):
        if not isinstance(y, I):
            y = I.rat(y)
        return I(self.lo+y.lo, self.hi+y.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, y):
        if not isinstance(y, I):
            y = I.rat(y)
        return self+-y

    def __rsub__(self, y):
        return -self+y

    def __mul__(self, y):
        if not isinstance(y, I):
            y = I.rat(y)
        products = (self.lo*y.lo, self.lo*y.hi, self.hi*y.lo, self.hi*y.hi)
        return I(min(products)//S, ceildiv(max(products), S))

    __rmul__ = __mul__

    def __truediv__(self, y):
        if not isinstance(y, I):
            y = I.rat(y)
        if y.lo <= 0:
            raise ValueError('Division requires a strictly positive denominator')
        pairs = ((a*S, b) for a in (self.lo, self.hi) for b in (y.lo, y.hi))
        bounds = [(a//b, ceildiv(a,b)) for a,b in pairs]
        return I(min(x[0] for x in bounds), max(x[1] for x in bounds))

    def sq(self):
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo*self.lo, self.hi*self.hi)
        upper = max(self.lo*self.lo, self.hi*self.hi)
        return I(lower//S, ceildiv(upper,S))

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('Square root of interval with negative lower endpoint')
        lo, hi = isqrt(self.lo*S), isqrt(self.hi*S)
        return I(lo, hi+(hi*hi < self.hi*S))

    def nonnegative(self):
        return I(max(0,self.lo), max(0,self.hi))


def arctan_bounds(denominator, terms=24):
    total = sum((Fraction((-1)**n, (2*n+1)*denominator**(2*n+1)) for n in range(terms)), Fraction())
    following = Fraction((-1)**terms, (2*terms+1)*denominator**(2*terms+1))
    return min(total,total+following), max(total,total+following)


def pi_bounds():
    a,b = arctan_bounds(5)
    c,d = arctan_bounds(239)
    lo,hi = 16*a-4*d, 16*b-4*c
    return I(lo.numerator*S//lo.denominator, ceildiv(hi.numerator*S,hi.denominator))


PI = pi_bounds()
HALF, EIGHTH = I.rat(1,2), I.rat(1,8)
THREE_HALVES = I.rat(3,2)


def cap_half_upper(pressure):
    # This choice affects efficiency only: a cap at any rational k >= 0 is
    # an admissible test in the minimization defining c_lambda.
    target = pressure.hi/S
    trial = min(target, (2*target)**(1/3)) if target else 0.
    for _ in range(12):
        trial -= (trial+.5*trial**3-target)/(1+1.5*trial**2)
    k = I(max(0, round(trial*S)))
    return PI*(pressure+(k.sq()/2-pressure*k)/(1+k.sq()).sqrt())


def setup(b, ratio):
    x = b*ratio
    pressure = b*(1-ratio)*(b.sq()*(1+ratio.sq())+4)/8
    negative_cbar = x.sq()*x/4+x+2*pressure
    positive_cbar = ratio*negative_cbar
    constant = x*negative_cbar
    j = (constant+pressure.sq()).sqrt()
    return x, pressure, positive_cbar, negative_cbar, j


def energy_lower(b, ratio, cells=48):
    total = 0
    for n in range(1,cells):
        t = I.rat(n,cells)+I(0,ceildiv(S,cells))
        s = 2*t/(1+t.sq())
        v = s.sq()
        qplus = 4+b.sq()*(v.sq()+(1-ratio)*v+1-ratio+ratio.sq())
        qminus = 4+b.sq()*(ratio.sq()*v.sq()-ratio*(1-ratio)*v+1-ratio+ratio.sq())
        common = 8*b.sq()*s.sq().sq()*s/(1+t.sq())
        plus = common/((ratio+v).sqrt()*qplus.sqrt())
        minus = common*ratio.sq()*ratio.sqrt()/((1+ratio*v).sqrt()*qminus.sqrt())
        total += (2*plus+minus).lo
    return total//cells


def root_bounds(amplitude, pressure, cbar, u2, sign, steps=15):
    cube = amplitude.sq()*amplitude

    def polynomial(v):
        z = I(v)
        return cube*z.sq().sq()+4*u2*(amplitude*z.sq()-sign*2*pressure*z-cbar)

    lo,hi = 0,S
    for _ in range(steps):
        mid = (lo+hi)//2
        if polynomial(mid).hi < 0:
            lo = mid
        else:
            hi = mid
    root_lo = lo
    lo,hi = 0,S
    for _ in range(steps):
        mid = (lo+hi)//2
        if polynomial(mid).lo > 0:
            hi = mid
        else:
            lo = mid
    return I(root_lo,hi)


def angle_lower(b, ratio, cells=48):
    x,pressure,cp,cn,j = setup(b,ratio)
    total = 0
    for n in range(1,cells):
        t = I.rat(n,cells)+I(0,ceildiv(S,cells))
        u2 = (2*t/(1+t.sq())).sq()
        for amp,cbar,sign,multiplicity in ((b,cp,1,2),(x,cn,-1,1)):
            v = root_bounds(amp,pressure,cbar,u2,sign)
            qbar = 2*cbar-amp*v.sq()+sign*3*pressure*v
            if qbar.hi <= 0:
                raise ValueError('Invalid positive Q upper bound')
            lower = I(4*j.lo)*I(v.lo)/(I(qbar.hi)*(1+t.sq()))
            total += multiplicity*lower.lo
    return total//cells


def certify_box(b,ratio,cells=48):
    _,pressure,_,_,_ = setup(b,ratio)
    if pressure.hi <= EIGHTH.lo:
        return 'small_pressure', 0
    if pressure.lo >= THREE_HALVES.hi:
        return 'large_pressure', 0
    cap = cap_half_upper(pressure).hi
    lower = (PI*I(b.lo).sq()/(2*(1+I(b.lo).sq()/2).sqrt())).lo
    if lower > cap:
        return 'elementary_energy', lower-cap
    try:
        lower = energy_lower(b,ratio,cells)
        if lower > cap:
            return 'integral_energy', lower-cap
        lower = angle_lower(b,ratio,cells)
        if lower > 2*PI.hi:
            return 'closure_angle', lower-2*PI.hi
    except ValueError:
        pass
    return None,0


def main():
    started = time.time()
    stack = [(0,2*S,0,S,0,'')]
    leaves,unresolved = [],[]
    visited = 0
    while stack:
        bl,bh,rl,rh,depth,path = stack.pop()
        visited += 1
        reason,margin = certify_box(I(bl,bh),I(rl,rh))
        if reason:
            leaves.append(dict(path=path,b=[bl,bh],ratio=[rl,rh],reason=reason,margin=margin))
        elif depth >= 18:
            unresolved.append(dict(path=path,b=[bl,bh],ratio=[rl,rh]))
        elif bh-bl >= 2*(rh-rl):
            mid = (bl+bh)//2
            stack.extend(((mid,bh,rl,rh,depth+1,path+'1'),(bl,mid,rl,rh,depth+1,path+'0')))
        else:
            mid = (rl+rh)//2
            stack.extend(((bl,bh,mid,rh,depth+1,path+'1'),(bl,bh,rl,mid,depth+1,path+'0')))
        if visited%100 == 0:
            print(json.dumps(dict(visited=visited,certified=len(leaves),unresolved=len(unresolved),pending=len(stack),seconds=round(time.time()-started,1))), flush=True)
    counts = {}
    for row in leaves:
        counts[row['reason']] = counts.get(row['reason'],0)+1
    report = dict(status='complete' if not unresolved else 'incomplete',
                  scale=S,pi=[PI.lo,PI.hi],cells=48,root_steps=15,
                  initial_rectangle=[0,2*S,0,S],visited=visited,
                  leaf_count=len(leaves),counts=counts,
                  seconds=time.time()-started,leaves=leaves,unresolved=unresolved)
    (ROOT/'outside_one_cycle_certificate.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('leaves','unresolved')},indent=2))


if __name__ == '__main__':
    main()
