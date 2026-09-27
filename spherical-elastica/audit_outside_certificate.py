"""Audit the full dyadic cover and replay each exact reported inequality.

This is a certificate replay using the proof-producing interval primitives,
not a second independent implementation of the geometric theorem. Arithmetic
enclosures are additionally checked against exact Fraction endpoint formulas.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
from random import Random
import hashlib
import json
import time
import certify_outside_one_cycle as c

ROOT = Path(__file__).parent


def arithmetic_audit():
    rng = Random(72019)
    checked = 0
    def contains(interval, lower, upper):
        nonlocal checked
        assert Fraction(interval.lo,c.S) <= lower <= upper <= Fraction(interval.hi,c.S)
        checked += 1
    for _ in range(2000):
        aa = sorted(rng.randrange(-4*c.S,4*c.S) for _ in range(2))
        bb = sorted(rng.randrange(1,4*c.S) for _ in range(2))
        a,b = c.I(*aa),c.I(*bb)
        x,y = [Fraction(z,c.S) for z in aa],[Fraction(z,c.S) for z in bb]
        contains(a+b,x[0]+y[0],x[1]+y[1])
        contains(a-b,x[0]-y[1],x[1]-y[0])
        products = [v*w for v in x for w in y]
        contains(a*b,min(products),max(products))
        quotients = [v/w for v in x for w in y]
        contains(a/b,min(quotients),max(quotients))
        squares = [v*v for v in x]
        contains(a.sq(),0 if x[0]<=0<=x[1] else min(squares),max(squares))
        root = b.sqrt()
        assert root.lo**2 <= b.lo*c.S and root.hi**2 >= b.hi*c.S
        checked += 1
    return checked


def main():
    started = time.time()
    source = ROOT/'outside_one_cycle_certificate_refined.json'
    report = json.loads(source.read_text(encoding='utf-8'))
    assert report['status']=='complete' and not report['unresolved']
    assert report['scale']==c.S and report['pi']==[c.PI.lo,c.PI.hi]
    rows = sorted(report['leaves'],key=lambda z:z['path'])
    assert len(rows)==report['leaf_count']
    max_depth = max(len(row['path']) for row in rows)
    assert sum(1 << (max_depth-len(row['path'])) for row in rows)==1 << max_depth
    previous = None
    for row in rows:
        path = row['path']
        assert previous is None or not path.startswith(previous)
        previous = path
        bl,bh,rl,rh = 0,2*c.S,0,c.S
        for bit in path:
            assert bit in '01'
            if bh-bl >= 2*(rh-rl):
                mid = (bl+bh)//2
                if bit=='0': bh=mid
                else: bl=mid
            else:
                mid = (rl+rh)//2
                if bit=='0': rh=mid
                else: rl=mid
        assert row['b']==[bl,bh] and row['ratio']==[rl,rh]
    print(json.dumps(dict(coverage='complete',leaf_count=len(rows),max_depth=max_depth)),flush=True)
    arithmetic_checks = arithmetic_audit()
    margins = {}
    counts = Counter()
    for index,row in enumerate(rows,1):
        b,ratio = c.I(*row['b']),c.I(*row['ratio'])
        _,pressure,_,_,_ = c.setup(b,ratio)
        reason = row['reason']
        if reason=='small_pressure':
            assert pressure.hi <= c.EIGHTH.lo
            margin = 0
        elif reason=='large_pressure':
            assert pressure.lo >= c.THREE_HALVES.hi
            margin = 0
        elif reason=='elementary_energy':
            lower = (c.PI*c.I(b.lo).sq()/(2*(1+c.I(b.lo).sq()/2).sqrt())).lo
            margin = lower-c.cap_half_upper(pressure).hi
        elif reason=='integral_energy':
            margin = c.energy_lower(b,ratio,row['cells'])-c.cap_half_upper(pressure).hi
        elif reason=='closure_angle':
            margin = c.angle_lower(b,ratio,row['cells'])-2*c.PI.hi
        else:
            raise AssertionError('Unknown certificate rule')
        assert margin == row['margin']
        if reason not in ('small_pressure','large_pressure'):
            assert margin>0
            margins[reason] = min(margins.get(reason,margin),margin)
        counts[reason] += 1
        if index%1000==0:
            print(json.dumps(dict(replayed=index,total=len(rows),seconds=round(time.time()-started,1))),flush=True)
    hashes = {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
              for name in ('certify_outside_one_cycle.py','refine_outside_certificate.py',source.name)}
    audit = dict(status='passed',coverage='complete prefix-free dyadic cover',
                 leaf_count=len(rows),max_depth=max_depth,counts=dict(counts),
                 arithmetic_checks=arithmetic_checks,min_margins_fixed_point=margins,
                 min_margins_display={k:v/c.S for k,v in margins.items()},
                 exact_scale=c.S,seconds=time.time()-started,sha256=hashes)
    (ROOT/'outside_one_cycle_certificate_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':
    main()
