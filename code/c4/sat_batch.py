from fractions import Fraction
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
import random, time, sys

def zhou_cols(n):
    r=2*n+1
    def A(i,j):
        s=-1 if (j+n)%2 else 1
        return int(Fraction(1-s*r,4)*i+Fraction(2-3*abs(i),2)*j)
    return [(A(-1,j), j, A(1,j)) for j in range(1,n+1)]

def deltas(rad):
    return [(a,b,c) for a in range(-rad,rad+1) for b in range(-rad,rad+1)
            for c in range(-rad,rad+1) if a+b+c==0]

def zhou_partition_sat(n, rad=2):
    cols=zhou_cols(n); DS=deltas(rad)
    pool=IDPool()
    options=[]   # (j, outs)
    by_col={j:[] for j in range(n)}
    by_val={x:[] for x in range(1,3*n+1)}
    for j,v in enumerate(cols):
        seen=set()
        for d in DS:
            outs=tuple(sorted(abs(3*v[k]+d[k]) for k in range(3)))
            raw=tuple(abs(3*v[k]+d[k]) for k in range(3))
            if any(o==0 or o>3*n for o in raw) or len(set(raw))<3: continue
            key=(j,raw)
            if key in seen: continue
            seen.add(key)
            oi=len(options); options.append((j,raw,d))
            var=pool.id(('o',oi))
            by_col[j].append(var)
            for x in raw: by_val[x].append(var)
    cnf=[]
    for j in range(n):
        cnf += CardEnc.equals(lits=by_col[j], bound=1, vpool=pool, encoding=EncType.pairwise).clauses
    for x in range(1,3*n+1):
        if not by_val[x]: return None
        cnf += CardEnc.equals(lits=by_val[x], bound=1, vpool=pool, encoding=EncType.pairwise).clauses
    with Cadical153(bootstrap_with=cnf) as s:
        if not s.solve(): return False
        model=set(l for l in s.get_model() if l>0)
    T=[None]*n
    for oi,(j,raw,d) in enumerate(options):
        if pool.id(('o',oi)) in model:
            T[j]=tuple(sorted(raw))
    return T

def signable(vals):
    S=sum(vals)
    if S%2: return False
    bits=1
    for v in vals: bits|=bits<<v
    return (bits>>(S//2))&1==1

def find_certificate(group, rng, tries=12000):
    m=len(group)
    MU=[(t[2],t[0],t[1]) for t in group]
    W =[(t[1]-t[0], t[1]+t[2], t[0]+t[2]) for t in group]
    for _ in range(tries):
        modes=[rng.randrange(3) for _ in range(m)]
        mu=[MU[i][modes[i]] for i in range(m)]
        w =[W[i][modes[i]] for i in range(m)]
        if sum(mu)%2:
            for i in rng.sample(range(m),m):
                hit=False
                for nm in (0,1,2):
                    if nm!=modes[i] and (MU[i][nm]-mu[i])%2:
                        modes[i]=nm; mu[i]=MU[i][nm]; w[i]=W[i][nm]; hit=True; break
                if hit: break
            if sum(mu)%2: continue
        if signable(mu) and signable(w): return modes
    return None

if __name__=="__main__":
    rng=random.Random(7)
    for m in [19,23,31,35,43,47,55,59,67,71,79,83]:
        n=3*m; t0=time.time()
        T=zhou_partition_sat(n, rad=2); note="|d|<=2"
        if T is False:
            T=zhou_partition_sat(n, rad=3); note="|d|<=3 (rad2 UNSAT!)"
        if not isinstance(T,list):
            print(f"m={m} n={n}: status={T} ({note})",flush=True); continue
        ok=sorted(x for t in T for x in t)==list(range(1,3*n+1))
        res=[find_certificate([T[j] for j in range(n) if j%3==r], rng) is not None for r in range(3)]
        print(f"m={m} n={n}: partition={ok} [{note}]  j-mod-3 certs={res}  ALL={all(res)}  ({time.time()-t0:.1f}s)",flush=True)
