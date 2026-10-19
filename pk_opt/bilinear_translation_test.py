import random, sys
sys.path.insert(0,'.')
from pk_opt_formulas import *
def nullspace_mod(rows,n):
    rows=[r[:] for r in rows]; piv=[]; r=0
    for c in range(n):
        pr=next((i for i in range(r,len(rows)) if rows[i][c]%p),None)
        if pr is None: continue
        rows[r],rows[pr]=rows[pr],rows[r]
        inv=pow(rows[r][c],-1,p); rows[r]=[x*inv%p for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]%p:
                f=rows[i][c]; rows[i]=[(x-f*y)%p for x,y in zip(rows[i],rows[r])]
        piv.append(c); r+=1
    free=[c for c in range(n) if c not in piv]; basis=[]
    for f in free:
        v=[0]*n; v[f]=1
        for i,c in enumerate(piv): v[c]=(-rows[i][f])%p
        basis.append(v)
    return basis
K=G                       # fixed addend (a table point)
# random points q on E and their translates q+K (D-centred group)
pts=[]; Q=dbl_ref(G); G2=dbl_ref(G)
while len(pts)<40:
    R=add_ref(Q,G2) if random.random()<.5 else dbl_ref(Q)
    if R==(0,0,0) or add_ref(R,K)==(0,0,0): Q=dbl_ref(G); continue
    Q=R; pts.append(Q)
rows=[]
for q in pts:
    q2=add_ref(q,K)
    rows.append([q[i]*q2[j]%p for i in range(3) for j in range(3)])
ns=nullspace_mod(rows,9)
print("dim of (1,1)-forms vanishing on graph of translation by K:",len(ns))
# build q' = r_i x r_j with r_i = M_i^T q
def mats(v): return [[v[3*i+j] for j in range(3)] for i in range(3)]
Mi,Mj=mats(ns[0]),mats(ns[1])
def rowvec(M,q): return [sum(q[i]*M[i][k] for i in range(3))%p for k in range(3)]
def cross(a,b): return ((a[1]*b[2]-a[2]*b[1])%p,(a[2]*b[0]-a[0]*b[2])%p,(a[0]*b[1]-a[1]*b[0])%p)
assert all(add_ref(q,K)!=(0,0,0) for q in pts)
ok=all(eq(cross(rowvec(Mi,q),rowvec(Mj,q)),add_ref(q,K)) for q in pts[:15])
print("translation by K == r_i(q) x r_j(q) (degree 2):",ok)
print("kernel basis entries small?", [max(min(x,p-x).bit_length() for x in v) for v in ns],"bits (max |coef|)")
