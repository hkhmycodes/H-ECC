import random
p=106839527430202782610735740077697524141755216002708025221579380756122101340723
u,v,w=5,1,7
G=(1,43885198696659819331868805857976446729580581201517994933321256281291666494433,52096118237682591038514212999531378342298973921021687132305677862559506645034)
cnt={'M':0,'S':0,'C':0}
def M(x,y): cnt['M']+=1; return x*y%p
def S(x): cnt['S']+=1; return x*x%p
def C(c,x): cnt['C']+=1; return c*x%p     # constant mult (u,v,w small)
def reset():
    for k in cnt: cnt[k]=0
def F(P):
    a,b,c=P
    return (u*a*(b*b-c*c)+v*b*(c*c-a*a)+w*c*(a*a-b*b))%p
def iota(P):
    a,b,c=P; return (b*c%p,c*a%p,a*b%p)
def eq(P,Q): return all((P[i]*Q[j]-P[j]*Q[i])%p==0 for i in range(3) for j in range(3))
# ---- reference
def grad(P):
    a,b,c=P
    return ((u*(b*b-c*c)+2*a*(w*c-v*b))%p,(v*(c*c-a*a)+2*b*(u*a-w*c))%p,(w*(a*a-b*b)+2*c*(v*b-u*a))%p)
def add_ref(P1,P2):
    g1=grad(P1);g2=grad(P2)
    s=sum(x*y for x,y in zip(P1,g2))%p; t=sum(x*y for x,y in zip(P2,g1))%p
    return iota(tuple((s*x-t*y)%p for x,y in zip(P1,P2)))
def dbl_ref(P):
    a,b,c=P; g=grad(P);B,Gm=g[1],g[2]; t=(0,Gm,(-B)%p)
    Ft=F(t); H=((u*a-w*c)*Gm*Gm-2*(v*c-w*b)*B*Gm+(v*b-u*a)*B*B)%p
    return iota(tuple((-Ft*x+H*y)%p for x,y in zip(P,t)))
# ---- optimised (v=1 so v-mults are free; u,w small constants)
def gradopt(P, need_alpha=True):
    a,b,c=P
    a2=S(a);b2=S(b);c2=S(c)
    ua=C(u,a); wc=C(w,c); vb=b   # v=1
    Fb=(c2-a2+M(2*b,(ua-wc)))%p                 # v*(c2-a2)+2b(ua-wc)
    Fg=(C(w,(a2-b2))+M(2*c,(vb-ua)))%p
    Fa=None
    if need_alpha: Fa=(C(u,(b2-c2))+M(2*a,(wc-vb)))%p
    return Fa,Fb,Fg,(ua,wc,vb)
def dbl_opt(P):
    a,b,c=P
    _,B,Gm,(ua,wc,vb)=gradopt(P,False)
    B2=S(B);G2=S(Gm);BG=M(B,Gm)
    Ft=M(BG,(B+C(w,Gm))%p)                       # BG(vB+wG)
    H=(M((ua-wc),G2)-M(2*(c-C(w,b)),BG)+M((vb-ua),B2))%p
    Rb=(M(-Ft,b)+M(H,Gm))%p; Rg=(M(-Ft,c)-M(H,B))%p; Ra=M(-Ft,a)
    return (M(Rb,Rg),M(Rg,Ra),M(Ra,Rb))
def add_opt(P1,P2,cache2=None):
    # cache2: precomputed gradient of P2 (fixed addend / table entry)
    Fa,Fb,Fg,_=gradopt(P1)
    g2=cache2 if cache2 else gradopt(P2)[:3]
    s=(M(P1[0],g2[0])+M(P1[1],g2[1])+M(P1[2],g2[2]))%p
    t=(M(P2[0],Fa)+M(P2[1],Fb)+M(P2[2],Fg))%p
    R=tuple((M(s,x)-M(t,y))%p for x,y in zip(P1,P2))
    return (M(R[1],R[2]),M(R[2],R[0]),M(R[0],R[1]))
if __name__=="__main__":
    P=G; ok=True
    for _ in range(5):
        Q=dbl_ref(P); reset(); Q2=dbl_opt(P); ok&=eq(Q,Q2); d=dict(cnt)
        R=add_ref(Q,P); reset(); R2=add_opt(Q,P); ok&=eq(R,R2); a_=dict(cnt)
        g2=gradopt(P)[:3]; reset(); R3=add_opt(Q,P,g2); ok&=eq(R,R3); ac=dict(cnt)
        P=R
    print("all agree:",ok); print("dbl_opt",d,"add_opt",a_,"add_opt cached addend",ac)
