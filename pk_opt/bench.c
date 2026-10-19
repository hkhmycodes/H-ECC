// Same 4x64 Montgomery field for pK and Jacobian (S implemented as M: bias in pK's favour).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
typedef unsigned __int128 u128; typedef uint64_t u64;
static const u64 P[4]={0}; 
static u64 PM[4], NINV, R2[4];
typedef struct{u64 l[4];}fe;
static inline int geq(const u64*a,const u64*b){for(int i=3;i>=0;i--){if(a[i]>b[i])return 1;if(a[i]<b[i])return 0;}return 1;}
static inline void fadd(fe*r,const fe*a,const fe*b){u64 c=0,t[4];for(int i=0;i<4;i++){u128 s=(u128)a->l[i]+b->l[i]+c;t[i]=(u64)s;c=s>>64;}
 u64 bw=0,d[4];for(int i=0;i<4;i++){u128 s=(u128)t[i]-PM[i]-bw;d[i]=(u64)s;bw=(s>>64)&1;}
 if(c||!bw)memcpy(r->l,d,32);else memcpy(r->l,t,32);}
static inline void fsub(fe*r,const fe*a,const fe*b){u64 bw=0,t[4];for(int i=0;i<4;i++){u128 s=(u128)a->l[i]-b->l[i]-bw;t[i]=(u64)s;bw=(s>>64)&1;}
 if(bw){u64 c=0;for(int i=0;i<4;i++){u128 s=(u128)t[i]+PM[i]+c;t[i]=(u64)s;c=s>>64;}}memcpy(r->l,t,32);}
static inline void fmul(fe*r,const fe*a,const fe*b){u64 t[6]={0};
 for(int i=0;i<4;i++){u64 c=0;for(int j=0;j<4;j++){u128 s=(u128)a->l[j]*b->l[i]+t[j]+c;t[j]=(u64)s;c=s>>64;}
  u128 s=(u128)t[4]+c;t[4]=(u64)s;t[5]=(u64)(s>>64);
  u64 m=t[0]*NINV;c=0;{u128 s2=(u128)m*PM[0]+t[0];c=s2>>64;}
  for(int j=1;j<4;j++){u128 s2=(u128)m*PM[j]+t[j]+c;t[j-1]=(u64)s2;c=s2>>64;}
  s=(u128)t[4]+c;t[3]=(u64)s;t[4]=t[5]+(u64)(s>>64);}
 u64 bw=0,d[4];for(int i=0;i<4;i++){u128 s=(u128)t[i]-PM[i]-bw;d[i]=(u64)s;bw=(s>>64)&1;}
 if(t[4]||!bw)memcpy(r->l,d,32);else memcpy(r->l,t,32);}
#define M(r,a,b) fmul(&r,&a,&b)
#define S(r,a) fmul(&r,&a,&a)
static inline fe dbl2(fe a){fe r;fadd(&r,&a,&a);return r;}
static inline fe mulc(fe x,int c){ // small constant by addition chain
 fe x2=dbl2(x),x4=dbl2(x2),x8=dbl2(x4),r;
 switch(c){case 5:fadd(&r,&x4,&x);return r;case 7:fsub(&r,&x8,&x);return r;case 3:fadd(&r,&x2,&x);return r;default:return x;}}
// ---- pK (u,v,w)=(5,1,7), doubling 15M+5S+5C ; addition (cached/mixed addend)
typedef struct{fe a,b,c;}pt;
static inline pt pk_dbl(pt p){
 fe a2,b2,c2;S(a2,p.a);S(b2,p.b);S(c2,p.c);
 fe ua=mulc(p.a,5),wc=mulc(p.c,7),vb=p.b,t1,t2,B,G,tmp;
 fsub(&t1,&ua,&wc); fe b2x=dbl2(p.b); M(tmp,b2x,t1); fsub(&B,&c2,&a2); fadd(&B,&B,&tmp);
 fsub(&t2,&vb,&ua); fe c2x=dbl2(p.c); M(tmp,c2x,t2); fsub(&G,&a2,&b2); G=mulc(G,7); fadd(&G,&G,&tmp);
 fe B2,G2,BG;S(B2,B);S(G2,G);M(BG,B,G);
 fe wG=mulc(G,7),s=B; fadd(&s,&s,&wG); fe Ft;M(Ft,BG,s);
 fe q=mulc(p.b,7); fsub(&q,&p.c,&q); q=dbl2(q);
 fe h1,h2,h3,H;M(h1,t1,G2);M(h2,q,BG);M(h3,t2,B2);fsub(&H,&h1,&h2);fadd(&H,&H,&h3);
 fe nF;fe z={{0}};fsub(&nF,&z,&Ft);
 fe Ra,Rb,Rg,x1,x2;M(Ra,nF,p.a);M(x1,nF,p.b);M(x2,H,G);fadd(&Rb,&x1,&x2);M(x1,nF,p.c);M(x2,H,B);fsub(&Rg,&x1,&x2);
 pt r;M(r.a,Rb,Rg);M(r.b,Rg,Ra);M(r.c,Ra,Rb);return r;}
typedef struct{fe a,b,c;}grad3;
static inline grad3 pk_grad(pt p){ // full gradient (3S+3M+4C)
 fe a2,b2,c2;S(a2,p.a);S(b2,p.b);S(c2,p.c);
 fe ua=mulc(p.a,5),wc=mulc(p.c,7),vb=p.b,t,x;grad3 g;
 fsub(&t,&wc,&vb);fe a2x=dbl2(p.a);M(x,a2x,t);fsub(&g.a,&b2,&c2);g.a=mulc(g.a,5);fadd(&g.a,&g.a,&x);
 fsub(&t,&ua,&wc);fe b2x=dbl2(p.b);M(x,b2x,t);fsub(&g.b,&c2,&a2);fadd(&g.b,&g.b,&x);
 fsub(&t,&vb,&ua);fe c2x=dbl2(p.c);M(x,c2x,t);fsub(&g.c,&a2,&b2);g.c=mulc(g.c,7);fadd(&g.c,&g.c,&x);return g;}
static inline pt pk_add_cached(pt p1,pt p2,grad3 g2,int mixed){ // g2 precomputed; mixed: p2.c==1
 grad3 g1=pk_grad(p1);fe s,t,x;
 M(s,p1.a,g2.a);M(x,p1.b,g2.b);fadd(&s,&s,&x);M(x,p1.c,g2.c);fadd(&s,&s,&x);
 M(t,p2.a,g1.a);M(x,p2.b,g1.b);fadd(&t,&t,&x);
 if(mixed)fadd(&t,&t,&g1.c);else{M(x,p2.c,g1.c);fadd(&t,&t,&x);}
 fe Ra,Rb,Rc,y;M(Ra,s,p1.a);M(y,t,p2.a);fsub(&Ra,&Ra,&y);M(Rb,s,p1.b);M(y,t,p2.b);fsub(&Rb,&Rb,&y);
 M(Rc,s,p1.c);if(mixed)fsub(&Rc,&Rc,&t);else{M(y,t,p2.c);fsub(&Rc,&Rc,&y);}
 pt r;M(r.a,Rb,Rc);M(r.b,Rc,Ra);M(r.c,Ra,Rb);return r;}

static inline pt pk_add_affine(pt p1,pt p2){ // addend affine (c=1), gradient NOT cached: 18M+5S
 grad3 g1=pk_grad(p1); fe a2,b2,x,t,u1;S(a2,p2.a);S(b2,p2.b); fe one=p2.c; grad3 g2;
 fsub(&t,&b2,&one);g2.a=mulc(t,5);fe w7=mulc(one,7);fsub(&t,&w7,&p2.b);fe a2x=dbl2(p2.a);M(x,a2x,t);fadd(&g2.a,&g2.a,&x);
 fsub(&t,&one,&a2);fe ua=mulc(p2.a,5);fsub(&u1,&ua,&w7);fe b2x=dbl2(p2.b);M(x,b2x,u1);fadd(&g2.b,&t,&x);
 fsub(&t,&a2,&b2);g2.c=mulc(t,7);fsub(&t,&p2.b,&ua);t=dbl2(t);fadd(&g2.c,&g2.c,&t);
 fe s,tt;M(s,p1.a,g2.a);M(x,p1.b,g2.b);fadd(&s,&s,&x);M(x,p1.c,g2.c);fadd(&s,&s,&x);
 M(tt,p2.a,g1.a);M(x,p2.b,g1.b);fadd(&tt,&tt,&x);fadd(&tt,&tt,&g1.c);
 fe Ra,Rb,Rc,y;M(Ra,s,p1.a);M(y,tt,p2.a);fsub(&Ra,&Ra,&y);M(Rb,s,p1.b);M(y,tt,p2.b);fsub(&Rb,&Rb,&y);
 M(Rc,s,p1.c);fsub(&Rc,&Rc,&tt);pt r;M(r.a,Rb,Rc);M(r.b,Rc,Ra);M(r.c,Ra,Rb);return r;}
// ---- Jacobian y^2=x^3-x (a=-1): the strongest simple baseline
typedef struct{fe X,Y,Z;}jp;
static inline jp j_dbl(jp p){fe XX,YY,YYYY,ZZ,ZZZZ,t,Sv,Mv,T;jp r;S(XX,p.X);S(YY,p.Y);S(YYYY,YY);S(ZZ,p.Z);
 fadd(&t,&p.X,&YY);S(t,t);fsub(&t,&t,&XX);fsub(&t,&t,&YYYY);Sv=dbl2(t);
 S(ZZZZ,ZZ);fe XX3=mulc(XX,3);fsub(&Mv,&XX3,&ZZZZ);S(T,Mv);fe S2=dbl2(Sv);fsub(&T,&T,&S2);
 r.X=T;fsub(&t,&Sv,&T);M(t,Mv,t);fe y8=dbl2(dbl2(dbl2(YYYY)));fsub(&r.Y,&t,&y8);
 fadd(&t,&p.Y,&p.Z);S(t,t);fsub(&t,&t,&YY);fsub(&r.Z,&t,&ZZ);return r;}
static inline jp j_add(jp p,jp q){fe Z1Z1,Z2Z2,U1,U2,S1,S2,H,I,J,rr,V,t;S(Z1Z1,p.Z);S(Z2Z2,q.Z);M(U1,p.X,Z2Z2);M(U2,q.X,Z1Z1);
 M(t,q.Z,Z2Z2);M(S1,p.Y,t);M(t,p.Z,Z1Z1);M(S2,q.Y,t);fsub(&H,&U2,&U1);fe H2=dbl2(H);S(I,H2);M(J,H,I);fsub(&rr,&S2,&S1);rr=dbl2(rr);M(V,U1,I);
 jp r;S(r.X,rr);fsub(&r.X,&r.X,&J);fe V2=dbl2(V);fsub(&r.X,&r.X,&V2);fsub(&t,&V,&r.X);M(t,rr,t);fe sj;M(sj,S1,J);sj=dbl2(sj);fsub(&r.Y,&t,&sj);
 fadd(&t,&p.Z,&q.Z);S(t,t);fsub(&t,&t,&Z1Z1);fsub(&t,&t,&Z2Z2);M(r.Z,t,H);return r;}
static inline jp j_madd(jp p,jp q){fe Z1Z1,U2,S2,H,HH,I,J,rr,V,t;S(Z1Z1,p.Z);M(U2,q.X,Z1Z1);M(t,p.Z,Z1Z1);M(S2,q.Y,t);fsub(&H,&U2,&p.X);S(HH,H);I=dbl2(dbl2(HH));M(J,H,I);
 fsub(&rr,&S2,&p.Y);rr=dbl2(rr);M(V,p.X,I);jp r;S(r.X,rr);fsub(&r.X,&r.X,&J);fe V2=dbl2(V);fsub(&r.X,&r.X,&V2);fsub(&t,&V,&r.X);M(t,rr,t);fe yj;M(yj,p.Y,J);yj=dbl2(yj);fsub(&r.Y,&t,&yj);
 fadd(&t,&p.Z,&H);S(t,t);fsub(&t,&t,&Z1Z1);fsub(&r.Z,&t,&HH);return r;}
static double now(){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec*1e9+ts.tv_nsec;}
static fe rnd(){fe r;for(int i=0;i<4;i++)r.l[i]=((u64)rand()<<33)^((u64)rand()<<11)^rand();r.l[3]%=PM[3];return r;}
int main(int argc,char**argv){
 // p from the paper (decimal) -> limbs via simple big-int parse
 const char*ps="106839527430202782610735740077697524141755216002708025221579380756122101340723";
 u64 l[4]={0};for(const char*c=ps;*c;c++){u128 carry=*c-'0';for(int i=0;i<4;i++){u128 x=(u128)l[i]*10+carry;l[i]=(u64)x;carry=x>>64;}}
 memcpy(PM,l,32);u64 inv=1;for(int i=0;i<7;i++)inv*=2-PM[0]*inv;NINV=-inv;
 int N=2000000;fe z={{0}};
 pt X={rnd(),rnd(),rnd()};pt Y={rnd(),rnd(),rnd()};grad3 gy=pk_grad(Y);
 jp A={rnd(),rnd(),rnd()},B={rnd(),rnd(),rnd()};
 double t0,t1;volatile u64 sink=0;
 t0=now();for(int i=0;i<N;i++)X=pk_dbl(X);t1=now();double pdbl=(t1-t0)/N;sink+=X.a.l[0];
 t0=now();for(int i=0;i<N;i++)X=pk_add_cached(X,Y,gy,0);t1=now();double padd=(t1-t0)/N;sink+=X.a.l[0];
 Y.c=(fe){{1,0,0,0}};gy=pk_grad(Y);
 t0=now();for(int i=0;i<N;i++)X=pk_add_cached(X,Y,gy,1);t1=now();double pmadd=(t1-t0)/N;sink+=X.a.l[0];
 t0=now();for(int i=0;i<N;i++)X=pk_add_affine(X,Y);t1=now();double pmaff=(t1-t0)/N;sink+=X.a.l[0];
 t0=now();for(int i=0;i<N;i++)A=j_dbl(A);t1=now();double jdbl=(t1-t0)/N;sink+=A.X.l[0];
 t0=now();for(int i=0;i<N;i++)A=j_add(A,B);t1=now();double jadd=(t1-t0)/N;sink+=A.X.l[0];
 t0=now();for(int i=0;i<N;i++)A=j_madd(A,B);t1=now();double jmadd=(t1-t0)/N;sink+=A.X.l[0];
 t0=now();fe m=rnd(),m2=rnd();for(int i=0;i<N*5;i++){fmul(&m,&m,&m2);}t1=now();double fm=(t1-t0)/(N*5);sink+=m.l[0];
 printf("field mul: %.1f ns\n",fm);
 printf("pK  dbl (15M+5S+5C):          %.0f ns  (%.1f M-eq)\n",pdbl,pdbl/fm);
 printf("pK  add, cached addend:       %.0f ns  (%.1f M-eq)\n",padd,padd/fm);
 printf("pK  add, cached+affine addend:%.0f ns  (%.1f M-eq)\n",pmadd,pmadd/fm);
 printf("pK  add, affine uncached:     %.0f ns  (%.1f M-eq)\n",pmaff,pmaff/fm);
 printf("Jac dbl  (a=-1, 1M+8S):       %.0f ns  (%.1f M-eq)\n",jdbl,jdbl/fm);
 printf("Jac add  (11M+5S):            %.0f ns  (%.1f M-eq)\n",jadd,jadd/fm);
 printf("Jac madd (7M+4S):             %.0f ns  (%.1f M-eq)\n",jmadd,jmadd/fm);
 // scalar-mult model: 254 dbl + 51 adds (w=5 fixed window, affine table)
 double pk=254*pdbl+51*pmadd, jc=254*jdbl+51*jmadd;
 printf("model 254-bit scalar mult (254 dbl + 51 table adds): pK %.1f us, Jacobian %.1f us, ratio %.2f\n",pk/1e3,jc/1e3,pk/jc);
 return (int)(sink&1);}
