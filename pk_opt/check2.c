#define main bench_main
#include "bench.c"
#undef main
#include "r2.h"
#include "aff.h"
static void toM(fe*r,const char*dec){u64 l[4]={0};for(const char*c=dec;*c;c++){u128 carry=*c-'0';for(int i=0;i<4;i++){u128 x=(u128)l[i]*10+carry;l[i]=(u64)x;carry=x>>64;}}
 fe a;memcpy(a.l,l,32);fe r2;memcpy(r2.l,R2L,32);fmul(r,&a,&r2);}
static void pr(const char*n,fe x){fe one={{1,0,0,0}},o;fmul(&o,&x,&one);printf("%s=%016llx%016llx%016llx%016llx\n",n,(unsigned long long)o.l[3],(unsigned long long)o.l[2],(unsigned long long)o.l[1],(unsigned long long)o.l[0]);}
int main(){const char*ps="106839527430202782610735740077697524141755216002708025221579380756122101340723";
 u64 l[4]={0};for(const char*c=ps;*c;c++){u128 carry=*c-'0';for(int i=0;i<4;i++){u128 x=(u128)l[i]*10+carry;l[i]=(u64)x;carry=x>>64;}}
 memcpy(PM,l,32);u64 inv=1;for(int i=0;i<7;i++)inv*=2-PM[0]*inv;NINV=-inv;
 pt G;toM(&G.a,"1");toM(&G.b,"43885198696659819331868805857976446729580581201517994933321256281291666494433");toM(&G.c,"52096118237682591038514212999531378342298973921021687132305677862559506645034");
 pt Ga;toM(&Ga.a,GA0);toM(&Ga.b,GA1);toM(&Ga.c,"1");
 pt D=pk_dbl(G);
 pt A1=pk_add_affine(D,Ga);pr("u0",A1.a);pr("u1",A1.b);pr("u2",A1.c);
 grad3 g=pk_grad(Ga);pt A2=pk_add_cached(D,Ga,g,1);pr("m0",A2.a);pr("m1",A2.b);pr("m2",A2.c);
 return 0;}
