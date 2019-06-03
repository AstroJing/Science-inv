#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<iostream>
#include<string>
#include<fstream>
//using namespace std
#define pi 3.14159265358979
#define G 6.673e-11
#define c 2.99792458e8
#define M 1.989e30
#define Cri 1e-7
#define Gc2 7.42471382405e-28
#define Gc4 8.26110825251e-45
#define Mscale 2.033931261665867e5

static double a[3000],b[3000],l=2,Ga[3000],Gb[3000];
static int j,scr,loopTOV=1;
char dysm[7]="/-\\/-\\";
FILE *inf,*outf;
double temp_rhoc,temp_M,temp_R,temp_freq,temp_dpt;
double temp_r,temp_m,temp_p,temp_E;
double *mfile,*pfile,*rhofile,*Bfile,*Bcor,*Wfile,*Wfile1,*Wfile2,*Vfile,*Vfile1,*Vfile2;
struct EoS
	{
		std::string Name;
		std::string tName;
	};
struct EoS list[50001];

double che(double e)
//linear interpolation
{
double le;
int i=0;
le=log10(e);
i=j-1;
for(;i>0;i--)
{
if(le>a[i-1])
  {
	return(pow(10, (b[i - 1] * (le - a[i]) - b[i] * (le - a[i - 1])) / (a[i - 1] - a[i])));
  }
}
return(0);
}

double ch(double p)
{
double lp;
double e,le;
int i=0;
lp=log10(p);
i=j-1;
for(;i>0;i--)
{
if(lp>b[i-1])
  {
  le=a[i-1]*(lp-b[i])/(b[i-1]-b[i])+a[i]*(lp-b[i-1])/(b[i]-b[i-1]);
  e=pow(10,le);
  return(e);
  }
}
return(0);
}

double fp(double r,double p,double e,double m)
{
	return(-G*(e+p/c/c)*(m+4*pi*r*r*r*p/c/c)/(r*r-2*r*m*Gc2));
}

double fm(double r,double e)
{
return(4*pi*r*r*e);
}

double Bf(double r,double p,double m,double B)
{
	return(2*Gc2/r/r*(m+4*pi*r*r*r*p/c/c)*B/(1-2*Gc2*m/r));
}

double DH1(double r,double m,double A,double p,double e,double H1,double H0,double K,double V)
{
	return(-1/r*(l+1+2*m*A/r+4*pi*r*r*A*(p-e))*H1+1/r*A*(H0+K-16*pi*(e+p)*V));
}

double DK(double r,double H0,double H1,double Dv,double K,double e,double p,double A,double W)
{
	return(1/r*H0+0.5*l*(l+1)/r*H1-((l+1)/r-0.5*Dv)*K-8*pi*(e+p)*sqrt(A)/r*W);
}

double DW(double r,double W,double A,double gamma,double p,double B,double X,double V,double H0,double K)
{
	return(-(l+1)/r*W+r*sqrt(A)*(1/gamma/p/sqrt(B)*X-l*(l+1)/r/r*V+0.5*H0+K));
}

double DX(double r,double X,double e,double p,double B,double Dv,double H0,double w,double H1,double K,double V,double A,double F,double W)
{
	return(-l/r*X+(e+p)*sqrt(B)*(0.5*(1/r-0.5*Dv)*H0+0.5*(r*w*w/B+0.5*l*(l+1)/r)*H1+0.5*(1.5*Dv-1/r)*K-0.5*l*(l+1)*Dv/r/r*V-1/r*(4*pi*(e+p)*sqrt(A)+w*w*sqrt(A)/B-0.5*r*r*F)*W));
}

double H0f(double r,double B,double X,double m,double p,double A,double H1,double K,double w)
{
double H0=8*pi*r*r*r/sqrt(B)*X-(0.5*l*(l+1)*(m+4*pi*r*r*r*p)-w*w*r*r*r/A/B)*H1+(0.5*(l+2)*(l-1)*r-w*w*r*r*r/B-1/r*A*(m+4*pi*r*r*r*p)*(3*m-r+4*pi*r*r*r*p))*K;
H0=H0/(3*m+0.5*(l+2)*(l-1)*r+4*pi*r*r*r*p);
return(H0);
}

double Vf(double r,double w,double e,double p,double B,double A,double Dp,double W,double H0,double X)
{
return(1/w/w/(e+p)*B*(1/sqrt(B)*X+1/r*Dp/sqrt(A)*W-0.5*(e+p)*H0));
}

double gf(double e)

{
double le, p, gamma;
int i=j-1;
e=e/G*c*c;
le=log10(e);
for(;i>0;i--)
  {
  if(le>a[i-1])
    {
    p=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
    p=pow(10,p);
    gamma=(e+p/c/c)/e*(b[i]-b[i-1])/(a[i]-a[i-1]);
    return(gamma);
    }
  }
return((pow(10,a[0])+pow(10,b[0])/c/c)/e*(b[i]-b[i-1])/(a[i]-a[i-1]));
}

double Ff(double r,double A,double e,double p,double m,double Dp)
{
double F;
F=8*pi/r/r*sqrt(A);
F=F*(e+3*p+r*Dp-(m/4/pi/r/r/r+p)*(4+A*(m/r-4*pi*r*r*e)));
return(F);
}

double interp_E2p(double cE)
{
double logp,logE=log10(cE),cp;
int n=j-1;
for(;n>0;n--)
{
if(logE>Ga[n-1])
{
logp=(logE-Ga[n-1])*(Gb[n]-Gb[n-1])/(Ga[n]-Ga[n-1])+Gb[n-1];
cp=pow(10,logp);
return(cp);
}
}
return(0);
}

double interp_p2E(double cp)
{
if(cp<0)
{loopTOV=0;return 0;}
int n=j-1;
double logp=log10(cp),logE,cE;
for(;n>0;n--)
{
if(logp>Gb[n-1])
{
logE=(logp-Gb[n-1])*(Ga[n]-Ga[n-1])/(Gb[n]-Gb[n-1])+Ga[n-1];
cE=pow(10,logE);
return(cE);
}
}
return(0);
}

void fEul(double fr,double fm,double fp,double fE,double h)
{
temp_r=fr+h;
temp_m=fm+4.0*pi*fr*fr*fE*h;
temp_p=h*(fE+fp)*(fm+4*pi*fr*fr*fr*fp)/(2*fm*fr-fr*fr)+fp;
temp_E=interp_p2E(temp_p);
return;
}

double getM(double pEc)
{
loopTOV=1;
if(pEc<5e16) return 0;
double fEc=pEc*6.6741e-11;
double h=3e-8,localRes,r;
double m,p=interp_E2p(fEc),E;
double m11,p11,m21,p21,E21,m22,p22,E22,m31,p31,E31,m61,p61,E61,m62,p62,E62;
double TR3m,TR3p,TR3E,TR2m,TR2p;
double res1,res2;
r=h;
E=fEc;
m=4.0/3*pi*r*r*r*E;
while(r<1e-4)
{
fEul(r,m,p,E,h);
m11=temp_m;
p11=temp_p;
fEul(r,m,p,E,0.5*h);
m21=temp_m;
p21=temp_p;
E21=temp_E;
fEul(r+0.5*h,m21,p21,E21,0.5*h);
m22=temp_m;
p22=temp_p;
E22=temp_E;
fEul(r,m,p,E,h/3.0);
m31=temp_m;
p31=temp_p;
E31=temp_E;
fEul(r,m,p,E,h/6.0);
m61=temp_m;
p61=temp_p;
E61=temp_E;
fEul(r+h/6,m61,p61,E61,h/6.0);
m62=temp_m;
p62=temp_p;
E62=temp_E;
if(loopTOV<1) break;
TR3m=m22+9*(m62-m31);
TR3p=p22+9*(p62-p31);
TR3E=E22+9*(E62-E31);
TR2m=2*m22-m11;
TR2p=2*p22-p11;

res1=(TR3m>TR2m)?((TR3m-TR2m)/TR3m):((TR2m-TR3m)/TR3m);
res2=(TR3p>TR2p)?((TR3p-TR2p)/TR3p):((TR2p-TR3p)/TR3p);
localRes=(res1>res2)?res1:res2;

if((localRes<5e-6) && (h<1.2e-7))
{
r=r+h;
m=TR3m;
p=TR3p;
E=TR3E;
h=h*2;
}else if((localRes>2e-4) && (h>7.5e-9)){
h=h*0.5;
}else{
r=r+h;
m=TR3m;
p=TR3p;
E=TR3E;
}
}
return(m*Mscale);
}

double M2rhoc(double fM)
{
	double e0=1e18,e1,de=1e11,m0,m1,objE;
	int pf;
	e1=e0+de;
	m0=getM(e0);
	m1=getM(e1);
	for(pf=0;pf<8;pf++)
	{
		e0=e0-(m1-fM)/(m1-m0)*de;
		if(e0>pow(10,a[j-1]))
		{
			m0=getM(pow(10,a[j-1]));
			objE=(fabs(m0-fM)>fabs(m1-fM))?e1:pow(10,a[j-1]);
			return objE;
		}
		e1=e0+de;
		m0=getM(e0);
		m1=getM(e1);
	}
	objE=e0;
return objE;
}

int fmode(double e0)
{
	double dr=0.5;
	double r,r0=1,R,RR,drx,rx;
	double ne,p,e,m,mR,A,B=1.0,BR,Bfactor,power,I=0,J,DDf,Df=0,f=1;
	double m1, m2, m3, m4;
	double p1, p2, p3, p4;
	double B1, B2, B3, B4;
	double H1,H0,K,W,X,F,V,Dv,gamma,Dp,N;
	double DH11,H01,H02,K1,K2,x,X1,X2,Xp1,Xp2,W1,W2,V1,V2,V01,V02,V0;
	double w,wcheck,o[2],wi,aR,bR,gR,hR,kR,n,Y1,Y2,Z,DZ,DDZ,VZ,Ar1,Ar2,Ai1,Ai2,ar,ai,Ar[2],Ai[2],Br[2],Bi[2];
	int t,q,wpf,rpf,pfEnd,n4;
	scr++;
	printf("\b%c",dysm[scr%6]);
	power=pow(10,b[0]);
	r=r0;
	ne=e0;
	e=e0;
	p=che(e);
	m=1.3333333*r*r*pi*e*r;
	wpf=-1;
	std::cout << p<<"testvalue p";
	for(;p>power;r=r+dr)
		{
			if(wpf>49995)
				{
					temp_rhoc=ne;
					temp_M=0;
					temp_R=0;
					temp_freq=0;
					temp_dpt=0;
					return 1;
				}
			wpf++;
			rhofile[wpf]=e*Gc2;   /*-- e,p,m in G=c=1 --*/
			pfile[wpf]=p*Gc4;
			Bfile[wpf]=B;
			mfile[wpf]=m*Gc2;
			A=1/(1-2*m*Gc2/r);
			p1=fp(r,p,e,m);
			m1=fm(r,e);
			/////////////////////////////////
			if((p+dr*p1/2)>power) e=ch(p+dr*p1/2); 
			else break;
			B1=Bf(r,p,m,B);
			p2=fp(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
			m2=fm(r+dr/2,e);
			if((p+dr*p2/2)>power) e=ch(p+dr*p2/2); else break;
			B2=Bf(r+dr/2,p+dr*p1/2,m+dr*m1/2,B+dr*B1/2);
			p3=fp(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
			m3=fm(r+dr/2,e);
			if((p+dr*p3)>power) e=ch(p+dr*p3); else break;
			B3=Bf(r+dr/2,p+dr*p2/2,m+dr*m2/2,B+dr*B2/2);
			p4=fp(r+dr,p+dr*p3,e,m+dr*m3);
			m4=fm(r+dr,e);
			B4=Bf(r+dr,p+dr*p3,m+dr*m3,B+dr*B3);
			J=-4*pi*(e+p/c/c)*Gc2*r*A;
			DDf=-(4/r*Df+J*Df+4/r*J*f);
			I=I-2.0/3*f*J/sqrt(A*B)*r*r*r*dr;
			p=p+dr*(p1+2*p2+2*p3+p4)/6;
			e=ch(p);
			m=m+dr*(m1+2*m2+2*m3+m4)/6;
			B=B+dr*(B1+2*B2+2*B3+B4)/6;
			f=f+Df*dr;
			Df=Df+DDf*dr;
		}
		pfEnd=wpf;
		R=r;
		mR=m*Gc2;
		BR=1-2*Gc2*m/r;
		Bfactor=BR/B;
		gamma=(b[1]-b[0])/(a[1]-a[0]);
		N=1/(gamma-1);
		std::cout << '\n' << p1 << '\n' << p2;
		RR=R-(N+1)*(p-dr*(p1+2*p2+2*p3+p4)/6)/(p1+2*p2+2*p3+p4)*6;
		//sovle TOV???
		wpf=-1;
		rpf=-1;
		scr++;printf("\b%c",dysm[scr%6]);
		for(n4=0;n4<pfEnd+1;n4++)
			{
				Bcor[n4]=Bfile[n4]*Bfactor;
			}
		I=I/sqrt(Bfactor);
		I=I/(f+2*I/r/r/r)/Gc2;
		I=m*sqrt(m/I)*Gc2;
		o[0]=(-0.0047+0.133*I+0.575*I*I)/mR-0.1e-5;
		o[1]=o[0]+0.2e-5;
		q=1;
		wcheck=0;
		for(t=0;;t++)
		{
			if(t>20)
			{
				temp_rhoc=ne;
				temp_M=0;
				temp_R=0;
				temp_freq=0;
				temp_dpt=0;
				return 2;
			}
			scr++;printf("\b%c",dysm[scr%6]);
			if(t==0) w=o[t];
			else w=o[q];
			e=rhofile[0];
			p=pfile[0];
			B=Bcor[0];
			W=1.0;
			K=(e+p);
			X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K);
			H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
			rpf=-1;
			wpf=-1;
			r=r0;
			while(rpf<pfEnd)
			{
				rpf++;
				p=pfile[rpf];
				e=rhofile[rpf];
				B=Bcor[rpf];
				m=mfile[rpf];
				Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
				Dv=-2*Dp/(e+p);
				A=1/(1-2*m/r);
				gamma=gf(e);
				H0=H0f(r,B,X,m,p,A,H1,K,w);
				V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
				if(r==r0)V01=V;
				if(fabs((wcheck-w)/w)<Cri)
				{
					wpf++;
					Wfile[wpf]=sqrt(1-2*m/r)*W;
					Vfile[wpf]=V;
				}
				F=Ff(r,A,e,p,m,Dp);
				DH11=DH1(r,m,A,p,e,H1,H0,K,V);
				K = K+DK(r,H0,H1,Dv,K,e,p,A,W)*dr;
				W = W + DW(r, W, A, gamma, p, B, X, V, H0, K)*dr;
				X = X + DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W)*dr;
				H1=H1+DH11*dr;
				r=r+dr;
			}
			wpf=-1;
			rpf=-1;
			X1=X;
			Xp1=DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W);
			K1=K;H01=H0f(r,B,X,m,p,A,H1,K,w);W1=W;V1=V;
			scr++;printf("\b%c",dysm[scr%6]);
			p=pfile[0];
			e=rhofile[0];
			B=Bcor[0];
			W=1.0;K=-(e+p);
			X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K);
			H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);
			r=r0;
			while(rpf<pfEnd)
			{
				rpf++;
				p=pfile[rpf];
				e=rhofile[rpf];
				B=Bcor[rpf];
				m=mfile[rpf];
				Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
				Dv=-2*Dp/(e+p);
				A=1/(1-2*m/r);
				gamma=gf(e);
				H0=H0f(r,B,X,m,p,A,H1,K,w);
				V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
				if(r==r0)V02=V;
				if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
				{
					wpf++;
					Wfile[wpf]=sqrt(1-2*m/r)*W;
					Vfile[wpf]=V;
				}
				F=Ff(r,A,e,p,m,Dp);
				DH11=DH1(r,m,A,p,e,H1,H0,K,V);
				K=K+DK(r,H0,H1,Dv,K,e,p,A,W)*dr;
				W=W+DW(r,W,A,gamma,p,B,X,V,H0,K)*dr;
				X=X+DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W)*dr;
				H1=H1+DH11*dr;
				r=r+dr;
			}
			wpf=-1;
			rpf=-1;
			X2=X;
			Xp2 = DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W);
			K2=K;H02=H0f(r,B,X,m,p,A,H1,K,w);W2=W;V2=V;
			x=-(X1-(RR-R)/(N+1)*Xp1)/(X2-(RR-R)/(N+1)*Xp2);
			H0=H01+x*H02; K=K1+x*K2; W=W1+x*W2; V=V1+x*V2;V0=V01+x*V02;
			if(fabs((wcheck-w)/w)<Cri)
			{
				r=r0;
				while(rpf<pfEnd)
				{
					rpf++;
					W1=Wfile1[rpf];
					W2=Wfile2[rpf];
					W=W1+x*W2;
					wpf++;
					Wfile[wpf]=W/(1+x);
					V1=Vfile1[rpf];
					V2=Vfile2[rpf];
					V=V1+x*V2;
					Vfile[wpf]=V/V0;
					r=r+dr;
				}
				break;
			}
			wcheck=w;
			n=0.5*(l-1)*(l+2);
			aR=-(n*R+3*mR)/(w*w*R*R-(n+1)*mR/R);
			bR=(n*R*(R-2*mR)-w*w*R*R*R*R+mR*(R-3*mR));
			bR=bR/(R-2*mR)/(w*w*R*R-(n+1)*mR/R);
			gR=n*(n+1)*R*R+3*n*mR*R+6*mR*mR;
			gR=gR/R/R/(n*R+3*mR);
			hR=-n*R*R+3*n*mR*R+3*mR*mR;
			hR=hR/(R-2*mR)/(n*R+3*mR);
			kR=-R*R/(R-2*mR);
			Y1=K;
			Y2=aR*H0+bR*K;
			Z=(kR*Y1-Y2)/(kR*gR-hR);
			DZ=(gR*Y2-hR*Y1)/(gR*kR-hR);
			if(w<2e-6)
				{
					temp_rhoc=ne;
					temp_M=0;
					temp_M=0;
					temp_R=0;
					temp_freq=0;
					temp_dpt=0;
					return 3;
				}
			for(r=R;r<25.0/w;r=r+dr)
			{
				drx=dr/(1-2*mR/r);
				VZ=(1-2*mR/r)/r/r/r/(n*r+3*mR)/(n*r+3*mR);
				VZ=VZ*(2.0*n*n*(n+1)*r*r*r+6.0*n*n*mR*r*r+18.0*n*mR*mR*r+18*mR*mR*mR);
				DDZ=(VZ-w*w)*Z;
				Z=Z+DZ*drx;
				DZ=DZ+DDZ*drx;
			}
			r=r-dr;
			rx=r+2*mR*log(r/2/mR-1);
			Ar1=2*cos(w*rx)-2*(n+1)/w/r*sin(w*rx)+1/w/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx));
			Ai1=2*sin(w*rx)+2*(n+1)/w/r*cos(w*rx)-1/w/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx));
			Ar2=-2*w*sin(w*rx)-2*(n+1)*cos(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx))+(1-2*mR/r)*2*(n+1)/w/r/r*sin(w*rx);
			Ai2=2*w*cos(w*rx)-2*(n+1)*sin(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx))-(1-2*mR/r)*2*(n+1)/w/r/r*cos(w*rx);
			ar=(Ai2*Z-Ai1*DZ)/(Ar1*Ai2-Ar2*Ai1);
			ai=-(Ar1*DZ-Ar2*Z)/(Ar1*Ai2-Ar2*Ai1);
			if(t==0)
			{
				Ar[t]=ar;
				Ai[t]=ai;
			}else{
				Ar[q]=ar;
				Ai[q]=ai;
				Br[0]=(o[0]*Ar[1]-o[1]*Ar[0])/(o[0]-o[1]);
				Br[1]=(Ar[0]-Ar[1])/(o[0]-o[1]);
				Bi[0]=(o[0]*Ai[1]-o[1]*Ai[0])/(o[0]-o[1]);
				Bi[1]=(Ai[0]-Ai[1])/(o[0]-o[1]);
				w=-(Br[0]*Br[1]+Bi[0]*Bi[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
				if (w<=o[0])
				{
					o[1]=o[0];o[0]=w;Ar[1]=Ar[0];Ai[1]=Ai[0];q=0;
				}else if(w>=o[1]){
					o[0]=o[1];o[1]=w;Ar[0]=Ar[1];Ai[0]=Ai[1];q=1;
				}else if((o[1]-w)>(w-o[0])){
					o[1]=w;q=1;
				}else{
					o[0]=w;q=0;
				}
			}
		}
		wi=(Br[0]*Bi[1]-Bi[0]*Br[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
		temp_rhoc=ne;
		temp_M=mR;
		temp_M=temp_M/(M*Gc2);
		temp_R=RR/1000;
		temp_freq=w*c/(2000*pi);
		temp_dpt=1/(wi*c);
		scr++;printf("\b%c",dysm[scr%6]);
	return 0;
}
void fPrint(void)
{
	time_t timep = time(NULL);
	struct tm *sysTime = (struct tm*)malloc(sizeof(struct tm));
	gmtime_s(sysTime, &timep);
	fprintf(outf, "\n[f-mode]\tTime: %d-%d-%d-%dh-%dmin-%ds\n", 1900 + sysTime->tm_year, 1 + sysTime->tm_mon, sysTime->tm_mday, 8 + sysTime->tm_hour, sysTime->tm_min, sysTime->tm_sec);
	fprintf(outf, "====================================\n");
	fprintf(outf, "%-10s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n", "EoS", "rhoc(kg.m-3)", "M(Msun)", "R(km)", "frequency(kHz)", "dampingTime(s)");
}
int main()
{
	int pf, n, len1, neos, mode;
	double srhoc[40], smass[40], sradius[40], sfreq[40], sdptime[40];
	char pathName[50];
	
	//char name[200][50];
	std::cout << "Advanced codes F-mode_v10 for Linux" << std::endl << std::endl;
	printf("================================================\n\n");
	//////////////////////////////////////////////////

	remove("temp_dir");
	system("dir/b EoS_lib >> temp_dir");

	//give attention 
	//fopen_s(&temp, "temp_dir", "r");

	std::ifstream eostable("temp_dir");

	if (eostable.is_open())
	{
		printf("Equation of State files detected:\n\n");
	}
	else
	{
		std::cout << "No Equation of State files detected, automatically quit.";
		return 0;
	}
	n = 0;
	char Eatthisshit[100];
	while (eostable.getline(Eatthisshit,100))
	{
		list[n].Name = Eatthisshit;
	//	printf("%d\t%s\n", n + 1, list[n].Name);
		std::cout <<n<<"   "<< list[n].Name<<'\n';
		n++;
	}
	len1 = n;
	for (neos = 0;neos<n; neos++)
	{
		list[neos].tName = list[neos].Name.substr(0, list[neos].Name.find(".txt", 0));
		sprintf_s(pathName, 50, "./EoS_lib/%s", list[neos].Name.c_str());//read eos table
		j = 0;
		if (fopen_s(&inf, pathName, "r"))
		{
			std::cout << '\n' << "Cannot load" << list[neos - 1].Name << '\n'; continue;
		}
		while (fscanf_s(inf, "%lf", &a[j]) == 1){ fscanf_s(inf, "%lf%*[^\n]", &b[j]); j++; }
		fclose(inf);
	}
	printf("\nCaution: please ensure all *.txt files above are available EoS.\n\nPress any key to start (Ctrl+C to quit).\n");
	////////////////////////////////////////////////////////////
	mfile = (double*)malloc(50000 * sizeof(double));
	pfile = (double*)malloc(50000 * sizeof(double));
	rhofile = (double*)malloc(50000 * sizeof(double));
	Bfile = (double*)malloc(50000 * sizeof(double));
	Bcor = (double*)malloc(50000 * sizeof(double));
	Wfile = (double*)malloc(50000 * sizeof(double));
	Wfile1 = (double*)malloc(50000 * sizeof(double));
	Wfile2 = (double*)malloc(50000 * sizeof(double));
	Vfile = (double*)malloc(50000 * sizeof(double));
	Vfile1 = (double*)malloc(50000 * sizeof(double));
	Vfile2 = (double*)malloc(50000 * sizeof(double));
	getchar();
	printf("\nPlease select a work mode:\n\t[0]\tCompute within a density range\n\t[1]\tSearch for a given mass\n\t[2]\tSearch for a given frequency\n\nEnter the mode number: ");
	std::cin >> mode;
	///////////////////////////////////////////////////////////////
	switch (mode)
	{
	case 0:
		{
			double strhoc, edrhoc, deltrho;//centre density
			printf("\nTips: 3e17<rhoc<3e18\n\nStart rhoc(kg.m-3)=");
			std::cin >> strhoc;
			while ((strhoc<3e17) || (strhoc>3e18))
			{
				printf("\nWrong input, Start at rhoc(kg.m-3)=");
				std::cin >> strhoc;
			}
			printf("\nEnd at rhoc(kg.m-3)=");
			std::cin >> edrhoc;
			while ((edrhoc<3e17) || (edrhoc>3e18))
			{
				printf("\nWrong input, End at rhoc(kg.m-3)=");
				std::cin >> edrhoc;
			}
			fopen_s(&outf, "Results_mod_0.dat", "a");
			fPrint();
			//main loop for case 0
			for (neos = 1; neos <= n; neos++)
			{
				 scr = 0;
				if (strhoc == edrhoc)
				{
					printf("\n%s is under computed...      %c", list[neos - 1].tName.c_str(), dysm[scr % 6]);
					fmode(strhoc);
					printf("\b \n%-8s\t%-8s\t%-15s\t%-18s\n", "M(Msun)", "R(km)", "frequency(kHz)", "dampingTime(s)");
					printf("%-8.2lf\t%-8.2lf\t%-15lf\t%-18lf\n", temp_M, temp_R, temp_freq, temp_dpt);
					fprintf(outf, "\n%-10s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n", list[neos - 1].tName.c_str(), strhoc, temp_M, temp_R, temp_freq, temp_dpt);
					continue;
				}
					if (edrhoc>strhoc)
					{
						edrhoc = (edrhoc<pow(10, a[j - 1])) ? edrhoc : pow(10, a[j - 1]);
					}
					else
					{
						strhoc = (strhoc<pow(10, a[j - 1])) ? strhoc : pow(10, a[j - 1]);
					}
				deltrho = 0.03125*(edrhoc - strhoc); //density interval 1/32
				printf("\n%s is under computed...      %c", list[neos - 1].tName.c_str(), dysm[scr % 6]);
				for (pf = 0; pf<33; pf++)
				{
					srhoc[pf] = strhoc + pf*deltrho;
					fmode(srhoc[pf]);
					smass[pf] = temp_M; sradius[pf] = temp_R; sfreq[pf] = temp_freq; sdptime[pf] = temp_dpt;
					scr++; printf("\b\b\b\b\b\b%3d%% %c", pf * 3, dysm[scr % 6]);
					fprintf(outf, "\n%-10s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n", list[neos - 1].tName.c_str(), srhoc[pf], smass[pf], sradius[pf], sfreq[pf], sdptime[pf]);
				}
				scr++; printf("\b\b\b\b\b\b100%% %c", dysm[scr % 6]);
				printf("\b ");
			}
			fclose(outf);
		}
		break;
	case 1:
		{
				  double objM, objrho;
				  printf("\nTips: 0.5<M<2.2\n\nM(Msun)=");
				  std::cin >> objM;
				  while ((objM<0.5) || (objM>2.2))
				  {
					  printf("\nWrong input, M(Msun)=");
					  std::cin >> objM;
				  }
				  fopen_s(&outf, "Results_mod_1.dat", "a");
				  fPrint();
				  for (neos = 1; neos <= len1; neos++)
				  {
					   scr = 0;				  
					  for (pf = 0; pf<j; pf++)
					  {
						  Ga[pf] = a[pf] - 10.175607290470733;
						  Gb[pf] = b[pf] - 27.129849799910058;
					  }
					  printf("\n%s is under computed...      %c", list[neos - 1].Name.c_str(), dysm[scr % 6]);
					  scr++; printf("\b\b\b\b\b\b  0%% %c", dysm[scr % 6]);
					  std::cout << objM<<std::endl;
					  objrho = M2rhoc(objM);
					  std::cout << objrho;
					  fmode(objrho);
					  scr++; printf("\b\b\b\b\b\b100%% %c", dysm[scr % 6]);
					  printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n", temp_M, temp_freq, temp_dpt);
					  n = pf;
					  fprintf(outf, "\n%-10s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n", list[neos - 1].tName.c_str(), objrho, temp_M, temp_R, temp_freq, temp_dpt);
				  }
				  fclose(outf);
		}
		break;
	case 2:
	{
			  double objF, lefrhoc, midrhoc, rigrhoc, lefF, midF, rigF;
			  printf("\nTips: 0.5<freq<3.6\n\nfreq(kHz)=");
			  fflush(stdin);
			  scanf_s("%lf", &objF);
			  while ((objF<0.5) || (objF>3.6))
			  {
				  printf("\nWrong input, freq(kHz)=");
				  fflush(stdin);
				  scanf_s("%lf", &objF);
			  }
			  fopen_s(&outf, "Results_mod_2.dat", "a");
			  fPrint();
			  for (neos = 1; neos <= len1; neos++)
			  {
				  scr = 0;
				  printf("\n%s is under computed...      %c", list[neos - 1].tName.c_str(), dysm[scr % 6]);
				  lefrhoc = 4.2e17;
				  rigrhoc = (pow(10, a[j - 1])<1.0e18) ? pow(10, a[j - 1]) : 1.0e18;
				  fmode(lefrhoc);
				  lefF = temp_freq;
				  fmode(rigrhoc);
				  rigF = temp_freq;
				  scr++; printf("\b\b\b\b\b\b  0%% %c", dysm[scr % 6]);
				  for (pf = 0; pf<32; pf++)
				  {
					  scr++; printf("\b\b\b\b\b\b%3d%% %c", pf * 3, dysm[scr % 6]);
					  if ((lefF<objF) && (objF<rigF))
					  {
						  midrhoc = 0.5*(lefrhoc + rigrhoc);
						  fmode(midrhoc);
						  midF = temp_freq;
						  if (midF<objF)
						  {
							  lefrhoc = midrhoc;
							  lefF = midF;
						  }
						  else{
							  rigrhoc = midrhoc;
							  rigF = midF;
						  }
					  }
					  else if (rigF<objF){
						  lefrhoc = rigrhoc;
						  lefF = rigF;
						  rigrhoc = 1.1*rigrhoc;
						  if (pow(10, a[j - 1])<rigrhoc)
						  {
							  rigrhoc = pow(10, a[j - 1]);
							  lefrhoc = rigrhoc;
							  break;
						  }
						  fmode(rigrhoc);
						  rigF = temp_freq;
						  if (rigF<lefF)
						  {
							  midrhoc = 0.5*(lefrhoc + rigrhoc);
							  fmode(midrhoc);
							  midF = temp_freq;
							  midrhoc = (midF>lefF) ? midrhoc : lefrhoc;
							  rigrhoc = midrhoc;
							  lefrhoc = midrhoc;
							  break;
						  }
					  }
					  else{
						  rigrhoc = lefrhoc;
						  rigF = lefF;
						  lefrhoc = 0.9*lefrhoc;
						  fmode(lefrhoc);
						  lefF = temp_freq;
						  if (temp_M<0.35)
						  {
							  rigrhoc = lefrhoc;
							  break;
						  }
					  }
				  }
				  midrhoc = 0.5*(lefrhoc + rigrhoc);
				  fmode(midrhoc);
				  scr++; 
				  printf("\b\b\b\b\b\b100%% %c", dysm[scr % 6]);
				  printf("\b f_%.2lf=%.3lf(kHz), dT=%.3lf(s)\n", temp_M, temp_freq, temp_dpt);
				  fprintf(outf, "\n%-10s\t%-15e\t%-10lf\t%-10lf\t%-15lf\t%-18lf\n", list[neos - 1].tName.c_str(), midrhoc, temp_M, temp_R, temp_freq, temp_dpt);
			  }
			  fclose(outf);
	}
	break;

	default:
		break;
	}

	free(mfile);
	free(pfile);
	free(rhofile);
	free(Bfile);
	free(Bcor);
	free(Wfile);
	free(Wfile1);
	free(Wfile2);
	free(Vfile);
	free(Vfile1);
	free(Vfile2);
	printf("Done!\n");
	return 0;
	
}
