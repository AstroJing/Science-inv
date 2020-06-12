#include"fmode.h"
#define pi 3.14159265358979
#define G 6.673e-11
#define c 2.99792458e8
#define M 1.989e30
#define Gc2 7.42471382405e-28
#define Gc4 8.26110825251e-45
#define Mscale 2.033931261665867e5
#define rho0   2.67e17
//a= rho   b=pres
static int loopTOV = 1;
double temp_r, temp_m, temp_p, temp_E;
double interp_E2p(double cE, double Ga[], double Gb[],int j)
{
	double logp, logE = log10(cE), cp;
	int n = j - 1;
	for (; n>0; n--)
	{
		if (logE>Ga[n - 1])
		{
			logp = (logE - Ga[n - 1])*(Gb[n] - Gb[n - 1]) / (Ga[n] - Ga[n - 1]) + Gb[n - 1];
			cp = pow(10, logp);
			return(cp);
		}
	}
	return(0);
}
double interp_p2E(double cp, double Ga[], double Gb[], int j)
{
	if (cp<0)
	{
		loopTOV = 0; return 0;
	}
	int n = j - 1;
	double logp = log10(cp), logE, cE;
	for (; n>0; n--)
	{
		if (logp>Gb[n - 1])
		{
			logE = (logp - Gb[n - 1])*(Ga[n] - Ga[n - 1]) / (Gb[n] - Gb[n - 1]) + Ga[n - 1];
			cE = pow(10, logE);
			return(cE);
		}
	}
	return(0);
}

void fEul(double fr, double fm, double fp, double fE, double h, double Ga[], double Gb[], int j)
{
	temp_r = fr + h;
	temp_m = fm + 4.0*pi*fr*fr*fE*h;
	temp_p = h*(fE + fp)*(fm + 4 * pi*fr*fr*fr*fp) / (2 * fm*fr - fr*fr) + fp;
	temp_E = interp_p2E(temp_p,Ga,Gb,j);
	return;
}
// use center density got the correspond Mass
double getM(double pEc, double Ga[], double Gb[], int j)
{
	loopTOV = 1;
	if (pEc<5e16) return 0;
	double fEc = pEc*6.6741e-11;
	double h = 3e-8, localRes, r;
	double m, p = interp_E2p(fEc, Ga, Gb, j), E;
	double m11, p11, m21, p21, E21, m22, p22, E22, m31, p31, E31, m61, p61, E61, m62, p62, E62;
	double TR3m, TR3p, TR3E, TR2m, TR2p;
	double res1, res2;
	r = h;
	E = fEc;
	m = 4.0 / 3 * pi*r*r*r*E;
	while (r<1e-4)
	{
		fEul(r, m, p, E, h, Ga, Gb, j);
		m11 = temp_m;
		p11 = temp_p;
		fEul(r, m, p, E, 0.5*h, Ga, Gb, j);
		m21 = temp_m;
		p21 = temp_p;
		E21 = temp_E;
		fEul(r + 0.5*h, m21, p21, E21, 0.5*h, Ga, Gb, j);
		m22 = temp_m;
		p22 = temp_p;
		E22 = temp_E;
		fEul(r, m, p, E, h / 3.0, Ga, Gb, j);
		m31 = temp_m;
		p31 = temp_p;
		E31 = temp_E;
		fEul(r, m, p, E, h / 6.0, Ga, Gb, j);
		m61 = temp_m;
		p61 = temp_p;
		E61 = temp_E;
		fEul(r + h / 6, m61, p61, E61, h / 6.0, Ga, Gb, j);
		m62 = temp_m;
		p62 = temp_p;
		E62 = temp_E;
		if (loopTOV<1) break;
		TR3m = m22 + 9 * (m62 - m31);
		TR3p = p22 + 9 * (p62 - p31);
		TR3E = E22 + 9 * (E62 - E31);
		TR2m = 2 * m22 - m11;
		TR2p = 2 * p22 - p11;

		res1 = (TR3m>TR2m) ? ((TR3m - TR2m) / TR3m) : ((TR2m - TR3m) / TR3m);
		res2 = (TR3p>TR2p) ? ((TR3p - TR2p) / TR3p) : ((TR2p - TR3p) / TR3p);
		localRes = (res1>res2) ? res1 : res2;

		if ((localRes<5e-6) && (h<1.2e-7))
		{
			r = r + h;
			m = TR3m;
			p = TR3p;
			E = TR3E;
			h = h * 2;
		}
		else if ((localRes>2e-4) && (h>7.5e-9)){
			h = h*0.5;
		}
		else{
			r = r + h;
			m = TR3m;
			p = TR3p;
			E = TR3E;
		}
	}
	return(m*Mscale);
}

double M2rho(double fM, double Ga[], double Gb[], int j)
{
	double e0 = 1e17, e1 = 4e18, temp, m0, m1;
	int pf = 0;
	m0 = getM(e0, Ga, Gb, j);
	m1 = getM(e1, Ga, Gb, j);
	if (fM<m0 || fM>m1)
	{
		std::cout << "EoS can not support this mass";
		return 0;  //erro process
	}
	else
	{
		while (pf<30)
		{
			temp = (e0 + e1) / 2;
			m0 = getM(temp, Ga, Gb, j);
		///	std::cout<<std::endl<<"     fm= "<<fM<<"    m0="<<m0<<std::endl;
			if (fM - m0<0.01&&m0-fM<0.01)
			{
				return temp;
			}
			if (fM > m0)
			{
				e0 = temp;
			}
			else
			{
				e1 = temp;
			}
			pf++;
		}
		std::cout << "can not fin the  proper center density of MASS=" << fM <<  "\n\n";
		return 0;
	}
}
