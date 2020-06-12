#include "fmode.h"
using namespace std;
static int loopTOV = 1;
double temp_r, temp_m, temp_p, temp_E;
int j;
double interp_E2p(double cE,double Ga[],double Gb[])
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
double interp_p2E(double cp,double Ga[],double Gb[])
{
//	cout<<"   j="<<j<<endl;
//	getchar();
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

void fEul(double fr, double fm, double fp, double fE, double h, double Ga[], double Gb[])
{
	temp_r = fr + h;
	temp_m = fm + 4.0*pi*fr*fr*fE*h;
	temp_p = h*(fE + fp)*(fm + 4 * pi*fr*fr*fr*fp) / (2 * fm*fr - fr*fr) + fp;
	temp_E = interp_p2E(temp_p,Ga,Gb);
	return;
}
double getM(double pEc, double Ga[], double Gb[])
{
	loopTOV = 1;
	if (pEc<5e16) return 0;
	double fEc = pEc*6.6741e-11;
	double h = 3e-8, localRes, r;
	double m, p = interp_E2p(fEc,Ga,Gb), E;
	double m11, p11, m21, p21, E21, m22, p22, E22, m31, p31, E31, m61, p61, E61, m62, p62, E62;
	double TR3m, TR3p, TR3E, TR2m, TR2p;
	double res1, res2;
	r = h;
	E = fEc;
	m = 4.0 / 3 * pi*r*r*r*E;
	while (r<1e-4)
	{
		fEul(r, m, p, E, h,  Ga,  Gb);
		m11 = temp_m;
		p11 = temp_p;
		fEul(r, m, p, E, 0.5*h, Ga, Gb);
		m21 = temp_m;
		p21 = temp_p;
		E21 = temp_E;
		fEul(r + 0.5*h, m21, p21, E21, 0.5*h, Ga, Gb);
		m22 = temp_m;
		p22 = temp_p;
		E22 = temp_E;
		fEul(r, m, p, E, h / 3.0, Ga, Gb);
		m31 = temp_m;
		p31 = temp_p;
		E31 = temp_E;
		fEul(r, m, p, E, h / 6.0, Ga, Gb);
		m61 = temp_m;
		p61 = temp_p;
		E61 = temp_E;
		fEul(r + h / 6, m61, p61, E61, h / 6.0, Ga, Gb);
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
double M2rhoc(double fM,double Ga[],double Gb[],int len)
{
	j = len;
	double e0 = 1e17, e1 = 1e18, mid, m0;
	int upline = 0;
	m0 = getM(e0,Ga,Gb);
	while (upline<1000)
	{
		mid = (e0 + e1) / 2;
		m0 = getM(mid,Ga,Gb);
		if (m0 == -2)
		{
			e0 = e0 + 0.01*e0;
			cout << "  m0   " << m0 << "     " << endl;
			if (e0>1e18)
				return 0;
		}
		if (fabs(fM - m0)<0.01)
		{
			return mid;
		}
		if (fM > m0)
		{
			e0 = mid;
		}
		else
		{
			e1 = mid;
		}
		upline++;
	}
	std::cout << "can not fin the  proper center density of MASS=1.4\n\n";
	return 0;
	//		}
}
