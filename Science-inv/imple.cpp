#include"Gen.h"
#define Esym0   31.6 //MeV
#define RHO0    0.16
#define E0RHO0  -15.9
#define xme     0.511
#define xmmu    105.658
#define H       1239.8522
#define r0m2    39.0625  // 1/(RHO0*RHO0)
using namespace std;
pair<double,double> GetTrho(double L,double Ksym,double Jsym,double J0,double K0)
{
	double gap = 0.001;
	double U=0.2;
	double Kmu;
	double RHOT;
	double delta;
	for (int i = 1; i < 10 / gap; i++)  // calculate delta
	{

		delta = solver( U,  L,  Ksym,  Jsym, J0,  K0,  1,&Kmu);
		if (Kmu>0)
		{
			RHOT = U*RHO0;
			break;
		}
		else
		{
			if (U>1.5)
			{
				RHOT = 0.16;
				delta = 1;
		//		rhoet = 1;
		//		rho_mu = 1;
				break;
			}
		}
		U = U + gap;
	}

	return pair<double,double> (RHOT,delta);
}
double solver(double U,double L, double Ksym, double Jsym, double J0, double K0,int c,double* point)
{
	double RHO, Esym, dEsym, ddEsym, E0, dE0, ddE0;
	double Kmu;
	double rho_p, rho_e, rho_mu, mu_e, mu_mu, pf_e, pf_mu, pf_p;
	double en, en1;
	double beta,delta;
	double pfe;
	/////////////need explanation///////////
	pf_e = 0;
	delta = 1;
	/////////////need explanation///////////
		RHO = RHO0*U;
		Esym = (Esym0 + L *(U - 1) / 3.0 + Ksym*(U - 1.0)*(U - 1.0) / 18.0 + Jsym*pow((U - 1.0), 3.0) / 162.0) / 197.3;
		////////////////delta calculate/////////////////////////
		en1 = 0;
		if (Esym>0)
		{
			for (beta = 0.0001; beta<1; beta += 0.0002)
			{
				rho_p = (1 - beta)*RHO / 2;
				pf_p = pow((3.0*pi*pi*rho_p), 0.3333333);
				mu_e = 4.0*Esym*197.3*beta;
				mu_mu = mu_e;
				if (mu_e - 0.511 < 0)
					mu_e = 0;
				else
				{
					pf_e = sqrt(mu_e*mu_e - 0.2611) / hbar;
					rho_e = pow(pf_e, 3)*pi*pi / 3.0;
				}
				if (mu_mu>105.66)
				{
					pf_mu = sqrt(mu_mu*mu_mu - 11164.0356)/ hbar;
					rho_mu = pow(pf_mu, 3)*pi*pi / 3.0;
				}
				else
				{
					pf_mu = 0;
				}
				en = pow(pf_p, 3.0) - pow(pf_e, 3.0) - pow(pf_mu, 3.0);
				if (en*en1<0)
				{
					delta = beta;
					break;
				}
				en1 = en;
			}
		}
		else
		{
			delta = 1;
		}
		dEsym = (L / 3.0 + Ksym*(U - 1.0) / 9.0 + Jsym*(U - 1.0)*(U - 1.0) / 54) / RHO0;
		ddEsym = r0m2*(Ksym / 9.0 + Jsym*(U - 1.0) / 27.0);
		E0 = E0RHO0 + K0 / 18.*((U - 1.0)*(U - 1.0)) + J0 / 162.*pow((U - 1.0), 3);
		dE0 = (K0 / 9.*(U - 1) + J0*(U - 1)*(U - 1) / 54) / RHO0;
		ddE0 = r0m2*(K0 / 9.0 + J0*(U - 1.0) / 27.0);
		if (c == 1)
		{
			Kmu = 2.0*RHO*dE0 + RHO*RHO * ddE0 + delta*delta*(RHO*RHO*ddEsym + 2.0*RHO*dEsym - 2.0*pow((RHO*dEsym), 2.0) / Esym / 197.3);
			*point = Kmu;
			return delta;
		}
		else if (c==0)
		{
			*point = E0;
			*(point + 1) = dE0;
			*(point + 2) = ddE0;
			*(point + 3) = Esym; 
			*(point + 4) = dEsym; 
			*(point + 5) = ddEsym;
			return delta;

		}


}
double PSolver(double rho,double delta,double Esym,double dE0,double dEsym,double* erhol)
{
	double pf_e,pf_mu,rho_e,mu_e, mu_mu;
	double rhoet=0;
	mu_e = 4.0*Esym*197.3*delta;
	mu_mu = mu_e;
	if (mu_e >= 0.511)

	{
		pf_e = sqrt(mu_e*mu_e - 1.022) / hbar;
		rho_e = pow(pf_e, 3)*pi*pi / 3.0;
		rhoet = rho_e;
		pf_mu = 0;
	}
	 if (mu_mu>105.66)
	{
		pf_mu = sqrt(mu_mu - 105.66)*(mu_mu + 105.66) / hbar;
		rhoet = pow(pf_mu, 3)*pi*pi / 3.0;
	}
	
	else
	{
		rhoet = 0;
		pf_mu = 0;
		pf_e = 0;
	}
//////////////////////////////////////////
	double t;
	double xlam = hbar / xme;
	t = xlam*pow(3.0*pi*pi*rhoet, 0.3333);
	double yita = xme / (8.0*pi*pi*xlam*xlam*xlam);
	double fai = t*sqrt(1.0 + t*t)*(2.0*t*t / 3.0 - 1.0) + log(t + sqrt(1.0 + t*t));//////////////need explaination///////////
	double pet = yita*fai;
	double xlammu = hbar / xmmu;
	double tmu = hbar*pf_mu / xmmu;
	double faimu = tmu*sqrt(1.0 + tmu*tmu)*(2.*tmu*tmu / 3.0 - 1.0) + log(tmu + sqrt(1.0 + tmu*tmu));
	double pmut = xmmu / (8 * pi*pi*xlammu*xlammu*xlammu)*faimu;    //yitamu*faimu
	//	double pbt = RHOT*RHOT*(dE0 + dEsym*delta*delta);
	double pt = (pet + pmut) + rho*rho*(dE0 + dEsym*delta*delta);
	double te = hbar*pf_e / xme;
	 *erhol = 0.125*(xme / (pi*pi*pow(xlam, 3))*te*(1.0 + pow(te, 2))*(2.0*pow(te, 2) + 1.0) - log(te + sqrt(1 + te*te))
		+ xmmu / (pi*pi*pow(xlammu, 3))*tmu*(1.0 + pow(tmu, 2))*(2.0*pow(tmu, 2) + 1.0) - log(tmu + sqrt(1 + tmu*tmu)));
	return pt;
}