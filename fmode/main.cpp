#include"fmode.h"
using namespace std;
static double a[3000], b[3000], l = 2;
//a= rho   b=pres
static int j, scr, loopTOV = 1;
char dysm[7] = "/-\\/-\\";
FILE *inf, *outf;
double temp_rhoc, temp_M, temp_R, temp_freq, temp_dpt;
double *mfile, *pfile, *rhofile, *Bfile, *Bcor, *Wfile, *Wfile1, *Wfile2, *Vfile, *Vfile1, *Vfile2;
double che(double e)
//linear interpolation
{
	double le;
	int i = 0;
	le = log10(e);
	i = j - 1;
	for (; i>0; i--)
	{
		if (le>a[i - 1])
		{
			return(pow(10, (b[i - 1] * (le - a[i]) - b[i] * (le - a[i - 1])) / (a[i - 1] - a[i])));
		}
	}
	return(0);
}

double ch(double p)
{
	double lp;
	double e, le;
	int i = 0;
	lp = log10(p);
	i = j - 1;
	for (; i>0; i--)
	{
		if (lp>b[i - 1])
		{
			le = a[i - 1] * (lp - b[i]) / (b[i - 1] - b[i]) + a[i] * (lp - b[i - 1]) / (b[i] - b[i - 1]);
			e = pow(10, le);
			return(e);
		}
	}
	return(0);
}

double fp(double r, double p, double e, double m)
{
	return(-G*(e + p / c / c)*(m + 4 * pi*r*r*r*p / c / c) / (r*r - 2 * r*m*Gc2));
}

double fm(double r, double e)
{
	return(4 * pi*r*r*e);
}

double Bf(double r, double p, double m, double B)
{
	return(2 * Gc2 / r / r*(m + 4 * pi*r*r*r*p / c / c)*B / (1 - 2 * Gc2*m / r));
}

double DH1(double r, double m, double A, double p, double e, double H1, double H0, double K, double V)
{
	return(-1 / r*(l + 1 + 2 * m*A / r + 4 * pi*r*r*A*(p - e))*H1 + 1 / r*A*(H0 + K - 16 * pi*(e + p)*V));
}

double DK(double r, double H0, double H1, double Dv, double K, double e, double p, double A, double W)
{
	return(1 / r*H0 + 0.5*l*(l + 1) / r*H1 - ((l + 1) / r - 0.5*Dv)*K - 8 * pi*(e + p)*sqrt(A) / r*W);
}

double DW(double r, double W, double A, double gamma, double p, double B, double X, double V, double H0, double K)
{
	return(-(l + 1) / r*W + r*sqrt(A)*(1 / gamma / p / sqrt(B)*X - l*(l + 1) / r / r*V + 0.5*H0 + K));
}

double DX(double r, double X, double e, double p, double B, double Dv, double H0, double w, double H1, double K, double V, double A, double F, double W)
{
	return(-l / r*X + (e + p)*sqrt(B)*(0.5*(1 / r - 0.5*Dv)*H0 + 0.5*(r*w*w / B + 0.5*l*(l + 1) / r)*H1 + 0.5*(1.5*Dv - 1 / r)*K - 0.5*l*(l + 1)*Dv / r / r*V - 1 / r*(4 * pi*(e + p)*sqrt(A) + w*w*sqrt(A) / B - 0.5*r*r*F)*W));
}

double H0f(double r, double B, double X, double m, double p, double A, double H1, double K, double w)
{
	double H0 = 8 * pi*r*r*r / sqrt(B)*X - (0.5*l*(l + 1)*(m + 4 * pi*r*r*r*p) - w*w*r*r*r / A / B)*H1 + (0.5*(l + 2)*(l - 1)*r - w*w*r*r*r / B - 1 / r*A*(m + 4 * pi*r*r*r*p)*(3 * m - r + 4 * pi*r*r*r*p))*K;
	H0 = H0 / (3 * m + 0.5*(l + 2)*(l - 1)*r + 4 * pi*r*r*r*p);
	return(H0);
}

double Vf(double r, double w, double e, double p, double B, double A, double Dp, double W, double H0, double X)
{
	return(1 / w / w / (e + p)*B*(1 / sqrt(B)*X + 1 / r*Dp / sqrt(A)*W - 0.5*(e + p)*H0));
}

double gf(double e)

{
	double le, p, gamma;
	int i = j - 1;
	e = e / G*c*c;
	le = log10(e);
	for (; i>0; i--)
	{
		if (le>a[i - 1])
		{
			p = b[i - 1] * (le - a[i]) / (a[i - 1] - a[i]) + b[i] * (le - a[i - 1]) / (a[i] - a[i - 1]);
			p = pow(10, p);
			gamma = (e + p / c / c) / e*(b[i] - b[i - 1]) / (a[i] - a[i - 1]);
			return(gamma);
		}
	}
	return((pow(10, a[0]) + pow(10, b[0]) / c / c) / e*(b[i] - b[i - 1]) / (a[i] - a[i - 1]));
}

double Ff(double r, double A, double e, double p, double m, double Dp)
{
	double F;
	F = 8 * pi / r / r*sqrt(A);
	F = F*(e + 3 * p + r*Dp - (m / 4 / pi / r / r / r + p)*(4 + A*(m / r - 4 * pi*r*r*e)));
	return(F);
}

double fmode(double e0)
{
	double tranM;
	double dr = 1.0;
	double r, r0 = 1, R, RR, drx, rx;
	double ne, p, e, m, mR, A, B = 1.0, BR, Bfactor, power, I = 0, J, DDf, Df = 0, f = 1;
	//e is  the density
	double m1, m2, m3, m4;
	double p1, p2, p3, p4;
	double B1, B2, B3, B4;
	double H1, H0, K, W, X, F, V, Dv, gamma, Dp, N;
	double DH11, H01, H02, K1, K2, x, X1, X2, Xp1, Xp2, W1, W2, V1, V2, V01, V02, V0;
	double w, wcheck, o[2], wi, aR, bR, gR, hR, kR, n, Y1, Y2, Z, DZ, DDZ, VZ, Ar1, Ar2, Ai1, Ai2, ar, ai, Ar[2], Ai[2], Br[2], Bi[2];
	int q, wpf, rpf, pfEnd;
	scr++;
	power = pow(10, b[0]);
	r = r0;
	ne = e0;
	e = e0;
	p = che(e);
	m = 1.3333333*r*r*pi*e*r;
	wpf = -1;
	int flag = 1;
	for (; p>power; r = r + dr)
	{
		if (wpf>49995)
		{
			temp_rhoc = ne;
			temp_M = 0;
			temp_R = 0;
			temp_freq = 0;
			temp_dpt = 0;
			return 1;
		}
		wpf++;
		rhofile[wpf] = e*Gc2;   /*-- e,p,m in G=c=1 --*/
		pfile[wpf] = p*Gc4;
		Bfile[wpf] = B;
		mfile[wpf] = m*Gc2;
		A = 1 / (1 - 2 * m*Gc2 / r);
		p1 = fp(r, p, e, m);
		m1 = 4 * pi*r*r*e;
		/////////////////////////////////
		if ((p + dr*p1 / 2)>power) e = ch(p + dr*p1 / 2);
		else break;
		B1 = Bf(r, p, m, B);
		p2 = fp(r + dr / 2, p + dr*p1 / 2, e, m + dr*m1 / 2);
		m2 = fm(r + dr / 2, e);
		if ((p + dr*p2 / 2)>power) e = ch(p + dr*p2 / 2); else break;
		B2 = Bf(r + dr / 2, p + dr*p1 / 2, m + dr*m1 / 2, B + dr*B1 / 2);
		p3 = fp(r + dr / 2, p + dr*p2 / 2, e, m + dr*m2 / 2);
		m3 = fm(r + dr / 2, e);
		if ((p + dr*p3)>power) e = ch(p + dr*p3); else break;
		B3 = Bf(r + dr / 2, p + dr*p2 / 2, m + dr*m2 / 2, B + dr*B2 / 2);
		p4 = fp(r + dr, p + dr*p3, e, m + dr*m3);
		m4 = fm(r + dr, e);
		B4 = Bf(r + dr, p + dr*p3, m + dr*m3, B + dr*B3);
		J = -4 * pi*(e + p / c / c)*Gc2*r*A;
		DDf = -(4 / r*Df + J*Df + 4 / r*J*f);
		I = I - 2.0 / 3 * f*J / sqrt(A*B)*r*r*r*dr;
		p = p + dr*(p1 + 2 * p2 + 2 * p3 + p4) / 6;
		e = ch(p);

		m = m + dr*(m1 + 2 * m2 + 2 * m3 + m4) / 6;
		B = B + dr*(B1 + 2 * B2 + 2 * B3 + B4) / 6;
		f = f + Df*dr;
		Df = Df + DDf*dr;
	}
	// got m and B
	pfEnd = wpf;
	R = r;
	mR = m*Gc2;
	BR = 1 - 2 * Gc2*m / r;
	Bfactor = BR / B;
	gamma = (b[1] - b[0]) / (a[1] - a[0]);
	N = 1 / (gamma - 1);
	//????????????????????????????
	RR = R - (N + 1)*(p - dr*(p1 + 2 * p2 + 2 * p3 + p4) / 6) / (p1 + 2 * p2 + 2 * p3 + p4) * 6;
	//sovle TOV???
	wpf = -1;
	rpf = -1;
	for (int n4 = 0; n4<pfEnd + 1; n4++)
	{
		Bcor[n4] = Bfile[n4] * Bfactor;
	}
	I = I / sqrt(Bfactor);
	I = I / (f + 2 * I / r / r / r) / Gc2;
	I = m*sqrt(m / I)*Gc2;
	o[0] = (-0.0047 + 0.133*I + 0.575*I*I) / mR - 0.1e-5;
	o[1] = o[0] + 0.2e-5;
	q = 1;
	wcheck = 0;
	for (int t = 0;; t++)
	{
		if (t>20)
		{
			temp_rhoc = ne;
			temp_M = 0;
			temp_R = 0;
			temp_freq = 0;
			temp_dpt = 0;
			return 2;
		}
		if (t == 0) w = o[0];
		else w = o[1];
		e = rhofile[0];
		p = pfile[0];
		B = Bcor[0];
		W = 1.0;
		K = (e + p);
		X = K*sqrt(B)*((4 * pi / 3 * (e + 3 * p) - w*w / B / l)*W + 0.5*K);
		H1 = (2 * l*K + 16 * pi*K*W) / l / (l + 1);
		rpf = -1;
		wpf = -1;
		r = r0;
		while (rpf<pfEnd)
		{
			rpf++;
			p = pfile[rpf];
			e = rhofile[rpf];
			B = Bcor[rpf];
			m = mfile[rpf];
			Dp = -(e + p)*(m + 4 * pi*r*r*r*p) / r / r / (1 - 2 * m / r);
			Dv = -2 * Dp / (e + p);
			A = 1 / (1 - 2 * m / r);
			gamma = gf(e);
			H0 = H0f(r, B, X, m, p, A, H1, K, w);
			V = Vf(r, w, e, p, B, A, Dp, W, H0, X);
			if (r == r0)V01 = V;
			if (fabs((wcheck - w) / w)<Cri)
			{
				wpf++;
				Wfile[wpf] = sqrt(1 - 2 * m / r)*W;
				Vfile[wpf] = V;
			}
			F = Ff(r, A, e, p, m, Dp);
			DH11 = DH1(r, m, A, p, e, H1, H0, K, V);
			K = K + DK(r, H0, H1, Dv, K, e, p, A, W)*dr;
			W = W + DW(r, W, A, gamma, p, B, X, V, H0, K)*dr;
			X = X + DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W)*dr;
			H1 = H1 + DH11*dr;
			r = r + dr;
		}
		wpf = -1;
		rpf = -1;
		X1 = X;
		Xp1 = DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W);
		K1 = K; H01 = H0f(r, B, X, m, p, A, H1, K, w); W1 = W; V1 = V;
		p = pfile[0];
		e = rhofile[0];
		B = Bcor[0];
		W = 1.0;
		K = -(e + p);
		X = (e + p)*sqrt(B)*((4 * pi / 3 * (e + 3 * p) - w*w / B / l)*W + 0.5*K);
		H1 = (2 * l*K + 16 * pi*(e + p)*W) / l / (l + 1);
		r = r0;
		while (rpf<pfEnd)
		{
			rpf++;
			p = pfile[rpf];
			e = rhofile[rpf];
			B = Bcor[rpf];
			m = mfile[rpf];
			Dp = -(e + p)*(m + 4 * pi*r*r*r*p) / r / r / (1 - 2 * m / r);
			Dv = -2 * Dp / (e + p);
			A = 1 / (1 - 2 * m / r);
			gamma = gf(e);
			H0 = H0f(r, B, X, m, p, A, H1, K, w);
			V = Vf(r, w, e, p, B, A, Dp, W, H0, X);
			if (r == r0)V02 = V;
			if ((wcheck - w) / w>-Cri && (wcheck - w) / w<Cri)
			{
				wpf++;
				Wfile[wpf] = sqrt(1 - 2 * m / r)*W;
				Vfile[wpf] = V;
			}
			F = Ff(r, A, e, p, m, Dp);
			DH11 = DH1(r, m, A, p, e, H1, H0, K, V);
			K = K + DK(r, H0, H1, Dv, K, e, p, A, W)*dr;
			W = W + DW(r, W, A, gamma, p, B, X, V, H0, K)*dr;
			X = X + DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W)*dr;
			H1 = H1 + DH11*dr;
			r = r + dr;
		}
		wpf = -1;
		rpf = -1;
		X2 = X;
		Xp2 = DX(r, X, e, p, B, Dv, H0, w, H1, K, V, A, F, W);
		K2 = K; H02 = H0f(r, B, X, m, p, A, H1, K, w); W2 = W; V2 = V;
		x = -(X1 - (RR - R) / (N + 1)*Xp1) / (X2 - (RR - R) / (N + 1)*Xp2);
		H0 = H01 + x*H02; K = K1 + x*K2; W = W1 + x*W2; V = V1 + x*V2; V0 = V01 + x*V02;
		if (fabs((wcheck - w) / w)<Cri)
		{
			r = r0;
			while (rpf<pfEnd)
			{
				rpf++;
				W1 = Wfile1[rpf];
				W2 = Wfile2[rpf];
				W = W1 + x*W2;
				wpf++;
				Wfile[wpf] = W / (1 + x);
				V1 = Vfile1[rpf];
				V2 = Vfile2[rpf];
				V = V1 + x*V2;
				Vfile[wpf] = V / V0;
				r = r + dr;
			}
			break;
		}
		wcheck = w;
		n = 0.5*(l - 1)*(l + 2);
		aR = -(n*R + 3 * mR) / (w*w*R*R - (n + 1)*mR / R);
		bR = (n*R*(R - 2 * mR) - w*w*R*R*R*R + mR*(R - 3 * mR));
		bR = bR / (R - 2 * mR) / (w*w*R*R - (n + 1)*mR / R);
		gR = n*(n + 1)*R*R + 3 * n*mR*R + 6 * mR*mR;
		gR = gR / R / R / (n*R + 3 * mR);
		hR = -n*R*R + 3 * n*mR*R + 3 * mR*mR;
		hR = hR / (R - 2 * mR) / (n*R + 3 * mR);
		kR = -R*R / (R - 2 * mR);
		Y1 = K;
		Y2 = aR*H0 + bR*K;
		Z = (kR*Y1 - Y2) / (kR*gR - hR);
		DZ = (gR*Y2 - hR*Y1) / (gR*kR - hR);
		if (w<2e-6)
		{
			temp_rhoc = ne;
			temp_M = 0;
			temp_M = 0;
			temp_R = 0;
			temp_freq = 0;
			temp_dpt = 0;
			return 3;
		}
		for (r = R; r<25.0 / w; r = r + dr)
		{
			drx = dr / (1 - 2 * mR / r);
			VZ = (1 - 2 * mR / r) / r / r / r / (n*r + 3 * mR) / (n*r + 3 * mR);
			VZ = VZ*(2.0*n*n*(n + 1)*r*r*r + 6.0*n*n*mR*r*r + 18.0*n*mR*mR*r + 18 * mR*mR*mR);
			DDZ = (VZ - w*w)*Z;
			Z = Z + DZ*drx;
			DZ = DZ + DDZ*drx;
		}
		r = r - dr;
		rx = r + 2 * mR*log(r / 2 / mR - 1);
		Ar1 = 2 * cos(w*rx) - 2 * (n + 1) / w / r*sin(w*rx) + 1 / w / w / r / r*(1.5*mR*w*(1 + 2 / n)*sin(w*rx) - n*(n + 1)*cos(w*rx));
		Ai1 = 2 * sin(w*rx) + 2 * (n + 1) / w / r*cos(w*rx) - 1 / w / w / r / r*(1.5*mR*w*(1 + 2 / n)*cos(w*rx) + n*(n + 1)*sin(w*rx));
		Ar2 = -2 * w*sin(w*rx) - 2 * (n + 1)*cos(w*rx) / r + 1 / w / r / r*(1.5*mR*w*(1 + 2 / n)*cos(w*rx) + n*(n + 1)*sin(w*rx)) + (1 - 2 * mR / r) * 2 * (n + 1) / w / r / r*sin(w*rx);
		Ai2 = 2 * w*cos(w*rx) - 2 * (n + 1)*sin(w*rx) / r + 1 / w / r / r*(1.5*mR*w*(1 + 2 / n)*sin(w*rx) - n*(n + 1)*cos(w*rx)) - (1 - 2 * mR / r) * 2 * (n + 1) / w / r / r*cos(w*rx);
		ar = (Ai2*Z - Ai1*DZ) / (Ar1*Ai2 - Ar2*Ai1);
		ai = -(Ar1*DZ - Ar2*Z) / (Ar1*Ai2 - Ar2*Ai1);
		if (t == 0)
		{
			Ar[t] = ar;
			Ai[t] = ai;
		}
		else{
			Ar[q] = ar;
			Ai[q] = ai;
			Br[0] = (o[0] * Ar[1] - o[1] * Ar[0]) / (o[0] - o[1]);
			Br[1] = (Ar[0] - Ar[1]) / (o[0] - o[1]);
			Bi[0] = (o[0] * Ai[1] - o[1] * Ai[0]) / (o[0] - o[1]);
			Bi[1] = (Ai[0] - Ai[1]) / (o[0] - o[1]);
			w = -(Br[0] * Br[1] + Bi[0] * Bi[1]) / (Br[1] * Br[1] + Bi[1] * Bi[1]);
			if (w <= o[0])
			{
				o[1] = o[0]; o[0] = w; Ar[1] = Ar[0]; Ai[1] = Ai[0]; q = 0;
			}
			else if (w >= o[1]){
				o[0] = o[1]; o[1] = w; Ar[0] = Ar[1]; Ai[0] = Ai[1]; q = 1;
			}
			else if ((o[1] - w)>(w - o[0])){
				o[1] = w; q = 1;
			}
			else{
				o[0] = w; q = 0;
			}
		}
	}
	wi = (Br[0] * Bi[1] - Bi[0] * Br[1]) / (Br[1] * Br[1] + Bi[1] * Bi[1]);
	temp_rhoc = ne;
	temp_M = mR / (M*Gc2);
	temp_R = RR / 1000;
	temp_freq = w*c / (2000 * pi);
	temp_dpt = 1 / (wi*c);
	return temp_M;
}
void fPrint(void)
{

	//	time_t timep = time(NULL);
	//	struct tm *sysTime = (struct tm*)malloc(sizeof(struct tm));
	//	gmtime_s(sysTime, &timep);
	//	fprintf(outf, "\n[f-mode]\tTime: %d-%d-%d-%dh-%dmin-%ds\n", 1900 + sysTime->tm_year, 1 + sysTime->tm_mon, sysTime->tm_mday, 8 + sysTime->tm_hour, sysTime->tm_min, sysTime->tm_sec);
	//	fprintf(outf, "====================================\n");
	//	fprintf(outf, "%-10s\t%-15s\t%-10s\t%-10s\t%-15s\t%-18s\n", "EoS", "rhoc(kg.m-3)", "M(Msun)", "R(km)", "frequency(kHz)", "dampingTime(s)");
}
int main(int argc, char *argv[])
{
	double *Ga, *Gb;
	int myrank;
	int numprocs;
	double *rres;
	int RankID = 0;
	double tranM;
	std::string eos = "e5585";
	std::string* Name = new string[1300000]();
	int pf, n, len1, neos, mode;
	double srhoc[40], smass[40], sradius[40], sfreq[40], sdptime[40];
	double objM, objrho;
	mode = 1;
	char pathName[50];
	double res[7];
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
	MPI_Status status;
	MPI::Init(argc, argv);
	myrank = MPI::COMM_WORLD.Get_rank();
	numprocs = MPI::COMM_WORLD.Get_size();
	len1=readEoS(Name,eos);
//	cout<<"here we are before loop";
//	getchar( );
	for (neos = myrank; neos <= len1; neos += numprocs)
	{
		sprintf(pathName, "%s/%s", eos.c_str(), Name[neos].c_str());//read eos table
		double KsymN = atol(Name[neos].substr(6, 4).c_str());
		double JsymN = atol(Name[neos].substr(10, 4).c_str());
		double J0N = atol(Name[neos].substr(14, 4).c_str());
		j = eoscan(a, b, pathName);
		Ga = (double *)malloc(j*sizeof(double));
		Gb = (double *)malloc(j*sizeof(double));
		for (pf = 0; pf<j; pf++)
		{
			Ga[pf] = a[pf] - 10.175607290470733;
			Gb[pf] = b[pf] - 27.129849799910058;
		}
	//	objrho = M2rhoc(objM, Ga, Gb,j);  //use M got the center density
	//	if (objrho == 0)
	//	{
	//		std::cout << "skip the EoS" << Name[neos].c_str() << "\n";
	//		continue;
	//	}
		FILE *outf;
		objrho = 5e17;
//		objrho = M2rho(1.4,Ga,Gb,j);
		if(objrho == 0) continue;
		outf = fopen("result_e5585.txt","a");
		while(objrho<1e18)
		{
			fmode(objrho);
			res[0] = temp_freq;
			res[1] = temp_dpt;
			res[2] = temp_M;
			res[3] = temp_R;
			res[4] = KsymN;
			res[5] = JsymN;
			res[6] = J0N;
			if(temp_M>2.5)
			{
				break;
			}
		//	rres = (double*)malloc((numprocs + 1) * 7 * sizeof(double));
//			std::cout << "\n mytank " << myrank << "    numprocs" << numprocs;
			if (a[0] != -1)
			fprintf(outf, "%e  %4.0lf %4.0lf %4.0lf %8.4lf %8.4lf %8.4lf %8.4lf\n",objrho,res[4], res[5], res[6], res[0], res[1], res[2], res[3]);
			objrho = objrho +1e16;		
//			cout<<"did we?";
//			getchar();
		}
		fclose(outf);
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
	MPI::Finalize();
	return 0;

}
