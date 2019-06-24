#define _CRT_SECURE_NO_DEPRECATE
#include"Gen.h"

#define XM 939.0   // Mass of nutron???
#define XME 0.511 // Mass of electron in MeV??
#define XMMU 105.658
#define Rho0 0.16
#define H 1239.8522
static long PtLen = 0, scr = 0;
double E[80], Pre[80], Den[80];
char dysm[6] = "\\/-\\", sT[50] = "EoS_lib";

int OEos(int n,double Ks, double Js, double J0, double peos[][280])
{
	char outName[80];
	int i=1;
	sprintf(outName, ".\\%s\\KsJsJ0%+04.0lf%+04.0lf%+04.0lf.txt", sT, Ks, Js, J0);
	std::ofstream outf;
	outf.open(outName);
	if (!outf) return 0;
	while (i < n)
	{
		outf  << peos[0][i] << "    "<< peos[1][i] << "    " << peos[2][i] << std::endl;
		i++;
	}

	return 1;
}
int paraEoS(double Ksym, double Jsym, double J0, double result[][280])
{

	int  XI = 0,n;
	double EN1 = 0, XL = 58.7, U, Rho,Erho,Ec=0;
	double RhoT,delta;
	double P,PC=-100;
	double dpde;
	double ArrayE[6];
	double ErhoL;
	std::pair<double, double>p1 = GetTrho(XL, Ksym, Jsym, -100, 230);
	RhoT = p1.first;
	delta = p1.second;
	P = PSolver(RhoT, delta, ArrayE[3], ArrayE[1], ArrayE[4], &ErhoL);  
	for ( n = 1; n<73; n++)
	{
		if (RhoT - Den[n]<0)
		{
			XI = n;
			break;
		}
		//output nv+bps eos
		result[0][n] = 15.0 + log10(E[n]*1.7802222);
		result[1][n] = 32.0 + log10(Pre[n]*1.6022);
		result[2][n] = 18.0 + log10(Den[n]*1.6749286);
	}
//	std::cout << "RhoT=" << RhoT << "  at  Ksym=" << Ksym << "  Jsym=" << Jsym << "   J0=-100" << J0 << "   delta=" << delta << std::endl;
	U = RhoT / Rho0;
	while (U<10) //density change from 0.025 times saturation density to 10 times.
	{
		Rho = U*Rho0;
		delta=solver(U, XL, Ksym, Jsym, J0, 230, 0,ArrayE);
		P=PSolver(Rho, delta, ArrayE[3], ArrayE[1], ArrayE[4],&ErhoL);
		/////////////////////Etot////////////////////
		Erho = Rho*(ArrayE[0] + ArrayE[3] * delta*delta*hbar) + Rho*XM + ErhoL;
		dpde = (P - PC) / (Erho - Ec);
		if (P - Pre[XI]>0.002){
			if (Rho - RhoT>0){
				if (P - PC>0){
					if (dpde<2.0)
					{
						result[0][n] = 15.0 + log10(Erho*1.7802222);
						result[1][n] = 32.0 + log10(P*1.6022);
						result[2][n] = 18.0 + log10(Rho*1.6749286);
						n = n + 1;
						PC = P;
						Ec = Erho;
					}
				}
				else{
					break;
				}

			}
		}
		U = U + 0.01;
	}
	return n;
}
int IniBPS()
{
	//BPS+NV EoS at low density
	for (int n = 0; n<1; n++)
	{
		E[0] = 0.0;
		E[1] = 1.044e4;
		E[2] = 2.622e4;
		E[3] = 6.587e4;
		E[4] = 1.654e5;
		E[5] = 4.156e5;
		E[6] = 1.044e6;
		E[7] = 2.622e6;
		E[8] = 6.588e6;
		E[9] = 8.293e6;
		E[10] = 1.655e7;
		E[11] = 3.302e7;
		E[12] = 6.589e7;
		E[13] = 1.315e8;
		E[14] = 2.624e8;
		E[15] = 3.304e8;
		E[16] = 5.237e8;
		E[17] = 8.301e8;
		E[18] = 1.045e9;
		E[19] = 1.316e9;
		E[20] = 1.657e9;
		E[21] = 2.626e9;
		E[22] = 4.164e9;
		E[23] = 6.601e9;
		E[24] = 8.312e9;
		E[25] = 1.046e10;
		E[26] = 1.318e10;
		E[27] = 1.659e10;
		E[28] = 2.09e10;
		E[29] = 2.631e10;
		E[30] = 3.313e10;
		E[31] = 4.172e10;
		E[32] = 5.254e10;
		E[33] = 6.617e10;
		E[34] = 8.332e10;
		E[35] = 1.049e11;
		E[36] = 1.322e11;
		E[37] = 1.664e11;
		E[38] = 1.844e11;
		E[39] = 2.096e11;
		E[40] = 2.64e11;
		E[41] = 3.325e11;
		E[42] = 4.188e11;
		E[43] = 4.299e11;
		E[44] = 6.69139e11;
		E[45] = 1.0043e12;
		E[46] = 1.6749e12;
		E[47] = 6.71e12;
		E[48] = 1.3432e13;
		E[49] = 2.0161e13;
		E[50] = 2.6897e13;
		E[51] = 3.3638e13;
		E[52] = 4.0384e13;
		E[53] = 4.71341e13;
		E[54] = 5.38889e13;
		E[55] = 6.06481e13;
		E[56] = 6.74101e13;
		E[57] = 7.4176e13;
		E[58] = 8.0946e13;
		E[59] = 8.7719e13;
		E[60] = 9.4496e13;
		E[61] = 1.0128e14;
		E[62] = 1.0806e14;
		E[63] = 1.1484e14;
		E[64] = 1.2163e14;
		E[65] = 1.2843e14;
		E[66] = 1.3522e14;
		E[67] = 1.5189e14;
		E[68] = 1.6893e14;
		E[69] = 1.8601e14;
		E[70] = 2.03139e14;
		E[71] = 2.2029e14;
		E[72] = 2.37491e14;
		E[73] = 2.5475e14;
		//-------------------------------Pressure IN dynEs cm^-2----------------------------
		Pre[0] = 0.0;
		Pre[1] = 9.744e18;
		Pre[2] = 4.968e19;
		Pre[3] = 2.431e20;
		Pre[4] = 1.151e21;
		Pre[5] = 5.266e21;
		Pre[6] = 2.318e22;
		Pre[7] = 9.755e22;
		Pre[8] = 3.911e23;
		Pre[9] = 5.259e23;
		Pre[10] = 1.435e24;
		Pre[11] = 3.833e24;
		Pre[12] = 1.006e25;
		Pre[13] = 2.604e25;
		Pre[14] = 6.676e25;
		Pre[15] = 8.738e25;
		Pre[16] = 1.629e26;
		Pre[17] = 3.029e26;
		Pre[18] = 4.129e26;
		Pre[19] = 5.036e26;
		Pre[20] = 6.86e26;
		Pre[21] = 1.272e27;
		Pre[22] = 2.356e27;
		Pre[23] = 4.362e27;
		Pre[24] = 5.662e27;
		Pre[25] = 7.702e27;
		Pre[26] = 1.048e28;
		Pre[27] = 1.425e28;
		Pre[28] = 1.938e28;
		Pre[29] = 2.503e28;
		Pre[30] = 3.404e28;
		Pre[31] = 4.628e28;
		Pre[32] = 5.949e28;
		Pre[33] = 8.089e28;
		Pre[34] = 1.1e29;
		Pre[35] = 1.495e29;
		Pre[36] = 2.033e29;
		Pre[37] = 2.597e29;
		Pre[38] = 2.892e29;
		Pre[39] = 3.29e29;
		Pre[40] = 4.473e29;
		Pre[41] = 5.816e29;
		Pre[42] = 7.538e29;
		Pre[43] = 7.805e29;
		Pre[44] = 8.7392e29;
		Pre[45] = 1.2608e30;
		Pre[46] = 1.8628e30;
		Pre[47] = 6.7577e30;
		Pre[48] = 1.7832e31;
		Pre[49] = 3.29451e31;
		Pre[50] = 5.1177e31;
		Pre[51] = 7.2036e31;
		Pre[52] = 9.52441e31;
		Pre[53] = 1.2065e32;
		Pre[54] = 1.4816e32;
		Pre[55] = 1.7775e32;
		Pre[56] = 2.09411e32;
		Pre[57] = 2.43159e32;
		Pre[58] = 2.7905e32;
		Pre[59] = 3.17111e32;
		Pre[60] = 3.5742e32;
		Pre[61] = 4.00029e32;
		Pre[62] = 4.45201e32;
		Pre[63] = 4.9246e32;
		Pre[64] = 5.42439e32;
		Pre[65] = 5.95039e32;
		Pre[66] = 6.5033e32;
		Pre[67] = 1.1652e33;
		Pre[68] = 1.4865e33;
		Pre[69] = 1.86999e33;
		Pre[70] = 2.2899e33;
		Pre[71] = 2.7806e33;
		Pre[72] = 3.33511e33;
		Pre[73] = 3.9504e33;

		//-------------------------------BARYON Density IN fm^-3----------------------------
		Den[0] = 0.0;
		Den[1] = 6.295e-12;
		Den[2] = 1.581e-11;
		Den[3] = 3.972e-11;
		Den[4] = 9.976e-11;
		Den[5] = 2.506e-10;
		Den[6] = 6.294e-10;
		Den[7] = 1.581e-9;
		Den[8] = 3.972e-9;
		Den[9] = 5.000e-9;
		Den[10] = 9.976e-9;
		Den[11] = 1.990e-8;
		Den[12] = 3.972e-8;
		Den[13] = 7.924e-8;
		Den[14] = 1.581e-7;
		Den[15] = 1.990e-7;
		Den[16] = 3.155e-7;
		Den[17] = 5.000e-7;
		Den[18] = 6.294e-7;
		Den[19] = 7.924e-7;
		Den[20] = 9.976e-7;
		Den[21] = 1.581e-6;
		Den[22] = 2.506e-6;
		Den[23] = 3.972e-6;
		Den[24] = 5.000e-6;
		Den[25] = 6.294e-6;
		Den[26] = 7.924e-6;
		Den[27] = 9.976e-6;
		Den[28] = 1.256e-5;
		Den[29] = 1.581e-5;
		Den[30] = 1.990e-5;
		Den[31] = 2.506e-5;
		Den[32] = 3.155e-5;
		Den[33] = 3.972e-5;
		Den[34] = 5.000e-5;
		Den[35] = 6.294e-5;
		Den[36] = 7.924e-5;
		Den[37] = 9.976e-5;
		Den[38] = 1.105e-4;
		Den[39] = 1.256e-4;
		Den[40] = 1.581e-4;
		Den[41] = 1.990e-4;
		Den[42] = 2.506e-4;
		Den[43] = 2.572e-4;
		Den[44] = 4.000e-4;
		Den[45] = 6.000e-4;
		Den[46] = 1.000e-3;
		Den[47] = 4.000e-3;
		Den[48] = 8.000e-3;
		Den[49] = 1.200e-2;
		Den[50] = 1.600e-2;
		Den[51] = 2.000e-2;
		Den[52] = 2.400e-2;
		Den[53] = 2.800e-2;
		Den[54] = 3.200e-2;
		Den[55] = 3.600e-2;
		Den[56] = 4.000e-2;
		Den[57] = 4.400e-2;
		Den[58] = 4.800e-2;
		Den[59] = 5.200e-2;
		Den[60] = 5.600e-2;
		Den[61] = 6.000e-2;
		Den[62] = 6.400e-2;
		Den[63] = 6.800e-2;
		Den[64] = 7.200e-2;
		Den[65] = 7.600e-2;
		Den[66] = 8.000e-2;
		Den[67] = 9.000e-2;
		Den[68] = 1.000e-1;
		Den[69] = 1.100e-1;
		Den[70] = 1.200e-1;
		Den[71] = 1.300e-1;
		Den[72] = 1.400e-1;
		Den[73] = 1.500e-1;
	}
	//BPS+NV EoS at low density
	// ...........change units...........
	for (int n = 1; n <= 73; n++)
	{
		Pre[n] = Pre[n] / 1.6022e33;
		E[n] = E[n] / 1.780222222e12;
		std::cout << 32+log10(Pre[n]*1.6022) << "    " << 18.0+log10(Den[n]*1.6749286) << '\n';	
	}
	//.............change units.........
	return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(void)
{
	double Ks, Js, J0, inKs0, inKs1, inJs0, inJs1, inJ00, inJ01, dKs, dJs, dJ0;
	double peos[3][280];
	long totalpf, pf;
	IniBPS();
	/*printf("\nTips:\t-400<Ksym<100\n\t-200<Jsym<800\n\t-400<J0<400\n\nAll in the unit of MeV.\n\n");
	printf("Set range Ksym, e.g. Ksym=-400:20:100\n(which means Ksym start from -400 to 100 with intervals 20)\n\n\tKsym= ");
	fflush(stdin);
	scanf("%lf:%lf:%lf", &inKs0, &dKs, &inKs1);
	////////////////////////////////////////////
	inKs0 = ((-400 <= inKs0) && (inKs0 <= 100)) ? inKs0 : -400;
	inKs1 = ((-400 <= inKs1) && (inKs1 <= 100)) ? inKs1 : 100;
	dKs = (dKs<1) ? 20 : dKs;
	printf("Set range Jsym, e.g. Jsym=-200:50:800\n(which means Jsym start from -200 to 800 with intervals 50)\n\n\tJsym= ");
	fflush(stdin);
	scanf("%lf:%lf:%lf", &inJs0, &dJs, &inJs1);
	inJs0 = ((-200 <= inJs0) && (inJs0 <= 800)) ? inJs0 : -200;
	inJs1 = ((-200 <= inJs1) && (inJs1 <= 800)) ? inJs1 : 800;
	dJs = (dJs<1) ? 50 : dJs;
	printf("Set range J0, e.g. J0=-400:50:400\n(which means J0 start from -400 to 400 with intervals 50)\n\n\tJ0= ");
	fflush(stdin);
	scanf("%lf:%lf:%lf", &inJ00, &dJ0, &inJ01);
	inJ00 = ((-400 <= inJ00) && (inJ00 <= 400)) ? inJ00 : -400;
	inJ01 = ((-400 <= inJ01) && (inJ01 <= 400)) ? inJ01 : 400;
	dJ0 = (dJ0<1) ? 50 : dJ0;
	*/
	inKs0 = 0;
	dKs = 20;
	inKs1 = 100;
	inJs0 = 0;
	dJs = 20;
	inJs1 = 100;	
	inJ00 = 0;
	dJ0 = 20;
	inJ01 = 100;
	///////////////////My dL///////////////////////////
	/*
	printf("Set range J0, e.g. J0=-400:50:400\n(which means J0 start from -400 to 400 with intervals 50)\n\n\tJ0= ");
	scanf("%lf:%lf:%lf", &inL0, &dL, &inL1);
	inJ00 = ((-400 <= inL0) && (inL0 <= 400)) ? inL0 : -400;
	inJ01 = ((-400 <= inL1) && (inL1 <= 400)) ? inL1 : 400;
	dJ0 = (dL<1) ? 50 : dL;
	*/
	//////////////////My dK0////////////////////////////
	///////////////////////////////////////////////////
	scr = 0;
	pf = 0;
	printf("Generating 	` EoS, please wait...      %c", dysm[scr % 6]);
	totalpf=int((inJ01-inJ00)/dJ0*(inJs1-inJs0)/dJs*(inKs1-inKs0)/dKs);//总任务量
	_mkdir(sT);
	for (J0 = inJ00; J0 <= inJ01; J0 = J0 + dJ0)
	{
		for (Js = inJs0; Js <= inJs1; Js = Js + dJs)
		{
			for (Ks = inKs0; Ks <= inKs1; Ks = Ks + dKs)
			{
				scr++; 
				printf("\b\b\b\b\b\b%3ld%% %c", pf * 100 / totalpf, dysm[scr % 6]);
				int n=paraEoS(Ks, Js, J0,peos);
				if (OEos(n, Ks, Js, J0, peos) == 0)
				{
					std::cout << "cannot open/create Eos file";
					return 0;
				}
//std::cout << "  Erho=" << peos[0][i] << "  P=" << peos[1][i] << "  rho=" << peos[2][i] << "  at  Ks=" << Ks << "  Jsym=" << Js << "   J0=" << J0 << std::endl;
				pf++;
			}
		}
	}
	printf("\b\b\b\b\b\b100%%  \n");
	system("pause");
	return 0;
}