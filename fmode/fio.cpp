#include"fmode.h"
int eoscan(double a[], double b[], std::string pathName)
{
	FILE *inf;
	int j = 0;
	double cccc;
	inf = fopen(pathName.c_str(), "r");
	if (inf == NULL)
	{
		std::cout << '\n' << "Cannot load" << pathName << '\n';
		return 0;

	}
	while (fscanf(inf, "%lf", a + j) == 1)
	{
		fscanf(inf, "%lf", b + j);
		      fscanf(inf,"%lf",&cccc);
		//      std::cout << a[j] << "       " << b[j] << "\n";
		j++;
	}
	fclose(inf);
	return j;
}
int errd(int n)
{
	switch (n)
	{
	case 0:
		return 0;
		break;
	default:
		break;
	}
}
int readEoS(std::string *name,std::string eos)
{
	int n = 0;
	std::string open = "ls "+ eos+"/> temp_dir";
	remove("temp_dir");
	system(open.c_str());
	std::ifstream eostable("temp_dir");
	if (!eostable.is_open())
	{
		std::cout << "No Equation of State files detected, automatically quit.\n";
		return 0;
	}
	char Eater[200];
	while (eostable.getline(Eater, 100))
	{
		name[n] = Eater;
		//	std::cout <<n<<"   "<< Name[n]<<'\n';
		n++;
	}
	//	}
	return n - 1;
}
