fmode:main.o fio.o m2rho.o
	mpicxx -o e5585f main.o fio.o m2rho.o
fmain.o: main.cpp fmode.h
	mpicxx -c main.cpp
fio.o: fio.cpp fmode.h
	mpicxx -c fio.cpp
m2rho.o: m2rho.cpp fmode.h
	mpicxx -c m2rho.cpp
clean:
	rm e5585f main.o fio.o m2rho.o
