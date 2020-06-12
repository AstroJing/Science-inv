pgenEoS : Gen.o imple.o
	mpicxx -o e5585 Gen.o imple.o
Gen.o:Gen.cpp Gen.h
	mpicxx -c Gen.cpp
imple.o: imple.cpp Gen.h
	mpicxx -c imple.cpp
clean:
	rm genEoS Gen.o imple.o
