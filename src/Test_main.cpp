// test-1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "CVct.h"
#include "CMtx.h"

// W poni¿szych plikach - impelementacja referencyjna
// #include "CVctRef.h"
// #include "CMtxRef.h"

#include <stdlib.h>
#include <stdio.h>


// ===================================================================
// FUNKCJE DO POMIARU CZASU
// ===================================================================

#ifdef _WIN32
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif
#include <time.h>
#include<math.h>

double mygettime(void) {
# ifdef _WIN32
    struct _timeb tb;
    _ftime(&tb);
    return (double)tb.time + (0.001 * (double)tb.millitm);
# else
    struct timeval tv;
    if(gettimeofday(&tv, 0) < 0) {
      perror("oops");
    }
    return (double)tv.tv_sec + (0.000001 * (double)tv.tv_usec);
# endif
}
 
// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

// Definiujemy szablon aby ³atwiej uruchamiaæ testy dla róznych implementacji
// klasy. Ró¿ne implementacje bêd¹ umieszczone w ró¿nych przestrzeniach nazw.
template<typename T>
double test()
{
	// Przyk³adowe testowe obliczenie macierzowe. Podobne obliczenia bêd¹ 
	// u¿ywane do oceny efektywnoœci implementacji w konkursie.
	const int SIZE = 100;
	const int ITER_CNT = 10;

	T A(SIZE, SIZE, true);
	T B(SIZE, SIZE, true);
	T W(1, 1);
	double t1 = mygettime();

	for (int i = 0; i < ITER_CNT; i++)
	{
		B = ((0.1 * i)*A + B * B) * 1.e-4;
		B = -B * ~(A + B);
	}
	W = (B - A);

	double exec_time = mygettime() - t1;

	W.display();

	return exec_time;
}


int _tmain(int argc, char* argv[])
{	
	double t_prog = test<MyAlgebra::CMtx>();
#if 0
	double t_ref = test<RefAlgebra::CMtx>();

	printf( "Czas wykonania referencyjny: %7.2lfs\n", t_ref );
	printf( "Czas wykonania testowany:    %7.2lfs\n", t_prog );
	printf("Wspó³czynnik przyspieszenia Q: %5.2lf", t_ref / t_prog);
#endif
}