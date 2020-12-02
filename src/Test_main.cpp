// Precompiled headers - don't use them
// #include "stdafx.h"

#include "CMtx.h"
#include "CVct.h"

// Reference implementation in those files
// #include "CVctRef.h"
// #include "CMtxRef.h"

#include <stdio.h>
#include <stdlib.h>

#include <iostream>

// ===================================================================
// FUNKCJE DO POMIARU CZASU
// ===================================================================

#ifdef _WIN32
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#include <math.h>
#include <time.h>

double mygettime(void) {
#ifdef _WIN32
  struct _timeb tb;
  _ftime(&tb);
  return (double)tb.time + (0.001 * (double)tb.millitm);
#else
  struct timeval tv;
  if (gettimeofday(&tv, 0) < 0) {
    perror("oops");
  }
  return (double)tv.tv_sec + (0.000001 * (double)tv.tv_usec);
#endif
}

// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

// Definiujemy szablon aby łatwiej uruchamiać testy dla roznych implementacji
// klasy. Rozne implementacje będą umieszczone w roznych przestrzeniach nazw.
template <typename T>
double test() {
  // Przykładowe testowe obliczenie macierzowe. Podobne obliczenia będą
  // uzywane do oceny efektywnosci implementacji w konkursie.
  const int SIZE = 100;
  const int ITER_CNT = 100;

  T A(SIZE, SIZE, true);
  // A.display();
  T B(SIZE, SIZE, true);
  // B.display();
  T W(1, 1);
  // W.display();
  double t1 = mygettime();

  for (int i = 0; i < ITER_CNT; i++) {
    B = ((0.1 * i) * A + B * B) * 1.e-4;
    B = -B * ~(A + B);
    // B.display();
    // std::cout << "\n";
  }
  W = (B - A);
  // W.display();

  double exec_time = mygettime() - t1;

  // W.display();

  return exec_time;
}

int main(int argc, char* argv[]) {
  double t_prog = test<MyAlgebra::CMtx>();
  printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);
  // MyAlgebra::CMtx matrix(3, 3, true);
  // MyAlgebra::CMtx second_matrix(3, 3, true);
  // matrix.display();
  // std::cout << '\n';
  // second_matrix.display();
  // std::cout << '\n';
  // MyAlgebra::CMtx tr_matrix = matrix * second_matrix;
  // tr_matrix.display();

#if 0
	double t_ref = test<RefAlgebra::CMtx>();

	printf( "Czas wykonania referencyjny: %7.2lfs\n", t_ref );
	printf( "Czas wykonania testowany:    %7.2lfs\n", t_prog );
	printf("Wsp�czynnik przyspieszenia Q: %5.2lf", t_ref / t_prog);
#endif
}