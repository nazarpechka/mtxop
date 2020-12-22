// Precompiled headers - don't use them
// #include "stdafx.h"

#include "matrix.h"
#include "vector.h"

// Reference implementation
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

#include "matrix_ref.h"
#include "vector_ref.h"

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

// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

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

// Definiujemy szablon aby łatwiej uruchamiać testy dla roznych implementacji
// klasy. Rozne implementacje będą umieszczone w roznych przestrzeniach nazw.
template <typename T>
double speedTest() {
  // Przykładowe testowe obliczenie macierzowe. Podobne obliczenia będą
  // uzywane do oceny efektywnosci implementacji w konkursie.
  const int SIZE = 100;
  const int ITER_CNT = 1000;

  T A(SIZE, SIZE, true);
  T B(SIZE, SIZE, true);
  T W(1, 1, false);
  double t1 = mygettime();

  for (int i = 0; i < ITER_CNT; i++) {
    B = ((0.1 * i) * A + B * B) * 1.e-4;
    B = -B * ~(A + B);
  }
  W = (B - A);

  double exec_time = mygettime() - t1;

  // W.display();

  return exec_time;
}

template <typename T>
void accuracyTest() {
  // Check if operations are performed correctly
  T first(5, 5, true);
  first.display();
  std::cout << '\n';
  T second(5, 5, true);
  second.display();
  std::cout << '\n';
  (first * second).display();

  return;
}

int main(void) {
  // accuracyTest<MyAlgebra::Matrix>();

#if 1
  const int kTestsAmount = 15;
  std::cout << "Matrix operations testing\n";
  double sum = 0;
  for (int i = 0; i < kTestsAmount; ++i) {
    sum += speedTest<MyAlgebra::Matrix>();
  }
  std::cout << "Average execution time: " << sum / kTestsAmount << '\n';
#endif

#if 0

  std::string command;
  while (command != "q") {
    std::cout << "> ";
    std::cin >> command;
    if (command == "test") {
      int size, iter;
      std::cout << "\tMatrices dimension = ";
      std::cin >> size;
      std::cout << "\tAmount of iterations = ";
      std::cin >> iter;

      double t_prog = test<MyAlgebra::Matrix>(size, iter);
      printf("\tExecution time:    %7.2lfs\n", t_prog);
    } else if (command == "q") {
      std::cout << "Goodbye!\n";
    } else {
      std::cout << "Unknown command!\n";
    }
  }
#endif

#if 0
	double t_ref = test<RefAlgebra::Matrix>();

	printf( "Czas wykonania referencyjny: %7.2lfs\n", t_ref );
	printf( "Czas wykonania testowany:    %7.2lfs\n", t_prog );
	printf("Wsp�czynnik przyspieszenia Q: %5.2lf", t_ref / t_prog);
#endif
}