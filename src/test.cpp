// Precompiled headers - don't use them
// #include "stdafx.h"

#include "matrix.h"
#include "vector.h"

// Reference implementation
#include <math.h>
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

bool operator==(RefAlgebra::Matrix &first, MyAlgebra::Matrix &second) {
  if (first.getRowCount() != second.getRowCount() ||
      first.getColCount() != second.getColCount())
    return false;

  for (int i = 0; i < first.getRowCount(); ++i) {
    for (int j = 0; j < first.getColCount(); ++j) {
      if (abs(first[i][j] - second[i][j]) > 10e-6) return false;
    }
  }

  return true;
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
  const int ROW_CNT = 500;
  const int COL_CNT = 500;

  srand(123);
  RefAlgebra::Matrix first_ref(ROW_CNT, COL_CNT, true);
  RefAlgebra::Matrix second_ref(ROW_CNT, COL_CNT, true);

  srand(123);
  T first(ROW_CNT, COL_CNT, true);
  T second(ROW_CNT, COL_CNT, true);

  std::cout << std::boolalpha;

  std::cout << "Equality test:\n";
  std::cout << (first_ref == first) << "\n";
  std::cout << (second_ref == second) << "\n";

  std::cout << "Multiplication test:\n";
  RefAlgebra::Matrix mult_ref(first_ref * second_ref);
  T mult(first * second);
  std::cout << (mult_ref == mult) << "\n";

  std::cout << "Addition test:\n";
  RefAlgebra::Matrix add_ref(first_ref + second_ref);
  T add(first + second);
  std::cout << (add_ref == add) << "\n";

  return;
}

int main(void) {
  // accuracyTest<MyAlgebra::Matrix>();

#if 1
  std::cout << "Matrix operations testing\n";
  const int TEST_AMOUNT = 25;

  double t_prog = 0;
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_prog += speedTest<MyAlgebra::Matrix>();
  }
  t_prog /= TEST_AMOUNT;
  printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);

  double t_ref = 0;
#if 0
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_ref += speedTest<RefAlgebra::Matrix>();
  }
  t_ref /= TEST_AMOUNT;
#else
  t_ref = 2.34;
#endif

  printf("Czas wykonania referencyjny: %7.2lfs\n", t_ref);

  printf("Wspolczynnik przyspieszenia Q: %5.2lf\n", t_ref / t_prog);

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
}