// Precompiled headers - don't use them
// #include "stdafx.h"

#define SELF_TEST 0
#define BENCHMARK 1

#include <stdlib.h>
#include <stdio.h>

#include <iomanip>
#include <iostream>
#include <immintrin.h>

// My implementation
#include "matrix.h"

// Reference implementation
#include "matrix_ref.h"

// ===================================================================
// FUNKCJE DO POMIARU CZASU
// ===================================================================

#ifdef _WIN32
#include <sys/timeb.h>
#else
#include <sys/time.h>
#endif

#include <time.h>
#include <cmath>

// ===================================================================
// FUNKCJA OCENY CZASU WYKONANIA
// ===================================================================

double mygettime() {
#ifdef _WIN32
  struct _timeb tb;
  _ftime(&tb);
  return (double)tb.time + (0.001 * (double)tb.millitm);
#else
  struct timeval tv;
  if (gettimeofday(&tv, nullptr) < 0) {
    perror("oops");
  }
  return static_cast<double>(tv.tv_sec) +
         (0.000001 * static_cast<double>(tv.tv_usec));
#endif
}

// ===================================================================
// Metody testujące poprawność działań
// ===================================================================

// For testing, compare different matrices
bool operator==(RefAlgebra::Matrix &first, MyAlgebra::Matrix &second) {
  if (first.getRowCount() != second.getRowCount() ||
      first.getColCount() != second.getColCount())
    return false;

  for (size_t i = 0; i < first.getRowCount(); ++i) {
    for (size_t j = 0; j < first.getColCount(); ++j) {
      if (std::abs(first[i][j] - second[i][j]) >
          RefAlgebra::Matrix::ALG_PRECISION)
        return false;
    }
  }

  return true;
}

template <typename T>
void selfTestGeneral() {
  const uint32_t rand_time = time(0);

  srand(rand_time);
  const size_t FIRST_ROW_CNT = rand() % 1000 + 1;
  const size_t FIRST_COL_CNT = rand() % 1000 + 1;
  const size_t SECOND_ROW_CNT = FIRST_COL_CNT;
  const size_t SECOND_COL_CNT = rand() % 1000 + 1;

  std::cout << "General self test started\n";
  std::cout << "Creating first_sq  = " << FIRST_ROW_CNT << "x" << FIRST_ROW_CNT
            << '\n';
  std::cout << "         second_sq = " << FIRST_ROW_CNT << "x" << FIRST_ROW_CNT
            << '\n';
  std::cout << "         first     = " << FIRST_ROW_CNT << "x" << FIRST_COL_CNT
            << '\n';
  std::cout << "         second    = " << SECOND_ROW_CNT << "x"
            << SECOND_COL_CNT << '\n';

  srand(rand_time);
  RefAlgebra::Matrix first_sq_ref(FIRST_ROW_CNT, FIRST_ROW_CNT, true);
  RefAlgebra::Matrix second_sq_ref(FIRST_ROW_CNT, FIRST_ROW_CNT, true);

  RefAlgebra::Matrix first_ref(FIRST_ROW_CNT, FIRST_COL_CNT, true);
  RefAlgebra::Matrix second_ref(SECOND_ROW_CNT, SECOND_COL_CNT, true);

  srand(rand_time);
  T first_sq(FIRST_ROW_CNT, FIRST_ROW_CNT, true);
  T second_sq(FIRST_ROW_CNT, FIRST_ROW_CNT, true);

  T first(FIRST_ROW_CNT, FIRST_COL_CNT, true);
  T second(SECOND_ROW_CNT, SECOND_COL_CNT, true);

  std::cout << std::boolalpha;

  std::cout << "Equality test:\n";
  std::cout << (first_sq_ref == first_sq) << "\n";
  std::cout << (second_sq_ref == second_sq) << "\n";

  std::cout << (first_ref == first) << "\n";
  std::cout << (second_ref == second) << "\n";

  std::cout << "Multiplication test:\n";

  RefAlgebra::Matrix mult_sq_ref(first_sq_ref * second_sq_ref);
  T mult_sq(first_sq * second_sq);
  std::cout << (mult_sq_ref == mult_sq) << "\n";

  RefAlgebra::Matrix mult_sq_ref_rev(second_sq_ref * first_sq_ref);
  T mult_sq_rev(second_sq * first_sq);
  std::cout << (mult_sq_ref_rev == mult_sq_rev) << "\n";

  RefAlgebra::Matrix mult_ref(first_ref * second_ref);
  T mult(first * second);
  std::cout << (mult_ref == mult) << "\n";

  std::cout << "Addition test:\n";
  RefAlgebra::Matrix add_ref(first_sq_ref + second_sq_ref);
  T add(first_sq + second_sq);
  std::cout << (add_ref == add) << "\n";

  RefAlgebra::Matrix add_ref_rev(second_sq_ref + first_sq_ref);
  T add_rev(second_sq + first_sq);
  std::cout << (add_ref_rev == add_rev) << "\n";
}

template <typename T>
void selfTestSpeedTest() {
  const uint32_t rand_time = time(0);
  std::cout << "\nspeedTest() test started\n";

  const int SIZE = 100;
  const int ITER_CNT = 1000;

  srand(rand_time);
  RefAlgebra::Matrix A_ref(SIZE, SIZE, true);
  RefAlgebra::Matrix B_ref(SIZE, SIZE, true);
  RefAlgebra::Matrix W_ref(1, 1, false);

  srand(rand_time);
  T A(SIZE, SIZE, true);
  T B(SIZE, SIZE, true);
  T W(1, 1, false);

  if (A_ref == A && B_ref == B && W_ref == W) {
    std::cout << "Input matrices are identical\n";
  } else {
    std::cout << "Input matrices are NOT identical, see differences:\n";
    std::cout << "A_ref == A: " << (A_ref == A) << '\n';
    std::cout << "B_ref == B: " << (B_ref == B) << '\n';
    std::cout << "W_ref == W: " << (W_ref == W) << '\n';
  }

  for (int i = 0; i < ITER_CNT; i++) {
    B_ref = ((0.1f * i) * A_ref + B_ref * B_ref) * 1.e-4f;
    B_ref = -B_ref * ~(A_ref + B_ref);
  }
  W_ref = (B_ref - A_ref);

  for (int i = 0; i < ITER_CNT; i++) {
    B = ((0.1f * i) * A + B * B) * 1.e-4f;
    B = -B * ~(A + B);
  }
  W = (B - A);

  if (W_ref == W) {
    std::cout << "Result matrices are identical\n";

  } else {
    std::cout << "Result matrices are NOT identical, see differences:\n";
    std::cout << "W_ref == W: " << (W_ref == W) << '\n';

    W_ref.display();
    W.display();
  }
}

template <typename T>
// Check if operations are performed correctly
void selfTest() {
  selfTestGeneral<T>();
  selfTestSpeedTest<T>();
}

// Definiujemy szablon aby łatwiej uruchamiać testy dla roznych implementacji
// klasy. Rozne implementacje będą umieszczone w roznych przestrzeniach nazw.
template <typename T>
double speedTest() {
  // Przykładowe testowe obliczenie macierzowe. Podobne obliczenia będą
  // uzywane do oceny efektywnosci implementacji w konkursie.
  srand(time(0));
  const int SIZE = 100;
  const int ITER_CNT = 1000;

  T A(SIZE, SIZE, true);
  T B(SIZE, SIZE, true);
  T W(1, 1, false);
  double t1 = mygettime();

  for (int i = 0; i < ITER_CNT; i++) {
    B = ((0.1 * i) * A + B * B) * 1.e-4f;
    B = -B * ~(A + B);
  }
  W = (B - A);

  double exec_time = mygettime() - t1;

  // W.display();

  return exec_time;
}

template <typename T>
double multSpeedTest() {
  // Przykładowe testowe obliczenie macierzowe na mnozenie.
  srand(time(0));
  const int SIZE = 1024;
  const int ITER_CNT = 10;

  T A(SIZE, SIZE, true);
  T B(SIZE, SIZE, true);
  T W(1, 1, false);
  double t1 = mygettime();

  for (int i = 0; i < ITER_CNT; i++) {
    W = A * B;
  }

  double exec_time = mygettime() - t1;

  // W.display();

  return exec_time;
}


int main() {

#if SELF_TEST
  selfTest<MyAlgebra::Matrix>();
#endif

#if BENCHMARK
  std::cout << "\nMatrix operations benchmark\n";
  const int TEST_AMOUNT = 15;

  double t_prog = 0;
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_prog += multSpeedTest<MyAlgebra::Matrix>();
    std::cout << "Test #" << i + 1 << std::endl;
  }
  std::cout << '\n';
  t_prog /= TEST_AMOUNT;
  printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);


  double t_ref = 0;
#if 0
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_ref += speedTest<RefAlgebra::Matrix>();
    std::cout << "Test #" << i + 1 << std::endl;
  }
  std::cout << '\n';
  t_ref /= TEST_AMOUNT;
#else
  t_ref = 2.23;
#endif

  printf("Czas wykonania referencyjny: %7.2lfs\n", t_ref);
  printf("Wspolczynnik przyspieszenia Q: %5.2lf\n", t_ref / t_prog);

#endif
}