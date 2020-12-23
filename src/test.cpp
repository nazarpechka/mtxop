// Precompiled headers - don't use them
// #include "stdafx.h"

#define USE_EIGEN 0

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iomanip>
#include <iostream>
#include <string>

#if USE_EIGEN
#include "Eigen/Dense"
#endif
#include "matrix.h"
#include "vector.h"

// Reference implementation
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

#include <time.h>

#include <cmath>

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
      if (std::abs(first[i][j] - second[i][j]) >
          RefAlgebra::Matrix::ALG_PRECISION)
        return false;
    }
  }

  return true;
}

#if USE_EIGEN
bool operator==(RefAlgebra::Matrix &first, Eigen::MatrixXd &second) {
  if (first.getRowCount() != second.rows() ||
      first.getColCount() != second.cols())
    return false;

  for (int i = 0; i < first.getRowCount(); ++i) {
    for (int j = 0; j < first.getColCount(); ++j) {
      if (std::abs(first[i][j] - second(i, j)) >
          RefAlgebra::Matrix::ALG_PRECISION)
        return false;
    }
  }

  return true;
}

#endif

// Definiujemy szablon aby łatwiej uruchamiać testy dla roznych implementacji
// klasy. Rozne implementacje będą umieszczone w roznych przestrzeniach nazw.
template <typename T>
double speedTest() {
  // Przykładowe testowe obliczenie macierzowe. Podobne obliczenia będą
  // uzywane do oceny efektywnosci implementacji w konkursie.
  srand(time(0));
  const int SIZE = 10;
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
double multSpeedTest() {
  // Przykładowe testowe obliczenie macierzowe na mnozenie.
  srand(time(NULL));
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

template <typename T>
// Check if operations are performed correctly
void selfTest() {
  uint32_t rand_time = time(0);
  srand(rand_time);
  const int FIRST_ROW_CNT = rand() % 1000 + 1;
  const int FIRST_COL_CNT = rand() % 1000 + 1;
  const int SECOND_ROW_CNT = FIRST_COL_CNT;
  const int SECOND_COL_CNT = rand() % 1000 + 1;

  std::cout << "Accuracy test started\n";
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

  std::cout << "speedTest() accuracy test:\n";
  const int SIZE = 10;
  const int ITER_CNT = 1000;

  srand(rand_time);
  RefAlgebra::Matrix A_ref(SIZE, SIZE, true);
  RefAlgebra::Matrix B_ref(SIZE, SIZE, true);
  RefAlgebra::Matrix W_ref(1, 1, false);

  srand(rand_time);
  T A(SIZE, SIZE, true);
  T B(SIZE, SIZE, true);
  T W(1, 1, false);

#if USE_EIGEN
  Eigen::MatrixXd A_e(SIZE, SIZE);
  Eigen::MatrixXd B_e(SIZE, SIZE);
  Eigen::MatrixXd W_e(1, 1);

  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      A_e(i, j) = A_ref[i][j];
      B_e(i, j) = B_ref[i][j];
    }
  }
#endif

  if (A_ref == A && B_ref == B && W_ref == W
#if USE_EIGEN
      && A_ref == A_e && B_ref == B_e && W_ref == W_e
#endif
  ) {
    std::cout << "Input matrices are identical\n";
  } else {
    std::cout << "Input matrices are NOT identical\n";
  }

  for (int i = 0; i < ITER_CNT; i++) {
    B_ref = ((0.1 * i) * A_ref + B_ref * B_ref) * 1.e-4;
    B_ref = -B_ref * ~(A_ref + B_ref);
  }
  W_ref = (B_ref - A_ref);

  for (int i = 0; i < ITER_CNT; i++) {
    B = ((0.1 * i) * A + B * B) * 1.e-4;
    B = -B * ~(A + B);
  }
  W = (B - A);

#if USE_EIGEN
  for (int i = 0; i < ITER_CNT; i++) {
    B_e = ((0.1 * i) * A_e + B_e * B_e) * 1.e-4;
    B_e = -B_e * (A_e + B_e).transpose();
  }
  W_e = (B_e - A_e);
#endif

  if (W_ref == W
#if USE_EIGEN
      && W_ref == W_e
#endif
  ) {
    std::cout << "Result matrices are identical\n";

  } else {
    std::cout << "Result matrices are NOT identical, see differences:\n";
#if USE_EIGEN
    std::cout << "W_ref == W_e: " << (W_ref == W_e) << '\n';
#endif
    std::cout << "W_ref == W: " << (W_ref == W) << '\n';

    W_ref.display();
    W.display();
#if USE_EIGEN
    std::cout << W_e << '\n';
#endif
  }

  return;
}

int main(void) {
  selfTest<MyAlgebra::Matrix>();

#if 0
  std::cout << "Matrix operations testing\n";
  const int TEST_AMOUNT = 1;

  double t_prog = 0;
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_prog += speedTest<MyAlgebra::Matrix>();
    std::cout << "Test #" << i + 1 << std::endl;
  }
  std::cout << '\n';
  t_prog /= TEST_AMOUNT;
  printf("Czas wykonania testowany:    %7.2lfs\n", t_prog);

  double t_ref = 0;
#if 1
  for (int i = 0; i < TEST_AMOUNT; ++i) {
    t_ref += speedTest<RefAlgebra::Matrix>();
    std::cout << "Test #" << i + 1 << std::endl;
  }
  std::cout << '\n';
  t_ref /= TEST_AMOUNT;
#else
  t_ref = 2.35;
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