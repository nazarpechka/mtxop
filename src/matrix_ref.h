#ifndef __MATRIX_REF_H__
#define __MATRIX_REF_H__

#include <stddef.h>

#include "vector_ref.h"

namespace RefAlgebra {
class Vector;

class Matrix {
 public:
  static const FPTYPE ALG_PRECISION;

  // =========================================================================
  // KONSTRUKTORY:
  // =========================================================================

  // Tworzy macierz z mozliwoscią losowej inicjalizacji
  // Removed default value for rand_init to avoid confusion with the next
  // constructor
  Matrix(size_t row_cnt, size_t col_cnt, bool rand_init);

  // Tworzy kwadratową macierz diagonalną
  Matrix(size_t row_cnt, FPTYPE diagonal);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other);

  ~Matrix();

  // =========================================================================
  // OPERATORY PRZYPISANIA:
  // =========================================================================

  const Matrix &operator=(const Matrix &other);

  // Operator przenoszący
  const Matrix &operator=(Matrix &&other);

  // Zamiana macierzy na macierz diagonalną
  const Matrix &operator=(const FPTYPE diagonal);

  // =========================================================================
  // OPERACJE ALGEBRAICZNE
  // =========================================================================

  Matrix operator+(const Matrix &other) const;

  Matrix operator-(const Matrix &other) const;

  // Minus unarny - zmiana znaku wszystkich wspołczynnikow macierzy
  Matrix operator-() const;

  Matrix operator*(const Matrix &other) const;

  // Mnozenie macierzy przez wektor, other musi być wektorem kolumnowym
  Vector operator*(const Vector &other) const;

  // Mnozenie macierzy przez stałą
  Matrix operator*(FPTYPE multiplier) const;

  // Transponowanie macierzy
  Matrix operator~() const;

  // Akceptuje tylko power >= -1:
  //    power = -1 - zwraca macierz odwroconą
  //    power = 0  - zwraca macierz jednostkową
  //    power = 1  - zwraca kopię macierzy
  //    power > 1  - zwraca iloczyn macierzy
  Matrix operator^(int power) const;

  // =========================================================================
  // INDEKSOWANIE MACIERZY
  // =========================================================================

  FPTYPE *operator[](int row_ind);

  int getRowCount() const;
  int getColCount() const;

  // Porownywanie macierzy z dokładnoscią do stałej ALG_PRECISION
  bool operator==(const Matrix &other) const;

  // Tylko do celow testowych - wypisuje macierz wierszami na stdout
  void display() const;

  // friend Matrix operator*( FPTYPE multiplier, const Matrix &other );

 private:
  int m_row_cnt;
  int m_col_cnt;
  FPTYPE **m_array;

  void copy(const Matrix &other);
  void move(Matrix &&other);
};

Matrix operator*(FPTYPE multiplier, const Matrix &other);
}  // namespace RefAlgebra

#endif  // __MATRIX_REF_H__