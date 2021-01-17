#ifndef __CMTX_H__
#define __CMTX_H__

#include <cstddef>
#include <string>

namespace MyAlgebra {

class Matrix {
 public:
  static const float ALG_PRECISION;

  // =========================================================================
  // KONSTRUKTORY:
  // =========================================================================

  // Tworzy macierz z mozliwoscią losowej inicjalizacji
  // Removed default value for rand_init to avoid confusion with the next
  // constructor
  Matrix(size_t row_cnt, size_t col_cnt, bool rand_init = false);

  // Tworzy kwadratową macierz diagonalną
  Matrix(size_t row_cnt, float diagonal);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other);

  ~Matrix();

  // =========================================================================
  // OPERATORY PRZYPISANIA:
  // =========================================================================

  const Matrix &operator=(const Matrix &other);

  // Operator przenoszący
  const Matrix &operator=(Matrix &&other) noexcept;

  // Zamiana macierzy na macierz diagonalną
  const Matrix &operator=(float diagonal);

  // =========================================================================
  // OPERACJE ALGEBRAICZNE
  // =========================================================================

  // Porownywanie macierzy z dokładnoscią do stałej ALG_PRECISION
  bool operator==(const Matrix &other) const;

  Matrix operator+(const Matrix &other) const;

  Matrix operator-(const Matrix &other) const;

  // Minus unarny - zmiana znaku wszystkich wspołczynnikow macierzy
  Matrix operator-() const;

  Matrix operator*(const Matrix &other);

  // Mnozenie macierzy przez stałą
  Matrix operator*(float multiplier) const;

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

  float *operator[](size_t row_ind);
  const float *operator[](size_t row_ind) const;

  size_t getRowCount() const;
  size_t getColCount() const;

  // Tylko do celow testowych - wypisuje macierz wierszami na stdout
  void display() const;

  static std::string authorName();

  // friend Matrix operator*( float multiplier, const Matrix &other );

 private:
  size_t m_row_cnt;
  size_t m_col_cnt;
  float *m_array;

  void multiply(const Matrix &res, const Matrix &other, size_t start,
                size_t end) const;

  void copy(const Matrix &other);
  void move(Matrix &&other);
};

Matrix operator*(float multiplier, const Matrix &other);
}  // namespace MyAlgebra

#endif  // __CMTX_H__