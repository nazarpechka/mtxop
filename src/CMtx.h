#ifndef __CMTX_H__
#define __CMTX_H__

#include <stddef.h>

#include "CVct.h"

namespace MyAlgebra {
class CVct;

class CMtx {
 public:
  static const FPTYPE ALG_PRECISION;

  // =========================================================================
  // KONSTRUKTORY:
  // =========================================================================

  // Tworzy macierz z mozliwoscią losowej inicjalizacji
  // Removed default value for rand_init to avoid confusion with the next
  // constructor
  CMtx(size_t row_cnt, size_t col_cnt, bool rand_init);

  // Tworzy kwadratową macierz diagonalną
  CMtx(size_t row_cnt, FPTYPE diagonal);

  CMtx(const CMtx &other);

  CMtx(CMtx &&other);

  ~CMtx();

  // =========================================================================
  // OPERATORY PRZYPISANIA:
  // =========================================================================

  const CMtx &operator=(const CMtx &other);

  // Zamiana macierzy na macierz diagonalną
  const CMtx &operator=(const FPTYPE diagonal);

  // Operator przenoszący
  const CMtx &operator=(CMtx &&other);

  // =========================================================================
  // OPERACJE ALGEBRAICZNE
  // =========================================================================

  // Minus unarny - zmiana znaku wszystkich wspołczynnikow macierzy
  CMtx operator-() const;

  CMtx operator-(const CMtx &other) const;

  // Transponowanie macierzy
  CMtx operator~() const;

  // Mnozenie macierzy przez wektor, other musi być wektorem kolumnowym
  CVct operator*(const CVct &other) const;

  CMtx operator*(const CMtx &other) const;

  // Mnozenie macierzy przez stałą
  CMtx operator*(FPTYPE multiplier) const;

  CMtx operator+(const CMtx &other) const;

  // Akceptuje tylko power >= -1:
  //    power = -1 - zwraca macierz odwroconą
  //    power = 0  - zwraca macierz jednostkową
  //    power = 1  - zwraca kopię macierzy
  //    power > 1  - zwraca iloczyn macierzy
  CMtx operator^(int power) const;

  // Wylicza determinant maciezy
  FPTYPE determinant() const;

  // =========================================================================
  // INDEKSOWANIE MACIERZY
  // =========================================================================

  FPTYPE *operator[](int row_ind);

  // Porownywanie macierzy z dokładnoscią do stałej ALG_PRECISION
  // TODO: Should it really be &&other?
  bool operator==(const CMtx &other) const;

  // Tylko do celow testowych - wypisuje macierz wierszami na stdout
  void display() const;

  // friend CMtx operator*( FPTYPE multiplier, const CMtx &other );

 private:
  int m_row_cnt;
  int m_col_cnt;
  FPTYPE *m_array;
};

CMtx operator*(FPTYPE multiplier, const CMtx &other);
}  // namespace MyAlgebra

#endif  // __CMTX_H__