#ifndef __CMTX_H__
#define __CMTX_H__

#include "CVct.h"

namespace MyAlgebra {
class CVct;

class CMtx {
 private:
  int row_cnt_;
  int col_cnt_;
  FPTYPE **row_ptr_;

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

  CMtx(const CMtx &rhs);

  // Jezeli potrzeba - nalezy zadeklarować i zaimplementować inne konstruktory
  ~CMtx();

  // =========================================================================
  // OPERATORY PRZYPISANIA:
  // =========================================================================

  const CMtx &operator=(const CMtx &rhs);

  // Zamiana macierzy na macierz diagonalną
  const CMtx &operator=(const FPTYPE diagonal);

  // Operator przenoszący
  //   const CMtx &operator=(CMtx &&rhs);

  // =========================================================================
  // OPERACJE ALGEBRAICZNE
  // =========================================================================

  // Minus unarny - zmiana znaku wszystkich wspołczynnikow macierzy
  CMtx operator-() const;

  CMtx operator-(const CMtx &rhs) const;

  // Transponowanie macierzy
  CMtx operator~() const;

  // Mnozenie macierzy przez wektor, rhs musi być wektorem kolumnowym
  CVct operator*(const CVct &rhs) const;

  CMtx operator*(const CMtx &rhs) const;

  // Mnozenie macierzy przez stałą
  CMtx operator*(FPTYPE multiplier) const;

  CMtx operator+(const CMtx &rhs) const;

  // Akceptuje tylko power >= -1:
  //    power = -1 - zwraca macierz odwroconą
  //    power = 0  - zwraca macierz jednostkową
  //    power = 1  - zwraca kopię macierzy
  //    power > 1  - zwraca iloczyn macierzy
  CMtx operator^(int power) const;

  // =========================================================================
  // INDEKSOWANIE MACIERZY
  // =========================================================================

  FPTYPE *operator[](int row_ind);

  // Porownywanie macierzy z dokładnoscią do stałej ALG_PRECISION
  bool operator==(const CMtx &&rhs) const;

  // Tylko do celow testowych - wypisuje macierz wierszami na stdout
  void display() const;

  // friend CMtx operator*( FPTYPE multiplier, const CMtx &rhs );
};

CMtx operator*(FPTYPE multiplier, const CMtx &rhs);
}  // namespace MyAlgebra

#endif  // __CMTX_H__