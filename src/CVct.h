#ifndef __CVCT_H__
#define __CVCT_H__

typedef double FPTYPE;

#include "CMtx.h"
// #include "CAlgError.h"

namespace MyAlgebra {
class CMtx;

class CVct {
 private:
  static const int kDefaultSize;
  FPTYPE *vector_;
  int size_;

 public:
  CVct(int size);
  ~CVct();

  const CVct &operator=(FPTYPE val);
  const CVct &operator=(const CVct &rhs);

  CVct operator-(const CVct &rhs);
  // Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
  CVct operator~();
  CVct operator*(const CMtx &rhs);
  CVct operator+(const CVct &rhs);

  FPTYPE &operator[](int ind) const;
};
}  // namespace MyAlgebra

#endif  // __CVCT_H__