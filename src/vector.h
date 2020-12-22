#ifndef __CVCT_H__
#define __CVCT_H__

typedef double FPTYPE;

#include <stddef.h>

#include "matrix.h"
// #include "CAlgError.h"

namespace MyAlgebra {
class Matrix;

class Vector {
 public:
  Vector(size_t size);
  ~Vector();

  const Vector &operator=(FPTYPE val);
  const Vector &operator=(const Vector &rhs);

  Vector operator-(const Vector &rhs);
  // Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
  Vector operator~();
  Vector operator*(const Matrix &rhs);
  Vector operator+(const Vector &rhs);

  FPTYPE &operator[](int ind) const;

 private:
  int size_;
  FPTYPE *vector_;
};
}  // namespace MyAlgebra

#endif  // __CVCT_H__