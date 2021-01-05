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

  Vector(const Vector &other);
  Vector(Vector &&other);

  ~Vector();

  const Vector &operator=(const Vector &other);
  const Vector &operator=(Vector &&other);
  const Vector &operator=(FPTYPE val);

  Vector operator+(const Vector &other);
  Vector operator-(const Vector &other);
  Vector operator*(const Matrix &other);
  // Transpozycja - zamiana wektora wierszowego na kolumnowy i odwrotnie
  Vector operator~();

  FPTYPE &operator[](int ind) const;

 private:
  int m_size;
  FPTYPE *m_array;

  void copy(const Vector &other);
  void move(Vector &&other);
};
}  // namespace MyAlgebra

#endif  // __CVCT_H__