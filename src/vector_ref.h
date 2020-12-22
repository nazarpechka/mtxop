#ifndef __VECTOR_REF_H__
#define __VECTOR_REF_H__

typedef double FPTYPE;

#include <stddef.h>

#include "matrix_ref.h"
// #include "CAlgError.h"

namespace RefAlgebra {
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
}  // namespace RefAlgebra

#endif  // __VECTOR_REF_H__