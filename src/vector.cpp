#include "vector.h"

#include <iostream>

namespace MyAlgebra {

Vector::Vector(size_t size) : size_(size), vector_(new FPTYPE[size_]) {}

Vector::~Vector() { delete[] vector_; }

const Vector &Vector::operator=(const Vector &rhs) {
  if (this == &rhs) return *this;

  delete[] vector_;

  size_ = rhs.size_;
  vector_ = new FPTYPE[size_];
  // TODO: Change shallow copy to deep copy?
  memcpy(vector_, rhs.vector_, size_ * sizeof(FPTYPE));

  return *this;
}

const Vector &Vector::operator=(FPTYPE val) {
  for (int i = 0; i < size_; ++i) {
    vector_[i] = val;
  }

  return *this;
}

Vector Vector::operator-(const Vector &rhs) {
  // TODO: Check size
  Vector res(size_);

  for (int i = 0; i < size_; ++i) {
    res.vector_[i] = vector_[i] - rhs.vector_[i];
  }

  return res;
}

Vector Vector::operator~() {}

Vector Vector::operator*(const Matrix &rhs) {}

Vector Vector::operator+(const Vector &rhs) {
  // TODO: Check size
  Vector res(size_);

  for (int i = 0; i < size_; ++i) {
    res.vector_[i] = vector_[i] + rhs.vector_[i];
  }

  return res;
}

FPTYPE &Vector::operator[](int ind) const { return vector_[ind]; }
}  // namespace MyAlgebra