#include "CVct.h"

#include <iostream>

namespace MyAlgebra {

const int CVct::kDefaultSize = 10;

CVct::CVct(int size) {
  if (size < 1) size = kDefaultSize;
  size_ = size;
  vector_ = new FPTYPE[size_];
}

CVct::~CVct() { delete[] vector_; }

const CVct &CVct::operator=(const CVct &rhs) {
  if (this == &rhs) return *this;

  delete[] vector_;

  size_ = rhs.size_;
  vector_ = new FPTYPE[size_];
  // TODO: Change shallow copy to deep copy?
  memcpy(vector_, rhs.vector_, size_ * sizeof(FPTYPE));

  return *this;
}

const CVct &CVct::operator=(FPTYPE val) {
  for (int i = 0; i < size_; ++i) {
    vector_[i] = val;
  }
}

CVct CVct::operator-(const CVct &rhs) {
  // TODO: Check size
  CVct res(size_);

  for (int i = 0; i < size_; ++i) {
    res.vector_[i] = vector_[i] - rhs.vector_[i];
  }

  return res;
}

CVct CVct::operator~() {}

CVct CVct::operator*(const CMtx &rhs) {}

CVct CVct::operator+(const CVct &rhs) {
  // TODO: Check size
  CVct res(size_);

  for (int i = 0; i < size_; ++i) {
    res.vector_[i] = vector_[i] + rhs.vector_[i];
  }

  return res;
}

FPTYPE &CVct::operator[](int ind) const { return vector_[ind]; }
}  // namespace MyAlgebra