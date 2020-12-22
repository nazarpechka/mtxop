#include "vector.h"

#include <iostream>

namespace MyAlgebra {

Vector::Vector(size_t size) : m_size(size), m_array(new FPTYPE[m_size]) {}

Vector::Vector(const Vector &other) { copy(other); }
Vector::Vector(Vector &&other) { move(std::move(other)); }

Vector::~Vector() { delete[] m_array; }

const Vector &Vector::operator=(const Vector &other) {
  if (this != &other) {
    delete[] m_array;

    copy(other);
  }

  return *this;
}

const Vector &Vector::operator=(Vector &&other) {
  if (this != &other) {
    delete[] m_array;

    move(std::move(other));
  }
}

const Vector &Vector::operator=(FPTYPE val) {
  for (int i = 0; i < m_size; ++i) {
    m_array[i] = val;
  }

  return *this;
}

Vector Vector::operator+(const Vector &other) {
  // TODO: Check size
  Vector res(m_size);

  for (int i = 0; i < m_size; ++i) {
    res.m_array[i] = m_array[i] + other.m_array[i];
  }

  return res;
}

Vector Vector::operator-(const Vector &other) {
  // TODO: Check size
  Vector res(m_size);

  for (int i = 0; i < m_size; ++i) {
    res.m_array[i] = m_array[i] - other.m_array[i];
  }

  return res;
}

Vector Vector::operator*(const Matrix &other) {}

Vector Vector::operator~() {}

FPTYPE &Vector::operator[](int ind) const { return m_array[ind]; }

void Vector::copy(const Vector &other) {
  m_size = other.m_size;
  m_array = new FPTYPE[m_size];

  memcpy(m_array, other.m_array, m_size * sizeof(FPTYPE));
}

void Vector::move(Vector &&other) {
  m_size = other.m_size;
  m_array = other.m_array;

  other.m_size = 0;
  other.m_array = nullptr;
}

}  // namespace MyAlgebra