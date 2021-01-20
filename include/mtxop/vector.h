//
// Created by Nazar Pechevystyi on 20.01.2021.
//

#ifndef MATRIX_OPERATIONS_VECTOR_H
#define MATRIX_OPERATIONS_VECTOR_H

#include <cstddef>

template <typename T>
class Matrix;

template <typename T>
class Vector {
 public:
  explicit Vector(size_t size, bool rand_init = false);

  Vector(const Vector &other);
  Vector(Vector &&other) noexcept;

  ~Vector();

  Vector &operator=(const Vector &other);
  Vector &operator=(Vector &&other) noexcept;
  Vector &operator=(T val);
  // Compare to another vector with precision CMP_PRECISION
  bool operator==(const Vector &other) const;

  Vector operator+(const Vector &other) const;
  Vector operator-(const Vector &other) const;
  Vector operator*(const Matrix<T> &other) const;
  Vector operator*(T multiplier) const;

  T &operator[](int ind);
  const T &operator[](int ind) const;

  [[nodiscard]] size_t size() const;
  [[nodiscard]] T get(size_t i) const;
  bool set(size_t i, T value);

  void display() const;

 private:
  size_t m_size;
  T *m_array;

  [[nodiscard]] T random_number() const;

  void copy(const Vector &other);
  void move(Vector &&other);
};

template <typename T>
Vector<T>::Vector(size_t size, bool rand_init) : m_size(size), m_array(new T[m_size]) {
  if (rand_init) {
    for (size_t i = 0; i < m_size; ++i) {
      m_array[i] = random_number();
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(T) * m_size);
  }
}

template <typename T>
Vector<T>::Vector(const Vector &other) { copy(other); }

template <typename T>
Vector<T>::Vector(Vector &&other) noexcept { move(std::move(other)); }

template <typename T>
Vector<T>::~Vector() { delete[] m_array; }

template <typename T>
Vector<T> &Vector<T>::operator=(const Vector &other) {
  if (this != &other) {
    delete[] m_array;

    copy(other);
  }

  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator=(Vector &&other) noexcept {
  if (this != &other) {
    delete[] m_array;

    move(std::move(other));
  }
}

template <typename T>
Vector<T> &Vector<T>::operator=(T val) {
  for (int i = 0; i < m_size; ++i) {
    m_array[i] = val;
  }

  return *this;
}

template <typename T>
bool Vector<T>::operator==(const Vector &other) const {
  if(this != &other) {
    if(m_size != other.m_size) return false;

    for(size_t i = 0; i < m_size; ++i) {
      if (std::abs(m_array[i] - other.m_array[i]) > Matrix<T>::CMP_PRECISION) return false;
    }
  }

  return true;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector &other) const {
  const size_t size = std::min(m_size, other.m_size);
  Vector res(size);

  for (size_t i = 0; i < size; ++i) {
    res.m_array[i] = m_array[i] + other.m_array[i];
  }

  return res;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector &other) const {
  const size_t size = std::min(m_size, other.m_size);
  Vector res(size);

  for (size_t i = 0; i < size; ++i) {
    res.m_array[i] = m_array[i] - other.m_array[i];
  }

  return res;
}

template <typename T>
Vector<T> Vector<T>::operator*(const Matrix<T> &other) const {
  //  x x x     x x      x x
  //         *  x x   =
  //            x x
  //            0 0
  //   1x3      4x2      1x2
  const size_t inner_dimension = std::min(m_size, other.rows());
  Vector res(other.columns());

  float sum;
  for (size_t col = 0; col < other.columns(); ++col) {
    sum = 0;

    for (size_t pos = 0; pos < inner_dimension; ++pos) {
      sum += m_array[pos] * other[pos][col];
    }

    res[col] = sum;
  }

  return res;
}

template <typename T>
Vector<T> Vector<T>::operator*(T multiplier) const {
  Vector res(m_size);

  for(size_t i = 0; i < m_size; ++i) {
    res.m_array[i] = m_array[i] * multiplier;
  }

  return res;
}

template <typename T>
T &Vector<T>::operator[](int ind) { return m_array[ind]; }

template <typename T>
const T &Vector<T>::operator[](int ind) const { return m_array[ind]; }

template <typename T>
size_t Vector<T>::size() const {
  return m_size;
}

template <typename T>
T Vector<T>::get(size_t i) const {
  if(i > (m_size - 1)) return -1;
  return m_array[i];
}

template <typename T>
bool Vector<T>::set(size_t i, T value) {
  if(i > (m_size - 1)) return false;
  m_array[i] = value;
  return true;
}

template <typename T>
void Vector<T>::display() const {
  std::cout << "| ";
  for (size_t i = 0; i < m_size; ++i) {
    std::cout << m_array[i] << " ";
  }
  std::cout << "|\n\n";
}

template <typename T>
T Vector<T>::random_number() const {
  return rand() % 10;
}

template <>
float Vector<float>::random_number() const {
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

template <>
double Vector<double>::random_number() const {
  return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}


template <typename T>
void Vector<T>::copy(const Vector &other) {
  m_size = other.m_size;
  m_array = new T[m_size];

  memcpy(m_array, other.m_array, m_size * sizeof(T));
}

template <typename T>
void Vector<T>::move(Vector &&other) {
  m_size = other.m_size;
  m_array = other.m_array;

  other.m_size = 0;
  other.m_array = nullptr;
}


#endif  // MATRIX_OPERATIONS_VECTOR_H
