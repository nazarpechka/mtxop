#ifndef MATRIX_OPERATIONS_MATRIX_H
#define MATRIX_OPERATIONS_MATRIX_H

#include <immintrin.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "vector.h"

#define TRANSPOSE_AND_SIMD 0
// If 1, use algorithm transpose second matrix
// and SIMD AVX2 instructions

template <typename T>
class Matrix {
 public:
  static constexpr float CMP_PRECISION = 10e-3f;

  // Matrix can be initialized with random values
  Matrix(size_t row_cnt, size_t col_cnt, bool rand_init = false);
  // Creates a square diagonal matrix
  Matrix(size_t dimension, T diagonal);

  Matrix(const Matrix &other);
  Matrix(Matrix &&other) noexcept;

  ~Matrix();

  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other) noexcept;
  Matrix &operator=(T diagonal);
  // Compare to another matrix with precision CMP_PRECISION
  bool operator==(const Matrix &other) const;

  Matrix operator+(const Matrix &other) const;
  Matrix operator-(const Matrix &other) const;
  Matrix operator*(const Matrix &other) const;
  Matrix operator*(const Vector<T> &other) const;
  Matrix operator*(T multiplier) const;

  // Change the sign of all matrix elements
  Matrix operator-() const;
  // Transpose the matrix
  Matrix operator~() const;

  T *operator[](size_t row_ind);
  const T *operator[](size_t row_ind) const;

  Vector<T> extractRow(size_t index) const;
  Vector<T> extractColumn(size_t index) const;

  static Matrix readFromFile(const std::string &filepath,
                             const char &delimiter = ' ');

  [[nodiscard]] size_t rows() const;
  [[nodiscard]] size_t columns() const;
  [[nodiscard]] T get(size_t i, size_t j) const;
  bool set(size_t i, size_t j, T value);

  // Only for testing - print the matrix
  void display() const;
  [[maybe_unused]] static std::string authorName();

 private:
  size_t m_row_cnt;
  size_t m_col_cnt;
  T *m_array;

  T random_number() const;

  void multiply(const Matrix &res, const Matrix &other, size_t start,
                size_t end, size_t inner_dim, size_t other_col_cnt) const;

  void copy(const Matrix &other);
  void move(Matrix &&other);
};

template <typename T>
Matrix<T> operator*(T multiplier, const Matrix<T> &other);

template <typename T>
Matrix<T>::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new T[m_row_cnt * m_col_cnt]) {
  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = random_number();
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);
  }
}

template <typename T>
Matrix<T>::Matrix(size_t dimension, T diagonal)
    : m_row_cnt(dimension),
      m_col_cnt(dimension),
      m_array(new T[dimension * dimension]) {
  memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &other) {
  copy(other);
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept {
  move(std::move(other));
}

template <typename T>
Matrix<T>::~Matrix() {
  delete[] m_array;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other) {
  if (this != &other) {
    delete[] m_array;
    copy(other);
  }
  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&other) noexcept {
  if (this != &other) {
    delete[] m_array;
    move(std::move(other));
  }
  return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(T diagonal) {
  if (m_row_cnt != m_col_cnt) {
    const size_t dimension = std::min(m_row_cnt, m_col_cnt);
    delete[] m_array;

    m_array = new T[dimension * dimension];
    m_row_cnt = m_col_cnt = dimension;
  }

  memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }

  return *this;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T> &other) const {
  if (this != &other) {
    if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
      return false;

    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      if (std::abs(m_array[i] - other.m_array[i]) > CMP_PRECISION) return false;
    }
  }

  return true;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
  const size_t row_cnt = std::min(m_row_cnt, other.m_row_cnt);
  const size_t col_cnt = std::min(m_col_cnt, other.m_col_cnt);
  Matrix<T> res(row_cnt, col_cnt, false);

  const size_t length = row_cnt * col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] + other.m_array[i];
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
  const size_t row_cnt = std::min(m_row_cnt, other.m_row_cnt);
  const size_t col_cnt = std::min(m_col_cnt, other.m_col_cnt);
  Matrix<T> res(row_cnt, col_cnt, false);

  const size_t length = row_cnt * col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] - other.m_array[i];
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
  // When multiplying matrices with incompatible dimensions (ie 3x3 * 4*4),
  // just cut a part of the bigger one from the multiplication
  // x x x      x x x x    x x x x
  // x x x   *  x x x x  = x x x x
  // x x x      x x x x    x x x x
  //            0 0 0 0
  //  3x3         4x4        3x4

  const size_t inner_dim = std::min(m_col_cnt, other.m_row_cnt);
  Matrix<T> res(m_row_cnt, other.m_col_cnt, false);

#if TRANSPOSE_AND_SIMD
  const Matrix transposed(~other);
#endif

  // Multithreading is only efficient for large matrices
  if (m_row_cnt > 95) {
    auto portion =
        static_cast<unsigned int>(ceil(static_cast<float>(m_row_cnt) / 4.0f));

#if TRANSPOSE_AND_SIMD
    std::thread t1([&] {
      multiply(res, transposed, 0, portion, inner_dim, other.m_col_cnt);
    });
    std::thread t2([&] {
      multiply(res, transposed, portion, 2 * portion, inner_dim,
               other.m_col_cnt);
    });
    std::thread t3([&] {
      multiply(res, transposed, 2 * portion, 3 * portion, inner_dim,
               other.m_col_cnt);
    });
    multiply(res, transposed, 3 * portion, m_row_cnt, inner_dim,
             other.m_col_cnt);
#else
    std::thread t1([&] {
      multiply(res, other, 0, portion, inner_dim, other.m_col_cnt);
    });
    std::thread t2([&] {
      multiply(res, other, portion, 2 * portion, inner_dim,
               other.m_col_cnt);
    });
    std::thread t3([&] {
      multiply(res, other, 2 * portion, 3 * portion, inner_dim,
               other.m_col_cnt);
    });
    multiply(res, other, 3 * portion, m_row_cnt, inner_dim,
             other.m_col_cnt);

#endif
    t1.join();
    t2.join();
    t3.join();
  } else {
#if TRANSPOSE_AND_SIMD
    multiply(res, transposed, 0, m_row_cnt, inner_dim, other.m_col_cnt);
#else
    multiply(res, other, 0, m_row_cnt, inner_dim, other.m_col_cnt);
#endif
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(T multiplier) const {
  Matrix<T> res(m_row_cnt, m_col_cnt, false);

  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] * multiplier;
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Vector<T> &other) const {
  //  x 0      x x x     x x x
  //  x 0   *         =  x x x
  //  x 0                x x x
  //  x 0                x x x
  //  4x2       1x3      4x3
  Matrix res(m_row_cnt, other.size());

  for (size_t row = 0; row < m_row_cnt; ++row) {
    for (size_t col = 0; col < other.size(); ++col) {
      res[row][col] = m_array[row * m_col_cnt] * other[col];
    }
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-() const {
  Matrix<T> res(m_row_cnt, m_col_cnt, false);

  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] * -1;
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator~() const {
  Matrix<T> res(m_col_cnt, m_row_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[j * m_row_cnt + i] = m_array[i * m_col_cnt + j];
    }
  }

  return res;
}

template <typename T>
T *Matrix<T>::operator[](size_t row_ind) {
  return &m_array[row_ind * m_col_cnt];
}

template <typename T>
const T *Matrix<T>::operator[](size_t row_ind) const {
  return &m_array[row_ind * m_col_cnt];
}

template <typename T>
Vector<T> Matrix<T>::extractRow(size_t index) const {
  if (index > (m_row_cnt - 1)) return Vector<T>(0);
  Vector<T> res(m_col_cnt);

  for (size_t j = 0; j < m_col_cnt; ++j) {
    res[j] = m_array[index][j];
  }

  return res;
}

template <typename T>
Vector<T> Matrix<T>::extractColumn(size_t index) const {
  if (index > (m_col_cnt - 1)) return Vector<T>(0);
  Vector<T> res(m_row_cnt);

  for (size_t j = 0; j < m_row_cnt; ++j) {
    res[j] = m_array[j * m_col_cnt + index];
  }

  return res;
}

template <typename T>
Matrix<T> Matrix<T>::readFromFile(const std::string &filepath,
                                  const char &delimiter) {
  std::ifstream file;
  file.open(filepath, std::ios::in);
  if (!file.is_open()) return Matrix<T>(0, 0, false);

  size_t rows = 0, columns;
  std::vector<size_t> ind_columns;
  std::string line, number;
  while (std::getline(file, line)) {
    rows++;
    std::stringstream row_stream(line);

    size_t column_len = 0;
    while (std::getline(row_stream, number, delimiter)) column_len++;
    ind_columns.push_back(column_len);
  }
  columns = *std::min_element(ind_columns.begin(), ind_columns.end());

  // We've read the dimensions, now read the numbers
  file.clear();
  file.seekg(0);

  Matrix<T> res(rows, columns);

  for (size_t i = 0; i < rows; ++i) {
    std::getline(file, line);
    std::replace(line.begin(), line.end(), ',', '.');
    std::stringstream row_stream(line);

    for (size_t j = 0;
         j < columns && std::getline(row_stream, number, delimiter); ++j) {
      res[i][j] = static_cast<T>(std::stod(number));
    }
  }

  file.close();
  return res;
}

template <typename T>
[[nodiscard]] size_t Matrix<T>::rows() const {
  return m_row_cnt;
}

template <typename T>
[[nodiscard]] size_t Matrix<T>::columns() const {
  return m_col_cnt;
}

template <typename T>
T Matrix<T>::get(size_t i, size_t j) const {
  if (i > (m_row_cnt - 1) || j > (m_col_cnt - 1)) return -1;
  return m_array[i * m_col_cnt + j];
}

template <typename T>
bool Matrix<T>::set(size_t i, size_t j, T value) {
  if (i > (m_row_cnt - 1) || j > (m_col_cnt - 1)) return false;
  m_array[i * m_col_cnt + j] = value;
  return true;
}

template <typename T>
void Matrix<T>::display() const {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    std::cout << "| ";
    for (size_t j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i * m_col_cnt + j] << "  ";
    }
    std::cout << "|\n";
  }
  std::cout << '\n';
}

template <typename T>
[[maybe_unused]] std::string Matrix<T>::authorName() {
  return "Nazar_Pechevystyi";
}

template <typename T>
T Matrix<T>::random_number() const {
  return rand() % 10;
}

template <>
float Matrix<float>::random_number() const {
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

template <>
double Matrix<double>::random_number() const {
  return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

template <typename T>
void Matrix<T>::multiply(const Matrix<T> &res, const Matrix<T> &other,
                         size_t start, size_t end, size_t inner_dim,
                         size_t other_col_cnt) const {
#if TRANSPOSE_AND_SIMD
  // Usual multiplication of A * ~B
  //    float sum = 0;
  //    for (size_t row = start; row < end; ++row) {
  //      for (size_t other_row = 0; other_row < other_col_cnt; ++other_row) {
  //        sum = 0;
  //        for (size_t pos = 0; pos < inner_dim; ++pos) {
  //          std::cout << "Dot product << this[" << row << "][" << pos << "]
  //          and other[" << other_row << "][" << pos << "]\n"; sum +=
  //          m_array[row * inner_dim + pos] *
  //                 other.m_array[other_row * inner_dim + pos];
  //        }
  //
  //        std::cout << "Writing to res[" << row << "][" << other_row << "]\n";
  //
  //        res.m_array[row * other_col_cnt + other_row] = sum;
  //      }
  //    }

  // Multiplication using SIMD AVX2 intrinsics
  //  if (inner_dim - col < 8) {
  ////              std::cout << "Finishing up with " << inner_dim - col << "
  /// floats\n";
  //    float sum = 0;
  //    for (size_t col_tmp = col; col_tmp < inner_dim; ++col_tmp) {
  //      sum += m_array[row * inner_dim + col_tmp] *
  //             other.m_array[other_row * inner_dim + col_tmp];
  //    }
  ////              std::cout << "Adding to res[" << row << "][" << other_row <<
  ///"], " << sum << "\n";
  //    res.m_array[row * other_col_cnt + other_row] += sum;
  //
  //  }

  __m256 _row_a, _row_b, _other_row_a, _other_row_b, _sum_a, _sum_b;

  for (size_t row = start; row < end; ++row) {
    for (size_t other_row = 0; other_row < other_col_cnt; ++other_row) {
      _sum_a = _mm256_setzero_ps();
      _sum_b = _mm256_setzero_ps();

      for (size_t col = 0; col < inner_dim; col += 16) {
        //        std::cout << "Grabbing 8 floats from this[" << row << "][" <<
        //        col
        //                  << ", this[" << row << "][" << col + 8
        //                  << "], 8 floats from other[" << other_row << "][" <<
        //                  col << "], and"
        //                  << " other[" << other_row << "][" << col + 8 <<
        //                  "]\n";
        _row_a = _mm256_load_ps(&m_array[row * inner_dim + col]);
        _row_b = _mm256_load_ps(&m_array[row * inner_dim + col + 8]);

        _other_row_a =
            _mm256_load_ps(&other.m_array[other_row * inner_dim + col]);
        _other_row_b =
            _mm256_load_ps(&other.m_array[other_row * inner_dim + col + 8]);
        _sum_a = _mm256_add_ps(_sum_a, _mm256_mul_ps(_row_a, _other_row_a));
        _sum_b = _mm256_add_ps(_sum_b, _mm256_mul_ps(_row_b, _other_row_b));
      }

      //      _mm256_storeu_ps(&res.m_array[row * other_col_cnt + other_row],
      //      _sum_res);

      //      std::cout << "Adding to res[" << row << "][" << other_row <<
      //      "]\n";
      for (size_t i = 0; i < 8; ++i) {
        res.m_array[row * other_col_cnt + other_row] += _sum_a[i];
      }
      for (size_t i = 0; i < 8; ++i) {
        res.m_array[row * other_col_cnt + other_row] += _sum_b[i];
      }
    }
  }

#else
  // Swapping loops algorithm
  // Swapping second and third loop to move over other
  // matrix columns, not rows (less cache misses)
  // 1. Iterate over first matrix rows
  // 2. Iterate over second matrix rows
  // 3. Iterate over columns

  for (size_t row = start; row < end; ++row) {
    for (size_t pos = 0; pos < inner_dim; ++pos) {
      for (size_t col = 0; col < other_col_cnt; ++col) {
        res.m_array[row * other_col_cnt + col] +=
            m_array[row * inner_dim + pos] *
            other.m_array[pos * other_col_cnt + col];
      }
    }
  }

#endif
}

template <typename T>
void Matrix<T>::copy(const Matrix<T> &other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  const size_t size = m_row_cnt * m_col_cnt;
  m_array = new T[size];

  memcpy(m_array, other.m_array, sizeof(T) * size);
}

template <typename T>
void Matrix<T>::move(Matrix<T> &&other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  m_array = other.m_array;

  other.m_row_cnt = 0;
  other.m_col_cnt = 0;
  other.m_array = nullptr;
}

template <typename T>
Matrix<T> operator*(T multiplier, const Matrix<T> &other) {
  return other * multiplier;
}

#endif  // MATRIX_OPERATIONS_MATRIX_H