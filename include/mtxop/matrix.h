#ifndef MATRIX_OPERATIONS_MATRIX_H
#define MATRIX_OPERATIONS_MATRIX_H

#include <cstddef>
#include <string>
#include <immintrin.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <thread>

#define TRANSPOSE_AND_SIMD 0
// If 1, use algorithm transpose second matrix
// and SIMD AVX2 instructions

template<typename T>
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
  // Change the sign of all matrix elements
  Matrix operator-() const;

  Matrix operator*(const Matrix &other);
  Matrix operator*(T multiplier) const;

  // Transpose the matrix
  Matrix operator~() const;

  float *operator[](size_t row_ind);
  const float *operator[](size_t row_ind) const;

  [[nodiscard]] size_t getRowCount() const;
  [[nodiscard]] size_t getColCount() const;

  // Only for testing - print the matrix
  void display() const;
  [[maybe_unused]] static std::string authorName();

 private:
  size_t m_row_cnt;
  size_t m_col_cnt;
  T *m_array;

  void multiply(const Matrix &res, const Matrix &other, size_t start,
                size_t end, size_t col_cnt, size_t other_col_cnt) const;

  void copy(const Matrix &other);
  void move(Matrix &&other);
};

template<typename T>
Matrix<T> operator*(T multiplier, const Matrix<T> &other);

template<typename T>
Matrix<T>::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new T[m_row_cnt * m_col_cnt]) {

  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = static_cast<T>(rand());
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);
  }
}

template<>
Matrix<int>::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new int[m_row_cnt * m_col_cnt]) {

  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = rand();
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(int) * m_row_cnt * m_col_cnt);
  }
}

template<>
Matrix<float>::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new float[m_row_cnt * m_col_cnt]) {

  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(float) * m_row_cnt * m_col_cnt);
  }
}

template<>
Matrix<double>::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new double[m_row_cnt * m_col_cnt]) {

  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(double) * m_row_cnt * m_col_cnt);
  }
}

template<typename T>
Matrix<T>::Matrix(size_t dimension, T diagonal)
    : m_row_cnt(dimension),
      m_col_cnt(dimension),
      m_array(new T[dimension * dimension]) {
  memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> &other) { copy(other); }

template<typename T>
Matrix<T>::Matrix(Matrix<T> &&other) noexcept { move(std::move(other)); }

template<typename T>
Matrix<T>::~Matrix() { delete[] m_array; }

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &other) {
  if (this != &other) {
    delete[] m_array;

    copy(other);
  }

  return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&other) noexcept {
  if (this != &other) {
    delete[] m_array;

    move(std::move(other));
  }

  return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(T diagonal) {
  memset(m_array, 0, sizeof(T) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }

  return *this;
}

template<typename T>
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

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return Matrix<T>(0, 0, false);

  Matrix<T> res(m_row_cnt, m_col_cnt, false);

  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] + other.m_array[i];
  }

  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return Matrix<T>(0, 0, false);

  Matrix<T> res(m_row_cnt, m_col_cnt, false);

  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] - other.m_array[i];
  }

  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const {
  Matrix<T> res(m_row_cnt, m_col_cnt, false);
  
  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] * -1;
  }

  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) {
  if (m_col_cnt != other.m_row_cnt) return Matrix<T>(0, 0, false);

  Matrix<T> res(m_row_cnt, other.m_col_cnt, false);

  // Multithreading is only efficient for large matrices
  if (m_row_cnt > 95) {
    auto portion =
        static_cast<unsigned int>(ceil(static_cast<float>(m_row_cnt) / 4.0f));

#if TRANSPOSE_AND_SIMD
    const Matrix transposed(~other);

    std::thread t1([&] {
      multiply(res, transposed, 0, portion, m_col_cnt, other.m_col_cnt);
    });
    std::thread t2([&] {
      multiply(res, transposed, portion, 2 * portion, m_col_cnt,
               other.m_col_cnt);
    });
    std::thread t3([&] {
      multiply(res, transposed, 2 * portion, 3 * portion, m_col_cnt,
               other.m_col_cnt);
    });
    multiply(res, transposed, 3 * portion, m_row_cnt, m_col_cnt,
             other.m_col_cnt);
#else
    std::thread t1(
        [&] { multiply(res, other, 0, portion, m_col_cnt, other.m_col_cnt); });
    std::thread t2([&] {
      multiply(res, other, portion, 2 * portion, m_col_cnt, other.m_col_cnt);
    });
    std::thread t3([&] {
      multiply(res, other, 2 * portion, 3 * portion, m_col_cnt,
               other.m_col_cnt);
    });
    multiply(res, other, 3 * portion, m_row_cnt, m_col_cnt, other.m_col_cnt);

#endif

    t1.join();
    t2.join();
    t3.join();

  } else {
#if TRANSPOSE_AND_SIMD
    const Matrix transposed(~other);
    multiply(res, transposed, 0, m_row_cnt, m_col_cnt, other.m_col_cnt);
#else
    multiply(res, other, 0, m_row_cnt, m_col_cnt, other.m_col_cnt);
#endif
  }

  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(T multiplier) const {
  Matrix<T> res(m_row_cnt, m_col_cnt, false);

  const size_t length = m_row_cnt * m_col_cnt;
  for (size_t i = 0; i < length; ++i) {
    res.m_array[i] = m_array[i] * multiplier;
  }

  return res;
}

template<typename T>
Matrix<T> Matrix<T>::operator~() const {
  Matrix<T> res(m_col_cnt, m_row_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[j * m_row_cnt + i] = m_array[i * m_col_cnt + j];
    }
  }

  return res;
}

template<typename T>
float *Matrix<T>::operator[](size_t row_ind) {
  return &m_array[row_ind * m_col_cnt];
}

template<typename T>
const float *Matrix<T>::operator[](size_t row_ind) const {
  return &m_array[row_ind * m_col_cnt];
}

template<typename T>
[[nodiscard]] size_t Matrix<T>::getRowCount() const { return m_row_cnt; }

template<typename T>
[[nodiscard]] size_t Matrix<T>::getColCount() const { return m_col_cnt; }

template<typename T>
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

template<typename T>
[[maybe_unused]] std::string Matrix<T>::authorName() { return "Nazar_Pechevystyi"; }

template<typename T>
void Matrix<T>::multiply(const Matrix<T> &res, const Matrix<T> &other, size_t start,
                      size_t end, size_t col_cnt, size_t other_col_cnt) const {
#if TRANSPOSE_AND_SIMD
  // Usual multiplication of A * ~B
  //    float sum = 0;
  //    for (size_t row = start; row < end; ++row) {
  //      for (size_t other_row = 0; other_row < other_col_cnt; ++other_row) {
  //        sum = 0;
  //        for (size_t pos = 0; pos < col_cnt; ++pos) {
  //          std::cout << "Dot product << this[" << row << "][" << pos << "]
  //          and other[" << other_row << "][" << pos << "]\n"; sum +=
  //          m_array[row * col_cnt + pos] *
  //                 other.m_array[other_row * col_cnt + pos];
  //        }
  //
  //        std::cout << "Writing to res[" << row << "][" << other_row << "]\n";
  //
  //        res.m_array[row * other_col_cnt + other_row] = sum;
  //      }
  //    }

  // Multiplication using SIMD AVX2 intrinsics
  //  if (col_cnt - col < 8) {
  ////              std::cout << "Finishing up with " << col_cnt - col << "
  /// floats\n";
  //    float sum = 0;
  //    for (size_t col_tmp = col; col_tmp < col_cnt; ++col_tmp) {
  //      sum += m_array[row * col_cnt + col_tmp] *
  //             other.m_array[other_row * col_cnt + col_tmp];
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

      for (size_t col = 0; col < col_cnt; col += 16) {
        //        std::cout << "Grabbing 8 floats from this[" << row << "][" <<
        //        col
        //                  << ", this[" << row << "][" << col + 8
        //                  << "], 8 floats from other[" << other_row << "][" <<
        //                  col << "], and"
        //                  << " other[" << other_row << "][" << col + 8 <<
        //                  "]\n";
        _row_a = _mm256_load_ps(&m_array[row * col_cnt + col]);
        _row_b = _mm256_load_ps(&m_array[row * col_cnt + col + 8]);

        _other_row_a =
            _mm256_load_ps(&other.m_array[other_row * col_cnt + col]);
        _other_row_b =
            _mm256_load_ps(&other.m_array[other_row * col_cnt + col + 8]);
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
    for (size_t pos = 0; pos < col_cnt; ++pos) {
      for (size_t col = 0; col < other_col_cnt; ++col) {
        res.m_array[row * other_col_cnt + col] +=
            m_array[row * col_cnt + pos] *
            other.m_array[pos * other_col_cnt + col];
      }
    }
  }

#endif
}

template<typename T>
void Matrix<T>::copy(const Matrix<T> &other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  const size_t size = m_row_cnt * m_col_cnt;
  m_array = new float[size];

  memcpy(m_array, other.m_array, sizeof(float) * size);
}

template<typename T>
void Matrix<T>::move(Matrix<T> &&other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  m_array = other.m_array;

  other.m_row_cnt = 0;
  other.m_col_cnt = 0;
  other.m_array = nullptr;
}

template<typename T>
Matrix<T> operator*(T multiplier, const Matrix<T> &other) {
  return other * multiplier;
}

#endif  // MATRIX_OPERATIONS_MATRIX_H