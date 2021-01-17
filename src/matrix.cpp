#include "matrix.h"

#include <cmath>
#include <iostream>
#include <thread>

#define MULTITHREADING 1
#define SWAP_LOOPS 1


namespace MyAlgebra {

const float Matrix::ALG_PRECISION = 10e-6f;

Matrix::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new float[m_row_cnt * m_col_cnt]) {
  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt; ++i) {
      for (size_t j = 0; j < m_col_cnt; ++j) {
        m_array[i * m_col_cnt + j] =
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
      }
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(float) * m_row_cnt * m_col_cnt);
  }
}

Matrix::Matrix(size_t row_cnt, float diagonal)
    : m_row_cnt(row_cnt),
      m_col_cnt(row_cnt),
      m_array(new float[m_row_cnt * m_row_cnt]) {
  memset(m_array, 0, sizeof(float) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }
}

Matrix::Matrix(const Matrix &other) { copy(other); }

Matrix::Matrix(Matrix &&other) { move(std::move(other)); }

Matrix::~Matrix() { delete[] m_array; }

const Matrix &Matrix::operator=(const Matrix &other) {
  if (this != &other) {
    delete[] m_array;

    copy(other);
  }

  return *this;
}

const Matrix &Matrix::operator=(Matrix &&other) noexcept {
  if (this != &other) {
    delete[] m_array;

    move(std::move(other));
  }

  return *this;
}

const Matrix &Matrix::operator=(const float diagonal) {
  memset(m_array, 0, sizeof(float) * m_row_cnt * m_col_cnt);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i * m_col_cnt + i] = diagonal;
  }

  return *this;
}

bool Matrix::operator==(const Matrix &other) const {
  if (this != &other) {
    if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
      return false;

    size_t tmp_idx;
    for (size_t i = 0; i < m_row_cnt; ++i) {
      for (size_t j = 0; j < m_col_cnt; ++j) {
        tmp_idx = i * m_col_cnt + j;
        if (std::abs(m_array[tmp_idx] - other.m_array[tmp_idx]) > ALG_PRECISION)
          return false;
      }
    }
  }

  return true;
}

Matrix Matrix::operator+(const Matrix &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return Matrix(0, 0, false);

  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
    res.m_array[i] = m_array[i] + other.m_array[i];
  }

  return res;
}

Matrix Matrix::operator-(const Matrix &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return Matrix(0, 0, false);

  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
    res.m_array[i] = m_array[i] - other.m_array[i];
  }

  return res;
}

Matrix Matrix::operator-() const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
    res.m_array[i] = m_array[i] * -1;
  }

  return res;
}

Matrix Matrix::operator*(const Matrix &other) {
  if (m_col_cnt != other.m_row_cnt) return Matrix(0, 0, false);

  Matrix res(m_row_cnt, other.m_col_cnt, false);

#if MULTITHREADING
  auto portion = static_cast<unsigned int>(ceil(m_row_cnt / 4.0f));

  std::thread t1([&] { multiply(res, other, 0, portion); });
  std::thread t2([&] { multiply(res, other, portion, 2 * portion); });
  std::thread t3([&] { multiply(res, other, 2 * portion, 3 * portion); });
  std::thread t4([&] { multiply(res, other, 3 * portion, m_row_cnt); });

  t1.join();
  t2.join();
  t3.join();
  t4.join();
#else

  multiply(res, other, 0, m_row_cnt);

#endif

  return res;
}

Matrix Matrix::operator*(float multiplier) const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
    res.m_array[i] = m_array[i] * multiplier;
  }

  return res;
}

Matrix Matrix::operator~() const {
  Matrix res(m_col_cnt, m_row_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[j * m_row_cnt + i] = m_array[i * m_col_cnt + j];
    }
  }

  return res;
}

Matrix Matrix::operator^(int power) const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  if (power == -1) {
  } else if (power == 0) {
    for (size_t i = 0; i < m_row_cnt; ++i) {
      for (size_t j = 0; j < m_col_cnt; ++j) {
        if (i != j)
          res.m_array[i * m_col_cnt + j] = 0;
        else
          res.m_array[i * m_col_cnt + j] = 1;
      }
    }
  } else if (power == 1) {
    res = *this;
  } else {
    res = *this;
    for (int i = 1; i < power; ++i) {
      res = res * (*this);
    }
  }

  return res;
}

float *Matrix::operator[](size_t row_ind) {
  return &m_array[row_ind * m_col_cnt];
}

const float *Matrix::operator[](size_t row_ind) const {
  return &m_array[row_ind * m_col_cnt];
}

size_t Matrix::getRowCount() const { return m_row_cnt; }

size_t Matrix::getColCount() const { return m_col_cnt; }

void Matrix::display() const {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i * m_col_cnt + j] << "  ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

void Matrix::multiply(const Matrix &res, const Matrix &other, size_t start,
                      size_t end) const {
  const size_t col_cnt = m_col_cnt;
  const size_t other_col_cnt = other.m_col_cnt;

#if SWAP_LOOPS
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

#else
  // Naive multiplication

  float sum = 0;
  for (int row = start; row < end; ++row) {
    for (int col = 0; col < other_col_cnt; ++col) {
      sum = 0;
      for (int pos = 0; pos < col_cnt; ++pos) {
        sum += m_array[row * col_cnt + pos] *
               other.m_array[pos * other_col_cnt + col];
      }

      res.m_array[row * other_col_cnt + col] = sum;
    }
  }
#endif
}

void Matrix::copy(const Matrix &other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  const size_t size = m_row_cnt * m_col_cnt;
  m_array = new float[size];

  memcpy(m_array, other.m_array, sizeof(float) * size);
}

void Matrix::move(Matrix &&other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  m_array = other.m_array;

  other.m_row_cnt = 0;
  other.m_col_cnt = 0;
  other.m_array = nullptr;
}

Matrix operator*(float multiplier, const Matrix &other) {
  return other * multiplier;
}

}  // namespace MyAlgebra