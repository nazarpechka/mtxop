#include "matrix_ref.h"

#include <cmath>
#include <iostream>

namespace RefAlgebra {

const float Matrix::ALG_PRECISION = 10e-6f;

Matrix::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt), m_col_cnt(col_cnt), m_array(new float *[m_row_cnt]) {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new float[m_col_cnt];
  }

  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt; ++i) {
      for (size_t j = 0; j < m_col_cnt; ++j) {
        m_array[i][j] =
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
      }
    }
  } else {
    // Initialize with zeros
    for (size_t i = 0; i < m_row_cnt; ++i) {
      memset(m_array[i], 0, sizeof(float) * m_col_cnt);
    }
  }
}

Matrix::Matrix(size_t row_cnt, float diagonal)
    : m_row_cnt(row_cnt), m_col_cnt(row_cnt), m_array(new float *[m_row_cnt]) {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new float[m_row_cnt];
    memset(m_array[i], 0, sizeof(float) * m_row_cnt);
    m_array[i][i] = diagonal;
  }
}

Matrix::Matrix(const Matrix &other) { copy(other); }

Matrix::Matrix(Matrix &&other) { move(std::move(other)); }

Matrix::~Matrix() {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    delete[] m_array[i];
  }
  delete[] m_array;
}

const Matrix &Matrix::operator=(const Matrix &other) {
  if (this != &other) {
    for (size_t i = 0; i < m_row_cnt; ++i) {
      delete[] m_array[i];
    }
    delete[] m_array;

    copy(other);
  }

  return *this;
}

const Matrix &Matrix::operator=(Matrix &&other) noexcept {
  if (this != &other) {
    for (size_t i = 0; i < m_row_cnt; ++i) {
      delete[] m_array[i];
    }
    delete[] m_array;

    move(std::move(other));
  }

  return *this;
}

const Matrix &Matrix::operator=(const float diagonal) {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    memset(m_array[i], 0, sizeof(float) * m_row_cnt);
    m_array[i][i] = diagonal;
  }

  return *this;
}

bool Matrix::operator==(const Matrix &other) const {
  if (this != &other) {
    if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
      return false;

    for (size_t i = 0; i < m_row_cnt; ++i) {
      for (size_t j = 0; j < m_col_cnt; ++j) {
        if (std::abs(m_array[i][j] - other.m_array[i][j]) > ALG_PRECISION)
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

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] + other.m_array[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator-(const Matrix &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return Matrix(0, 0, false);

  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] - other.m_array[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator-() const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] * -1;
    }
  }

  return res;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (m_col_cnt != other.m_row_cnt) return Matrix(0, 0, false);

  Matrix res(m_row_cnt, other.m_col_cnt, false);

  float sum;
  for (size_t row = 0; row < m_row_cnt; ++row) {
    for (size_t col = 0; col < other.m_col_cnt; ++col) {
      sum = 0;

      for (size_t pos = 0; pos < m_col_cnt; ++pos) {
        sum += m_array[row][pos] * other.m_array[pos][col];
      }

      res.m_array[row][col] = sum;
    }
  }

  return res;
}

Matrix Matrix::operator*(float multiplier) const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] * multiplier;
    }
  }

  return res;
}

Matrix Matrix::operator~() const {
  Matrix res(m_col_cnt, m_row_cnt, false);

  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      res.m_array[j][i] = m_array[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator^(int power) const {}

float *Matrix::operator[](size_t row_ind) { return m_array[row_ind]; }

const float *Matrix::operator[](size_t row_ind) const { return m_array[row_ind]; }

size_t Matrix::getRowCount() const { return m_row_cnt; }

size_t Matrix::getColCount() const { return m_col_cnt; }

void Matrix::display() const {
  for (size_t i = 0; i < m_row_cnt; ++i) {
    for (size_t j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i][j] << "  ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

void Matrix::copy(const Matrix &other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  m_array = new float *[m_row_cnt];

  for (size_t i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new float[m_col_cnt];
    memcpy(m_array[i], other.m_array[i], sizeof(float) * m_col_cnt);
  }
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

}  // namespace RefAlgebra
