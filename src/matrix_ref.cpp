#include "matrix_ref.h"

#include <cstdlib>
#include <iostream>

namespace RefAlgebra {

const FPTYPE Matrix::ALG_PRECISION = 10e-6;

Matrix::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt), m_col_cnt(col_cnt), m_array(new FPTYPE *[m_row_cnt]) {
  for (int i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new FPTYPE[m_col_cnt];
  }

  if (rand_init) {
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        m_array[i][j] = rand();
      }
    }
  }
}

Matrix::Matrix(size_t row_cnt, FPTYPE diagonal)
    : m_row_cnt(row_cnt), m_col_cnt(row_cnt), m_array(new FPTYPE *[m_row_cnt]) {
  for (int i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new FPTYPE[m_row_cnt];

    for (int j = 0; j < m_row_cnt; ++j) {
      if (i != j)
        m_array[i][j] = 0;
      else
        m_array[i][j] = diagonal;
    }
  }
}

Matrix::Matrix(const Matrix &other) { copy(other); }

Matrix::Matrix(Matrix &&other) { move(std::move(other)); }

Matrix::~Matrix() {
  for (int i = 0; i < m_row_cnt; ++i) {
    delete[] m_array[i];
  }
  delete[] m_array;
}

const Matrix &Matrix::operator=(const Matrix &other) {
  if (this != &other) {
    for (int i = 0; i < m_row_cnt; ++i) {
      delete[] m_array[i];
    }
    delete[] m_array;

    copy(other);
  }

  return *this;
}

const Matrix &Matrix::operator=(Matrix &&other) {
  if (this != &other) {
    for (int i = 0; i < m_row_cnt; ++i) {
      delete[] m_array[i];
    }
    delete[] m_array;

    move(std::move(other));
  }

  return *this;
}

const Matrix &Matrix::operator=(const FPTYPE diagonal) {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      if (i != j)
        m_array[i][j] = 0;
      else
        m_array[i][j] = diagonal;
    }
  }

  return *this;
}

Matrix Matrix::operator+(const Matrix &other) const {
  // TODO: Check size, should be the same

  Matrix res(m_row_cnt, m_col_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] + other.m_array[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator-(const Matrix &other) const {
  // TODO: Check size, should be the same

  Matrix res(m_row_cnt, m_col_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] - other.m_array[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator-() const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      // TODO: Change for other types
      res.m_array[i][j] = m_array[i][j] * -1;
    }
  }

  return res;
}

Matrix Matrix::operator*(const Matrix &other) const {
  // TODO: Check size
  Matrix res(m_row_cnt, other.m_col_cnt, false);

  for (int l_row = 0; l_row < m_row_cnt; ++l_row) {
    for (int r_col = 0; r_col < other.m_col_cnt; ++r_col) {
      FPTYPE sum = 0;
      int iter = 0;
      while (iter < m_row_cnt) {
        sum += m_array[l_row][iter] * other.m_array[iter][r_col];
        iter++;
      }

      res.m_array[l_row][r_col] = sum;
    }
  }

  return res;
}

Vector Matrix::operator*(const Vector &other) const {}

Matrix Matrix::operator*(FPTYPE multiplier) const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[i][j] * multiplier;
    }
  }

  return res;
}

Matrix Matrix::operator~() const {
  Matrix res(m_col_cnt, m_row_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i][j] = m_array[j][i];
    }
  }

  return res;
}

Matrix Matrix::operator^(int power) const {}

FPTYPE *Matrix::operator[](int row_ind) { return m_array[row_ind]; }

int Matrix::getRowCount() const { return m_row_cnt; }

int Matrix::getColCount() const { return m_col_cnt; }

bool Matrix::operator==(const Matrix &other) const {
  if (this != &other) {
    if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
      return false;

    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        if (m_array[i][j] - other.m_array[i][j] > ALG_PRECISION) return false;
      }
    }
  }

  return true;
}

void Matrix::display() const {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i][j] << " ";
    }
    std::cout << '\n';
  }
}

void Matrix::copy(const Matrix &other) {
  m_row_cnt = other.m_row_cnt;
  m_col_cnt = other.m_col_cnt;
  m_array = new FPTYPE *[m_row_cnt];

  for (int i = 0; i < m_row_cnt; ++i) {
    m_array[i] = new FPTYPE[m_col_cnt];
    memcpy(m_array[i], other.m_array[i], sizeof(FPTYPE) * m_col_cnt);
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

Matrix operator*(FPTYPE multiplier, const Matrix &other) {
  return other * multiplier;
}

}  // namespace RefAlgebra
