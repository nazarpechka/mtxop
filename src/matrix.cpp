#include "matrix.h"

#include <math.h>

#include <cstdlib>
#include <iostream>
#include <thread>

#define MULTITHREAD 1
#define CHANGE_LOOPS_ORDER 1

namespace MyAlgebra {

const FPTYPE Matrix::ALG_PRECISION = 10e-6;

Matrix::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new FPTYPE[m_row_cnt * m_col_cnt]) {
  if (rand_init) {
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        m_array[i * m_col_cnt + j] =
            static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      }
    }
  } else {
    // Initialize with zeros
    memset(m_array, 0, sizeof(FPTYPE) * m_row_cnt * m_col_cnt);
  }
}

Matrix::Matrix(size_t row_cnt, FPTYPE diagonal)
    : m_row_cnt(row_cnt),
      m_col_cnt(row_cnt),
      m_array(new FPTYPE[m_row_cnt * m_row_cnt]) {
  memset(m_array, 0, sizeof(FPTYPE) * m_row_cnt * m_col_cnt);

  for (int i = 0; i < m_row_cnt; ++i) {
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

const Matrix &Matrix::operator=(Matrix &&other) {
  if (this != &other) {
    delete[] m_array;

    move(std::move(other));
  }

  return *this;
}

const Matrix &Matrix::operator=(const FPTYPE diagonal) {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      if (i != j)
        m_array[i * m_col_cnt + j] = 0;
      else
        m_array[i * m_col_cnt + j] = diagonal;
    }
  }

  return *this;
}

bool Matrix::operator==(const Matrix &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return false;

  if (this != &other) {
    int tmp_idx;
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        tmp_idx = i * m_col_cnt + j;
        if (abs(m_array[tmp_idx] - other.m_array[tmp_idx]) > ALG_PRECISION)
          return false;
      }
    }
  }

  return true;
}

Matrix Matrix::operator+(const Matrix &other) const {
  // TODO: Check size, should be the same

  Matrix res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] + other.m_array[tmp_idx];
    }
  }

  return res;
}

Matrix Matrix::operator-(const Matrix &other) const {
  // TODO: Check size, should be the same

  Matrix res(m_row_cnt, m_col_cnt, false);

  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] - other.m_array[tmp_idx];
    }
  }

  return res;
}

Matrix Matrix::operator-() const {
  Matrix res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] * -1;
    }
  }

  return res;
}

Matrix Matrix::operator*(const Matrix &other) const {
  if (m_col_cnt != other.m_row_cnt) {
    std::cout << "YOU VIOLATED THE LAW!\n";
    return Matrix(0, 0, false);
  }

  Matrix res(m_row_cnt, other.m_col_cnt, false);

#if MULTITHREAD
  int portion = ceil(m_row_cnt / 4.0f);
  // std::cout << "portion = " << portion << '\n';

  // std::cout << "Thread 1 start = " << 0 << ", end = " << portion << '\n';
  // std::cout << "Thread 2 start = " << portion << ", end = " << 2 * portion
  //           << '\n';
  // std::cout << "Thread 3 start = " << 2 * portion << ", end = " << 3 *
  // portion
  //           << '\n';
  // std::cout << "Thread 4 start = " << 3 * portion << ", end = " << m_row_cnt
  //           << '\n';

  std::thread t1([&] { multiplyThreaded(res, other, 0, portion); });
  std::thread t2([&] { multiplyThreaded(res, other, portion, 2 * portion); });
  std::thread t3(
      [&] { multiplyThreaded(res, other, 2 * portion, 3 * portion); });
  multiplyThreaded(res, other, 3 * portion, m_row_cnt);
  t1.join();
  t2.join();
  t3.join();

#else
  FPTYPE sum = 0;
  const int row_cnt = m_row_cnt;
  const int col_cnt = m_col_cnt;
  const int other_col_cnt = other.m_col_cnt;

#if CHANGE_LOOPS_ORDER

  for (int row = 0; row < row_cnt; ++row) {
    for (int pos = 0; pos < col_cnt; ++pos) {
      sum = 0;
      for (int col = 0; col < other_col_cnt; ++col) {
        res.m_array[row * other_col_cnt + col] +=
            m_array[row * col_cnt + pos] *
            other.m_array[pos * other_col_cnt + col];
      }
    }
  }

#else
  for (int row = 0; row < row_cnt; ++row) {
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

#endif

  return res;
}

Vector Matrix::operator*(const Vector &other) const {}

Matrix Matrix::operator*(FPTYPE multiplier) const {
  Matrix res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] * multiplier;
    }
  }

  return res;
}

Matrix Matrix::operator~() const {
  Matrix res(m_col_cnt, m_row_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i * m_col_cnt + j] = m_array[j * m_col_cnt + i];
    }
  }

  return res;
}

Matrix Matrix::operator^(int power) const {
  Matrix res(m_row_cnt, m_col_cnt, false);

  if (power == -1) {
  } else if (power == 0) {
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
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

FPTYPE *Matrix::operator[](int row_ind) {
  return &m_array[row_ind * m_col_cnt];
}

int Matrix::getRowCount() const { return m_row_cnt; }

int Matrix::getColCount() const { return m_col_cnt; }

void Matrix::display() const {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i * m_col_cnt + j] << "  ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

void Matrix::multiplyThreaded(const Matrix &res, const Matrix &other, int start,
                              int end) const {
  FPTYPE sum = 0;
  const int col_cnt = m_col_cnt;
  const int other_col_cnt = other.m_col_cnt;
#if CHANGE_LOOPS_ORDER
  for (int row = start; row < end; ++row) {
    for (int pos = 0; pos < col_cnt; ++pos) {
      sum = 0;
      for (int col = 0; col < other_col_cnt; ++col) {
        res.m_array[row * other_col_cnt + col] +=
            m_array[row * col_cnt + pos] *
            other.m_array[pos * other_col_cnt + col];
      }
    }
  }

#else
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

  uint32_t size = m_row_cnt * m_col_cnt;
  m_array = new FPTYPE[size];
  memcpy(m_array, other.m_array, sizeof(FPTYPE) * size);
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

}  // namespace MyAlgebra