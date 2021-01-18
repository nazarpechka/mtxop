#include "matrix.h"

#include <immintrin.h>

#include <cmath>
#include <iostream>
#include <thread>

#define SWAP_LOOPS 0

namespace MyAlgebra {

const float Matrix::ALG_PRECISION = 10e-3f;

Matrix::Matrix(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new float[m_row_cnt * m_col_cnt]) {
  if (rand_init) {
    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      m_array[i] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      //      m_array[i] = rand() % 10;
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

    for (size_t i = 0; i < m_row_cnt * m_col_cnt; ++i) {
      if (std::abs(m_array[i] - other.m_array[i]) > ALG_PRECISION) return false;
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

  const Matrix transposed(~other);

  // Multithreading is only efficient for large matrices
  if (m_row_cnt > 127) {
    auto portion =
        static_cast<unsigned int>(ceil(static_cast<float>(m_row_cnt) / 4.0f));

    //    std::thread t1([&] { multiply(res, other, 0, portion); });
    //    std::thread t2([&] { multiply(res, other, portion, 2 * portion); });
    //    std::thread t3([&] { multiply(res, other, 2 * portion, 3 * portion);
    //    }); std::thread t4([&] { multiply(res, other, 3 * portion, m_row_cnt);
    //    });

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
    std::thread t4([&] {
      multiply(res, transposed, 3 * portion, m_row_cnt, m_col_cnt,
               other.m_col_cnt);
    });

    t1.join();
    t2.join();
    t3.join();
    t4.join();
  } else {
    multiply(res, transposed, 0, m_row_cnt, m_col_cnt, other.m_col_cnt);
  }

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

std::string Matrix::authorName() { return "Nazar_Pechevystyi"; }

void Matrix::multiply(const Matrix &res, const Matrix &other, size_t start,
                      size_t end, size_t col_cnt, size_t other_col_cnt) const {
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

  //  __m256 _row, _other_row, _mult_res, _res_row, _add_res;
  //
  //  for (size_t row = start; row < end; ++row) {
  //    for (size_t pos = 0; pos < col_cnt; ++pos) {
  //      _row = _mm256_set1_ps(m_array[row * col_cnt + pos]);
  //      for (size_t col = 0; col < other_col_cnt; col += 8) {
  //
  //        if (other_col_cnt - col < 8) {
  //          for (size_t tmp_col = col; tmp_col < other_col_cnt; ++tmp_col) {
  //            res.m_array[row * other_col_cnt + tmp_col] +=
  //                m_array[row * col_cnt + pos] *
  //                other.m_array[pos * other_col_cnt + tmp_col];
  //          }
  //
  //        } else {
  //          _other_row =
  //              _mm256_load_ps(&other.m_array[pos * other_col_cnt + col]);
  //
  //          _mult_res = _mm256_mul_ps(_row, _other_row);
  //
  //          _res_row = _mm256_load_ps(&res.m_array[row * other_col_cnt +
  //          col]);
  //
  //          _add_res = _mm256_add_ps(_res_row, _mult_res);
  //
  ////          _mm256_store_ps(&res.m_array[row * other_col_cnt + col],
  ///_add_res);
  //          for (size_t i = 0; i < 8; ++i) {
  //            res.m_array[row * other_col_cnt + col + i] = _add_res[i];
  //          }
  //        }
  //      }
  //    }
  //  }

#else
  // Naive multiplication

//    float sum = 0;
//    for (size_t row = start; row < end; ++row) {
//      for (size_t other_row = 0; other_row < other_col_cnt; ++other_row) {
//        sum = 0;
//        for (size_t pos = 0; pos < col_cnt; ++pos) {
//          std::cout << "Dot product << this[" << row << "][" << pos << "] and other[" << other_row << "][" << pos << "]\n";
//          sum += m_array[row * col_cnt + pos] *
//                 other.m_array[other_row * col_cnt + pos];
//        }
//
//        std::cout << "Writing to res[" << row << "][" << other_row << "]\n";
//
//        res.m_array[row * other_col_cnt + other_row] = sum;
//      }
//    }

  __m256 _row, _other_row, _mult_res, _sum_res;

  for (size_t row = start; row < end; ++row) {
    for (size_t other_row = 0; other_row < other_col_cnt; ++other_row) {
      _sum_res = _mm256_setzero_ps();

      for (size_t col = 0; col < col_cnt; col += 8) {
          if (col_cnt - col < 8) {
//              std::cout << "Finishing up with " << col_cnt - col << " floats\n";
              float sum = 0;
              for (size_t col_tmp = col; col_tmp < col_cnt; ++col_tmp) {
                  sum += m_array[row * col_cnt + col_tmp] *
                         other.m_array[other_row * col_cnt + col_tmp];
              }
//              std::cout << "Adding to res[" << row << "][" << other_row << "], " << sum << "\n";
              res.m_array[row * other_col_cnt + other_row] += sum;

          } else {
//            std::cout << "Grabbing 8 floats from this[" << row << "][" << col << "] and other[" << other_row << "][" << col << "]\n";
            _row = _mm256_loadu_ps(&m_array[row * col_cnt + col]);

            _other_row =
                _mm256_loadu_ps(&other.m_array[other_row * col_cnt + col]);
            _mult_res = _mm256_mul_ps(_row, _other_row);
            _sum_res = _mm256_add_ps(_sum_res, _mult_res);
          }
      }

      //        _mm256_store_ps(&res.m_array[row * other_col_cnt + col],
      //        _sum_res);
//      std::cout << "Adding to res[" << row << "][" << other_row << "]\n";
      for (size_t i = 0; i < 8; ++i) {
        res.m_array[row * other_col_cnt + other_row] += _sum_res[i];
      }
    }
  }

  // Naive
  //  float sum = 0;
  //  for (size_t row = start; row < end; ++row) {
  //    for (size_t col = 0; col < other_col_cnt; ++col) {
  //      sum = 0;
  //      for (size_t pos = 0; pos < col_cnt; ++pos) {
  //        sum += m_array[row * col_cnt + pos] *
  //               other.m_array[pos * other_col_cnt + col];
  //      }
  //
  //      res.m_array[row * other_col_cnt + col] = sum;
  //    }
  //  }

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