#include "CMtx.h"

#include <cstdlib>
#include <iostream>

namespace MyAlgebra {

const FPTYPE CMtx::ALG_PRECISION = 10e-6;

CMtx::CMtx(size_t row_cnt, size_t col_cnt, bool rand_init)
    : m_row_cnt(row_cnt),
      m_col_cnt(col_cnt),
      m_array(new FPTYPE[m_row_cnt * m_col_cnt]) {
  if (rand_init) {
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        m_array[i * m_col_cnt + j] = rand() % 15;
      }
    }
  }
}

CMtx::CMtx(size_t row_cnt, FPTYPE diagonal)
    : m_row_cnt(row_cnt),
      m_col_cnt(row_cnt),
      m_array(new FPTYPE[m_row_cnt * m_row_cnt]) {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_row_cnt; ++j) {
      if (i != j)
        m_array[i * m_col_cnt + j] = 0;
      else
        m_array[i * m_col_cnt + j] = diagonal;
    }
  }
}

CMtx::CMtx(const CMtx &other)
    : m_row_cnt(other.m_row_cnt),
      m_col_cnt(other.m_col_cnt),
      m_array(new FPTYPE[m_row_cnt * m_col_cnt]) {
  // TODO: Change shallow copy to deep copy?
  memcpy(m_array, other.m_array, sizeof(FPTYPE) * m_row_cnt * m_col_cnt);
}

CMtx::CMtx(CMtx &&other)
    : m_row_cnt(other.m_row_cnt),
      m_col_cnt(other.m_col_cnt),
      m_array(other.m_array) {
  other.m_row_cnt = 0;
  other.m_col_cnt = 0;
  other.m_array = nullptr;
}

CMtx::~CMtx() { delete[] m_array; }

const CMtx &CMtx::operator=(const CMtx &other) {
  if (this != &other) {
    delete[] m_array;

    m_row_cnt = other.m_row_cnt;
    m_col_cnt = other.m_col_cnt;

    size_t tmp_size = m_row_cnt * m_col_cnt;

    m_array = new FPTYPE[tmp_size];
    // TODO: Change shallow copy to deep copy?
    memcpy(m_array, other.m_array, sizeof(FPTYPE) * tmp_size);
  }

  return *this;
}

const CMtx &CMtx::operator=(const FPTYPE diagonal) {
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

const CMtx &CMtx::operator=(CMtx &&other) {
  if (this != &other) {
    delete[] m_array;

    m_row_cnt = other.m_row_cnt;
    m_col_cnt = other.m_col_cnt;
    m_array = other.m_array;

    other.m_row_cnt = 0;
    other.m_col_cnt = 0;
    other.m_array = nullptr;
  }

  return *this;
}

CMtx CMtx::operator-() const {
  CMtx res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      // TODO: Change for other types
      res.m_array[tmp_idx] = m_array[tmp_idx] * -1;
    }
  }

  return res;
}

CMtx CMtx::operator-(const CMtx &other) const {
  // TODO: Check size, should be the same

  CMtx res(m_row_cnt, m_col_cnt, false);

  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] - other.m_array[tmp_idx];
    }
  }

  return res;
}

CMtx CMtx::operator~() const {
  CMtx res(m_col_cnt, m_row_cnt, false);

  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      res.m_array[i * m_col_cnt + j] = m_array[j * m_col_cnt + i];
    }
  }

  return res;
}

FPTYPE CMtx::determinant() const {
  if (m_row_cnt != m_col_cnt) return 0;
}

CVct CMtx::operator*(const CVct &other) const {}

CMtx CMtx::operator*(const CMtx &other) const {
  // TODO: Check size
  CMtx res(m_row_cnt, other.m_col_cnt, false);

  for (int l_row = 0; l_row < m_row_cnt; ++l_row) {
    for (int r_col = 0; r_col < other.m_col_cnt; ++r_col) {
      FPTYPE sum = 0;
      int iter = 0;
      while (iter < m_row_cnt) {
        sum += m_array[l_row * m_col_cnt + iter] *
               other.m_array[iter * other.m_col_cnt + r_col];
        iter++;
      }

      res.m_array[l_row * other.m_col_cnt + r_col] = sum;
    }
  }

  return res;
}

CMtx CMtx::operator*(FPTYPE multiplier) const {
  CMtx res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] * multiplier;
    }
  }

  return res;
}

CMtx CMtx::operator+(const CMtx &other) const {
  // TODO: Check size, should be the same

  CMtx res(m_row_cnt, m_col_cnt, false);
  int tmp_idx;
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      tmp_idx = i * m_col_cnt + j;
      res.m_array[tmp_idx] = m_array[tmp_idx] + other.m_array[tmp_idx];
    }
  }

  return res;
}

CMtx CMtx::operator^(int power) const {
  CMtx res(m_row_cnt, m_col_cnt, false);

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

FPTYPE *CMtx::operator[](int row_ind) { return &m_array[row_ind * m_col_cnt]; }

bool CMtx::operator==(const CMtx &other) const {
  if (m_row_cnt != other.m_row_cnt || m_col_cnt != other.m_col_cnt)
    return false;

  if (this != &other) {
    int tmp_idx;
    for (int i = 0; i < m_row_cnt; ++i) {
      for (int j = 0; j < m_col_cnt; ++j) {
        tmp_idx = i * m_col_cnt + j;
        if (m_array[tmp_idx] - other.m_array[tmp_idx] > ALG_PRECISION)
          return false;
      }
    }
  }

  return true;
}

void CMtx::display() const {
  for (int i = 0; i < m_row_cnt; ++i) {
    for (int j = 0; j < m_col_cnt; ++j) {
      std::cout << m_array[i * m_col_cnt + j] << " ";
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

CMtx operator*(FPTYPE multiplier, const CMtx &other) {
  return other * multiplier;
}

}  // namespace MyAlgebra