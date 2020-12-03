#include "CMtx.h"

#include <cstdlib>
#include <iostream>

namespace MyAlgebra {

CMtx::CMtx(size_t row_cnt, size_t col_cnt, bool rand_init)
    : row_cnt_(row_cnt), col_cnt_(col_cnt), row_ptr_(new FPTYPE *[row_cnt_]) {
  for (int i = 0; i < row_cnt_; ++i) {
    row_ptr_[i] = new FPTYPE[col_cnt_];
  }

  if (rand_init) {
    for (int i = 0; i < row_cnt_; ++i) {
      for (int j = 0; j < col_cnt_; ++j) {
        row_ptr_[i][j] = rand() % 100;
      }
    }
  }
}

CMtx::CMtx(size_t row_cnt, FPTYPE diagonal)
    : row_cnt_(row_cnt), col_cnt_(row_cnt), row_ptr_(new FPTYPE *[row_cnt_]) {
  for (int i = 0; i < row_cnt_; ++i) {
    row_ptr_[i] = new FPTYPE[row_cnt_];

    for (int j = 0; j < row_cnt_; ++j) {
      if (i != j)
        row_ptr_[i][j] = 0;
      else
        row_ptr_[i][j] = diagonal;
    }
  }
}

CMtx::CMtx(const CMtx &rhs)
    : row_cnt_(rhs.row_cnt_),
      col_cnt_(rhs.col_cnt_),
      row_ptr_(new FPTYPE *[row_cnt_]) {
  for (int i = 0; i < row_cnt_; ++i) {
    row_ptr_[i] = new FPTYPE[col_cnt_];
    // TODO: Change shallow copy to deep copy?
    memcpy(row_ptr_[i], rhs.row_ptr_[i], sizeof(FPTYPE) * col_cnt_);
  }
}

CMtx::~CMtx() {
  for (int i = 0; i < row_cnt_; ++i) {
    delete[] row_ptr_[i];
  }
  delete[] row_ptr_;
}

const CMtx &CMtx::operator=(const CMtx &rhs) {
  if (this == &rhs) return *this;

  for (int i = 0; i < row_cnt_; ++i) {
    delete[] row_ptr_[i];
  }
  delete[] row_ptr_;

  row_cnt_ = rhs.row_cnt_;
  col_cnt_ = rhs.col_cnt_;

  row_ptr_ = new FPTYPE *[row_cnt_];
  for (int i = 0; i < row_cnt_; ++i) {
    row_ptr_[i] = new FPTYPE[col_cnt_];
    // TODO: Change shallow copy to deep copy?
    memcpy(row_ptr_[i], rhs.row_ptr_[i], sizeof(FPTYPE) * col_cnt_);
  }

  return *this;
}

const CMtx &CMtx::operator=(const FPTYPE diagonal) {
  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      if (i != j)
        row_ptr_[i][j] = 0;
      else
        row_ptr_[i][j] = diagonal;
    }
  }

  return *this;
}

// const CMtx &CMtx::operator=(CMtx &&rhs) {}

CMtx CMtx::operator-() const {
  CMtx res(row_cnt_, col_cnt_, false);

  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      // TODO: Change for other types
      res.row_ptr_[i][j] = row_ptr_[i][j] * -1;
    }
  }

  return res;
}

CMtx CMtx::operator-(const CMtx &rhs) const {
  // TODO: Check size, should be the same

  CMtx res(row_cnt_, col_cnt_, false);

  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      res.row_ptr_[i][j] = row_ptr_[i][j] - rhs.row_ptr_[i][j];
    }
  }

  return res;
}

CMtx CMtx::operator~() const {
  CMtx res(col_cnt_, row_cnt_, false);

  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      res.row_ptr_[i][j] = row_ptr_[j][i];
    }
  }

  return res;
}

CVct CMtx::operator*(const CVct &rhs) const {}

CMtx CMtx::operator*(const CMtx &rhs) const {
  // TODO: Check size
  CMtx res(row_cnt_, rhs.col_cnt_, false);

  for (int l_row = 0; l_row < row_cnt_; ++l_row) {
    for (int r_col = 0; r_col < rhs.col_cnt_; ++r_col) {
      FPTYPE sum = 0;
      int iter = 0;
      while (iter < row_cnt_) {
        sum += row_ptr_[l_row][iter] * rhs.row_ptr_[iter][r_col];
        iter++;
      }

      res.row_ptr_[l_row][r_col] = sum;
    }
  }

  return res;
}

CMtx CMtx::operator*(FPTYPE multiplier) const {
  CMtx res(row_cnt_, col_cnt_, false);

  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      res.row_ptr_[i][j] = row_ptr_[i][j] * multiplier;
    }
  }

  return res;
}

CMtx CMtx::operator+(const CMtx &rhs) const {
  // TODO: Check size, should be the same

  CMtx res(row_cnt_, col_cnt_, false);

  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      res.row_ptr_[i][j] = row_ptr_[i][j] + rhs.row_ptr_[i][j];
    }
  }

  return res;
}

CMtx CMtx::operator^(int power) const {}

FPTYPE *CMtx::operator[](int row_ind) { return row_ptr_[row_ind]; }

bool CMtx::operator==(const CMtx &&rhs) const {}

void CMtx::display() const {
  for (int i = 0; i < row_cnt_; ++i) {
    for (int j = 0; j < col_cnt_; ++j) {
      std::cout << row_ptr_[i][j] << " ";
    }
    std::cout << '\n';
  }
}

CMtx operator*(FPTYPE multiplier, const CMtx &rhs) { return rhs * multiplier; }

}  // namespace MyAlgebra