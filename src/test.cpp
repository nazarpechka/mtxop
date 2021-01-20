
// My implementation
#include "matrix.h"


int main() {
  Matrix double_mat = Matrix<double>::readFromFile("../doc/mat_float_double_example.txt");
  Matrix int_mat = Matrix<int>::readFromFile("../doc/mat_int_example.txt");

  std::cout << "Double example:\n";
  double_mat.display();

  std::cout << "Int example:\n";
  int_mat.display();

  std::cout << "Second column:\n";
  int_mat.extractColumn(1).display();
}