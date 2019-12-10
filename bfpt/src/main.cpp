#include <iostream>
#include<extensions/adaptors.hpp>
int main() {
  std::cout << "a" << std::endl;
  int x = 10;
  x = 30;

  #ifdef NDEBUG
  x=10;
  #endif
}