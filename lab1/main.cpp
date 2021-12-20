#include <iostream>
#include "modular.hpp"
#include "polynomial.hpp"

int main(int argc, const char * argv[]) {
//    write menu
    
//    auto res = inf_secure::primes_n(2000);
//    std::for_each(res.cbegin(), res.cend(), [](auto val){ std::cout << val << " "; });
    using namespace inf_secure;
    Poly1D_type a { Poly1D_type{5,6,7,8},  Poly1D_type{5,-6} };
    Poly1D_type b { Poly1D_type{4}, Poly1D_type{3,3}, Poly1D_type{1}, Poly1D_type{0} };
    a -= b;
    std::cout << a.deg() << std::endl;
    std::cout << a << std::endl;
    
    Poly1D_type j1 = {-4, 0, -2, 1}, j2 = {-3, 1};
    std::cout << j1 << std::endl << j2 << std::endl;
    auto [qu, re] = j1 / j2;
    std::cout << qu << std::endl << re << std::endl;
    std::cout << j1 * j2 << std::endl;
    
    return 0;
}

