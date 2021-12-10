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
    
    return 0;
}

