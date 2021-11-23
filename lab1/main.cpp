#include <iostream>
#include <algorithm>
#include "primes.hpp"
#include "modular.hpp"

int main(int argc, const char * argv[]) {
//    write menu
    
//    auto res = inf_secure::primes_n(2000);
//    std::for_each(res.cbegin(), res.cend(), [](auto val){ std::cout << val << " "; });
    static_assert(inf_secure::fast_pow_mod(249, 321, 499) == 447);
        
    static_assert(inf_secure::sum_mod(1723345, 2124945, 11) == 6);
    static_assert(inf_secure::sub_mod(1723345, 2124945, 11) == 10);
    static_assert(inf_secure::mul_mod(1723345, 2124945, 11) == 6);

    auto [d, x, y] = inf_secure::ext_gcd(26, 11);
    std::cout << d << " " << x << " " << y << " " << inf_secure::mod(26, 11);
    
    static_assert(inf_secure::div_mod(145, 11, 26) == 25); //145 * 19 mod 26
    return 0;
}

