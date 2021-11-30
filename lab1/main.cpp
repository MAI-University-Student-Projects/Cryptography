#include <iostream>
#include <algorithm>
#include "primes.hpp"
#include "modular.hpp"

int main(int argc, const char * argv[]) {
//    write menu
    
//    auto res = inf_secure::primes_n(2000);
//    std::for_each(res.cbegin(), res.cend(), [](auto val){ std::cout << val << " "; });

    static_assert(fast_pow_mod(249, 321, 499) == 447);
    constexpr Zn_type<11> a{1723345l};
    constexpr Zn_type<11> b{2124945l};
    static_assert(a + b == 6);
    static_assert(a - b == 10);
    static_assert(a * b == 6);
    
    constexpr Zn_type<26> f{11};
    constexpr Zn_type<26> f1{145};
    static_assert(f.invert().get_remainder() == 19);
    static_assert(f1 / f == 25);

    return 0;
}

