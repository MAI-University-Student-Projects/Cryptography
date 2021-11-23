#ifndef __PRIMES_HPP__
#define __PRIMES_HPP__

#include <type_traits>
#include <vector>
#include <cmath>

namespace inf_secure {

template<typename T>
    requires std::is_integral_v<T> && (!std::is_same_v<T, bool>)
std::vector<T> primes_n(T num) {
    num = (num > 2) ? num : 3;
    std::vector<bool> sieve(num, true);
    sieve[0] = sieve[1] = false;

    int sq_num = std::sqrt(num);

    for (size_t i = 2; i < sq_num + 1; ++i) {
        if(sieve[i]) {
            for(size_t j = i * i; j < num + 1; j += i)
                sieve[j] = false;
        }
    }

    std::vector<T> res;
    res.reserve(num / log(num)); // estimation for amount of primes before n
    for(size_t i = 0; i < sieve.size(); ++i) {
        if(sieve[i]) res.push_back(i);
    }
    return res;
}

}

#endif
