#ifndef __MODULAR_HPP__
#define __MODULAR_HPP__

#include <type_traits>
#include <limits>
#include <tuple>

namespace inf_secure {

template<typename T>
concept Integer = std::is_integral_v<T> && (!std::is_same_v<T, bool>);

template<Integer T>
constexpr T mod(T a, T b) { 
    T m = a % b;
    return (m > 0) ? m : (b > 0) ? (m + b) : (m - b); 
}

// right-to-left
template<Integer T>
constexpr T fast_pow_mod(T val, T ord, T m) {
    T res = 1;
    while(ord > 0) {
        if(ord & 1) res = mod(res * val, m);
        val = mod(val * val, m);
        ord = ord >> 1;
    }
    return res;
}

template<Integer T>
constexpr std::tuple<T, std::make_signed_t<T>, std::make_signed_t<T>> ext_gcd(T a, T b)
{
    using T_ = std::make_signed_t<T>;
    if constexpr(std::is_unsigned_v<T>) {
        if(a / b > std::numeric_limits<T_>::max())
            throw "overflow threat";
    }
    
    T r = a, r_ = b; 
    T_ x = 1, x_ = 0, y = 0, y_ = 1;
    while(r_ > 0) {
        x = x - (r / r_) * x_; std::swap(x, x_);
        y = y - (r / r_) * y_; std::swap(y, y_);
        r = r % r_; std::swap(r, r_);
    }
    
    return {r, x, y};
}

template<Integer T>
constexpr T pow_mod(T val, T ord, T m) {
    if(ord < 0) {
        auto [d, x, y] = ext_gcd(m, val);
        return (d == 1) ? fast_pow_mod(y, -ord, m) : throw "uninvertable";
    }
    return fast_pow_mod(val, ord, m);
}

template<Integer T>
constexpr T sum_mod(T lhs, T rhs, T m) { return mod(mod(lhs, m) + mod(rhs, m), m); }

template<Integer T>
constexpr T sub_mod(T lhs, T rhs, T m) { return mod(mod(lhs, m) - mod(rhs, m), m); }

template<Integer T>
constexpr T mul_mod(T lhs, T rhs, T m) { return mod(mod(lhs, m) * mod(rhs, m), m); }

template<Integer T>
constexpr T div_mod(T lhs, T rhs, T m) { 
    T rhs_1 = pow_mod(rhs, -1, m);
    return mul_mod(lhs, rhs_1, m);
}

}

#endif
