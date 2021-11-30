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
constexpr T fast_pow_mod(T val, T power, const T& m) {
    T res = 1;
    while(power > 0) {
        if (power & 1) res = mod(res * val, m);
        val = mod(val * val, m);
        power = power >> 1;
    }
    return res;
}

template<size_t Modulus_>
    requires (Modulus_ > 0)
class Zn_type final {

    size_t value_ = 0;

    constexpr size_t extend_gcd_modulus_(size_t a, size_t b, long long& x, long long& y) const {
        if(a == 0) {
            x = 0;
            y = 1;
            return b;
        }
        long long x_tmp, y_tmp;
        size_t gcd = extend_gcd_modulus_(b % a, a, x_tmp, y_tmp);
        x = y_tmp - (b / a) * x_tmp;
        y = x_tmp;
        return gcd;
    }

public:
    constexpr Zn_type() = default;

    template<Integer T>
    explicit constexpr Zn_type(const T& val) : value_(val < 0 ? Modulus_ == 1 ? 0 : Modulus_ - ((-val) % Modulus_) : val % Modulus_) {}

    constexpr Zn_type invert() const {
        long long x, y;
        return extend_gcd_modulus_(value_, Modulus_, x, y) == 1 ? Zn_type{x} : throw "uninvertable"; 
    }

    constexpr Zn_type pow(long long ord) const {
        return ord < 0 ? Zn_type{ fast_pow_mod(value_, static_cast<size_t>(-ord), Modulus_) } 
                        : Zn_type{ fast_pow_mod(value_, static_cast<size_t>(ord), Modulus_) };
    }

    constexpr Zn_type& operator+=(const Zn_type& oth) { 
        if((value_ += oth.value_) >= Modulus_) value_ -= Modulus_;
        return *this;
    }

    constexpr Zn_type& operator-=(const Zn_type& oth) {
        //avoiding possible overflow with value_<unsigned> -= oth.value
        if(value_ < oth.value_) value_ += Modulus_; 
        value_ -= oth.value_;
        return *this;
    }

    constexpr Zn_type& operator*=(const Zn_type& oth) {
        value_ *= oth.value_;
        value_ %= Modulus_;
        return *this;
    }

    constexpr size_t get_remainder() const { return value_; }
};

template<size_t Mod_>
constexpr size_t operator+(Zn_type<Mod_> lhs, const Zn_type<Mod_>& rhs) { lhs += rhs; return lhs.get_remainder(); }

template<size_t Mod_>
constexpr size_t operator-(Zn_type<Mod_> lhs, const Zn_type<Mod_>& rhs) { lhs -= rhs; return lhs.get_remainder(); }

template<size_t Mod_>
constexpr size_t operator*(Zn_type<Mod_> lhs, const Zn_type<Mod_>& rhs) { lhs *= rhs; return lhs.get_remainder(); }

template<size_t Mod_>
constexpr size_t operator/(Zn_type<Mod_> lhs, const Zn_type<Mod_>& rhs) { lhs *= rhs.invert(); return lhs.get_remainder(); }


#endif
