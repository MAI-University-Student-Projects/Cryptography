#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <vector>
#include <type_traits>
#include <numeric>

namespace inf_secure {

template<class _Coef_type>
class Poly1D_type {
    std::vector<_Coef_type> coefs_;
    using ccoef_it = typename std::vector<_Coef_type>::const_iterator;

    template<class U>
    struct initials { inline static const U zero = U{0}; inline static const U one = U{1}; };

    template<class T>
    struct initials<Poly1D_type<T>> { 
        inline static const Poly1D_type<T> zero = { initials<T>::zero }; 
        inline static const Poly1D_type<T> one = { initials<T>::one }; 
    };

    void cut_senior_zero_coefs() {
        for(auto it = coefs_.crbegin(); it != coefs_.crend() - 1; ++it){ 
            if(*it == initials<_Coef_type>::zero) coefs_.pop_back();
            else break;
        } 
    }

    Poly1D_type() = default; 
    Poly1D_type(std::vector<_Coef_type> cf_vec) : coefs_(std::move(cf_vec)) {}
    
public:
    explicit Poly1D_type(size_t deg) : coefs_(deg + 1, initials<_Coef_type>::one) { } 
    Poly1D_type(std::initializer_list<_Coef_type> cf_l) : coefs_(cf_l) { this->cut_senior_zero_coefs(); }
    
    Poly1D_type& operator+=(const Poly1D_type& oth) {
        if(oth.coefs_.size() > coefs_.size()) coefs_.resize(oth.coefs_.size(), initials<_Coef_type>::zero);
        std::transform(oth.coefs_.cbegin(), oth.coefs_.cend(), coefs_.cbegin(), coefs_.begin(), std::plus<>());
        this->cut_senior_zero_coefs();
        return *this;
    }
    
    Poly1D_type& operator-=(const Poly1D_type& oth) {
        if(oth.coefs_.size() > coefs_.size()) coefs_.resize(oth.coefs_.size(), initials<_Coef_type>::zero);
        std::transform(oth.coefs_.cbegin(), oth.coefs_.cend(), coefs_.cbegin(), coefs_.begin(), [](auto& oth_, auto& th_) { return th_ - oth_; });
        this->cut_senior_zero_coefs();
        return *this;
    }

    Poly1D_type& operator*=(const Poly1D_type& oth) {
        *this = *this * oth;
        return *this;
    }

    int deg() const { return this->is_zero() ? -1 : coefs_.size() - 1; } 

    ccoef_it begin() const { return coefs_.cbegin(); }
    ccoef_it end() const { return coefs_.cend(); }

    _Coef_type& operator[](const size_t i) { return coefs_[i]; }
    const _Coef_type& operator[](const size_t i) const { return coefs_[i]; }

    bool is_zero() const { return *this == initials<Poly1D_type>::zero; }

    static Poly1D_type zero_polynom() { return initials<Poly1D_type>::zero; }
    
    friend bool operator==(const Poly1D_type& lhs, const Poly1D_type& rhs) = default;

    friend Poly1D_type operator*(const Poly1D_type& lhs, const Poly1D_type& rhs) {
        if(lhs.deg() * rhs.deg() < 0) return initials<Poly1D_type>::zero;
        std::vector<_Coef_type> res_cf(lhs.coefs_.size() + rhs.coefs_.size() - 1, initials<_Coef_type>::zero); //if coef_type is polynom of smkindof degree
        for(size_t i = 0; i < lhs.coefs_.size(); ++i) {
            for(size_t j = 0; j < rhs.coefs_.size(); ++j)
                res_cf[i + j] += lhs[i] * rhs[j];
        }
        return Poly1D_type(std::move(res_cf));
    }    

    friend auto operator/(const Poly1D_type& lhs, const Poly1D_type& rhs) {
        Poly1D_type quotient = Poly1D_type::zero_polynom();
        Poly1D_type remainder = lhs;
        _Coef_type tmp_cf;
        while(!remainder.is_zero() && remainder.coefs_.size() >= rhs.coefs_.size()) {
            tmp_cf = remainder.coefs_.back() / rhs.coefs_.back();
            auto tmp_size = remainder.coefs_.size() - rhs.coefs_.size();
            Poly1D_type tmp_poly(std::vector<_Coef_type>(tmp_size + 1, initials<_Coef_type>::zero));
            tmp_poly.coefs_.back() = tmp_cf;
            quotient += tmp_poly;
            remainder -= tmp_poly * rhs;
        }

        return std::make_pair(quotient, remainder);
    }   

};

template<class T>
Poly1D_type<T> operator+(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) { return lhs += rhs; }

template<class T>
Poly1D_type<T> operator-(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) { return lhs -= rhs; }

template<class T>
std::ostream& operator<<(std::ostream& os, const Poly1D_type<T>& poly) {
    for(size_t degr = 0; auto& coef : poly) {
        degr == 0 ? os << "(" << coef << ")" : os << " + (" << coef << ")x^" << degr;
        ++degr;
    }
    return os;
}

}

#endif
