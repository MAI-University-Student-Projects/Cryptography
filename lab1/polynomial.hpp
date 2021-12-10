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

    int deg() const { return this->is_zero() ? -1 : coefs_.size() - 1; } 

    ccoef_it begin() const { return coefs_.cbegin(); }
    ccoef_it end() const { return coefs_.cend(); }

    _Coef_type& operator[](const size_t i) { return i > this->deg() ? coefs_[coefs_.size() - 1] : coefs_[i]; }
    const _Coef_type& operator[](const size_t i) const { return i > this->deg() ? coefs_[coefs_.size() - 1] : coefs_[i]; }

    bool is_zero() const { return *this == initials<Poly1D_type>::zero; }

    static Poly1D_type zero_polynom() { return initials<Poly1D_type>::zero; }
    
    friend bool operator==(const Poly1D_type& lhs, const Poly1D_type& rhs) = default;

};

template<class T>
Poly1D_type<T> operator+(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) {
    lhs += rhs;
    return lhs;
}

template<class T>
Poly1D_type<T> operator-(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) {
    lhs -= rhs;
    return lhs;
}

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
