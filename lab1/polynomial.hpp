#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

template<class _Coef_type>
class Poly1D_type {
    using coef_map = typename std::map<int, _Coef_type, std::greater<int> >;
    using ccoef_it = typename coef_map::const_iterator;

    coef_map coefs_;

    template<class _U>
    struct initials {
        inline static const _U zero = _U{0};
        inline static const _U one = _U{1};
    };

    template<class _T>
    struct initials< Poly1D_type<_T> > {
        inline static const Poly1D_type<_T> zero = { {-1, initials<_T>::zero} };
        inline static const Poly1D_type<_T> one = { {-1, initials<_T>::zero}, {0, initials<_T>::one} };
    };

    void cut_zero_coefs() {
        for(auto it = coefs_.cbegin(), end_it = std::prev(coefs_.cend()); it != end_it; ){ 
            if(it->second == initials<_Coef_type>::zero) it = coefs_.erase(it);
            else ++it;
        } 
    }

public:
    explicit Poly1D_type(size_t deg) {
        auto [it, succs] = coefs_.emplace(-1, initials<_Coef_type>::zero);
        for(int i = 0; i < deg + 1; ++i) { it = coefs_.emplace_hint(it, i, initials<_Coef_type>::one); }
    }

// key < -1!
    Poly1D_type(std::initializer_list< std::pair<const int, _Coef_type> > cf_l) : coefs_(cf_l) {
        coefs_.emplace_hint(coefs_.end(), -1, initials<_Coef_type>::zero);
        this->cut_zero_coefs(); 
    }

    Poly1D_type& operator+=(const Poly1D_type& oth) {
        for(const auto& [key_oth, val_oth] : oth) {
            auto [it_tmp, success] = coefs_.insert({key_oth, val_oth});
            if(!success) coefs_.insert_or_assign(it_tmp, key_oth, it_tmp->second + val_oth);
        }
        this->cut_zero_coefs(); 
        return *this;
    }

    Poly1D_type& operator-=(const Poly1D_type& oth) {
        for(const auto& [key_oth, val_oth] : oth) {
            auto [it_tmp, success] = coefs_.insert({key_oth, initials<_Coef_type>::zero - val_oth});
            if(!success) coefs_.insert_or_assign(it_tmp, key_oth, it_tmp->second - val_oth);
        }
        this->cut_zero_coefs();
        return *this;
    }

    Poly1D_type& operator*=(const Poly1D_type& oth) {
        coef_map new_coefs_;
        for(auto it = coefs_.cbegin(), end_it = std::prev(coefs_.cend()); it != end_it; ++it) {
            for(const auto& [key_oth, val_oth] : oth) {
                auto [it_tmp, success] = new_coefs_.insert({it->first+key_oth, it->second * val_oth});
                if(!success) new_coefs_.insert_or_assign(it_tmp, it->first+key_oth, it_tmp->second + it->second * val_oth);
            }
        }
        coefs_ = std::move(new_coefs_);
        this->cut_zero_coefs();
        return *this;
    }

    Poly1D_type& operator*=(const _Coef_type& val) {
        for(auto it = coefs_.cbegin(), end_it = std::prev(coefs_.cend()); it != end_it; ++it)
            coefs_.insert_or_assign(it, it->first, it->second * val);
        return *this;
    }

    int deg() const { return coefs_.cbegin()->first; } 

    friend bool operator==(const Poly1D_type& lhs, const Poly1D_type& rhs) = default;
    
    static Poly1D_type constant_polynom() { return initials<Poly1D_type>::one; }
    static Poly1D_type zero_polynom() { return initials<Poly1D_type>::zero; }
    bool is_zero() const { return coefs_.cbegin()->first == -1; }

    ccoef_it begin() const { return coefs_.cbegin(); }
    ccoef_it end() const { return coefs_.cend(); }

};

template<class T>
std::ostream& operator<<(std::ostream& os, const Poly1D_type<T>& poly) {
    os << "[" << poly.begin()->second << "]x^" << poly.begin()->first;
    std::for_each(std::next(poly.begin()), std::prev(poly.end()), [&os](const auto& val){ os << " + [" << val.second << "]x^" << val.first; });
    return os;
}

template<class T>
Poly1D_type<T> operator+(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) { return lhs += rhs; }

template<class T>
Poly1D_type<T> operator-(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) { return lhs -= rhs; }

template<class T>
Poly1D_type<T> operator*(Poly1D_type<T> lhs, const Poly1D_type<T>& rhs) { return lhs *= rhs; }

template<class T>
Poly1D_type<T> operator*(Poly1D_type<T> lhs, const T& rhs) { return lhs *= rhs; }

template<class T>
auto operator/(const Poly1D_type<T>& lhs, const Poly1D_type<T>& rhs) {
    Poly1D_type quotient = Poly1D_type<T>::zero_polynom();
    Poly1D_type remainder = lhs;
    T tmp_cf;
    while(!remainder.is_zero() && remainder.deg() >= rhs.deg()) {
        Poly1D_type<T> tmp_poly = { {remainder.deg() - rhs.deg(), remainder.begin()->second / rhs.begin()->second } };
        quotient += tmp_poly;
        remainder -= tmp_poly * rhs;
    }
    return std::make_pair(quotient, remainder);    
}

