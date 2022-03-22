#ifndef __VPOLYNOMIAL_HPP__
#define __VPOLYNOMIAL_HPP__

#include <vector>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <iostream>

namespace inf_secure {

template<class _Coef_type>
class VectorPoly1D_type {
    using coef_vector = typename std::vector<_Coef_type>;
    using ccoef_it = typename coef_vector::const_iterator;
    
    coef_vector coefs_;
    
    template<class U>
    struct initials {
        inline static const U zero = U{0};
        inline static const U one = U{1};
    };

    template<class T>
    struct initials<VectorPoly1D_type<T>> {
        inline static const VectorPoly1D_type<T> zero = { initials<T>::zero };
        inline static const VectorPoly1D_type<T> one = { initials<T>::one };
    };
    
    size_t inline static const sz_from_degree(size_t deg) { return deg + 2; }
    size_t inline static const dg_from_size(size_t size) { return size - 2; }

    void cut_senior_zero_coefs() {
        for(auto it = coefs_.crbegin(), end_it = coefs_.crend() - 1; it != end_it; ++it){
            if(*it == initials<_Coef_type>::zero) coefs_.pop_back();
            else break;
        }
    }
    
//    template<std::random_access_iterator InputIt, std::random_access_iterator OutputIt> //c++20 iterator concepts
    template<class InputIt, class OutputIt>
    static void karatsuba_impl(
        InputIt first1, InputIt last1,
        InputIt first2, InputIt last2, OutputIt first_out,
        size_t limit) {
        
        size_t n = std::distance(first1, last1);
        if(n <= limit) {
            for(size_t i = 0; i < n; ++i) {
                for(size_t j = 0; j < n; ++j)
                    *(first_out + i + j) += *(first1 + i) * *(first2 + j);
            }
        }
        else {
            auto k = n / 2; auto half1 = first1 + k, half2 = first2 + k;

            coef_vector a, b; a.reserve(k+1), b.reserve(k+1);
            std::transform(first1, half1, half1, std::back_inserter(a), std::plus<>());
            std::transform(first2, half2, half2, std::back_inserter(b), std::plus<>());

            coef_vector t(n, initials<_Coef_type>::zero);
            karatsuba_impl(a.cbegin(), a.cend(), b.cbegin(), b.cend(), t.begin(), limit);
            karatsuba_impl(first1, half1, first2, half2, first_out, limit); //p1
            karatsuba_impl(first1, half1, first2, half2, first_out + n, limit); //p2

            auto out_it_1 = first_out + k; //p1
            auto out_it_2 = first_out + n; //p2
            

            for(size_t i = 0; i < k; ++i) {
                auto kn1 = *(out_it_1 + i) + t[i] - *(first_out + i) - *(out_it_2 + i);
                auto nnk2 = *(out_it_2 + i) + t[k + i] - *(out_it_1 + i) - *(out_it_2 + k + i); //trail zero needed
                *(out_it_1 + i) = kn1;
                *(out_it_2 + i) = nnk2;
            }
        }
    }
    
    VectorPoly1D_type() = default;
    
public:
    explicit VectorPoly1D_type(size_t deg) {
        coefs_.reserve(sz_from_degree(deg));
        coefs_.emplace_back(initials<_Coef_type>::zero);
        std::fill_n(std::back_inserter(coefs_), sz_from_degree(deg) - 1, initials<_Coef_type>::one);
    }
    
    VectorPoly1D_type(std::initializer_list<_Coef_type> cf_l) {
        coefs_.reserve(cf_l.size() + 1);
        coefs_.emplace_back(initials<_Coef_type>::zero);
        coefs_.insert(coefs_.cend(), cf_l.begin(), cf_l.end());
        this->cut_senior_zero_coefs();
    }
    
//reminder: coefs[0] <-> zero_polynom tag
    VectorPoly1D_type(coef_vector cf_vec) {
        coefs_.reserve(cf_vec.size() + 1);
        coefs_.emplace_back(initials<_Coef_type>::zero);
        std::move(cf_vec.begin(), cf_vec.end(), std::back_inserter(coefs_));
    }
    
    VectorPoly1D_type& operator+=(const VectorPoly1D_type& oth) {
        if(oth.coefs_.size() > coefs_.size()) coefs_.resize(oth.coefs_.size(), initials<_Coef_type>::zero);
        std::transform(oth.coefs_.cbegin(), oth.coefs_.cend(), coefs_.cbegin(), coefs_.begin(), std::plus<>());
        this->cut_senior_zero_coefs();
        return *this;
    }

    VectorPoly1D_type& operator-=(const VectorPoly1D_type& oth) {
        if(oth.coefs_.size() > coefs_.size()) coefs_.resize(oth.coefs_.size(), initials<_Coef_type>::zero);
        std::transform(oth.coefs_.cbegin(), oth.coefs_.cend(), coefs_.cbegin(), coefs_.begin(), [](auto& oth_, auto& th_) { return th_ - oth_; });
        this->cut_senior_zero_coefs();
        return *this;
    }

//reminder: coefs[0] <-> zero_polynom tag
    VectorPoly1D_type& operator*=(const VectorPoly1D_type& oth) {
        coef_vector res_cf(coefs_.size() + oth.coefs_.size() - 2, initials<_Coef_type>::zero);
        for(size_t i = 0; i < coefs_.size() - 1; ++i) {
            for(size_t j = 0; j < oth.coefs_.size() - 1; ++j)
                res_cf[i + j + 1] += coefs_[i + 1] * oth.coefs_[j + 1];
        }
        coefs_ = std::move(res_cf);
        return *this;
    }

    ccoef_it begin() const { return coefs_.cbegin() + 1; }
    ccoef_it end() const { return coefs_.cend(); }
    
    int deg() const { return this->is_zero() ? -1 : dg_from_size(coefs_.size()); }
    bool is_zero() const { return coefs_.size() == 1; }
    static VectorPoly1D_type zero_polynom() { return initials<VectorPoly1D_type>::zero; }
    
    friend bool operator==(const VectorPoly1D_type& lhs, const VectorPoly1D_type& rhs) = default;
    
//    reminder: coefs[0] <-> zero_polynom tag
    friend auto operator/(const VectorPoly1D_type& lhs, const VectorPoly1D_type& rhs) {
        VectorPoly1D_type quotient = VectorPoly1D_type::zero_polynom();
        VectorPoly1D_type remainder = lhs;
        _Coef_type tmp_cf;
        while(!remainder.is_zero() && remainder.coefs_.size() >= rhs.coefs_.size()) {
            tmp_cf = remainder.coefs_.back() / rhs.coefs_.back();
            auto tmp_size = remainder.coefs_.size() - rhs.coefs_.size();
            VectorPoly1D_type tmp_poly(coef_vector(tmp_size + 1, initials<_Coef_type>::zero));
            tmp_poly.coefs_.back() = tmp_cf;
            quotient += tmp_poly;
            remainder -= tmp_poly * rhs;
        }

        return std::make_pair(quotient, remainder);
    }

//    reminder: coefs[0] <-> zero_polynom tag
    VectorPoly1D_type karatsuba(const VectorPoly1D_type& oth, size_t limit = 64) {
        // assert(coefs_.size() == oth.coefs_.size())
        coef_vector res_cf(coefs_.size() + oth.coefs_.size() - 2, initials<_Coef_type>::zero);
        VectorPoly1D_type::karatsuba_impl(this->begin(), this->end(), oth.begin(), oth.end(), res_cf.begin(), limit);
        res_cf.pop_back(); //trail zero appearing
        return VectorPoly1D_type(std::move(res_cf));
    }

};

template<class T>
VectorPoly1D_type<T> operator+(VectorPoly1D_type<T> lhs, const VectorPoly1D_type<T>& rhs) { return lhs += rhs; }

template<class T>
VectorPoly1D_type<T> operator-(VectorPoly1D_type<T> lhs, const VectorPoly1D_type<T>& rhs) { return lhs -= rhs; }

template<class T>
VectorPoly1D_type<T> operator*(VectorPoly1D_type<T> lhs, const VectorPoly1D_type<T>& rhs) { return lhs *= rhs; }

template<class T>
std::ostream& operator<<(std::ostream& os, const VectorPoly1D_type<T>& poly) {
    os << "(" << *poly.begin() << ")";
    std::for_each(poly.begin() + 1, poly.end(), [&os, deg=1](const auto& val) mutable { os << " + (" << val << ")x^" << deg; ++deg; });
    return os;
}

}

#endif
