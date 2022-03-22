#include <iostream>
#include "modular.hpp"
#include "vector_polynomial.hpp"
#include "polynomial.hpp"

int main(int argc, const char * argv[]) {
//    write menu
    
//    auto res = inf_secure::primes_n(2000);
//    std::for_each(res.cbegin(), res.cend(), [](auto val){ std::cout << val << " "; });
    using namespace inf_secure;
//
    std::cout << "======VECTOR|POLYNOM=====" << std::endl;
    VectorPoly1D_type<int> oo(7);
    std::cout << oo << std::endl;
    VectorPoly1D_type a { VectorPoly1D_type{5,6,7,8},  VectorPoly1D_type{5,-6} };
    VectorPoly1D_type b { VectorPoly1D_type{4}, VectorPoly1D_type{3,3}, VectorPoly1D_type{1}, VectorPoly1D_type{0} };
    a -= b;
    std::cout << a.deg() << std::endl;
    std::cout << a << std::endl;

    VectorPoly1D_type j1 = {-4, 0, -2, 1}, j2 = {-3, 1};
    std::cout << j1 << std::endl << j2 << std::endl;
    auto [qu, re] = j1 / j2;
    std::cout << qu << std::endl << re << std::endl;
    std::cout << qu * j2 + re << std::endl;

    std::cout << "======|POLYNOM|=====" << std::endl;
    Poly1D_type<Poly1D_type<int>> pol({
        {0, Poly1D_type<int>{{0, 5},{1, 6},{2, 7},{3, 8}}},
        {1, Poly1D_type<int>{{0, 5},{1, -6}}}
    });
    Poly1D_type<Poly1D_type<int>> pol_b = {
        {0, Poly1D_type<int>{{0, 4}}},
        {1, Poly1D_type<int>{{0,3},{1,3}}},
        {2, Poly1D_type<int>{{0,1}}},
        {3, Poly1D_type<int>{{0,0}}}
    };
    pol -= pol_b;
    std::cout << pol.deg() << std::endl;
    std::cout << pol << std::endl;

    Poly1D_type<int> j11 = {{0, -4}, {2, -2}, {3, 1} }, j21 = {{0, -3}, {1, 1}};
    std::cout << j11 << std::endl << j21 << std::endl;
    auto [qu1, re1] = j11 / j21;
    std::cout << qu1 << std::endl << re1 << std::endl;
    std::cout << j21 * qu1 + re1 << std::endl;
    
    VectorPoly1D_type<int> j_1(7), j_2(7);
    std::cout << "Vanila multipl:" << j_1 * j_2 << std::endl;
    std::cout << "karatsuba:" << j_1.karatsuba(j_2, 2) << std::endl;
    
    return 0;
}

