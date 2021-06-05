#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <exception>

void sieve_eratosthenes(size_t m) {
    std::vector<bool> is_prime(m, true);
    std::ofstream f_out("eratosth_output.txt", std::ios::out);
//    p * q = M, p <= q -> max p: p * p = M -> p <= sqrt(M)
    for(int i = 2; i * i < m; ++i) {
        if(!is_prime[i])
            continue;
        for(int j = i * i; j < m; j += i)
            is_prime[j] = false;
    }
    
    for(int i = 0; i < is_prime.size(); ++i) {
        if(is_prime[i])
            f_out << i << ' ';
    }
}

int pow_mod(int a, int n, int m) {
    if(n < 1 || m < 1)
        throw std::invalid_argument("non-natural argument");
    int res = 1;
    while(n > 0) {
//        проверку на четность можно реализовать по последнему биту двоичной записи (2^0 ==? 1)
        if(n & 1) {
//            при нечетной степени и при первой
            res = res * a % m;
            n--;
        }
        else {
            a = a * a % m;
//            >> ~ /2 (см перемещ '1' вправо в двоич. записи); соотв << ~ *2
            n = n >> 1;
        }
    }
    return res;
}

int main(int argc, const char * argv[]) {
//    write menu
//    sieve_eratosthenes(100000);
    std::cout << pow_mod(5, 3, 11) << std::endl;
    return 0;
}

