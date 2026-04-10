#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using namespace boost::multiprecision;

// ---------------- 模运算 ----------------
cpp_int modexp(cpp_int base, cpp_int exp, const cpp_int& mod) {
    cpp_int res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return res;
}

cpp_int modinv(const cpp_int& a, const cpp_int& p) {
    return modexp(a, p - 2, p);
}

// ---------------- Legendre ----------------
int legendre(const cpp_int& a, const cpp_int& p) {
    if (a % p == 0) return 0;
    cpp_int res = modexp(a, (p - 1) / 2, p);
    if (res == 1) return 1;
    if (res == p - 1) return -1;
    return 0;
}

// ---------------- 256-bit 随机 ----------------
cpp_int rand_mod(const cpp_int& p) {
    static random_device rd;
    static mt19937_64 gen(rd());

    cpp_int r = 0;
    for (int i = 0; i < 4; i++) { // 4×64 = 256
        r <<= 64;
        r += gen();
    }
    return r % p;
}

// ---------------- 快速 Lucas V_n ----------------
cpp_int lucasV(cpp_int P, cpp_int n, const cpp_int& mod) {
    cpp_int V = 2;
    cpp_int Q = 1;

    cpp_int Vk = 2;
    cpp_int Vk1 = P;

    // 二进制展开
    vector<int> bits;
    while (n > 0) {
        bits.push_back((n & 1) ? 1 : 0);
        n >>= 1;
    }
    reverse(bits.begin(), bits.end());

    for (size_t i = 0; i <= bits.size()-2; i++) {
        // doubling
        cpp_int V2k = (Vk * Vk - 2) % mod;
        if (V2k < 0) V2k += mod;

        cpp_int V2k1 = (Vk * Vk1 - P) % mod;
        if (V2k1 < 0) V2k1 += mod;

        if (bits[i] == 0) {
            Vk = V2k;
            Vk1 = V2k1;
        }
        else {
            cpp_int tmp = (Vk1 * Vk1 - 2) % mod;
            if (tmp < 0) tmp += mod;

            Vk = V2k1;
            Vk1 = tmp;
        }
    }
    if (bits[bits.size()-1] == 0)
        return (Vk*Vk-2) % mod;
    else
        return (Vk*Vk1-P) % mod;
}

// ---------------- Müller sqrt ----------------
cpp_int muller_sqrt(cpp_int Q, const cpp_int& q) {
    if (Q == 0) return 0;
    if (Q == 4) return 2;

    cpp_int t;

    if (legendre(Q - 4, q) == -1) {
        t = 1;
    }
    else {
        while (true) {
            t = rand_mod(q);
            if (t == 0) continue;

            cpp_int val = (Q * t % q * t % q - 4) % q;
            if (val < 0) val += q;

            if (legendre(val, q) == -1) break;
        }
    }

    cpp_int P = (Q * t % q * t % q - 2) % q;
    if (P < 0) P += q;

    cpp_int k = (q - 1) / 4;

    cpp_int V = lucasV(P, k, q);

    cpp_int a = (V * modinv(t, q)) % q;

    return a;
}

// ---------------- 测试 ----------------
int main() {
    cpp_int q = (cpp_int(1) << 224) - (cpp_int(1) << 96) + 1; // 256-bit prime
    cpp_int Q = q-1;


    cpp_int root = muller_sqrt(Q, q);

    cout << "Q = " << Q << endl;
    cout << "sqrt = " << root << endl;
    cout << "check = " << (root * root % q) << endl;

    return 0;
}