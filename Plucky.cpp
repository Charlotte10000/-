#include <iostream>
#include <vector>
#include <optional>
#include <random>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;
using BigInt = cpp_int;

using std::optional;
using std::nullopt;
using std::cout;
using std::endl;
using std::vector;

// --- 模运算 ---
BigInt mod_add(const BigInt& a, const BigInt& b, const BigInt& q) {
    return (a + b) % q;
}
BigInt mod_sub(const BigInt& a, const BigInt& b, const BigInt& q) {
    BigInt res = (a - b) % q;
    if (res < 0) res += q;
    return res;
}
BigInt mod_mul(const BigInt& a, const BigInt& b, const BigInt& q) {
    return (a * b) % q;
}

// 快速幂
BigInt mod_pow(BigInt base, BigInt exp, const BigInt& q) {
    BigInt res = 1;
    base %= q;
    while (exp > 0) {
        if ((exp & 1) != 0) res = mod_mul(res, base, q);
        base = mod_mul(base, base, q);
        exp >>= 1;
    }
    return res;
}

// 逆元（q 为素数）
BigInt mod_inv(const BigInt& a, const BigInt& q) {
    return mod_pow(a, q - 2, q);
}

// --- 随机数生成 ---
BigInt random_bigint(const BigInt& q) {
    static std::random_device rd;
    static std::mt19937_64 gen(rd());

    BigInt r = 0;
    for (int i = 0; i < 4; i++) { // 4*64 = 256 bits
        r <<= 64;
        r += gen();
    }
    return r % (q - 1) + 1;
}

// --- Plucky 算法 ---
optional<BigInt> plucky_sqrt(const BigInt& x, const BigInt& q) {
    if (x == 0) return BigInt(0);

    while (true) {
        BigInt a = random_bigint(q);

        if (mod_mul(a, a, q) == mod_sub(0, x, q)) continue;

        BigInt denom = mod_add(mod_mul(a, a, q), x, q);
        BigInt inv = mod_inv(denom, q);

        BigInt P = mod_mul(2,
            mod_mul(mod_sub(mod_mul(a, a, q), x, q), inv, q), q);

        BigInt V = 2, W = P;

        // q 的二进制展开
        vector<int> b;
        BigInt temp = q;
        while (temp > 0) {
            b.push_back((temp & 1) ? 1 : 0);
            temp >>= 1;
        }
        std::reverse(b.begin(), b.end());

        for (size_t i = 0; i < b.size() - 1; i++) {
            if (b[i] == 0) {
                BigInt newV = mod_sub(mod_mul(V, V, q), 2, q);
                BigInt newW = mod_sub(mod_mul(V, W, q), P, q);
                V = newV; W = newW;
            }
            else {
                BigInt newV = mod_sub(mod_mul(V, W, q), P, q);
                BigInt newW = mod_sub(mod_mul(W, W, q), 2, q);
                V = newV; W = newW;
            }

            if (V == 0) {
                return mod_mul(
                    mod_mul(W, mod_add(mod_mul(a, a, q), x, q), q),
                    mod_inv((4 * a) % q, q),
                    q
                );
            }
        }

        if (V == 2 || V == q - 2) continue;
        return nullopt;
    }
}

int main() {
    BigInt q = BigInt((cpp_int(1) << 224) - (cpp_int(1) << 96) + 1); // 2^224 - 2^96 + 1
    BigInt x = q-1;

    auto y = plucky_sqrt(x, q);

    if (y.has_value()) {
        cout << "sqrt = " << y.value() << endl;
        cout << "other root = " << q - y.value() << endl;
    }
    else {
        cout << "not a quadratic residue" << endl;
    }
}