#include <iostream>
#include <vector>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using namespace boost::multiprecision;

// ---------------- 基本运算 ----------------
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

cpp_int mod_sub(const cpp_int& a, const cpp_int& b, const cpp_int& p) {
    cpp_int r = (a - b) % p;
    if (r < 0) r += p;
    return r;
}

// ---------------- 随机 ----------------
cpp_int rand_mod(const cpp_int& p) {
    static random_device rd;
    static mt19937_64 gen(rd());

    cpp_int r = 0;
    for (int i = 0; i < 4; i++) {
        r <<= 64;
        r += gen();
    }
    return r % p;
}

// ---------------- 二进制展开 ----------------
vector<int> to_binary(const cpp_int& n) {
    vector<int> bits;
    cpp_int t = n;
    while (t > 0) {
        bits.push_back((t & 1) ? 1 : 0);
        t >>= 1;
    }
    reverse(bits.begin(), bits.end());
    return bits;
}

// ---------------- KKK sqrt ----------------
cpp_int KKK_sqrt(const cpp_int& x, const cpp_int& q) {

    while (true) {
        // STEP 1: 随机 a
        cpp_int a = rand_mod(q);
        if (a == 0) continue;

        // N(theta) = a^2 + x
        cpp_int aa = (a * a) % q;
        cpp_int norm = (aa + x) % q;
        if (norm == 0) continue;

        // b = 2a / (a^2 + x)
        cpp_int inv_norm = modinv(norm, q);
        cpp_int b = (2 * a % q) * inv_norm % q;
        a = ((a * a %q - x)% q * inv_norm) % q;
        // 初始化
        cpp_int V = 2 % q;
        cpp_int W = (2 * a) % q;

        // n = (q-1)/2^m
        cpp_int temp = q - 1;
        int m = 0;
        while (temp % 2 == 0) {
            temp /= 2;
            m++;
        }
        cpp_int n = temp;

        vector<int> bits = to_binary(n);

        // 主循环
        for (int j = 0; j < (int)bits.size(); j++) {
            cpp_int V0 = V;
            cpp_int W0 = W;

            if (bits[j] == 0) {
                // W = V W - 2a
                cpp_int t1 = (V0 * W0) % q;
                W = mod_sub(t1, 2 * a % q, q);

                // V = V^2 - 2
                cpp_int t2 = (V0 * V0) % q;
                V = mod_sub(t2, 2, q);
            }
            else {
                // V = V W - 2a
                cpp_int t1 = (V0 * W0) % q;
                V = mod_sub(t1, 2 * a % q, q);

                // W = W^2 - 2
                cpp_int t2 = (W0 * W0) % q;
                W = mod_sub(t2, 2, q);
            }
        }

        // 判定
        if (V == 2 || V == q - 2) continue;

        // 情况1：V = 0
        if (V == 0) {
            cpp_int numerator = (2 * b % q) * x % q;
            cpp_int invW = modinv(W, q);
            return numerator * invW % q;
        }

        // 情况2：循环平方
        while (V != 0) {
            cpp_int V0 = V;
            cpp_int W0 = W;

            cpp_int t1 = (V0 * W0) % q;
            W = mod_sub(t1, 2 * a % q, q);

            cpp_int t2 = (V0 * V0) % q;
            V = mod_sub(t2, 2, q);

            if (V == 0) {
                cpp_int numerator = ( V0 % q) * b % q * x % q;

                cpp_int denom = (a * V0) % q;
                denom = mod_sub(denom, W0, q);

                cpp_int inv_denom = modinv(denom, q);

                return numerator * inv_denom % q;
            }
        }
    }
}

// ---------------- 测试 ----------------
int main() {
    cpp_int q = (cpp_int(1) << 224) - (cpp_int(1) << 96) + 1;
    cpp_int Q = q-1;


    cpp_int root = KKK_sqrt(Q, q);

    cout << "root = " << root << " or " << (q-root)%q << endl;
    cout << "check = " << (root * root % q) << endl;

    return 0;
}