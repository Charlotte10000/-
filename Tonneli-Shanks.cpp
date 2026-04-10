#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using namespace boost::multiprecision;

// ---------------- 模幂 ----------------
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

// ---------------- 模逆 ----------------
cpp_int modinv(const cpp_int& a, const cpp_int& p) {
    return modexp(a, p - 2, p);
}

// ---------------- Legendre symbol ----------------
int legendre(const cpp_int& a, const cpp_int& p) {
    if (a % p == 0) return 0;
    cpp_int r = modexp(a, (p - 1) / 2, p);
    if (r == 1) return 1;
    if (r == p - 1) return -1;
    return 0;
}

// ---------------- 256-bit 随机 ----------------
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

// ---------------- Tonelli-Shanks ----------------
cpp_int tonelli_shanks(const cpp_int& a, const cpp_int& q) {

    if (a == 0) return 0;

    // 1. 判断是否为二次剩余
    if (legendre(a, q) != 1) {
        throw runtime_error("not a quadratic residue");
    }

    // 2. 分解 q-1 = 2^m * n
    cpp_int temp = q - 1;
    int m = 0;
    while ((temp & 1) == 0) {
        temp >>= 1;
        m++;
    }
    cpp_int n = temp;

    // 3. 找非二次剩余 z
    cpp_int z;
    while (true) {
        z = rand_mod(q);
        if (z == 0) continue;
        if (legendre(z, q) == -1) break;
    }

    // 4. 初始化
    cpp_int x = modexp(z, n, q);
    cpp_int t = modexp(a, n, q);
    cpp_int R = modexp(a, (n + 1) / 2, q);
    int k = m;

    // 5. 主循环
    while (t != 1) {
        // 找最小 i 使得 t^{2^i} = 1
        int i = 0;
        cpp_int tmp = t;
        while (tmp != 1) {
            tmp = (tmp * tmp) % q;
            i++;
            if (i == k) {
                throw runtime_error("failure in TS");
            }
        }

        // b = x^{2^{k-i-1}}
        cpp_int b = x;
        for (int j = 0; j < k - i - 1; j++) {
            b = (b * b) % q;
        }

        // 更新
        R = (R * b) % q;

        cpp_int bb = (b * b) % q;
        t = (t * bb) % q;

        x = bb;
        k = i;
    }

    return R;
}

// ---------------- 测试 ----------------
int main() {
    cpp_int q = (cpp_int(1) << 224) - (cpp_int(1) << 96) + 1;
    cpp_int a = q-1;



    cpp_int root = tonelli_shanks(a, q);

    cout << "root = " << root << endl;
    cout << "check = " << (root * root % q) << endl;

    return 0;
}