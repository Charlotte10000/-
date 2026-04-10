import random
from sympy import isprime

class Zmod:
    def __init__(self, q):
        self.q = q

    def cardinality(self):
        return self.q

    def __call__(self, x):
        return x % self.q

    def random_element(self):
        return random.randrange(0, self.q)
def is_square(a, q):
    if a % q == 0:
        return True
    return pow(a, (q - 1) // 2, q) == 1
class SQRT:

    # The initialization method computes all relevant tables. In a
    # practical implementation, this would be done once before compilation,
    # and the resulting table values would be hardcoded as read-only data.
    def __init__(self, K, w, leaf_w=None):
        # Get the field order q, and split it into q = m*2^n + 1, with m odd.
        self.K = K

        q = K.cardinality()
        self.q = q
        assert (q & 1) == 1
        assert q >= 3
        m = q - 1
        n = 0
        while (m & 1) == 0:
            m >>= 1
            n += 1
        self.m = m
        self.n = n

        # Find g, a primitive 2^n-th root of 1; we find a non-square and
        # raise if to the power m.
        while True:
            g = K.random_element()
            if is_square(g, q):
                continue
            g = pow(g, m, q)

            # ✅ 核心检查：是否真的是 primitive 2^n root
            if pow(g, 1 << (n - 1), q) == q - 1:
                break

        # Compute gpp[i] = g^i for i = 0 to n-1. We check that gpp[n-1] = -1
        # (this confirms that g is indeed a primitive 2^n-th root).
        gpp = []
        gpp.append(g)
        t = g
        for j in range(1, n):
            t = (t*t) % self.q
            gpp.append(t)
        assert gpp[n - 1] == q-1
        self.gpp = gpp

        # Precompute powers of g for our window size:
        #   gw[i][j] = g^(j*2^(i*w))  for 0 <= j < 2^w, and 0 <= i*w < n
        # We ensure that the window is not larger than n.
        if w > n:
            w = n
        assert w >= 1
        gw = []
        i = 0
        while i < n:
            gwt = []
            t = 1
            gwt.append(t)
            for j in range(1, 1 << w):
                t = (t*gpp[i]) %self.q
                gwt.append(t)
            gw.append(gwt)
            i += w
        self.w = w
        self.gw = gw

        if leaf_w is None:
            leaf_w = w
        rll = {}
        rll[0] = 0
        s = pow(g,  (1 << (n - leaf_w)),q)
        t = 1
        rll[t] = 0
        if leaf_w < n:
            u = pow(g , ((1 << n) - (1 << (n - leaf_w - 1))),q)
        else:
            u = pow(g , ((1 << n) - (1 << (n - leaf_w))),q)
        v = 1
        fhl = []
        fhl.append(v)
        for i in range(1, 1 << leaf_w):
            t = (t*s)%self.q
            if leaf_w < n or (i & 1) == 1:
                v = (v*u)%self.q
            rll[t] = i
            fhl.append(v)
        self.rll = rll
        self.fhl = fhl
        self.leaf_w = leaf_w

        # We keep track of the cost of the last operation (in multiplications
        # and squarings, not counting the initial exponentiation).
        self.costM = 0
        self.costS = 0

    def __call__(self, x):
        # Make sure the value is a field element.
        x = self.K(x)

        # v = x^((m-1)/2)
        # w = x*v
        # h = w*v = x^m
        v = (x ** ((self.m - 1) >> 1))%self.q
        w = (x * v)%self.q
        h = (w * v)%self.q

        # We initialize the cost accounting. We do not include the cost
        # of computing the initial exponentation (v).
        self.costM = 2
        self.costS = 0

        e, d = self.solve_dlp_pow2(0, h, True)

        y = (w * d) %self.q
        self.costM += 1

        # If x is not a square then we replace the value with None.
        y = self.SELECT(y, None, x != 0 and (e & 1) != 0)
        return y

    # Return a1 if ctl == True, or a0 if ctl == False.
    def SELECT(self, a0, a1, ctl):
        # CT: in a constant-time implementation, ctl is secret and this
        # should use a constant-time selection.
        if ctl:
            return a1
        else:
            return a0

    def GPOW(self, i, e, elen):

        e &= (1 << elen) - 1
        w = self.w
        wm = (1 << w) - 1
        ri = i % w
        i = i // w
        if ri != 0:
            e <<= ri
            elen += ri
        t = self.gw[i][e & wm]
        while True:
            elen -= w
            if elen <= 0:
                break
            e >>= w
            i += 1
            t = (t*self.gw[i][e & wm]) % self.q
            self.costM += 1
        return t

    def GPOW_cost(self, i, elen):
        if elen == 0:
            return 0
        w = self.w
        ri = i % w
        if ri != 0:
            elen += ri
        return ((elen + w - 1) // w) - 1

    def solve_dlp_pow2(self, i, h, ret_d, helper=None):
        # Base is gpp[i] and has order exactly 2^lb.
        n = self.n
        w = self.w
        leaf_w = self.leaf_w
        lb = n - i

        if lb <= leaf_w:

            e = self.rll.get(h)
            if e is None:
                e = 1
            else:
                e >>= (leaf_w - lb)
            if ret_d:
                d = self.fhl[e << (leaf_w - lb)]
            else:
                d = None
            return (e, d)

        # Split the order.
        lb0 = lb >> 1
        lb1 = lb - lb0

        hlp = None
        nlb1 = lb1 - (lb1 >> 1)
        if helper is None:

            cost_hlp = 1 + self.GPOW_cost(i + nlb1, lb0)
            do_hlp = True
            if lb1 <= leaf_w or cost_hlp >= nlb1:
                do_hlp = False

            h0 = h
            for j in range(0, lb1):
                if do_hlp and j == nlb1:
                    hlp = h0
                h0 =(h0 * h0)%self.q
            self.costS += lb1
        else:
            h0 = helper
        e0, _ = self.solve_dlp_pow2(i + lb1, h0, False)

        if ret_d:
            if i == 0:
                f = self.GPOW(0, ((1 << lb0) - e0 + 1) >> 1, lb0 - 1)
                f = self.SELECT(f, self.gpp[lb0 - 1], e0 == 0)
            else:
                f = self.GPOW(i - 1, (1 << lb0) - e0, lb0)
                f = self.SELECT(f, self.gpp[i - 1 + lb0], e0 == 0)
            h1 = (f ** 2)%self.q
            self.costS += 1

        else:
            h1 = self.GPOW(i, (1 << lb0) - e0, lb0)
            h1 = self.SELECT(h1, self.gpp[i + lb0], e0 == 0)

        h1 =(h1 * h)%self.q
        self.costM += 1
        if not (hlp is None):
            hlp1 = self.GPOW(i + nlb1, (1 << lb0) - e0, lb0)
            hlp1 = self.SELECT(hlp1, self.gpp[i + nlb1 + lb0], e0 == 0)
            hlp =(hlp* hlp1)%self.q
            self.costM += 1
        e1, d1 = self.solve_dlp_pow2(i + lb0, h1, ret_d, helper=hlp)

        if ret_d:
            # We here have 1/d1^2 = (b^(2^lb0))^e1 = h1 = h*f^2.
            # Thus: 1/(d1*f)^2 = h
            #
            # If i == 0 and e0 is odd, then we used h1 = h*g*f^2, and we thus
            # have: 1/(d1*f)^2 = h*g = g^(e+1)
            # In that case, d1*f = g^(-(e+1)/2)
            d = (d1 * f)%self.q
            self.costM += 1
        else:
            d = None

        # Since we used h1 = h*b^(2^lb0 - e0) instead of h1 = h/b^e0, the
        # obtained e1 must be decremented as a corrective action.
        e1 = (e1 - 1) & ((1 << lb1) - 1)

        return e0 + (e1 << lb0), d

    def last_cost(self):
        return self.costS, self.costM

    # Run a self-test and report the cost (S,M) of each square root.
    def self_test(self):
        assert self(0) == 0

        for i in range(0, 10):
            x = (self.K.random_element() ** 2)%self.q
            y = self(x)
            assert y is not None
            assert (y * y) % self.q == x
            if (y * y) % self.q == x:
                x = (x*self.gpp[0])%self.q
            assert self(x) is None
        return self.last_cost()


# Run tests and get costs for various degrees and window sizes.
def mkstats():
    nn = [32,64, 96, 128,256,512]

    for n in nn:
        print(f'n = {n}')
        q = (1 << n) + 1

        while not isprime(q):
            q += 1 << (n + 1)


        K = Zmod(q)

        for w in [2,4,8,16]:
            S = SQRT(K, w)
            S.self_test()
            costS, costM = S.last_cost()
            print(f'   w = {w:2d}: {costS:4d} S + {costM:4d} M   total: {costS + costM:4d}')

def minhc(K, w, R=None):
    q = K.cardinality()
    assert (q & 1) == 1
    assert q >= 3
    assert q.is_prime()
    m = q - 1
    n = 0
    while (m & 1) == 0:
        m >>= 1
        n += 1
    assert w <= n
    while True:
        g = K.random_element()
        if not (g.is_square()):
            break
    g = pow(g , m ,q)
    g = pow(g , (1 << (n - w)) ,q)
    rr = []
    t_init = 1
    if not (R is None):
        t_init = (t_init*R)%q
    t = t_init
    rr.append(int(t))
    for i in range(1, 1 << w):
        t = (t*g)%q
        rr.append(int(t))
    assert rr[1 << (w - 1)] == int(-t_init)
    t = (t*g)%q
    assert t == t_init

    qlen = len(q.bits())
    min_i = 0
    min_j = len(q.bits())
    for i in range(0, max(qlen - w, 0)):
        for j in range(w, min_j):
            tt = []
            jm = (1 << j) - 1
            for k in range(0, 1 << w):
                tt.append((rr[k] >> i) & jm)
            tt.sort()
            uu = True
            for k in range(1, 1 << w):
                if tt[k - 1] == tt[k]:
                    uu = False
                    break
            if uu:
                min_i = i
                min_j = j
                break
    return min_i, min_j

mkstats()