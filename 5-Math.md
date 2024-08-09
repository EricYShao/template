# Math
## Binary exponentiation
```cpp
ll power(ll a, ll b){
  ll res = 1;
  for (; b; a = a * a % MOD, b >>= 1){
    if (b & 1) res = res * a % MOD;
  }
  return res;
}
```

## Matrix Exponentiation: $O(n^3 \log{b})$
```cpp
const int N = 100, MOD = 1e9 + 7;

struct matrix{
  ll m[N][N];
  int n;
  matrix(){
    n = N;
    memset(m, 0, sizeof(m));
  };
  matrix(int n_){
    n = n_;
    memset(m, 0, sizeof(m));
  };
  matrix(int n_, ll val){
    n = n_;
    memset(m, 0, sizeof(m));
    for (int i = 0; i < n; i++) m[i][i] = val;
  };

  matrix operator* (matrix oth){
    matrix res(n);
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        for (int k = 0; k < n; k++){
          res.m[i][j] = (res.m[i][j] + m[i][k] * oth.m[k][j]) % MOD;
        }
      }
    }
    return res;
  }
};

matrix power(matrix a, ll b){
  matrix res(a.n, 1);
  for (; b; a = a * a, b >>= 1){
    if (b & 1) res = res * a;
  }
  return res;
}
```

## Extended Euclidean Algorithm
+ $O(\max(\log a, \log b))$
+ Finds solution $(x, y)$ to $ax + by = \gcd(a, b)$
+ Can find all solutions given $(x_0, y_0): \forall k, a(x_0 + kb/g) + b(y_0 - ka/g) = \gcd(a, b)$.
```cpp
ll euclid(ll a, ll b, ll &x, ll &y) {
  if (!b) return x = 1, y = 0, a;
  ll d = euclid(b, a % b, y, x);
  return y -= a/b * x, d;
}
```

## CRT
+ $crt(a, m, b, n)$ computes $x$ such that $x\equiv a \pmod m$, $x\equiv b \pmod n$
+ If $|a| < m$ and $|b| < n$, $x$ will obey $0 \le x < \text{lcm}(m, n)$.
+ Assumes $mn < 2^{62}$.
+ $O(\max(\log m, \log n))$
```cpp
ll crt(ll a, ll m, ll b, ll n) {
  if (n > m) swap(a, b), swap(m, n);
  ll x, y, g = euclid(m, n, x, y);
  assert((a - b) % g == 0); // else no solution
  // can replace assert with whatever needed
  x = (b - a) % n * x % n / g * m + a;
  return x < 0 ? x + m*n/g : x;
}
```

## Linear Sieve

+ Mobius Function

```cpp
vector<int> prime;
bool is_composite[MAX_N];
int mu[MAX_N];

void sieve(int n){
  fill(is_composite, is_composite + n, 0);
  mu[1] = 1;
  for (int i = 2; i < n; i++){
    if (!is_composite[i]){
      prime.push_back(i);
      mu[i] = -1; //i is prime
      }
  for (int j = 0; j < prime.size() && i * prime[j] < n; j++){
    is_composite[i * prime[j]] = true;
    if (i % prime[j] == 0){
      mu[i * prime[j]] = 0; //prime[j] divides i
      break;
      } else {
      mu[i * prime[j]] = -mu[i]; //prime[j] does not divide i
      }
    }
  }
}
```



+ Euler's Totient Function

```cpp
vector<int> prime;
bool is_composite[MAX_N];
int phi[MAX_N];

void sieve(int n){
  fill(is_composite, is_composite + n, 0);
  phi[1] = 1;
  for (int i = 2; i < n; i++){
    if (!is_composite[i]){
      prime.push_back (i);
      phi[i] = i - 1; //i is prime
      }
  for (int j = 0; j < prime.size () && i * prime[j] < n; j++){
    is_composite[i * prime[j]] = true;
    if (i % prime[j] == 0){
      phi[i * prime[j]] = phi[i] * prime[j]; //prime[j] divides i
      break;
      } else {
      phi[i * prime[j]] = phi[i] * phi[prime[j]]; //prime[j] does not divide i
      }
    }
  }
}
```

## Gaussian Elimination

```cpp
bool is_0(Z v) { return v.x == 0; }
Z abs(Z v) { return v; }
bool is_0(double v) { return abs(v) < 1e-9; }

// 1 => unique solution, 0 => no solution, -1 => multiple solutions
template <typename T>
int gaussian_elimination(vector<vector<T>> &a, int limit) {
  if (a.empty() || a[0].empty()) return -1;
  int h = (int)a.size(), w = (int)a[0].size(), r = 0;
  for (int c = 0; c < limit; c++) {
    int id = -1;
    for (int i = r; i < h; i++) {
      if (!is_0(a[i][c]) && (id == -1 || abs(a[id][c]) < abs(a[i][c]))) {
        id = i;
      }
    }
    if (id == -1) continue;
    if (id > r) {
      swap(a[r], a[id]);
      for (int j = c; j < w; j++) a[id][j] = -a[id][j];
    }
    vector<int> nonzero;
    for (int j = c; j < w; j++) {
      if (!is_0(a[r][j])) nonzero.push_back(j);
    }
    T inv_a = 1 / a[r][c];
    for (int i = r + 1; i < h; i++) {
      if (is_0(a[i][c])) continue;
      T coeff = -a[i][c] * inv_a;
      for (int j : nonzero) a[i][j] += coeff * a[r][j];
    }
    ++r;
  }
  for (int row = h - 1; row >= 0; row--) {
    for (int c = 0; c < limit; c++) {
      if (!is_0(a[row][c])) {
        T inv_a = 1 / a[row][c];
        for (int i = row - 1; i >= 0; i--) {
          if (is_0(a[i][c])) continue;
          T coeff = -a[i][c] * inv_a;
          for (int j = c; j < w; j++) a[i][j] += coeff * a[row][j];
        }
        break;
      }
    }
  } // not-free variables: only it on its line
  for(int i = r; i < h; i++) if(!is_0(a[i][limit])) return 0;   
  return (r == limit) ? 1 : -1;
}

template <typename T>
pair<int,vector<T>> solve_linear(vector<vector<T>> a, const vector<T> &b, int w) {
  int h = (int)a.size();
  for (int i = 0; i < h; i++) a[i].push_back(b[i]);
  int sol = gaussian_elimination(a, w);
  if(!sol) return {0, vector<T>()};
  vector<T> x(w, 0);
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      if (!is_0(a[i][j])) {
        x[j] = a[i][w] / a[i][j];
        break;
      }
    }
  }
  return {sol, x};
}
```

## is_prime

+ (Miller–Rabin primality test)

```cpp
typedef __int128_t i128;

i128 power(i128 a, i128 b, i128 MOD = 1, i128 res = 1) {
  for (; b; b /= 2, (a *= a) %= MOD)
    if (b & 1) (res *= a) %= MOD;
  return res;
}

bool is_prime(ll n) {
  if (n < 2) return false;
  static constexpr int A[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
  int s = __builtin_ctzll(n - 1);
  ll d = (n - 1) >> s;
  for (auto a : A) {
    if (a == n) return true;
    ll x = (ll)power(a, d, n);
    if (x == 1 || x == n - 1) continue;
    bool ok = false;
    for (int i = 0; i < s - 1; ++i) {
      x = ll((i128)x * x % n);  // potential overflow!
      if (x == n - 1) {
        ok = true;
        break;
      }
    }
    if (!ok) return false;
  }
  return true;
}
```

```cpp
typedef __int128_t i128;

ll pollard_rho(ll x) {
  ll s = 0, t = 0, c = rng() % (x - 1) + 1;
  ll stp = 0, goal = 1, val = 1;
  for (goal = 1;; goal *= 2, s = t, val = 1) {
    for (stp = 1; stp <= goal; ++stp) {
      t = ll(((i128)t * t + c) % x);
      val = ll((i128)val * abs(t - s) % x);
      if ((stp % 127) == 0) {
        ll d = gcd(val, x);
        if (d > 1) return d;
      }
    }
    ll d = gcd(val, x);
    if (d > 1) return d;
  }
}

ll get_max_factor(ll _x) {
  ll max_factor = 0;
  function<void(ll)> fac = [&](ll x) {
    if (x <= max_factor || x < 2) return;
    if (is_prime(x)) {
      max_factor = max_factor > x ? max_factor : x;
      return;
    }
    ll p = x;
    while (p >= x) p = pollard_rho(x);
    while ((x % p) == 0) x /= p;
    fac(x), fac(p);
  };
  fac(_x);
  return max_factor;
}
```

## Berlekamp-Massey

+ Recovers any $n$-order linear recurrence relation from the first $2n$ terms of the sequence.
+ Input $s$ is the sequence to be analyzed.
+ Output $c$ is the shortest sequence $c_1, ..., c_n$, such that
$$s_m=\sum_{i=1}^{n} c_i \cdot s_{m-i} \text{, for all } m \ge n.$$
+ Be careful since $c$ is returned in 0-based indexation.
+ Complexity: $O(N^2)$

```cpp
vector<ll> berlekamp_massey(vector<ll> s) {
  int n = sz(s), l = 0, m = 1;
  vector<ll> b(n), c(n);
  ll ldd = b[0] = c[0] = 1;
  for (int i = 0; i < n; i++, m++) {
    ll d = s[i];
    for (int j = 1; j <= l; j++) d = (d + c[j] * s[i - j]) % MOD;
    if (d == 0) continue;
    vector<ll> temp = c;
    ll coef = d * power(ldd, MOD - 2) % MOD;
    for (int j = m; j < n; j++){
      c[j] = (c[j] + MOD - coef * b[j - m]) % MOD;
      if (c[j] < 0) c[j] += MOD;
    }
    if (2 * l <= i) {
      l = i + 1 - l;
      b = temp;
      ldd = d;
      m = 0;
    }
  }
  c.resize(l + 1);
  c.erase(c.begin());
  for (ll &x : c)
    x = (MOD - x) % MOD;
  return c;
}
```

## Calculating k-th term of a linear recurrence
+ Given the first $n$ terms $s_0, s_1, ..., s_{n-1}$ and the sequence $c_1, c_2, ..., c_n$ such that
$$s_m=\sum_{i=1}^{n} c_i \cdot s_{m-i} \text{, for all } m \ge n,$$
the function calc_kth computes $s_k$.
+ Complexity: $O(n^2 \log{k})$
```cpp
vector<ll> poly_mult_mod(vector<ll> p, vector<ll> q, vector<ll>& c){
  vector<ll> ans(sz(p) + sz(q) - 1);
  for (int i = 0; i < sz(p); i++){
    for (int j = 0; j < sz(q); j++){
      ans[i + j] = (ans[i + j] + p[i] * q[j]) % MOD;
    }
  }
  int n = sz(ans), m = sz(c);
  for (int i = n - 1; i >= m; i--){
    for (int j = 0; j < m; j++){
      ans[i - 1 - j] = (ans[i - 1 - j] + c[j] * ans[i]) % MOD;
    }
  }
  ans.resize(m);
  return ans;
}

ll calc_kth(vector<ll> s, vector<ll> c, ll k){
  assert(sz(s) >= sz(c)); // size of s can be greater than c, but not less
  if (k < sz(s)) return s[k];
  vector<ll> res{1};
  for (vector<ll> poly = {0, 1}; k; poly = poly_mult_mod(poly, poly, c), k >>= 1){
    if (k & 1) res = poly_mult_mod(res, poly, c);
  }
  ll ans = 0;
  for (int i = 0; i < min(sz(res), sz(c)); i++) ans = (ans + s[i] * res[i]) % MOD;
  return ans;
}
```

## Partition Function
+ Returns number of partitions of $n$ in $O(n^{1.5})$
```cpp
int partition(int n) {
  int dp[n + 1];
  dp[0] = 1;
  for (int i = 1; i <= n; i++) {
    dp[i] = 0;
    for (int j = 1, r = 1; i - (3 * j * j - j) / 2 >= 0; ++j, r *= -1) {
      dp[i] += dp[i - (3 * j * j - j) / 2] * r;
      if (i - (3 * j * j + j) / 2 >= 0) dp[i] += dp[i - (3 * j * j + j) / 2] * r;
    }
  }
  return dp[n];
}
```

## NTT

```cpp
void ntt(vector<ll>& a, int f) {
  int n = int(a.size());
  vector<ll> w(n);
  vector<int> rev(n);
  for (int i = 0; i < n; i++) rev[i] = (rev[i / 2] / 2) | ((i & 1) * (n / 2));
  for (int i = 0; i < n; i++) {
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  }
  ll wn = power(f ? (MOD + 1) / 3 : 3, (MOD - 1) / n);
  w[0] = 1;
  for (int i = 1; i < n; i++) w[i] = w[i - 1] * wn % MOD;
  for (int mid = 1; mid < n; mid *= 2) {
    for (int i = 0; i < n; i += 2 * mid) {
      for (int j = 0; j < mid; j++) {
        ll x = a[i + j], y = a[i + j + mid] * w[n / (2 * mid) * j] % MOD;
        a[i + j] = (x + y) % MOD, a[i + j + mid] = (x + MOD - y) % MOD;
      }
    }
  }
  if (f) {
    ll iv = power(n, MOD - 2);
    for (auto& x : a) x = x * iv % MOD;
  }
}
vector<ll> mul(vector<ll> a, vector<ll> b) {
  int n = 1, m = (int)a.size() + (int)b.size() - 1;
  while (n < m) n *= 2;
  a.resize(n), b.resize(n);
  ntt(a, 0), ntt(b, 0); // if squaring, you can save one NTT here
  for (int i = 0; i < n; i++) a[i] = a[i] * b[i] % MOD;
  ntt(a, 1);
  a.resize(m);
  return a;
}
```

## FFT

```cpp
const ld PI = acosl(-1);
auto mul = [&](const vector<ld>& aa, const vector<ld>& bb) {
  int n = (int)aa.size(), m = (int)bb.size(), bit = 1;
  while ((1 << bit) < n + m - 1) bit++;
  int len = 1 << bit;
  vector<complex<ld>> a(len), b(len);
  vector<int> rev(len);
  for (int i = 0; i < n; i++) a[i].real(aa[i]);
  for (int i = 0; i < m; i++) b[i].real(bb[i]);
  for (int i = 0; i < len; i++) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (bit - 1));
  auto fft = [&](vector<complex<ld>>& p, int inv) {
    for (int i = 0; i < len; i++)
      if (i < rev[i]) swap(p[i], p[rev[i]]);
    for (int mid = 1; mid < len; mid *= 2) {
      auto w1 = complex<ld>(cos(PI / mid), (inv ? -1 : 1) * sin(PI / mid));
      for (int i = 0; i < len; i += mid * 2) {
        auto wk = complex<ld>(1, 0);
        for (int j = 0; j < mid; j++, wk = wk * w1) {
          auto x = p[i + j], y = wk * p[i + j + mid];
          p[i + j] = x + y, p[i + j + mid] = x - y;
        }
      }
    }
    if (inv == 1) {
      for (int i = 0; i < len; i++) p[i].real(p[i].real() / len);
    }
  };
  fft(a, 0), fft(b, 0);
  for (int i = 0; i < len; i++) a[i] = a[i] * b[i];
  fft(a, 1);
  a.resize(n + m - 1);
  vector<ld> res(n + m - 1);
  for (int i = 0; i < n + m - 1; i++) res[i] = a[i].real();
  return res;
};
```

## MIT's FFT/NTT, Polynomial mod/log/exp Template
+ For integers rounding works if $(|a| + |b|)\max(a, b) < \mathtt{\sim} 10^9$, or in theory maybe $10^6$
+ $\frac{1}{P(x)}$ in $O(n \log n)$, $e^{P(x)}$ in $O(n \log n)$, $\ln (P(x))$ in $O(n \log n)$, $P(x)^k$ in $O(n \log n)$, Evaluates $P(x_1), \cdots, P(x_n)$ in $O(n \log^2 n)$, Lagrange Interpolation in $O(n \log^2 n)$
```cpp
// use #define FFT 1 to use FFT instead of NTT (default)
// Examples:
// poly a(n+1); // constructs degree n poly
// a[0].v = 10; // assigns constant term a_0 = 10
// poly b = exp(a);
// poly is vector<num>
// for NTT, num stores just one int named v
// for FFT, num stores two doubles named x (real), y (imag)

#define sz(x) ((int)x.size())
#define rep(i, j, k) for (int i = int(j); i < int(k); i++)
#define trav(a, x) for (auto &a : x)
#define per(i, a, b) for (int i = (b)-1; i >= (a); --i)
using ll = long long;
using vi = vector<int>;

namespace fft {
#if FFT
// FFT
using dbl = double;
struct num { 
  dbl x, y;
  num(dbl x_ = 0, dbl y_ = 0): x(x_), y(y_) {}
};
inline num operator+(num a, num b) {
  return num(a.x + b.x, a.y + b.y);
}
inline num operator-(num a, num b) {
  return num(a.x - b.x, a.y - b.y);
}
inline num operator*(num a, num b) {
  return num(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
inline num conj(num a) { return num(a.x, -a.y); }
inline num inv(num a) {
  dbl n = (a.x * a.x + a.y * a.y);
  return num(a.x / n, -a.y / n);
}

#else
// NTT
const int mod = 998244353, g = 3;
// For p < 2^30 there is also (5 << 25, 3), (7 << 26, 3),
// (479 << 21, 3) and (483 << 21, 5). Last two are > 10^9.
struct num { 
  int v;
  num(ll v_ = 0): v(int(v_ % mod)) {
    if (v < 0) v += mod;
  }
  explicit operator int() const { return v; }
};
inline num operator+(num a, num b) { return num(a.v + b.v); }
inline num operator-(num a, num b) {
  return num(a.v + mod - b.v);
}
inline num operator*(num a, num b) {
  return num(1ll * a.v * b.v);
}
inline num pow(num a, int b) {
  num r = 1;
  do {
    if (b & 1) r = r * a;
    a = a * a;
  } while (b >>= 1);
  return r;
}
inline num inv(num a) { return pow(a, mod - 2); }

#endif
using vn = vector<num>;
vi rev({0, 1});
vn rt(2, num(1)), fa, fb;
inline void init(int n) { 
  if (n <= sz(rt)) return;
  rev.resize(n);
  rep(i, 0, n) rev[i] = (rev[i >> 1] | ((i & 1) * n)) >> 1;
  rt.reserve(n);
  for (int k = sz(rt); k < n; k *= 2) {
    rt.resize(2 * k);
#if FFT
    double a = M_PI / k;
    num z(cos(a), sin(a)); // FFT
#else
    num z = pow(num(g), (mod - 1) / (2 * k)); // NTT
#endif
    rep(i, k / 2, k) rt[2 * i] = rt[i],
                            rt[2 * i + 1] = rt[i] * z;
  }
} 
inline void fft(vector<num>& a, int n) { 
  init(n);
  int s = __builtin_ctz(sz(rev) / n);
  rep(i, 0, n) if (i < rev[i] >> s) swap(a[i], a[rev[i] >> s]);
  for (int k = 1; k < n; k *= 2)
    for (int i = 0; i < n; i += 2 * k) rep(j, 0, k) {
        num t = rt[j + k] * a[i + j + k];
        a[i + j + k] = a[i + j] - t;
        a[i + j] = a[i + j] + t;
      }
} 
// Complex/NTT
vn multiply(vn a, vn b) { 
  int s = sz(a) + sz(b) - 1;
  if (s <= 0) return {};
  int L = s > 1 ? 32 - __builtin_clz(s - 1) : 0, n = 1 << L;
  a.resize(n), b.resize(n);
  fft(a, n);
  fft(b, n);
  num d = inv(num(n));
  rep(i, 0, n) a[i] = a[i] * b[i] * d;
  reverse(a.begin() + 1, a.end());
  fft(a, n);
  a.resize(s);
  return a;
} 
// Complex/NTT power-series inverse
// Doubles b as b[:n] = (2 - a[:n] * b[:n/2]) * b[:n/2]
vn inverse(const vn& a) { 
  if (a.empty()) return {};
  vn b({inv(a[0])});
  b.reserve(2 * a.size());
  while (sz(b) < sz(a)) {
    int n = 2 * sz(b);
    b.resize(2 * n, 0);
    if (sz(fa) < 2 * n) fa.resize(2 * n);
    fill(fa.begin(), fa.begin() + 2 * n, 0);
    copy(a.begin(), a.begin() + min(n, sz(a)), fa.begin());
    fft(b, 2 * n);
    fft(fa, 2 * n);
    num d = inv(num(2 * n));
    rep(i, 0, 2 * n) b[i] = b[i] * (2 - fa[i] * b[i]) * d;
    reverse(b.begin() + 1, b.end());
    fft(b, 2 * n);
    b.resize(n);
  }
  b.resize(a.size());
  return b;
} 
#if FFT
// Double multiply (num = complex)
using vd = vector<double>;
vd multiply(const vd& a, const vd& b) { 
  int s = sz(a) + sz(b) - 1;
  if (s <= 0) return {};
  int L = s > 1 ? 32 - __builtin_clz(s - 1) : 0, n = 1 << L;
  if (sz(fa) < n) fa.resize(n);
  if (sz(fb) < n) fb.resize(n);
  fill(fa.begin(), fa.begin() + n, 0);
  rep(i, 0, sz(a)) fa[i].x = a[i];
  rep(i, 0, sz(b)) fa[i].y = b[i];
  fft(fa, n);
  trav(x, fa) x = x * x;
  rep(i, 0, n) fb[i] = fa[(n - i) & (n - 1)] - conj(fa[i]);
  fft(fb, n);
  vd r(s);
  rep(i, 0, s) r[i] = fb[i].y / (4 * n);
  return r;
} 
// Integer multiply mod m (num = complex) 
vi multiply_mod(const vi& a, const vi& b, int m) {
  int s = sz(a) + sz(b) - 1;
  if (s <= 0) return {};
  int L = s > 1 ? 32 - __builtin_clz(s - 1) : 0, n = 1 << L;
  if (sz(fa) < n) fa.resize(n);
  if (sz(fb) < n) fb.resize(n);
  rep(i, 0, sz(a)) fa[i] =
    num(a[i] & ((1 << 15) - 1), a[i] >> 15);
  fill(fa.begin() + sz(a), fa.begin() + n, 0);
  rep(i, 0, sz(b)) fb[i] =
    num(b[i] & ((1 << 15) - 1), b[i] >> 15);
  fill(fb.begin() + sz(b), fb.begin() + n, 0);
  fft(fa, n);
  fft(fb, n);
  double r0 = 0.5 / n; // 1/2n
  rep(i, 0, n / 2 + 1) {
    int j = (n - i) & (n - 1);
    num g0 = (fb[i] + conj(fb[j])) * r0;
    num g1 = (fb[i] - conj(fb[j])) * r0;
    swap(g1.x, g1.y);
    g1.y *= -1;
    if (j != i) {
      swap(fa[j], fa[i]);
      fb[j] = fa[j] * g1;
      fa[j] = fa[j] * g0;
    }
    fb[i] = fa[i] * conj(g1);
    fa[i] = fa[i] * conj(g0);
  }
  fft(fa, n);
  fft(fb, n);
  vi r(s);
  rep(i, 0, s) r[i] =
    int((ll(fa[i].x + 0.5) + (ll(fa[i].y + 0.5) % m << 15) +
          (ll(fb[i].x + 0.5) % m << 15) +
          (ll(fb[i].y + 0.5) % m << 30)) %
      m);
  return r;
} 
#endif
} // namespace fft
// For multiply_mod, use num = modnum, poly = vector<num>
using fft::num;
using poly = fft::vn;
using fft::multiply;
using fft::inverse;

poly& operator+=(poly& a, const poly& b) {
  if (sz(a) < sz(b)) a.resize(b.size());
  rep(i, 0, sz(b)) a[i] = a[i] + b[i];
  return a;
}
poly operator+(const poly& a, const poly& b) {
  poly r = a;
  r += b;
  return r;
}
poly& operator-=(poly& a, const poly& b) {
  if (sz(a) < sz(b)) a.resize(b.size());
  rep(i, 0, sz(b)) a[i] = a[i] - b[i];
  return a;
}
poly operator-(const poly& a, const poly& b) {
  poly r = a;
  r -= b;
  return r;
}
poly operator*(const poly& a, const poly& b) {
  return multiply(a, b);
}
poly& operator*=(poly& a, const poly& b) { return a = a * b; }

poly& operator*=(poly& a, const num& b) { // Optional
  trav(x, a) x = x * b;
  return a;
}
poly operator*(const poly& a, const num& b) {
  poly r = a;
  r *= b;
  return r;
}
// Polynomial floor division; no leading 0's please
poly operator/(poly a, poly b) { 
  if (sz(a) < sz(b)) return {};
  int s = sz(a) - sz(b) + 1;
  reverse(a.begin(), a.end());
  reverse(b.begin(), b.end());
  a.resize(s);
  b.resize(s);
  a = a * inverse(move(b));
  a.resize(s);
  reverse(a.begin(), a.end());
  return a;
} 
poly& operator/=(poly& a, const poly& b) { return a = a / b; }
poly& operator%=(poly& a, const poly& b) { 
  if (sz(a) >= sz(b)) {
    poly c = (a / b) * b;
    a.resize(sz(b) - 1);
    rep(i, 0, sz(a)) a[i] = a[i] - c[i];
  }
  return a;
} 
poly operator%(const poly& a, const poly& b) {
  poly r = a;
  r %= b;
  return r;
}
// Log/exp/pow
poly deriv(const poly& a) { 
  if (a.empty()) return {};
  poly b(sz(a) - 1);
  rep(i, 1, sz(a)) b[i - 1] = a[i] * i;
  return b;
} 
poly integ(const poly& a) { 
  poly b(sz(a) + 1);
  b[1] = 1; // mod p
  rep(i, 2, sz(b)) b[i] =
    b[fft::mod % i] * (-fft::mod / i); // mod p
  rep(i, 1, sz(b)) b[i] = a[i - 1] * b[i]; // mod p
  //rep(i,1,sz(b)) b[i]=a[i-1]*inv(num(i)); // else
  return b;
} 
poly log(const poly& a) { // MUST have a[0] == 1 
  poly b = integ(deriv(a) * inverse(a));
  b.resize(a.size());
  return b;
} 
poly exp(const poly& a) { // MUST have a[0] == 0 
  poly b(1, num(1));
  if (a.empty()) return b;
  while (sz(b) < sz(a)) {
    int n = min(sz(b) * 2, sz(a));
    b.resize(n);
    poly v = poly(a.begin(), a.begin() + n) - log(b);
    v[0] = v[0] + num(1);
    b *= v;
    b.resize(n);
  }
  return b;
} 
poly pow(const poly& a, int m) { // m >= 0 
  poly b(a.size());
  if (!m) {
    b[0] = 1;
    return b;
  }
  int p = 0;
  while (p < sz(a) && a[p].v == 0) ++p;
  if (1ll * m * p >= sz(a)) return b;
  num mu = pow(a[p], m), di = inv(a[p]);
  poly c(sz(a) - m * p);
  rep(i, 0, sz(c)) c[i] = a[i + p] * di;
  c = log(c);
  trav(v, c) v = v * m;
  c = exp(c);
  rep(i, 0, sz(c)) b[i + m * p] = c[i] * mu;
  return b;
} 
// Multipoint evaluation/interpolation

vector<num> eval(const poly& a, const vector<num>& x) {
  int n = sz(x);
  if (!n) return {};
  vector<poly> up(2 * n);
  rep(i, 0, n) up[i + n] = poly({0 - x[i], 1});
  per(i, 1, n) up[i] = up[2 * i] * up[2 * i + 1];
  vector<poly> down(2 * n);
  down[1] = a % up[1];
  rep(i, 2, 2 * n) down[i] = down[i / 2] % up[i];
  vector<num> y(n);
  rep(i, 0, n) y[i] = down[i + n][0];
  return y;
} 

poly interp(const vector<num>& x, const vector<num>& y) {
  int n = sz(x);
  assert(n);
  vector<poly> up(n * 2);
  rep(i, 0, n) up[i + n] = poly({0 - x[i], 1});
  per(i, 1, n) up[i] = up[2 * i] * up[2 * i + 1];
  vector<num> a = eval(deriv(up[1]), x);
  vector<poly> down(2 * n);
  rep(i, 0, n) down[i + n] = poly({y[i] * inv(a[i])});
  per(i, 1, n) down[i] =
    down[i * 2] * up[i * 2 + 1] + down[i * 2 + 1] * up[i * 2];
  return down[1];
}
```

## Simplex method for linear programs
+ Maximize $c^T x$ subject to $Ax \leq b$, $x \geq 0$.
+ Returns $-\infty$ if there is no solution, $+\infty$ if there are arbitrarily good solutions, or the maximum value of $c^T x$ otherwise. The (arbitrary) input vector is set to an optimal $x$ (or in the unbounded case, an arbitrary solution fulfilling the constraints). Numerical stability is not guaranteed. For better performance, define variables such that $x = 0$ is viable.
+ Complexity: $O(NM \cdot pivots)$. $O(2^n)$ in general (very hard to achieve).
```cpp
typedef double T; // might be much slower with long doubles
typedef vector<T> vd;
typedef vector<vd> vvd;
const T eps = 1e-8, inf = 1/.0;
#define MP make_pair
#define ltj(X) if(s == -1 || MP(X[j],N[j]) < MP(X[s],N[s])) s=j
#define rep(i, a, b) for(int i = a; i < (b); ++i)

struct LPSolver {
  int m, n;
  vector<int> N,B;
  vvd D;
  LPSolver(const vvd& A, const vd& b, const vd& c) : m(sz(b)), n(sz(c)), N(n+1), B(m), D(m+2, vd(n+2)){
    rep(i,0,m) rep(j,0,n) D[i][j] = A[i][j];
    rep(i,0,m) { B[i] = n+i; D[i][n] = -1; D[i][n+1] = b[i];} rep(j,0,n) { N[j] = j; D[m][j] = -c[j]; }
    N[n] = -1; D[m+1][n] = 1;
  };
  void pivot(int r, int s){
    T *a = D[r].data(), inv = 1 / a[s];
    rep(i,0,m+2) if (i != r && abs(D[i][s]) > eps) {
      T *b = D[i].data(), inv2 = b[s] * inv;
      rep(j,0,n+2) b[j] -= a[j] * inv2;
      b[s] = a[s] * inv2;
    }
    rep(j,0,n+2) if (j != s) D[r][j] *= inv;
    rep(i,0,m+2) if (i != r) D[i][s] *= -inv;
    D[r][s] = inv;
    swap(B[r], N[s]);
  }
  bool simplex(int phase){
    int x = m + phase - 1;
    for (;;) {
      int s = -1;
      rep(j,0,n+1) if (N[j] != -phase) ltj(D[x]); if (D[x][s] >= -eps) return true;
      int r = -1;
      rep(i,0,m) {
        if (D[i][s] <= eps) continue;
        if (r == -1 || MP(D[i][n+1] / D[i][s], B[i]) < MP(D[r][n+1] / D[r][s], B[r])) r = i;
      }
      if (r == -1) return false;
      pivot(r, s);
    }
  }
  T solve(vd &x){
    int r = 0;
    rep(i,1,m) if (D[i][n+1] < D[r][n+1]) r = i;
    if (D[r][n+1] < -eps) {
      pivot(r, n);
      if (!simplex(2) || D[m+1][n+1] < -eps) return -inf;
      rep(i,0,m) if (B[i] == -1) {
        int s = 0;
        rep(j,1,n+1) ltj(D[i]);
        pivot(i, s);
      }
    }
    bool ok = simplex(1); x = vd(n);
    rep(i,0,m) if (B[i] < n) x[B[i]] = D[i][n+1];
    return ok ? D[m][n+1] : inf;
  }
};
```
