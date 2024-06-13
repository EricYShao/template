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

## Extended Euclidean Algorithm

```cpp
// gives (x, y) for ax + by = g
// solutions given (x0, y0): a(x0 + kb/g) + b(y0 - ka/g) = g
int gcd(int a, int b, int& x, int& y) {
  x = 1, y = 0; int sum1 = a;
  int x2 = 0, y2 = 1, sum2 = b;
  while (sum2) {
    int q = sum1 / sum2;
    tie(x, x2) = make_tuple(x2, x - q * x2);
    tie(y, y2) = make_tuple(y2, y - q * y2);
    tie(sum1, sum2) = make_tuple(sum2, sum1 - q * sum2);
  }
  return sum1;
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
