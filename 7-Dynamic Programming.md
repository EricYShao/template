# Dynamic Programming
## Sum over Subset DP
+ Computes $f[A] = \sum_{B \subseteq A}{a[B]}$.
+ Complexity: $O(2^n \cdot n)$.
```cpp
for (int i = 0; i < (1 << n); i++) f[i] = a[i];
for (int i = 0; i < n; i++) for (int mask = 0; mask < (1 << n); mask++) if ((mask >> i) & 1){
  f[mask] += f[mask ^ (1 << i)];
}
```
## Divide and Conquer DP
+ Helps to compute 2D DP of the form:
+ $dp[i][j] = \min\limits_{0 \leq k \leq j-1}{\left(dp[i - 1][k] + cost(k + 1, j)\right)}$
+ **Necessary condition:** let $opt(i, j)$ be the optimal $k$ for the state $(i, j)$. Then, $opt(i, j) \leq opt(i, j + 1)$.
+ **Sufficient condition:** $cost(a, d) + cost(b, c) \ge cost(a, c) + cost(b, d)$ where $a < b < c < d$.
+ Complexity: $O(M \cdot N \cdot \log N)$ for computing $dp[M][N]$.
```cpp
vector<ll> dp_old(N), dp_new(N);

void rec(int l, int r, int optl, int optr){
  if (l > r) return;
  int mid = (l + r) / 2;
  pair<ll, int> best = {INF, optl};
  for (int i = optl; i <= min(mid - 1, optr); i++){ // If k can be j, change to "i <= min(mid, optr)".
    ll cur = dp_old[i] + cost(i + 1, mid);
    if (cur < best.fi) best = {cur, i};
  }
  dp_new[mid] = best.fi;

  rec(l, mid - 1, optl, best.se);
  rec(mid + 1, r, best.se, optr);
}

// Computes the DP "by layers"
fill(all(dp_old), INF);
dp_old[0] = 0;
while (layers--){
   rec(0, n, 0, n);
   dp_old = dp_new;
 }
```
## Knuth's DP Optimization
+ Computes DP of the form
+ $dp[i][j] = \min\limits_{i \leq k \leq j-1}{\left(dp[i][k] + dp[k+1][j] + cost(i, j)\right)}$
+ **Necessary Condition:** $opt(i, j-1) \leq opt(i, j) \leq opt(i+1, j)$
+ **Sufficient Condition:** For $a \leq b \leq c \leq d$, $cost(b, c) \le cost(a, d)$ AND $cost(a, d) + cost(b, c) \ge cost(a, c) + cost(b, d)$
+ Complexity: $O(n^2)$
```cpp
int N;
int dp[N][N], opt[N][N];
auto C = [&](int i, int j) {
  // Implement cost function C.
};
for (int i = 0; i < N; i++) {
  opt[i][i] = i;
  // Initialize dp[i][i] according to the problem
}
for (int i = N-2; i >= 0; i--) {
  for (int j = i+1; j < N; j++) {
    int mn = INT_MAX;
    int cost = C(i, j);
    for (int k = opt[i][j-1]; k <= min(j-1, opt[i+1][j]); k++) {
      if (mn >= dp[i][k] + dp[k+1][j] + cost) {
        opt[i][j] = k; 
        mn = dp[i][k] + dp[k+1][j] + cost; 
      }
    }
    dp[i][j] = mn; 
  }
}
```
