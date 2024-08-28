# Templates
## Ken’s template

```cpp
#include <bits/stdc++.h>
using namespace std;
#define all(v) (v).begin(), (v).end()
typedef long long ll;
typedef long double ld;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;
#define pb push_back
#define sz(x) (int)(x).size()
#define fi first
#define se second
#define forn(i, n) for (int i = 0; i < int(n); i++)
#define endl '\n'
```

## Kevin’s template

```cpp
// paste Ken's Template, minus last line
const char nl = '\n';
ll k, n, m, u, v, w, x, y, z;
string s;

bool multiTest = 1;
void solve(int tt){
}
 
int main(){
  ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
  cout<<fixed<< setprecision(14);
 
  int t = 1;
  if (multiTest) cin >> t;
  forn(ii, t) solve(ii);
}
```

## Kevin’s Template Extended

+ to type after the start of the contest

```cpp
typedef pair<double, double> pdd;
const ld PI = acosl(-1);
const ll mod7 = 1e9 + 7;
const ll mod9 = 998244353;
const ll INF = 2*1024*1024*1023;
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T> using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
vi d4x = {1, 0, -1, 0};
vi d4y = {0, 1, 0, -1};
vi d8x = {1, 0, -1, 0, 1, 1, -1, -1};
vi d8y = {0, 1, 0, -1, 1, -1, 1, -1};
mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
```
