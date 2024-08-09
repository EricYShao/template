# Data Structures
## Fenwick Tree
```cpp
ll sum(int r) {
  ll ret = 0;
  for (; r >= 0; r = (r & r + 1) - 1) ret += bit[r];
  return ret;
}
void add(int idx, ll delta) {
  for (; idx < n; idx |= idx + 1) bit[idx] += delta;
}
```
## Lazy Propagation SegTree
```cpp
// Clear: clear() or build()
const int N = 2e5 + 10; // Change the constant!
template<typename T>
struct LazySegTree{
  T t[4 * N];
  T lazy[4 * N];
  int n;

  // Change these functions, default return, and lazy mark.
  T default_return = 0, lazy_mark = numeric_limits<T>::min();
  // Lazy mark is how the algorithm will identify that no propagation is needed.
  function<T(T, T)> f = [&] (T a, T b){
    return a + b;
  };
  // f_on_seg calculates the function f, knowing the lazy value on segment,
  // segment's size and the previous value.
  // The default is segment modification for RSQ. For increments change to:
  //     return cur_seg_val + seg_size * lazy_val;
  // For RMQ.   Modification: return lazy_val;   Increments: return cur_seg_val + lazy_val;
  function<T(T, int, T)> f_on_seg = [&] (T cur_seg_val, int seg_size, T lazy_val){
    return seg_size * lazy_val;
  };
  // upd_lazy updates the value to be propagated to child segments.
  // Default: modification. For increments change to:
  //     lazy[v] = (lazy[v] == lazy_mark? val : lazy[v] + val);
  function<void(int, T)> upd_lazy = [&] (int v, T val){
    lazy[v] = val;
  };
  // Tip: for "get element on single index" queries, use max() on segment: no overflows.

  LazySegTree(int n_) : n(n_) {
    clear(n);
  }

  void build(int v, int tl, int tr, vector<T>& a){
    if (tl == tr) {
      t[v] = a[tl];
      return;
    }
    int tm = (tl + tr) / 2;
    // left child: [tl, tm]
    // right child: [tm + 1, tr]
    build(2 * v + 1, tl, tm, a);
    build(2 * v + 2, tm + 1, tr, a);
    t[v] = f(t[2 * v + 1], t[2 * v + 2]);
  }

  LazySegTree(vector<T>& a){
    build(a);
  }

  void push(int v, int tl, int tr){
    if (lazy[v] == lazy_mark) return;
    int tm = (tl + tr) / 2;
    t[2 * v + 1] = f_on_seg(t[2 * v + 1], tm - tl + 1, lazy[v]);
    t[2 * v + 2] = f_on_seg(t[2 * v + 2], tr - tm, lazy[v]);
    upd_lazy(2 * v + 1, lazy[v]), upd_lazy(2 * v + 2, lazy[v]);
    lazy[v] = lazy_mark;
  }
  
  void modify(int v, int tl, int tr, int l, int r, T val){
    if (l > r) return;
    if (tl == l && tr == r){
      t[v] = f_on_seg(t[v], tr - tl + 1, val);
      upd_lazy(v, val);
      return;
    }
    push(v, tl, tr);
    int tm = (tl + tr) / 2;
    modify(2 * v + 1, tl, tm, l, min(r, tm), val);
    modify(2 * v + 2, tm + 1, tr, max(l, tm + 1), r, val);
    t[v] = f(t[2 * v + 1], t[2 * v + 2]);
  }
  
  T query(int v, int tl, int tr, int l, int r) {
    if (l > r) return default_return;
    if (tl == l && tr == r) return t[v];
    push(v, tl, tr);
    int tm = (tl + tr) / 2;
    return f(
      query(2 * v + 1, tl, tm, l, min(r, tm)),
      query(2 * v + 2, tm + 1, tr, max(l, tm + 1), r)
    );
  }
  
  void modify(int l, int r, T val){
    modify(0, 0, n - 1, l, r, val);
  }

  T query(int l, int r){
    return query(0, 0, n - 1, l, r);
  }
  
  T get(int pos){
    return query(pos, pos);
  }
  
  // Change clear() function to t.clear() if using unordered_map for SegTree!!!
  void clear(int n_){
    n = n_;
    for (int i = 0; i < 4 * n; i++) t[i] = 0, lazy[i] = lazy_mark;
  }

  void build(vector<T>& a){
    n = sz(a);
    clear(n);
    build(0, 0, n - 1, a);
  }
};
```

## Sparse Table

```cpp
const int N = 2e5 + 10, LOG = 20; // Change the constant!
template<typename T>
struct SparseTable{
int lg[N];
T st[N][LOG];
int n;

// Change this function
function<T(T, T)> f = [&] (T a, T b){
  return min(a, b);
};

void build(vector<T>& a){
  n = sz(a);
  lg[1] = 0;
  for (int i = 2; i <= n; i++) lg[i] = lg[i / 2] + 1;

  for (int k = 0; k < LOG; k++){
    for (int i = 0; i < n; i++){
      if (!k) st[i][k] = a[i];
      else st[i][k] = f(st[i][k - 1], st[min(n - 1, i + (1 << (k - 1)))][k - 1]);
    }
  }
}

T query(int l, int r){
  int sz = r - l + 1;
  return f(st[l][lg[sz]], st[r - (1 << lg[sz]) + 1][lg[sz]]);
}
};
```

## Suffix Array and LCP array
+ (uses SparseTable above)

```cpp
struct SuffixArray{
  vector<int> p, c, h;
  SparseTable<int> st;
  /*
  In the end, array c gives the position of each suffix in p
  using 1-based indexation!
  */

  SuffixArray() {}

  SuffixArray(string s){
    buildArray(s);
    buildLCP(s);
    buildSparse();
  }

  void buildArray(string s){
    int n = sz(s) + 1;
    p.resize(n), c.resize(n);
    for (int i = 0; i < n; i++) p[i] = i;
    sort(all(p), [&] (int a, int b){return s[a] < s[b];});
    c[p[0]] = 0;
    for (int i = 1; i < n; i++){
      c[p[i]] = c[p[i - 1]] + (s[p[i]] != s[p[i - 1]]);
    }
    vector<int> p2(n), c2(n);
    // w is half-length of each string.
    for (int w = 1; w < n; w <<= 1){
      for (int i = 0; i < n; i++){
        p2[i] = (p[i] - w + n) % n;
      }
      vector<int> cnt(n);
      for (auto i : c) cnt[i]++;
      for (int i = 1; i < n; i++) cnt[i] += cnt[i - 1];
      for (int i = n - 1; i >= 0; i--){
        p[--cnt[c[p2[i]]]] = p2[i];
      }
      c2[p[0]] = 0;
      for (int i = 1; i < n; i++){
        c2[p[i]] = c2[p[i - 1]] + 
        (c[p[i]] != c[p[i - 1]] ||
        c[(p[i] + w) % n] != c[(p[i - 1] + w) % n]);
      }
      c.swap(c2);
    }
    p.erase(p.begin());
  }

  void buildLCP(string s){
    // The algorithm assumes that suffix array is already built on the same string.
    int n = sz(s);
    h.resize(n - 1);
    int k = 0;
    for (int i = 0; i < n; i++){
      if (c[i] == n){
        k = 0;
        continue;
      }
      int j = p[c[i]];
      while (i + k < n && j + k < n && s[i + k] == s[j + k]) k++;
      h[c[i] - 1] = k;
      if (k) k--;
    }
    /*
    Then an RMQ Sparse Table can be built on array h
    to calculate LCP of 2 non-consecutive suffixes.
    */
  }

  void buildSparse(){
    st.build(h);
  }

  // l and r must be in 0-BASED INDEXATION
  int lcp(int l, int r){
    l = c[l] - 1, r = c[r] - 1;
    if (l > r) swap(l, r);
    return st.query(l, r - 1);
  }
};
```

## Aho Corasick Trie
+ For each node in the trie, the suffix link points to the longest proper suffix of the represented string. The terminal-link tree has square-root height (can be constructed by DFS).

```cpp
const int S = 26;

// Function converting char to int.
int ctoi(char c){
  return c - 'a';
}

// To add terminal links, use DFS
struct Node{
  vector<int> nxt;
  int link;
  bool terminal;

  Node() {
    nxt.assign(S, -1), link = 0, terminal = 0;
  }
};

vector<Node> trie(1);

// add_string returns the terminal vertex.
int add_string(string& s){
  int v = 0;
  for (auto c : s){
    int cur = ctoi(c);
    if (trie[v].nxt[cur] == -1){
      trie[v].nxt[cur] = sz(trie);
      trie.emplace_back();
    }
    v = trie[v].nxt[cur];
  }
  trie[v].terminal = 1;
  return v;
}

/*
Suffix links are compressed.
This means that:
  If vertex v has a child by letter x, then:
    trie[v].nxt[x] points to that child.
  If vertex v doesn't have such child, then:
    trie[v].nxt[x] points to the suffix link of that child
    if we would actually have it.
*/
void add_links(){
  queue<int> q;
  q.push(0);
  while (!q.empty()){
    auto v = q.front();
    int u = trie[v].link;
    q.pop();
    for (int i = 0; i < S; i++){
      int& ch = trie[v].nxt[i];
      if (ch == -1){
        ch = v? trie[u].nxt[i] : 0;
      }
      else{
        trie[ch].link = v? trie[u].nxt[i] : 0;
        q.push(ch);
      }
    }
  }
}

bool is_terminal(int v){
  return trie[v].terminal;
}

int get_link(int v){
  return trie[v].link;
}

int go(int v, char c){
  return trie[v].nxt[ctoi(c)];
}
```

## Convex Hull Trick
+ Allows to insert a linear function to the hull in Ө(1) and get the minimum/maximum value of the stored function at a point in O(log n).
+ NOTE: The lines must be added in the order of decreasing/increasing gradients. CAREFULLY CHECK THE SETUP BEFORE USING!
+ IMPORTANT: THE DEFAULT VERSION SURELY WORKS. IF MODIFIED VERSIONS DON'T WORK, TRY TRANSFORMING THEM TO THE DEFAULT ONE BY CHANGING SIGNS.

```cpp
struct line{
  ll k, b;
  ll f(ll x){
    return k * x + b;
  };
};
 
vector<line> hull;
 
void add_line(line nl){
  if (!hull.empty() && hull.back().k == nl.k){
    nl.b = min(nl.b, hull.back().b); // Default: minimum. For maximum change "min" to "max".
    hull.pop_back();
  }
  while (sz(hull) > 1){
    auto& l1 = hull.end()[-2], l2 = hull.back();
    if ((nl.b - l1.b) * (l2.k - nl.k) >= (nl.b - l2.b) * (l1.k - nl.k)) hull.pop_back(); // Default: decreasing gradient k. For increasing k change the sign to <=.
    else break;
  }
  hull.pb(nl);
}
 
ll get(ll x){
  int l = 0, r = sz(hull);
  while (r - l > 1){
    int mid = (l + r) / 2;
    if (hull[mid - 1].f(x) >= hull[mid].f(x)) l = mid; // Default: minimum. For maximum change the sign to <=.
    else r = mid;
  }
  return hull[l].f(x);
}
```

## Li-Chao Segment Tree

+ allows to add linear functions in any order and query minimum/maximum value of those at a point, all in O(log n).
+ Clear: clear()

```cpp
const ll INF = 1e18; // Change the constant!
struct LiChaoTree{
  struct line{
    ll k, b;
    line(){
      k = b = 0;
    };
    line(ll k_, ll b_){
      k = k_, b = b_;
    };
    ll f(ll x){
      return k * x + b;
    };
  };
  int n;
  bool minimum, on_points;
  vector<ll> pts;
  vector<line> t;

  void clear(){
    for (auto& l : t) l.k = 0, l.b = minimum? INF : -INF;
  }
  
  LiChaoTree(int n_, bool min_){ // This is a default constructor for numbers in range [0, n - 1].
    n = n_, minimum = min_, on_points = false;
    t.resize(4 * n);
    clear();
  };

  LiChaoTree(vector<ll> pts_, bool min_){ // This constructor will build LCT on the set of points you pass. The points may be in any order and contain duplicates.
    pts = pts_, minimum = min_;
    sort(all(pts));
    pts.erase(unique(all(pts)), pts.end());
    on_points = true;
    n = sz(pts);
    t.resize(4 * n);
    clear();
  };

  void add_line(int v, int l, int r, line nl){
    // Adding on segment [l, r)
    int m = (l + r) / 2;
    ll lval = on_points? pts[l] : l, mval = on_points? pts[m] : m;
    if ((minimum && nl.f(mval) < t[v].f(mval)) || (!minimum && nl.f(mval) > t[v].f(mval))) swap(t[v], nl);
    if (r - l == 1) return;
    if ((minimum && nl.f(lval) < t[v].f(lval)) || (!minimum && nl.f(lval) > t[v].f(lval))) add_line(2 * v + 1, l, m, nl);
    else add_line(2 * v + 2, m, r, nl);
  }

  ll get(int v, int l, int r, int x){
    int m = (l + r) / 2;
    if (r - l == 1) return t[v].f(on_points? pts[x] : x);
    else{
      if (minimum) return min(t[v].f(on_points? pts[x] : x), x < m? get(2 * v + 1, l, m, x) : get(2 * v + 2, m, r, x));
      else return max(t[v].f(on_points? pts[x] : x), x < m? get(2 * v + 1, l, m, x) : get(2 * v + 2, m, r, x));
    }
  }

  void add_line(ll k, ll b){
    add_line(0, 0, n, line(k, b));
  }

  ll get(ll x){
    return get(0, 0, n, on_points? lower_bound(all(pts), x) - pts.begin() : x);
  }; // Always pass the actual value of x, even if LCT is on points.
};
```

## Persistent Segment Tree

+ for RSQ

```cpp
struct Node {
  ll val;
  Node *l, *r;

  Node(ll x) : val(x), l(nullptr), r(nullptr) {}
  Node(Node *ll, Node *rr) {
    l = ll, r = rr;
    val = 0;
    if (l) val += l->val;
    if (r) val += r->val;
  }
  Node(Node *cp) : val(cp->val), l(cp->l), r(cp->r) {}
};
const int N = 2e5 + 20;
ll a[N];
Node *roots[N];
int n, cnt = 1;
Node *build(int l = 1, int r = n) {
  if (l == r) return new Node(a[l]);
  int mid = (l + r) / 2;
  return new Node(build(l, mid), build(mid + 1, r));
}
Node *update(Node *node, int val, int pos, int l = 1, int r = n) {
  if (l == r) return new Node(val);
  int mid = (l + r) / 2;
  if (pos > mid)
    return new Node(node->l, update(node->r, val, pos, mid + 1, r));
  else return new Node(update(node->l, val, pos, l, mid), node->r);
}
ll query(Node *node, int a, int b, int l = 1, int r = n) {
  if (l > b || r < a) return 0;
  if (l >= a && r <= b) return node->val;
  int mid = (l + r) / 2;
  return query(node->l, a, b, l, mid) + query(node->r, a, b, mid + 1, r);
}
```
