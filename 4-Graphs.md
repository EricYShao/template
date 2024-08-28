# Graphs
## Kuhn’s algorithm for bipartite matching

```cpp
/*
The graph is split into 2 halves of n1 and n2 vertices.
Complexity: O(n1 * m). Usually runs much faster. MUCH FASTER!!!
*/
const int N = 305;

vi g[N]; // Stores edges from left half to right.
bool used[N]; // Stores if vertex from left half is used.
int mt[N]; // For every vertex in right half, stores to which vertex in left half it's matched (-1 if not matched).

bool try_dfs(int v){
  if (used[v]) return false;
  used[v] = 1;
  for (auto u : g[v]){
    if (mt[u] == -1 || try_dfs(mt[u])){
      mt[u] = v;
      return true;
    }
  }
  return false;
}

int main(){
// ......
  for (int i = 1; i <= n2; i++) mt[i] = -1;
  for (int i = 1; i <= n1; i++) used[i] = 0;
  for (int i = 1; i <= n1; i++){
    if (try_dfs(i)){
      for (int j = 1; j <= n1; j++) used[j] = 0;
    }
  }
  vector<pair<int, int>> ans;
  for (int i = 1; i <= n2; i++){
    if (mt[i] != -1) ans.pb({mt[i], i});
  }
}

// Finding maximal independent set: size = # of nodes - # of edges in matching.
// To construct: launch Kuhn-like DFS from unmatched nodes in the left half.
// Independent set = visited nodes in left half + unvisited in right half.
// Finding minimal vertex cover: complement of maximal independent set.
```


## Hungarian algorithm for Assignment Problem
+ Given a 1-indexed $(n \times m)$ matrix $A$, select a number in each row such that each column has at most 1 number selected, and the sum of the selected numbers is minimized.

```cpp
int INF = 1e9; // constant greater than any number in the matrix
vi u(n+1), v(m+1), p(m+1), way(m+1);
for (int i=1; i<=n; ++i) {
  p[0] = i;
  int j0 = 0;
  vi minv (m+1, INF);
  vector<bool> used (m+1, false);
  do {
    used[j0] = true;
    int i0 = p[j0],  delta = INF,  j1;
    for (int j=1; j<=m; ++j)
      if (!used[j]) {
        int cur = A[i0][j]-u[i0]-v[j];
        if (cur < minv[j])
          minv[j] = cur,  way[j] = j0;
        if (minv[j] < delta)
          delta = minv[j],  j1 = j;
      }
    for (int j=0; j<=m; ++j)
      if (used[j])
        u[p[j]] += delta,  v[j] -= delta;
      else
        minv[j] -= delta;
    j0 = j1;
  } while (p[j0] != 0);
  do {
    int j1 = way[j0];
    p[j0] = p[j1];
    j0 = j1;
  } while (j0);
}
vi ans (n+1); // ans[i] stores the column selected for row i
for (int j=1; j<=m; ++j)
  ans[p[j]] = j;
int cost = -v[0]; // the total cost of the matching
```

## Dijkstra’s Algorithm

```cpp
priority_queue<pair<ll, ll>, vector<pair<ll, ll>>, greater<pair<ll, ll>>> q;
dist[start] = 0;
q.push({0, start});
while (!q.empty()){
  auto [d, v] = q.top();
  q.pop();
  if (d != dist[v]) continue;
  for (auto [u, w] : g[v]){
    if (dist[u] > dist[v] + w){
      dist[u] = dist[v] + w;
      q.push({dist[u], u});
    }
  }
}
```

## Eulerian Cycle DFS

```cpp
void dfs(int v){
  while (!g[v].empty()){
    int u = g[v].back();
    g[v].pop_back();
    dfs(u);
    ans.pb(v);
  }
}
```

## SCC and 2-SAT

```cpp
void scc(vector<vi>& g, int* idx) {
  int n = g.size(), ct = 0;
  int out[n];
  vi ginv[n];
  memset(out, -1, sizeof out);
  memset(idx, -1, n * sizeof(int));
  function<void(int)> dfs = [&](int cur) {
    out[cur] = INT_MAX;
    for(int v : g[cur]) {
      ginv[v].push_back(cur);
      if(out[v] == -1) dfs(v);
    }
    ct++; out[cur] = ct;
  };
  vi order;
  for(int i = 0; i < n; i++) {
    order.push_back(i);
    if(out[i] == -1) dfs(i);
  }
  sort(order.begin(), order.end(), [&](int& u, int& v) {
    return out[u] > out[v];
  });
  ct = 0;
  stack<int> s;
  auto dfs2 = [&](int start) {
    s.push(start);
    while(!s.empty()) {
      int cur = s.top();
      s.pop();
      idx[cur] = ct;
      for(int v : ginv[cur])
        if(idx[v] == -1) s.push(v);
    }
  };
  for(int v : order) {
    if(idx[v] == -1) {
      dfs2(v);
      ct++;
    }
  }
}

// 0 => impossible, 1 => possible
pair<int,vi> sat2(int n, vector<pii>& clauses) {
  vi ans(n);
  vector<vi> g(2*n + 1);
  for(auto [x, y] : clauses) {
    x = x < 0 ? -x + n : x;
    y = y < 0 ? -y + n : y;
    int nx = x <= n ? x + n : x - n;
    int ny = y <= n ? y + n : y - n;
    g[nx].push_back(y);
    g[ny].push_back(x);
  }
  int idx[2*n + 1];
  scc(g, idx);
  for(int i = 1; i <= n; i++) {
    if(idx[i] == idx[i + n]) return {0, {}};
    ans[i - 1] = idx[i + n] < idx[i];
  }
  return {1, ans};
}
```

## Finding Bridges

```cpp
/*
Bridges.
Results are stored in a map "is_bridge".
For each connected component, call "dfs(starting vertex, starting vertex)".
*/
const int N = 2e5 + 10; // Careful with the constant!
 
vi g[N];
int tin[N], fup[N], timer;
map<pair<int, int>, bool> is_bridge;
 
void dfs(int v, int p){
  tin[v] = ++timer;
  fup[v] = tin[v];
  for (auto u : g[v]){
    if (!tin[u]){
      dfs(u, v);
      if (fup[u] > tin[v]){
        is_bridge[{u, v}] = is_bridge[{v, u}] = true;
      }
      fup[v] = min(fup[v], fup[u]);
    }
    else{
      if (u != p) fup[v] = min(fup[v], tin[u]);
    }
  }
}
```


## Virtual Tree

```cpp
// order stores the nodes in the queried set
sort(all(order), [&] (int u, int v){return tin[u] < tin[v];});
int m = sz(order);
for (int i = 1; i < m; i++){
  order.pb(lca(order[i], order[i - 1]));
}
sort(all(order), [&] (int u, int v){return tin[u] < tin[v];});
order.erase(unique(all(order)), order.end());
vi stk{order[0]};
for (int i = 1; i < sz(order); i++){
  int v = order[i];
  while (tout[stk.back()] < tout[v]) stk.pop_back();
  int u = stk.back();
  vg[u].pb({v, dep[v] - dep[u]});
  stk.pb(v);
}
```

## HLD on Edges DFS

```cpp
void dfs1(int v, int p, int d){
  par[v] = p;
  for (auto e : g[v]){
    if (e.fi == p){
      g[v].erase(find(all(g[v]), e));
      break;
    }
  }
  dep[v] = d;
  sz[v] = 1;
  for (auto [u, c] : g[v]){
    dfs1(u, v, d + 1);
    sz[v] += sz[u];
  }
  if (!g[v].empty()) iter_swap(g[v].begin(), max_element(all(g[v]), comp));
}
void dfs2(int v, int rt, int c){
  pos[v] = sz(a);
  a.pb(c);
  root[v] = rt;
  for (int i = 0; i < sz(g[v]); i++){
    auto [u, c] = g[v][i];
    if (!i) dfs2(u, rt, c);
    else dfs2(u, u, c);
  }
}
int getans(int u, int v){
  int res = 0;
  for (; root[u] != root[v]; v = par[root[v]]){
    if (dep[root[u]] > dep[root[v]]) swap(u, v);
    res = max(res, rmq(0, 0, n - 1, pos[root[v]], pos[v]));
  }
  if (pos[u] > pos[v]) swap(u, v);
  return max(res, rmq(0, 0, n - 1, pos[u] + 1, pos[v]));
}
```

## Centroid Decomposition

```cpp
vector<char> res(n), seen(n), sz(n);
function<int(int, int)> get_size = [&](int node, int fa) {
  sz[node] = 1;
  for (auto& ne : g[node]) {
    if (ne == fa || seen[ne]) continue;
    sz[node] += get_size(ne, node);
  }
  return sz[node];
};
function<int(int, int, int)> find_centroid = [&](int node, int fa, int t) {
  for (auto& ne : g[node])
    if (ne != fa && !seen[ne] && sz[ne] > t / 2) return find_centroid(ne, node, t);
  return node;
};
function<void(int, char)> solve = [&](int node, char cur) {
  get_size(node, -1); auto c = find_centroid(node, -1, sz[node]);
  seen[c] = 1, res[c] = cur;
  for (auto& ne : g[c]) {
    if (seen[ne]) continue;
    solve(ne, char(cur + 1)); // we can pass c here to build tree
  }
};
```

## Biconnected Components and Block-Cut Tree
+ Biconnected components are the ones that have no articulation points.
+ They are defined by edge sets that are "bounded" by articulation points in the original graph.
+ The corresponding vertex sets are stored in $comps$.
+ Block-Cut tree is constructed by creating a fictive node for each component, and attaching edges to its members.
+ Articulation points in the original graph are the non-leaf non-fictive nodes in the BC tree.
+ Complexity: $O(n)$.
```cpp
// Usage: pass in adjacency list in 0-based indexation.
// Return: adjacency list of block-cut tree (nodes 0...n-1 represent original nodes, the rest are component nodes).
vector<vi> biconnected_components(vector<vi> g) {
	int n = sz(g);
	vector<vi> comps;
	vi stk, num(n), low(n);
  int timer = 0;
	// Finds the biconnected components
	function<void(int, int)> dfs = [&](int v, int p) {
		num[v] = low[v] = ++timer;
		stk.pb(v);
		for (int son : g[v]) {
			if (son == p) continue;
			if (num[son]) low[v] = min(low[v], num[son]);
      else{
				dfs(son, v);
				low[v] = min(low[v], low[son]);
				if (low[son] >= num[v]){
					comps.pb({v});
					while (comps.back().back() != son){
						comps.back().pb(stk.back());
						stk.pop_back();
					}
				}
			}
		}
	};
	dfs(0, -1);
	// Build the block-cut tree
	auto build_tree = [&]() {
		vector<vi> t(n);
		for (auto &comp : comps){
			t.push_back({});
			for (int u : comp){
				t.back().pb(u);
        t[u].pb(sz(t) - 1);
      }
		}
		return t;
	};
	return build_tree();
}
```
