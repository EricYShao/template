# Flows
## O(N^2 * M), on unit networks O(N^(1/2) * M)

```cpp
struct FlowEdge {
    int v, u;
    long long cap, flow = 0;
    FlowEdge(int v, int u, long long cap) : v(v), u(u), cap(cap) {}
};
struct Dinic {
    const long long flow_inf = 1e18;
    vector<FlowEdge> edges;
    vector<vector<int>> adj;
    int n, m = 0;
    int s, t;
    vector<int> level, ptr;
    queue<int> q;
    Dinic(int n, int s, int t) : n(n), s(s), t(t) {
        adj.resize(n);
        level.resize(n);
        ptr.resize(n);
    }
    void add_edge(int v, int u, long long cap) {
        edges.emplace_back(v, u, cap);
        edges.emplace_back(u, v, 0);
        adj[v].push_back(m);
        adj[u].push_back(m + 1);
        m += 2;
    }
    bool bfs() {
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (int id : adj[v]) {
                if (edges[id].cap - edges[id].flow < 1)
                    continue;
                if (level[edges[id].u] != -1)
                    continue;
                level[edges[id].u] = level[v] + 1;
                q.push(edges[id].u);
            }
        }
        return level[t] != -1;
    }
    long long dfs(int v, long long pushed) {
        if (pushed == 0)
            return 0;
        if (v == t)
            return pushed;
        for (int& cid = ptr[v]; cid < (int)adj[v].size(); cid++) {
            int id = adj[v][cid];
            int u = edges[id].u;
            if (level[v] + 1 != level[u] || edges[id].cap - edges[id].flow < 1)
                continue;
            long long tr = dfs(u, min(pushed, edges[id].cap - edges[id].flow));
            if (tr == 0)
                continue;
            edges[id].flow += tr;
            edges[id ^ 1].flow -= tr;
            return tr;
        }
        return 0;
    }
    long long flow() {
        long long f = 0;
        while (true) {
            fill(level.begin(), level.end(), -1);
            level[s] = 0;
            q.push(s);
            if (!bfs())
                break;
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, flow_inf)) {
                f += pushed;
            }
        }
        return f;
    }
};
// To recover flow through original edges: iterate over even indices in edges.
```
## MCMF – maximize flow, then minimize its cost. O(F*m*n).

```cpp
#include <ext/pb_ds/priority_queue.hpp>
template <typename T, typename C>
class MCMF {
 public:
   static constexpr T eps = (T) 1e-9;

   struct edge {
     int from;
     int to;
     T c;
     T f;
     C cost;
   };

   int n;
   vector<vector<int>> g;
   vector<edge> edges;
   vector<C> d;
   vector<C> pot;
   __gnu_pbds::priority_queue<pair<C, int>> q;
   vector<typename decltype(q)::point_iterator> its;
   vector<int> pe;
   const C INF_C = numeric_limits<C>::max() / 2;

   explicit MCMF(int n_) : n(n_), g(n), d(n), pot(n, 0), its(n), pe(n) {}

   int add(int from, int to, T forward_cap, C edge_cost, T backward_cap = 0) {
     assert(0 <= from && from < n && 0 <= to && to < n);
     assert(forward_cap >= 0 && backward_cap >= 0);
     int id = static_cast<int>(edges.size());
     g[from].push_back(id);
     edges.push_back({from, to, forward_cap, 0, edge_cost});
     g[to].push_back(id + 1);
     edges.push_back({to, from, backward_cap, 0, -edge_cost});
     return id;
   }

   void expath(int st) {
     fill(d.begin(), d.end(), INF_C);
     q.clear();
     fill(its.begin(), its.end(), q.end());
     its[st] = q.push({pot[st], st});
     d[st] = 0;
     while (!q.empty()) {
       int i = q.top().second;
       q.pop();
       its[i] = q.end();
       for (int id : g[i]) {
         const edge &e = edges[id];
         int j = e.to;
         if (e.c - e.f > eps && d[i] + e.cost < d[j]) {
           d[j] = d[i] + e.cost;
           pe[j] = id;
           if (its[j] == q.end()) {
             its[j] = q.push({pot[j] - d[j], j});
           } else {
             q.modify(its[j], {pot[j] - d[j], j});
           }
         }
       }
     }
     swap(d, pot);
   }

   pair<T, C> max_flow(int st, int fin) {
     T flow = 0;
     C cost = 0;
     bool ok = true;
     for (auto& e : edges) {
       if (e.c - e.f > eps && e.cost + pot[e.from] - pot[e.to] < 0) {
         ok = false;
         break;
       }
     }
     if (ok) {
       expath(st);
     } else {
       vector<int> deg(n, 0);
       for (int i = 0; i < n; i++) {
         for (int eid : g[i]) {
           auto& e = edges[eid];
           if (e.c - e.f > eps) {
             deg[e.to] += 1;
           }
         }
       }
       vector<int> que;
       for (int i = 0; i < n; i++) {
         if (deg[i] == 0) {
           que.push_back(i);
         }
       }
       for (int b = 0; b < (int) que.size(); b++) {
         for (int eid : g[que[b]]) {
           auto& e = edges[eid];
           if (e.c - e.f > eps) {
             deg[e.to] -= 1;
             if (deg[e.to] == 0) {
               que.push_back(e.to);
             }
           }
         }
       }
       fill(pot.begin(), pot.end(), INF_C);
       pot[st] = 0;
       if (static_cast<int>(que.size()) == n) {
         for (int v : que) {
           if (pot[v] < INF_C) {
             for (int eid : g[v]) {
               auto& e = edges[eid];
               if (e.c - e.f > eps) {
                 if (pot[v] + e.cost < pot[e.to]) {
                   pot[e.to] = pot[v] + e.cost;
                   pe[e.to] = eid;
                 }
               }
             }
           }
         }
       } else {
         que.assign(1, st);
         vector<bool> in_queue(n, false);
         in_queue[st] = true;
         for (int b = 0; b < (int) que.size(); b++) {
           int i = que[b];
           in_queue[i] = false;
           for (int id : g[i]) {
             const edge &e = edges[id];
             if (e.c - e.f > eps && pot[i] + e.cost < pot[e.to]) {
               pot[e.to] = pot[i] + e.cost;
               pe[e.to] = id;
               if (!in_queue[e.to]) {
                 que.push_back(e.to);
                 in_queue[e.to] = true;
               }
             }
           }
         }
       }
     }
     while (pot[fin] < INF_C) {
       T push = numeric_limits<T>::max();
       int v = fin;
       while (v != st) {
         const edge &e = edges[pe[v]];
         push = min(push, e.c - e.f);
         v = e.from;
       }
       v = fin;
       while (v != st) {
         edge &e = edges[pe[v]];
         e.f += push;
         edge &back = edges[pe[v] ^ 1];
         back.f -= push;
         v = e.from;
       }
       flow += push;
       cost += push * pot[fin];
       expath(st);
     }
     return {flow, cost};
   }
};

// Examples: MCMF<int, int> g(n); g.add(u,v,c,w,0); g.max_flow(s,t).
// To recover flow through original edges: iterate over even indices in edges.
```
