# Strings

```cpp
vector<int> prefix_function(string s){
  int n = sz(s);
  vector<int> pi(n);
  for (int i = 1; i < n; i++){
    int k = pi[i - 1];
    while (k > 0 && s[i] != s[k]){
      k = pi[k - 1];
    }
    pi[i] = k + (s[i] == s[k]);
  }
  return pi;
}
// Returns the positions of the first character
vector<int> kmp(string s, string k){
  string st = k + "#" + s;
  vector<int> res;
  auto pi = prefix_function(st);
  for (int i = 0; i < sz(st); i++){
    if (pi[i] == sz(k)){
      res.pb(i - 2 * sz(k));
    }
  }
  return res;
}
vector<int> z_function(string s){
  int n = sz(s);
  vector<int> z(n);
  int l = 0, r = 0;
  for (int i = 1; i < n; i++){
    if (r >= i) z[i] = min(z[i - l], r - i + 1);
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]){
      z[i]++;
    }
    if (i + z[i] - 1 > r){
      l = i, r = i + z[i] - 1;
    }
  }
  return z;
}
```
## Manacher’s algorithm
```cpp
/*
Finds longest palindromes centered at each index
even[i] = d --> [i - d, i + d - 1] is a max-palindrome
odd[i] = d  --> [i - d, i + d] is a max-palindrome
*/
pair<vector<int>, vector<int>> manacher(string s) {
  vector<char> t{'^', '#'};
  for (char c : s) t.push_back(c), t.push_back('#');
  t.push_back('$');
  int n = t.size(), r = 0, c = 0;
  vector<int> p(n, 0);
  for (int i = 1; i < n - 1; i++) {
    if (i < r + c) p[i] = min(p[2 * c - i], r + c - i);
    while (t[i + p[i] + 1] == t[i - p[i] - 1]) p[i]++;
    if (i + p[i] > r + c) r = p[i], c = i;
  }
  vector<int> even(sz(s)), odd(sz(s));
  for (int i = 0; i < sz(s); i++){
    even[i] = p[2 * i + 1] / 2, odd[i] = p[2 * i + 2] / 2;
  }
  return {even, odd};
}
```
## Aho-Corasick Trie
+ Given a set of strings, constructs a trie with suffix links.
+ For a particular node, $link$ points to the longest proper suffix of this node that's contained in the trie.
+ $nxt$ encodes suffix links in a compressed format:
  + If vertex $v$ has a child by letter $x$, then $trie[v].nxt[x]$ points to that child.
  + If vertex $v$ doesn't have such child, then $trie[v].nxt[x]$ points to the suffix link of that child if we would actually have it.
+ **Facts:** suffix link graph can be seen as a tree; terminal link tree has height $O(\sqrt{N})$, where $N$ is the sum of strings' lengths.
+ **Usage:** add all strings, then call $add \textunderscore links()$.
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
## Suffix Automaton
+ Given a string $S$, constructs a DAG that is an automaton of all suffixes of $S$.
+ The automaton has $\le 2n$ nodes and $\le 3n$ edges.
+ Properties (let all paths start at node 0):
    - Every path represents a unique substring of $S$.
    - A path ends at a terminal node iff it represents a suffix of $S$.
    - All paths ending at a fixed node $v$ have the same set of right endpoints of their occurences in $S$.
    - Let $endpos(v)$ represent this set. Then, $link(v) := u$ such that $endpos(v) \subset endpos(u)$ and $|endpos(u)|$ is smallest possible. $link(0):=-1$. Links form a tree.
    - Let $len(v)$ be the longest path ending at $v$. All paths ending at $v$ have distinct lengths: every length from interval $[len(link(v))+1, len(v)]$.
+ One of the main applications is dealing with **distinct** substrings. Such problems can be solved with DFS and DP.
+ Complexity: $O(|S| \cdot \log{|\Sigma|)}$. Perhaps replace map with vector if $|\Sigma|$ is small.
```cpp
const int MAXLEN = 1e5 + 20;
 
struct suffix_automaton{
  struct state {
    int len, link;
    bool terminal = 0, used = 0;
    map<char, int> next;
  };
 
  state st[MAXLEN * 2];
  int sz = 0, last;
 
  suffix_automaton(){
    st[0].len = 0;
    st[0].link = -1;
    sz++;
    last = 0;
  };
 
  void extend(char c) {
    int cur = sz++;
    st[cur].len = st[last].len + 1;
    int p = last;
    while (p != -1 && !st[p].next.count(c)) {
      st[p].next[c] = cur;
      p = st[p].link;
    }
    if (p == -1) {
      st[cur].link = 0;
    } else {
      int q = st[p].next[c];
      if (st[p].len + 1 == st[q].len) {
        st[cur].link = q;
      } else {
        int clone = sz++;
        st[clone].len = st[p].len + 1;
        st[clone].next = st[q].next;
        st[clone].link = st[q].link;
        while (p != -1 && st[p].next[c] == q) {
            st[p].next[c] = clone;
            p = st[p].link;
        }
        st[q].link = st[cur].link = clone;
      }
    }
    last = cur;
  }
 
  void mark_terminal(){
    int cur = last;
    while (cur) st[cur].terminal = 1, cur = st[cur].link;
  }
};
/*
Usage:
suffix_automaton sa;
for (int i = 0; i < sz(str); i++) sa.extend(str[i]);
sa.mark_terminal();
*/
```
