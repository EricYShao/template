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
