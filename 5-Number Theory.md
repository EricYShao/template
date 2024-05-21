# Number Theory
## EXTENDED EUCLIDEAN ALGORITHM

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

```cpp
vector<int> prime;
bool is_composite[MAX_N];
int phi[MAX_N];

void sieve (int n){
  fill(is_composite, is_composite + n, false);
  phi[1] = 1;
  for (int i = 2; i < n; ++i) {
    if (!is_composite[i]) {
      prime.push_back (i);
      phi[i] = i - 1; //i is prime
      }
  for (int j = 0; j < prime.size () && i * prime[j] < n; ++j) {
    is_composite[i * prime[j]] = true;
    if (i % prime[j] == 0) {
      phi[i * prime[j]] = phi[i] * prime[j]; //prime[j] divides i
      break;
      } else {
      phi[i * prime[j]] = phi[i] * phi[prime[j]]; //prime[j] does not divide i
      }
    }
  }
}
```
