# Geometry

```cpp
template<typename T>
struct TPoint{
  T x, y;
  int id;
  static constexpr T eps = static_cast<T>(1e-9);
  TPoint() : x(0), y(0), id(-1) {}
  TPoint(const T& x_, const T& y_) : x(x_), y(y_), id(-1) {}
  TPoint(const T& x_, const T& y_, const int id_) : x(x_), y(y_), id(id_) {}
 
  TPoint operator + (const TPoint& rhs) const {
    return TPoint(x + rhs.x, y + rhs.y);
  }
  TPoint operator - (const TPoint& rhs) const {
    return TPoint(x - rhs.x, y - rhs.y);
  }
  TPoint operator * (const T& rhs) const {
    return TPoint(x * rhs, y * rhs);
  }
  TPoint operator / (const T& rhs) const {
    return TPoint(x / rhs, y / rhs);
  }
  TPoint ort() const {
    return TPoint(-y, x);
  }
  T abs2() const {
    return x * x + y * y;
  }
};
template<typename T>
bool operator< (TPoint<T>& A, TPoint<T>& B){
  return make_pair(A.x, A.y) < make_pair(B.x, B.y);
}
template<typename T>
bool operator== (TPoint<T>& A, TPoint<T>& B){
  return abs(A.x - B.x) <= TPoint<T>::eps && abs(A.y - B.y) <= TPoint<T>::eps;
}
template<typename T>
struct TLine{
  T a, b, c;
  TLine() : a(0), b(0), c(0) {}
  TLine(const T& a_, const T& b_, const T& c_) : a(a_), b(b_), c(c_) {}
  TLine(const TPoint<T>& p1, const TPoint<T>& p2){
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = -a * p1.x - b * p1.y;
  }
};
template<typename T>
T det(const T& a11, const T& a12, const T& a21, const T& a22){
  return a11 * a22 - a12 * a21;
}
template<typename T>
T sq(const T& a){
  return a * a;
}
template<typename T>
T smul(const TPoint<T>& a, const TPoint<T>& b){
  return a.x * b.x + a.y * b.y;
}
template<typename T>
T vmul(const TPoint<T>& a, const TPoint<T>& b){
  return det(a.x, a.y, b.x, b.y);
}
template<typename T>
bool parallel(const TLine<T>& l1, const TLine<T>& l2){
  return abs(vmul(TPoint<T>(l1.a, l1.b), TPoint<T>(l2.a, l2.b))) <= TPoint<T>::eps;
}
template<typename T>
bool equivalent(const TLine<T>& l1, const TLine<T>& l2){
  return parallel(l1, l2) &&
  abs(det(l1.b, l1.c, l2.b, l2.c)) <= TPoint<T>::eps &&
  abs(det(l1.a, l1.c, l2.a, l2.c)) <= TPoint<T>::eps;
}
template<typename T>
TPoint<T> intersection(const TLine<T>& l1, const TLine<T>& l2){
  return TPoint<T>(
    det(-l1.c, l1.b, -l2.c, l2.b) / det(l1.a, l1.b, l2.a, l2.b),
    det(l1.a, -l1.c, l2.a, -l2.c) / det(l1.a, l1.b, l2.a, l2.b)
  );
}
template<typename T>
int sign(const T& x){
  if (abs(x) <= TPoint<T>::eps) return 0;
  return x > 0? +1 : -1;
}
template<typename T>
T area(const vector<TPoint<T>>& pts){
  int n = sz(pts);
  T ans = 0;
  for (int i = 0; i < n; i++){
    ans += vmul(pts[i], pts[(i + 1) % n]);
  }
  return abs(ans) / 2;
}
template<typename T>
T dist_pp(const TPoint<T>& a, const TPoint<T>& b){
  return sqrt(sq(a.x - b.x) + sq(a.y - b.y));
}
template<typename T>
TLine<T> perp_line(const TLine<T>& l, const TPoint<T>& p){
  T na = -l.b, nb = l.a, nc = - na * p.x - nb * p.y;
  return TLine<T>(na, nb, nc);
}
template<typename T>
TPoint<T> projection(const TPoint<T>& p, const TLine<T>& l){
  return intersection(l, perp_line(l, p));
}
template<typename T>
T dist_pl(const TPoint<T>& p, const TLine<T>& l){
  return dist_pp(p, projection(p, l));
}
template<typename T>
struct TRay{
  TLine<T> l;
  TPoint<T> start, dirvec;
  TRay() : l(), start(), dirvec() {}
  TRay(const TPoint<T>& p1, const TPoint<T>& p2){
    l = TLine<T>(p1, p2);
    start = p1, dirvec = p2 - p1;
  }
};
template<typename T>
bool is_on_line(const TPoint<T>& p, const TLine<T>& l){
  return abs(l.a * p.x + l.b * p.y + l.c) <= TPoint<T>::eps;
}
template<typename T>
bool is_on_ray(const TPoint<T>& p, const TRay<T>& r){
  if (is_on_line(p, r.l)){
    return sign(smul(r.dirvec, TPoint<T>(p - r.start))) != -1;
  }
  else return false;
}
template<typename T>
bool is_on_seg(const TPoint<T>& P, const TPoint<T>& A, const TPoint<T>& B){
  return is_on_ray(P, TRay<T>(A, B)) && is_on_ray(P, TRay<T>(B, A));
}
template<typename T>
T dist_pr(const TPoint<T>& P, const TRay<T>& R){
  auto H = projection(P, R.l);
  return is_on_ray(H, R)? dist_pp(P, H) : dist_pp(P, R.start);
}
template<typename T>
T dist_ps(const TPoint<T>& P, const TPoint<T>& A, const TPoint<T>& B){
  auto H = projection(P, TLine<T>(A, B));
  if (is_on_seg(H, A, B)) return dist_pp(P, H);
  else return min(dist_pp(P, A), dist_pp(P, B));
}
template<typename T>
bool acw(const TPoint<T>& A, const TPoint<T>& B){
  T mul = vmul(A, B);
  return mul > 0 || abs(mul) <= TPoint<T>::eps;
}
template<typename T>
bool cw(const TPoint<T>& A, const TPoint<T>& B){
  T mul = vmul(A, B);
  return mul < 0 || abs(mul) <= TPoint<T>::eps;
}
template<typename T>
vector<TPoint<T>> convex_hull(vector<TPoint<T>> pts){
  sort(all(pts));
  pts.erase(unique(all(pts)), pts.end());
  vector<TPoint<T>> up, down;
  for (auto p : pts){
    while (sz(up) > 1 && acw(up.end()[-1] - up.end()[-2], p - up.end()[-2])) up.pop_back();
    while (sz(down) > 1 && cw(down.end()[-1] - down.end()[-2], p - down.end()[-2])) down.pop_back();
    up.pb(p), down.pb(p);
  }
  for (int i = sz(up) - 2; i >= 1; i--) down.pb(up[i]);
  return down;
}

template<typename T>
bool in_triangle(TPoint<T>& P, TPoint<T>& A, TPoint<T>& B, TPoint<T>& C){
  if (is_on_seg(P, A, B) || is_on_seg(P, B, C) || is_on_seg(P, C, A)) return true;
  return cw(P - A, B - A) == cw(P - B, C - B) && 
  cw(P - A, B - A) == cw(P - C, A - C);
}
template<typename T>
void prep_convex_poly(vector<TPoint<T>>& pts){
  rotate(pts.begin(), min_element(all(pts)), pts.end());
}
// 0 - Outside, 1 - Exclusively Inside, 2 - On the Border
template<typename T>
int in_convex_poly(TPoint<T>& p, vector<TPoint<T>>& pts){
  int n = sz(pts);
  if (!n) return 0;
  if (n <= 2) return is_on_seg(p, pts[0], pts.back());
  int l = 1, r = n - 1;
  while (r - l > 1){
    int mid = (l + r) / 2;
    if (acw(pts[mid] - pts[0], p - pts[0])) l = mid;
    else r = mid;
  }
  if (!in_triangle(p, pts[0], pts[l], pts[l + 1])) return 0;
  if (is_on_seg(p, pts[l], pts[l + 1]) ||
    is_on_seg(p, pts[0], pts.back()) ||
    is_on_seg(p, pts[0], pts[1])
  ) return 2;
  return 1;
}
// 0 - Outside, 1 - Exclusively Inside, 2 - On the Border
template<typename T>
int in_simple_poly(TPoint<T> p, vector<TPoint<T>>& pts){
  int n = sz(pts);
  bool res = 0;
  for (int i = 0; i < n; i++){
    auto a = pts[i], b = pts[(i + 1) % n];
    if (is_on_seg(p, a, b)) return 2;
    if (((a.y > p.y) - (b.y > p.y)) * vmul(b - p, a - p) > TPoint<T>::eps){
      res ^= 1;
    }
  }
  return res;
}
template<typename T>
void minkowski_rotate(vector<TPoint<T>>& P){
  int pos = 0;
  for (int i = 1; i < sz(P); i++){
    if (abs(P[i].y - P[pos].y) <= TPoint<T>::eps){
      if (P[i].x < P[pos].x) pos = i;
    }
    else if (P[i].y < P[pos].y) pos = i;
  }
  rotate(P.begin(), P.begin() + pos, P.end());
}
// P and Q are strictly convex, points given in counterclockwise order
template<typename T>
vector<TPoint<T>> minkowski_sum(vector<TPoint<T>> P, vector<TPoint<T>> Q){
  minkowski_rotate(P);
  minkowski_rotate(Q);
  P.pb(P[0]);
  Q.pb(Q[0]);
  vector<TPoint<T>> ans;
  int i = 0, j = 0;
  while (i < sz(P) - 1 || j < sz(Q) - 1){
    ans.pb(P[i] + Q[j]);
    T curmul;
    if (i == sz(P) - 1) curmul = -1;
    else if (j == sz(Q) - 1) curmul = +1;
    else curmul = vmul(P[i + 1] - P[i], Q[j + 1] - Q[j]);
    if (abs(curmul) < TPoint<T>::eps || curmul > 0) i++;
    if (abs(curmul) < TPoint<T>::eps || curmul < 0) j++;
  }
  return ans;
}
using Point = TPoint<ll>; using Line = TLine<ll>; using Ray = TRay<ll>; const ld PI = acos(-1);
```
