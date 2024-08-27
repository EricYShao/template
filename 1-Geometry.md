# Geometry

## Point and vector basics
```cpp
const ld EPS = 1e-9;

struct point{
  ld x, y;
  point() : x(0), y(0) {}
  point(ld x_, ld y_) : x(x_), y(y_) {}
 
  point operator+ (point rhs) const{
    return point(x + rhs.x, y + rhs.y); }
  point operator- (point rhs) const{
    return point(x - rhs.x, y - rhs.y); }
  point operator* (ld rhs) const{
    return point(x * rhs, y * rhs); }
  point operator/ (ld rhs) const{
    return point(x / rhs, y / rhs); }
  point ort() const{
    return point(-y, x); }
  ld abs2() const{
    return x * x + y * y; }
  ld len() const{
    return sqrtl(abs2()); }
  point unit() const{
    return point(x, y) / len(); }
  point rotate(ld a) const{
    return point(x * cosl(a) - y * sinl(a), x * sinl(a) + y * cosl(a));
  }
  friend ostream& operator<<(ostream& os, point p){
    return os << "(" << p.x << "," << p.y << ")";
  }

  bool operator< (point rhs) const{
    return make_pair(x, y) < make_pair(rhs.x, rhs.y);
  }
  bool operator== (point rhs) const{
    return abs(x - rhs.x) < EPS && abs(y - rhs.y) < EPS;
  }
};

ld sq(ld a){
  return a * a; }
ld dot(point a, point b){
  return a.x * b.x + a.y * b.y; }
ld cross(point a, point b){
  return a.x * b.y - a.y * b.x; }
ld dist(point a, point b){
  return (a - b).len(); }
bool acw(point a, point b){
  return cross(a, b) > -EPS; }
bool cw(point a, point b){
  return cross(a, b) < EPS; }
int sgn(ld x){
  return (x > EPS) - (x < EPS); } // for integer: EPS = 0
int half(point p) { return p.y != 0 ? sgn(p.y) : sgn(p.x); } // +1: [0, pi), -1: [pi, 2*pi)
bool angle_comp(point a, point b) { int A = half(a), B = half(b);
  return A == B ? cross(a, b) > 0 : A > B; }
```
## Line basics
```cpp
struct line{
  ld a, b, c;
  line() : a(0), b(0), c(0) {}
  line(ld a_, ld b_, ld c_) : a(a_), b(b_), c(c_) {}
  line(point p1, point p2){
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = -a * p1.x - b * p1.y;
  }
};

ld det(ld a11, ld a12, ld a21, ld a22){
  return a11 * a22 - a12 * a21;
}
bool parallel(line l1, line l2){
  return abs(cross(point(l1.a, l1.b), point(l2.a, l2.b))) < EPS;
}
bool operator==(line l1, line l2){
  return parallel(l1, l2) &&
  abs(det(l1.b, l1.c, l2.b, l2.c)) < EPS &&
  abs(det(l1.a, l1.c, l2.a, l2.c)) < EPS;
}
```
# Line and segment intersections
```cpp
// {p, 0} - unique intersection, {p, 1} - infinite, {p, 2} - none
pair<point, int> line_inter(line l1, line l2){
  if (parallel(l1, l2)){
    return {point(), l1 == l2? 1 : 2};
  }
  return {point(
    det(-l1.c, l1.b, -l2.c, l2.b) / det(l1.a, l1.b, l2.a, l2.b),
    det(l1.a, -l1.c, l2.a, -l2.c) / det(l1.a, l1.b, l2.a, l2.b)
  ), 0};
}


// Checks if p lies on ab
bool is_on_seg(point p, point a, point b){
  return abs(cross(p - a, p - b)) < EPS && dot(p - a, p - b) < EPS;
}

/*
If a unique intersection point between the line segments going from a to b and from c to d exists then it is returned.
If no intersection point exists an empty vector is returned.
If infinitely many exist a vector with 2 elements is returned, containing the endpoints of the common line segment.
*/
vector<point> segment_inter(point a, point b, point c, point d) {
  auto oa = cross(d - c, a - c), ob = cross(d - c, b - c), oc = cross(b - a, c - a), od = cross(b - a, d - a);
  if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0) return {(a * ob - b * oa) / (ob - oa)};
  set<point> s;
  if (is_on_seg(a, c, d)) s.insert(a);
  if (is_on_seg(b, c, d)) s.insert(b);
  if (is_on_seg(c, a, b)) s.insert(c);
  if (is_on_seg(d, a, b)) s.insert(d);
  return {all(s)};
}
```
## Distances from a point to line and segment
```cpp
// Distance from p to line ab
ld line_dist(point p, point a, point b){
  return cross(b - a, p - a) / (b - a).len();
}

// Distance from p to segment ab
ld segment_dist(point p, point a, point b){
  if (a == b) return (p - a).len();
  auto d = (a - b).abs2(), t = min(d, max((ld)0, dot(p - a, b - a)));
  return ((p - a) * d - (b - a) * t).len() / d;
}
```
## Polygon area and Centroid
```cpp
pair<point,ld> cenArea(const vector<point>& v) { assert(sz(v) >= 3);
  point cen(0, 0); ld area = 0; 
  forn(i,sz(v)) {
    int j = (i+1)%sz(v); ld a = cross(v[i],v[j]);
    cen = cen + a*(v[i]+v[j]); area += a; }
  return {cen/area/(ld)3,area/2}; // area is SIGNED
}
```
## Convex hull
+ Complexity: $O(n \log n)$.
```cpp
vector<point> convex_hull(vector<point> pts){
  sort(all(pts));
  pts.erase(unique(all(pts)), pts.end());
  vector<point> up, down;
  for (auto p : pts){
    while (sz(up) > 1 && acw(up.end()[-1] - up.end()[-2], p - up.end()[-2])) up.pop_back();
    while (sz(down) > 1 && cw(down.end()[-1] - down.end()[-2], p - down.end()[-2])) down.pop_back();
    up.pb(p), down.pb(p);
  }
  for (int i = sz(up) - 2; i >= 1; i--) down.pb(up[i]);
  return down;
}
```
## Point location in a convex polygon
+ Complexity: $O(n)$ precalculation and $O(\log n)$ query.
```cpp
void prep_convex_poly(vector<point>& pts){
  rotate(pts.begin(), min_element(all(pts)), pts.end());
}

// 0 - Outside, 1 - Exclusively Inside, 2 - On the Border
int in_convex_poly(point p, vector<point>& pts){
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
```
## Point location in a simple polygon
+ Complexity: $O(n)$.
```cpp
// 0 - Outside, 1 - Exclusively Inside, 2 - On the Border
int in_simple_poly(point p, vector<point>& pts){
  int n = sz(pts);
  bool res = 0;
  for (int i = 0; i < n; i++){
    auto a = pts[i], b = pts[(i + 1) % n];
    if (is_on_seg(p, a, b)) return 2;
    if (((a.y > p.y) - (b.y > p.y)) * cross(b - p, a - p) > EPS){
      res ^= 1;
    }
  }
  return res;
}
```
## Minkowski Sum
+ For two convex polygons $P$ and $Q$, returns the set of points $(p + q)$, where $p \in P, q \in Q$.
+ This set is also a convex polygon.
+ Complexity: $O(n)$.
```cpp
void minkowski_rotate(vector<point>& P){
  int pos = 0;
  for (int i = 1; i < sz(P); i++){
    if (abs(P[i].y - P[pos].y) <= EPS){
      if (P[i].x < P[pos].x) pos = i;
    }
    else if (P[i].y < P[pos].y) pos = i;
  }
  rotate(P.begin(), P.begin() + pos, P.end());
}
// P and Q are strictly convex, points given in counterclockwise order.
vector<point> minkowski_sum(vector<point> P, vector<point> Q){
  minkowski_rotate(P);
  minkowski_rotate(Q);
  P.pb(P[0]);
  Q.pb(Q[0]);
  vector<point> ans;
  int i = 0, j = 0;
  while (i < sz(P) - 1 || j < sz(Q) - 1){
    ans.pb(P[i] + Q[j]);
    ld curmul;
    if (i == sz(P) - 1) curmul = -1;
    else if (j == sz(Q) - 1) curmul = +1;
    else curmul = cross(P[i + 1] - P[i], Q[j + 1] - Q[j]);
    if (abs(curmul) < EPS || curmul > 0) i++;
    if (abs(curmul) < EPS || curmul < 0) j++;
  }
  return ans;
}
```
## Half-plane intersection
+ Given $N$ half-plane conditions in the form of a ray, computes the vertices of their intersection polygon.
+ Complexity: $O(N \log{N})$.
+ A ray is defined by a point $p$ and direction vector $dp$. The half-plane is to the **left** of the direction vector.
```cpp
// Extra functions needed: point operations, dot, cross
const ld EPS = 1e-9;

int sgn(ld a){
  return (a > EPS) - (a < -EPS);
}
int half(point p){
  return p.y != 0? sgn(p.y) : -sgn(p.x);
}
bool angle_comp(point a, point b){
  int A = half(a), B = half(b);
  return A == B? cross(a, b) > 0 : A < B;
}
struct ray{
  point p, dp; // origin, direction
  ray(point p_, point dp_){
    p = p_, dp = dp_;
  }
  point isect(ray l){
    return p + dp * (cross(l.dp, l.p - p) / cross(l.dp, dp));
  }
  bool operator<(ray l){
    return angle_comp(dp, l.dp);
  }
};
vector<point> half_plane_isect(vector<ray> rays, ld DX = 1e9, ld DY = 1e9){
  // constrain the area to [0, DX] x [0, DY]
  rays.pb({point(0, 0), point(1, 0)});
  rays.pb({point(DX, 0), point(0, 1)});
  rays.pb({point(DX, DY), point(-1, 0)});
  rays.pb({point(0, DY), point(0, -1)});
  sort(all(rays));
  {
    vector<ray> nrays;
    for (auto t : rays){
      if (nrays.empty() || cross(nrays.back().dp, t.dp) > EPS){
        nrays.pb(t);
        continue;
      }
      if (cross(t.dp, t.p - nrays.back().p) > 0) nrays.back() = t;
    }
    swap(rays, nrays);
  }
  auto bad = [&] (ray a, ray b, ray c){
    point p1 = a.isect(b), p2 = b.isect(c);
    if (dot(p2 - p1, b.dp) <= EPS){
      if (cross(a.dp, c.dp) <= 0) return 2;
      return 1;
    }
    return 0;
  };
  #define reduce(t) \
    while (sz(poly) > 1){ \
      int b = bad(poly[sz(poly) - 2], poly.back(), t); \
      if (b == 2) return {}; \
      if (b == 1) poly.pop_back(); \
      else break; \
    }
  deque<ray> poly;
  for (auto t : rays){
    reduce(t);
    poly.pb(t);
  }
  for (;; poly.pop_front()){
    reduce(poly[0]);
    if (!bad(poly.back(), poly[0], poly[1])) break;
  }
  assert(sz(poly) >= 3); // expect nonzero area
  vector<point> poly_points;
  for (int i = 0; i < sz(poly); i++){
    poly_points.pb(poly[i].isect(poly[(i + 1) % sz(poly)]));
  }
  return poly_points;
}
```
## Circles
+ Finds minimum enclosing circle of vector of points in expected $O(N)$
```cpp
// necessary point functions
ld sq(ld a) { return a*a; }
point operator+(const point& l, const point& r) { 
  return point(l.x+r.x,l.y+r.y); }
point operator*(const point& l, const ld& r) { 
  return point(l.x*r,l.y*r); }
point operator*(const ld& l, const point& r) { return r*l; }
ld abs2(const point& p) { return sq(p.x)+sq(p.y); }
ld abs(const point& p) { return sqrt(abs2(p)); }
point conj(const point& p) { return point(p.x,-p.y); }
point operator-(const point& l, const point& r) { 
  return point(l.x-r.x,l.y-r.y); }
point operator*(const point& l, const point& r) { 
   return point(l.x*r.x-l.y*r.y,l.y*r.x+l.x*r.y); }
point operator/(const point& l, const ld& r) { 
   return point(l.x/r,l.y/r); }
point operator/(const point& l, const point& r) { 
   return l*conj(r)/abs2(r); }

// circle code
using circ = pair<point,ld>;

circ ccCenter(point a, point b, point c) { 
  b = b-a; c = c-a;
  point res = b*c*(conj(c)-conj(b))/(b*conj(c)-conj(b)*c);
  return {a+res,abs(res)};
}

circ mec(vector<point> ps) {
  // expected O(N)
  shuffle(all(ps), rng);
  point o = ps[0]; ld r = 0, EPS = 1+1e-8;
  forn(i,sz(ps)) if (abs(o-ps[i]) > r*EPS) {
    o = ps[i], r = 0; // point is on MEC
    forn(j,i) if (abs(o-ps[j]) > r*EPS) {
      o = (ps[i]+ps[j])/2, r = abs(o-ps[i]);
      forn(k,j) if (abs(o-ps[k]) > r*EPS) 
        tie(o,r) = ccCenter(ps[i],ps[j],ps[k]);
    }
  }
  return {o,r};
}
```
+ Circle tangents, external and internal
```cpp
point unit(const point& p) { return p * (1/abs(p)); }

point tangent(point p, circ c, int t = 0) {
  c.se = abs(c.se); // abs needed because internal calls y.s < 0
  if (c.se == 0) return c.fi;
  ld d = abs(p-c.fi);
  point a = pow(c.se/d,2)*(p-c.fi)+c.fi;
  point b = sqrt(d*d-c.se*c.se)/d*c.se*unit(p-c.fi)*point(0,1); 
  return t == 0 ? a+b : a-b;
}
vector<pair<point,point>> external(circ a, circ b) { 
  vector<pair<point,point>> v; 
  if (a.se == b.se) {
    point tmp = unit(a.fi-b.fi)*a.se*point(0, 1);
    v.emplace_back(a.fi+tmp,b.fi+tmp);
    v.emplace_back(a.fi-tmp,b.fi-tmp);
  } else {
    point p = (b.se*a.fi-a.se*b.fi)/(b.se-a.se);
    forn(i,2) v.emplace_back(tangent(p,a,i),tangent(p,b,i));
  }
  return v;
}
vector<pair<point,point>> internal(circ a, circ b) { 
  return external({a.fi,-a.se},b); }
```
