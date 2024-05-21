# Miscellaneous

## ORDERED SET

```cpp
#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp> 
using namespace __gnu_pbds; 
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;
```cpp

## MEASURING EXECUTION TIME

```cpp
ld tic = clock();
// execute algo…
ld tac = clock();
// Time in milliseconds
cerr << (tac - tic) / CLOCKS_PER_SEC * 1000 << endl;
// No need to comment out the print because it’s done to cerr.
```

## SETTING FIXED D.P. PRECISION

```cpp
cout << setprecision(d) << fixed;
// Each number is rounded to d digits after the decimal point, and truncated.
```
