# Miscellaneous

## Ordered Set

```cpp
#include <ext/pb_ds/assoc_container.hpp> 
#include <ext/pb_ds/tree_policy.hpp> 
using namespace __gnu_pbds; 
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;
```

## Measuring Execution Time

```cpp
ld tic = clock();
// execute algo…
ld tac = clock();
// Time in milliseconds
cerr << (tac - tic) / CLOCKS_PER_SEC * 1000 << endl;
// No need to comment out the print because it’s done to cerr.
```

## Setting Fixed D.P. Precision

```cpp
cout << setprecision(d) << fixed;
// Each number is rounded to d digits after the decimal point, and truncated.
```
## Debugging Flags
+ Converts segfaults into Wrong Answers. Similarly one can catch SIGABRT (assertion failures) and SIGFPE (zero divisions)
```cpp
signal(SIGSEGV, [](int) { Exit(0); });
```

## Common Bugs and General Advice

+ Check overflow, array bounds
+ Check variable overloading
+ Check special cases (n=1?)
+ Do something instead of nothing, stay organized
+ Write stuff down!
+ Don't get stuck on one approach!

## MIT troubleshoot.txt
+ General:
+ Write down most of your thoughts, even if you're not sure whether they're useful.
+ Stay organized and don't leave papers all over the place! \vspace{5pt}

+ Pre-submit:
+ Write a few simple test cases if sample is not enough.
+ Are time limits close? If so, generate max cases.
+ Is the memory usage fine?
+ Could anything overflow?
+ Remove debug output.
+ Make sure to submit the right file. \vspace{5pt}

+ Wrong answer:
+ Print your solution! Print debug output as well.
+ Read the full problem statement again.
+ Have you understood the problem correctly?
+ Are you sure your algorithm works? 
+ Try writing a slow (but correct) solution.
+ Can your algorithm handle the whole range of input?
+ Did you consider corner cases (ex. n=1)?
+ Is your output format correct? (including whitespace)
+ Are you clearing all data structures between test cases?
+ Any uninitialized variables?
+ Any undefined behavior (array out of bounds)?
+ Any overflows or NaNs (or shifting ll by >=64 bits)?
+ Confusing N and M, i and j, etc.?
+ Confusing ++i and i++?
+ Return vs continue vs break?
+ Are you sure the STL functions you use work as you think?
+ Add some assertions, maybe resubmit.
+ Create some test cases to run your algorithm on.
+ Go through the algorithm for a simple case.
+ Go through this list again.
+ Explain your algorithm to a teammate.
+ Ask a teammate to look at your code.
+ Go for a small walk, e.g. to the toilet.
+ Rewrite your solution from the start or let a teammate do it. \vspace{5pt}

+ Geometry:
+ Work with ints if possible.
+ Correctly account for numbers close to (but not) zero.
+ Related: for functions like acos make sure absolute val of input is not (slightly) greater than one.
+ Correctly deal with vertices that are collinear, concyclic, coplanar (in 3D), etc.
+ Subtracting a point from every other (but not itself)? \vspace{5pt}

+ Runtime error:
+ Have you tested all corner cases locally?
+ Any uninitialized variables?
+ Are you reading or writing outside the range of any vector?
+ Any assertions that might fail?
+ Any possible division by 0? (mod 0 for example)
+ Any possible infinite recursion?
+ Invalidated pointers or iterators?
+ Are you using too much memory?
+ Debug with resubmits (e.g. remapped signals, see Various). \vspace{5pt}

+ Time limit exceeded:
+ Do you have any possible infinite loops?
+ What's your complexity? Large TL does not mean that something simple (like NlogN) isn't intended.
+ Are you copying a lot of unnecessary data? (References)
+ Avoid vector, map. (use arrays/unordered_map)
+ How big is the input and output? (consider FastIO)
+ What do your teammates think about your algorithm?
+ Calling count() on multiset? \vspace{5pt}

+ Memory limit exceeded:
+ What is the max amount of memory your algorithm should need?
+ Are you clearing all data structures between test cases?
