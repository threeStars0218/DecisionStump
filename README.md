# DecisionStump
A Simple C++ code for decision stump using STL libralies.

# Time Complexity
Assume that we are given a set of $m$-labeled examples:
$\{(x_1, y_1), \ldots (x_m, y_m)\}
\subseteq \mathbb R^d \times \{-1, +1\}$,
Then, this program
- $O(md)$ pre-processing time
- $O(md)$ time to compute max edge stump.

# HOW TO USE
Here is an example to use decision stump:
In this example, we consider binary classification data.
we use 4 training examples and each example is 3-dimensional real vector.

```cpp
#include "decisionstump.hpp"
#include <vector>

using namespace std;

int main() {
    vector<vector<double> > dat = { {  1,  2,  3 },
                                    {  4,  5,  6 },
                                    {  7,  8,  9 },
                                    { 10, 11, 12 }};
    vector<int> lab = { 1, -1, -1, 1 };
    vector<double> distribution = {0.2, 0.3, 0.1, 0.4};
    decisionstump dstump = decisionstump(dat, lab);
    std::function<int(double)> func = dstump.stump(distribution);

    return 0;
}
```

The constructer needs 3 argument:
- A data matrix
(each row vector is a data, thus # of row is equals to $m$).
- A label vector(length must be equals to $m$).

When you want a max edge hypothesis at distribution `d`,
just call `stump` method as follows:
```cpp
auto classifier = dstump.stump(d);
```
after you receive a function `classifier = dstump.stump(d)`,
you can use this function for prediction
(This is just a decision stump,
so the accuracy is not good generally).
```cpp
vector<double> new_data = {1, 4, 2};
int prediction = classifier(new_data);
```

This function returns a function of type `std::function<int(double)>`.

Note that the distribution vector must be non-negative, and its 1-norm is equals to 1.