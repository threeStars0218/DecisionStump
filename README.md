# DecisionStump
A C++ code for decision stump

# HOW TO USE
Here is an example to use decision stump:
In this example, we consider binary classification data.
we use 4 training examples and each example is 3-dimensional real vector.

```cpp
#include <vector>
#include "decisionstump.hpp"

int main() {
    std::vector< std::vector< double > > dat = { {  1,  2,  3 }
                                               , {  4,  5,  6 }
                                               , {  7,  8,  9 }
                                               , { 10, 11, 12 }
                                               };
    std::vector< int > lab = { 1, -1, -1, 1 };
    bool one_side = false;
    decisionstump dstump = decisionstump( dat, lab, one_side );

    return 0;
}
```
The constructer needs 3 argument:
- A data matrix.
- A label vector.
- A boolean that determine one-side hypothesis or not.

