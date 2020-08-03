#include "decisionstump.hpp"
#include <iostream>

int main() {
    std::vector<std::vector<double> >
        dat = {{1}, {-5}, {3}, {2}, {2}, {-9}, {4}, {-2}};
    std::vector<int> lab = {-1, -1,  1,  1,  1, -1,  1, 1};
    std::vector<double> dist(8, 1.0/8.0);

    decisionstump dstump = decisionstump(dat, lab, false);
    size_t num = dstump.get_number_of_stump();
    std::cout << "dat: ";
    std::cout << "\nlab: ";
    for(int i=0; i<8; ++i)
        std::cout << '(' << dat[i][0] << ", " << lab[i] << ')'
                  << (i == 7 ? " \n" : ", ");
    std::cout << "# of stumps: " << num << ",\n";
    std::function<int(std::vector<double>)> clf = dstump.stump(dist);
    double max_edge = 0.0;
    std::cout << "data, edge vector:\n";
    for(int i=0; i<8; ++i) {
        max_edge += dist[i] * lab[i] * clf(dat[i]);
        std::cout << '(' << dat[i][0] << ", " << lab[i] << "), "
                  << lab[i] * clf(dat[i]) << std::endl;
    }
    std::cout << "max edge: " << max_edge << std::endl;
    return 0;
}

