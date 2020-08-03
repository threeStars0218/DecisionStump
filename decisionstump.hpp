#include "decisionstump_1d.hpp"
#include <functional>
#include <iostream>
#include <vector>

#ifndef DECISIONSTUMP_HPP
#define DECISIONSTUMP_HPP

class decisionstump {
private:
    std::vector<std::vector<double> > data;
    std::vector<int> label;
    std::vector<decisionstump_1d> dstumps;
    bool one_side;
    size_t m; // サンプル数
    size_t d; // 各サンプルの次元数
    size_t number_of_stump;
public:
    decisionstump();
    decisionstump(std::vector<std::vector<double> >, std::vector<int>);
    decisionstump(std::vector<std::vector<double> >, std::vector<int>, bool);
    std::function<int(std::vector<double>)> stump(std::vector<double>);
    // std::function<int(double)> stump_gumbel(std::vector<double>, double);
    std::vector<double> edge_vector(const std::vector< double >&);
    // std::vector<double> edge_vector_gumbel(const std::vector<double>&, double);
    std::vector<std::vector<double> >
        all_edge_vector_at(size_t);
    size_t get_number_of_stump();
    std::vector<double> get_threshoulds_at(size_t);
};

inline
std::vector<std::vector<double> >
decisionstump::
all_edge_vector_at(size_t idx)
{
    return this->dstumps[idx].all_edge_vector(this->data, this->label, idx);
}

inline
size_t
decisionstump::
get_number_of_stump()
{
    return this->number_of_stump;
}

inline
std::vector<double>
decisionstump::
get_threshoulds_at(size_t idx)
{
    return this->dstumps[idx].get_threshoulds();
}
#endif
