#include "decisionstump.hpp"

int pred(std::function<int(double)> f,
         int idx,
         std::vector<double> data)
{
    return f(data[idx]);
}

decisionstump::
decisionstump() {};

decisionstump::
decisionstump(std::vector<std::vector<double> > data_,
              std::vector<int> label_)
    : data(data_)
    , label(label_)
    , one_side(false)
    , m(data_.size())
    , d(data_[0].size())
    , number_of_stump(0)
{
    std::cout << "pre-processing decisionstump.." << std::flush;
    for(int i=0; i<this->d; ++i) {
        std::vector<double> tmp(this->m);
        for(int j=0; j<this->m; ++j)
            tmp[j] = this->data[j][i];
        dstumps.push_back(decisionstump_1d(tmp, label, one_side));
    }
    for(auto &each_stump : this->dstumps)
        this->number_of_stump += each_stump.get_number_of_stump();

    std::cout << "done!" << std::endl;
}

decisionstump::
decisionstump(std::vector<std::vector<double> > data_,
              std::vector<int> label_,
              bool one_side_)
    : data(data_)
    , label(label_)
    , one_side(one_side_)
    , m(data_.size())
    , d(data_[0].size())
    , number_of_stump(0)
{
    std::cout << "pre-processing decisionstump.." << std::flush;
    for(int i=0; i<this->d; ++i) {
        std::vector<double> tmp(this->m);
        for(int j=0; j<this->m; ++j)
            tmp[j] = this->data[j][i];
        dstumps.push_back(decisionstump_1d(tmp, label, one_side));
    }
    for(auto &each_stump : this->dstumps)
        this->number_of_stump += each_stump.get_number_of_stump();

    std::cout << "done!" << std::endl;
}

std::function<int(std::vector<double>)>
decisionstump::
stump(std::vector<double> dist)
{
    auto   p = this->dstumps[0].stump(dist);
    auto   max_edge_f = p.second;
    double max_edge = p.first;
    int    max_edge_idx = 0;

    for (int i=1; i<this->d; ++i) {
        auto tmp = dstumps[i].stump(dist);
        if (tmp.first > max_edge) {
            max_edge   = tmp.first;
            max_edge_f = tmp.second;
            max_edge_idx = i;
        }
    }
    return std::bind(pred, max_edge_f, max_edge_idx, std::placeholders::_1);
}

std::vector<double>
decisionstump::
edge_vector(const std::vector<double> &dist)
{
    std::vector<double> ev(this->m);
    std::function<int(std::vector<double>)> clf = this->stump(dist);
    for (int k=0; k<this->m; ++k)
        ev[k] = this->label[k] * clf(this->data[k]);
    return ev;
}

