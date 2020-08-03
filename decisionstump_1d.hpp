#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <vector>
#include <functional>

#ifndef DECISIONSTUMP_1D_HPP
#define DECISIONSTUMP_1D_HPP

class decisionstump_1d {
private:
    std::vector<size_t> sorted_index;
    std::vector<double> sorted_data;
    std::vector<double> sorted_label;
    std::vector<std::vector<double> > ev_mat;
    std::vector<std::function<int(double)> > classifiers;
    std::vector<double> threshoulds;
    size_t m;
    size_t number_of_stump;
    bool   one_side;
    double eps;
public:
    decisionstump_1d();
    decisionstump_1d(std::vector<double>, std::vector<int>, bool);
    std::pair<double, std::function<int(double)> >
            stump_one_side(std::vector<double>);
    std::pair<double, std::function<int(double)> >
            stump(std::vector<double>);
    //std::pair<double, std::function<int(double)> >
    //        stump_one_side_gumbel(std::vector<double>, double);
    //std::pair<double, std::pair<double, bool>>
    //        stump_gumbel(std::vector<double>, double);
    double naive_calc_edge(std::vector<double>, std::function<int(double)>);
    std::vector<std::vector<double> >
            all_edge_vector(const std::vector< std::vector< double > >&,
                            const std::vector< int >&,
                            const size_t& );
    int h( double, double, bool );
    std::vector< double > get_threshoulds();
    size_t get_number_of_stump();
};

inline
std::vector<double>
decisionstump_1d::
get_threshoulds(void)
{
    return this->threshoulds;
}

inline
size_t
decisionstump_1d::
get_number_of_stump(void)
{
    return this->number_of_stump;
}


#endif

