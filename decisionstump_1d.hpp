#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <vector>

#ifndef DECISIONSTUMP_1D_HPP
#define DECISIONSTUMP_1D_HPP

class decisionstump_1d {
private:
    std::vector< size_t > sorted_index;
    std::vector< double > sorted_data;
    std::vector< double > sorted_label;
    std::vector< double > threshoulds;
    size_t  m;
    bool    one_side;
    double eps;
public:
    decisionstump_1d();
    decisionstump_1d( std::vector< double >, std::vector< int >, bool );
    std::pair< double, std::pair< double, bool > > stump_one_side( std::vector< double > );
    std::pair< double, std::pair< double, bool > > stump( std::vector< double > );
    double naive_calc_edge( double, std::vector< double >, bool );
    std::vector< std::vector< double > > all_edge_vector( const std::vector< std::vector< double > >&
                                                        , const std::vector< int >&
                                                        , const size_t& );
    int h( double, double, bool );
    std::vector< double > get_threshoulds();
};

inline std::vector< double > decisionstump_1d::get_threshoulds( void ) {
    return this->threshoulds;
}
#endif

