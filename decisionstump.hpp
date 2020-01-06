#include <iostream>
#include <vector>
#include "decisionstump_1d.hpp"

#ifndef DECISIONSTUMP_HPP
#define DECISIONSTUMP_HPP

class decisionstump {
private:
    std::vector< std::vector< double > > data;
    std::vector< int > label;
    std::vector< decisionstump_1d > dstumps;
    bool one_side;
    size_t m; // サンプル数
    size_t d; // 各サンプルの次元数
    size_t number_of_stump;
public:
    decisionstump();
    decisionstump( std::vector< std::vector< double > >, std::vector< int >, bool );
    std::pair< size_t, std::pair< double, bool > > stump( std::vector< double > );
    std::vector< double > edge_vector( std::vector< double > );
    std::vector< std::vector< double > > all_edge_vector_at( size_t );
    // std::pair< std::vector< double >, std::vector< std::vector< double > > > all_edge_vector_with_edge_at( size_t );
    int h( double, double, bool );
    size_t get_number_of_stump();
    std::vector< double > get_threshoulds_at( size_t );
};

inline std::vector< std::vector< double > > decisionstump::all_edge_vector_at( size_t idx ) {
    return this->dstumps[idx].all_edge_vector( this->data, this->label, idx );
}
inline size_t decisionstump::get_number_of_stump() {
    return this->number_of_stump;
}

inline std::vector< double > decisionstump::get_threshoulds_at( size_t idx ) {
    return this->dstumps[idx].get_threshoulds();
}
#endif
