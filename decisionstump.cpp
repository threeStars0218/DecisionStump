#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
// #include "decisionstump_1d.hpp"
#include "decisionstump.hpp"

#define MATRIX  std::vector< std::vector< double > >
#define VECTOR  std::vector< double >
#define INDICES std::vector< size_t >
#define DAT     std::vector< std::vector< double > >
#define DAT_1D  std::vector< double >
#define LAB     std::vector< int >
#define DIST    std::vector< double >

#define SENSE        bool   // true -> h(x) = 1 if x ≥ θ else 0
#define EDGE         double // エッジ
#define THR          double // 閾値
#define COMP         int    // どの次元で比較したか
#define INDEX        size_t
#define HYPOTHESIS_1D std::pair< THR, SENSE > // 1次元データに対する仮説. (閾値, 不等号の向き) の二つ組で表現.
// #define HYPOTHESIS    std::pair< INDEX, std::pair< EDGE, HYPOTHESIS_1D > > // 仮説.
#define HYPOTHESIS    std::pair< INDEX, HYPOTHESIS_1D > // 仮説.

decisionstump::decisionstump() {};
decisionstump::decisionstump( DAT data_
                            , LAB label_
                            , bool one_side_
                            )
    : data( data_ )
    , label( label_ )
    , one_side( one_side_ )
    , m( data_.size() )
    , d( data_[0].size() )
    , number_of_stump( data_[0].size() * (data_.size() + 1) )
{
    std::cout << "preprocessing decisionstump.." << std::flush;
    for (int i=0; i<this->d; ++i) {
        DAT_1D tmp( this->m );
        for (int j=0; j<this->m; ++j) tmp[j] = this->data[j][i];
        dstumps.push_back( decisionstump_1d( tmp, label, one_side ) );
    }
    std::cout << "done!" << std::endl;
}

// HYPOTHESIS_1D std::pair< THR, SENSE > // 1次元データに対する仮説. (閾値, 不等号の向き) の二つ組で表現.
// HYPOTHESIS    std::pair< INDEX, HYPOTHESIS_1D > // 仮説.
HYPOTHESIS decisionstump::stump( DIST dist ) {
    std::pair< EDGE, HYPOTHESIS_1D > edge_and_hyp = dstumps[0].stump( dist );
    INDEX max_edge_index = 0;

    for (int i=1; i<this->d; ++i) {
        std::pair< EDGE, HYPOTHESIS_1D > tmp = dstumps[i].stump( dist );
        if (tmp.first > edge_and_hyp.first) {
            edge_and_hyp = tmp;
            max_edge_index = i;
        }
    }
    return std::make_pair( max_edge_index, edge_and_hyp.second );
//     return std::make_pair(max_edge_index, hypothesis);
}

// i番目のデータを第k成分でstump
int decisionstump::h( THR theta, double x, SENSE sns ) {
    int prediction;
    if (sns) {
        prediction = (theta < x) ? 1 : -1;
    } else {
        prediction = (theta > x) ? 1 : -1;
    }
    return prediction;
}

// HYPOTHESIS_1D std::pair< THR, SENSE > // 1次元データに対する仮説. (閾値, 不等号の向き) の二つ組で表現.
// HYPOTHESIS    std::pair< INDEX, HYPOTHESIS_1D > // 仮説.
VECTOR decisionstump::edge_vector( DIST dist ) {
    VECTOR     u( this->m, 0.0 );
    std::pair< INDEX, HYPOTHESIS_1D > idx_and_hyp = this->stump( dist );
    INDEX idx = idx_and_hyp.first;
    HYPOTHESIS_1D hh = idx_and_hyp.second;
    for (int k=0; k<this->m; ++k) {
        u[k] = (double) this->label[k] * this->dstumps[idx].h(hh.first, this->data[k][idx], hh.second);
    }
    double gamma = 0.0;
    for (int i=0; i<this->m; ++i) gamma += dist[i] * u[i];
    return u;
}

// MATRIX decisionstump::all_edge_vector_at( size_t idx ) {
//     std::cout << "idx = " << idx << std::endl;
//     MATRIX eg_mat = this->dstumps[idx].all_edge_vector( this->data, this->label, idx );
//     return eg_mat;
// }
// std::pair< VECTOR, MATRIX > all_edge_vector_with_edge_at( size_t ) {
// }
