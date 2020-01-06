#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <vector>

#include "decisionstump_1d.hpp"

using MATRIX = std::vector< std::vector< double > >;
using VECTOR = std::vector< double >;
using INDICES = std::vector< size_t >;
using DAT = std::vector< std::vector< double > >;
using DAT_1D = std::vector< double >;
using LAB = std::vector< int >;
using DIST = std::vector< double >;

using SENSE = bool; // true -> h(x) = 1 if x ≥ θ else 0
using EDGE = double; // エッジ
using THR = double; // 閾値
using COMP = int; // どの次元で比較したか
using INDEX = size_t;
using HYPOTHESIS_1D = std::pair< THR, SENSE >; // 1次元データに対する仮説. (閾値, 不等号の向き) の二つ組で表現.;
using HYPOTHESIS = std::pair< COMP, HYPOTHESIS_1D >;// 仮説.

decisionstump_1d::decisionstump_1d() {};
decisionstump_1d::decisionstump_1d( DAT_1D data_
                                  , LAB    label_
                                  , bool   one_side_
                                  )
    : sorted_index( label_.size() )
    , sorted_data( label_.size() )
    , sorted_label( label_.size() )
    , m( label_.size() )
    , one_side( one_side_ )
    , eps( 0.5 )
{
    std::iota( this->sorted_index.begin(), this->sorted_index.end(), 0 );
    std::sort( this->sorted_index.begin()
             , this->sorted_index.end()
             , [&data_] (size_t i, size_t j) {return data_[i] < data_[j];} );
    for (int i=0; i<this->m; ++i) {
        this->sorted_data[i]  = data_[ this->sorted_index[i] ];
        this->sorted_label[i] = label_[ this->sorted_index[i] ];
    }

    this->threshoulds.push_back( this->sorted_data[0] - this->eps );
    for (int i=0; i<this->m-1; ++i) {
        while (i < this->m-1 && this->sorted_data[i] == this->sorted_data[i+1]) ++i;
        if (i == this->m-1) break;
        double thr = (this->sorted_data[i] + this->sorted_data[i+1]) / 2.0;
        this->threshoulds.push_back( thr );
    }
    this->threshoulds.push_back( this->sorted_data[m-1] + this->eps );
}

std::pair< EDGE, HYPOTHESIS_1D > decisionstump_1d::stump_one_side( DIST dist ) {
    THR  theta = this->sorted_data[0] - this->eps;
    EDGE edge  = this->naive_calc_edge( theta, dist, true );
    EDGE edge_tmp = edge;
    HYPOTHESIS_1D hypothesis_1d;
    for (int k=0; k<this->m; ++k) {
        int idx = k;
        while (idx < this->m && this->sorted_data[idx] == this->sorted_data[idx+1]) ++idx;
        THR  theta_tmp =  (this->sorted_data[idx] + this->sorted_data[idx+1]) / 2.0;
        if (idx == this->m-1) theta_tmp = this->sorted_data[m-1] + this->eps;
        for (int s=k; s<idx+1; ++s) {
            edge_tmp += 2*dist[this->sorted_index[s]]*this->sorted_label[s]*this->h(theta_tmp, this->sorted_data[s], true);
        }
        k = idx;
        if ( edge_tmp > edge ) {
            edge  = edge_tmp;
            theta = theta_tmp;
        }
    }

    hypothesis_1d.first = theta; hypothesis_1d.second = true;
    return std::make_pair( edge, hypothesis_1d );
}
std::pair< EDGE, HYPOTHESIS_1D > decisionstump_1d::stump( DIST dist ) {
    if (this->one_side) return this->stump_one_side( dist );
    THR  theta  = this->sorted_data[0] - this->eps;
    EDGE edge_1 = this->naive_calc_edge( theta, dist, true );
    EDGE edge_2 = this->naive_calc_edge( theta, dist, false );
    EDGE edge   = std::max(edge_1, edge_2);
    SENSE sns   = (edge_1 >= edge_2) ? true : false;
    HYPOTHESIS_1D hypothesis_1d;
    for (int k=0; k<this->m; ++k) {
        int idx = k;
        while (idx < this->m && this->sorted_data[idx] == this->sorted_data[idx+1]) ++idx;
        THR  theta_tmp =  (this->sorted_data[idx] + this->sorted_data[idx+1]) / 2.0;
        if (idx == this->m-1) theta_tmp = this->sorted_data[m-1] + this->eps;
        // THR  theta_tmp = (this->sorted_data[k] + this->sorted_data[k+1]) / 2.0;

        for (int s=k; s<idx+1; ++s) {
            edge_1 += 2*dist[this->sorted_index[s]]*this->sorted_label[s]*this->h(theta_tmp, this->sorted_data[s], true);
            edge_2 += 2*dist[this->sorted_index[s]]*this->sorted_label[s]*this->h(theta_tmp, this->sorted_data[s], false);
        }
        k = idx;
        if ( std::max(edge_1, edge_2) >= edge ) {
            edge  = std::max(edge_1, edge_2);
            theta = theta_tmp;
            sns   = (edge_1 >= edge_2) ? true : false;
        }
    }

    hypothesis_1d.first = theta; hypothesis_1d.second = sns;
    return std::make_pair( edge, hypothesis_1d );
}

int decisionstump_1d::h( THR theta, double val, SENSE sns ) {
    int prediction;
    if (sns) {
        prediction = (theta < val) ? 1 : -1;
    } else {
        prediction = (theta > val) ? 1 : -1;
    }
    return prediction;
}


EDGE decisionstump_1d::naive_calc_edge( THR theta, DIST dist, SENSE sns ) {
    EDGE edge = 0.0;
    for (int k=0; k<this->m; ++k) {
        double yh = this->sorted_label[k] * this->h( theta, this->sorted_data[k], sns );
        edge += dist[ this->sorted_index[k] ] * yh;
    }
    return edge;
}

MATRIX decisionstump_1d::all_edge_vector( const DAT &dat, const LAB &lab, const size_t &idx ) {
    MATRIX ev_mat;
    THR theta;
    VECTOR v(this->m), w(this->m);
    for (const double theta : this->threshoulds) {
        for (int i=0; i<this->m; ++i) {
            v[i] = lab[i] * this->h( theta, dat[i][idx], true );
            if (!this->one_side) w[i] = lab[i] * this->h( theta, dat[i][idx], false );
        }
        ev_mat.push_back(v);
        if (!one_side) ev_mat.push_back(w);
    }
    return ev_mat;
}

// double decisionstump_1d::flip_prediction( THR theta, EDGE edge, INDEX idx, SENSE sns, DIST dist ) {
//     double yh = this->sorted_label[idx] * this->h( theta, this->sorted_data[idx], sns );
//     return edge + 2 * dist[ this->sorted_index[idx] ] * yh;
// }

