#include "decisionstump_1d.hpp"

int
predict(double threshould,
        bool sense,
        double data)
{
    if (sense)
        return (data > threshould ? 1 : -1);
    else
        return (data < threshould ? 1 : -1);
}

decisionstump_1d::
decisionstump_1d() {};

decisionstump_1d::
decisionstump_1d(std::vector<double> data_,
                 std::vector<int> label_,
                 bool one_side_)
    : sorted_index(label_.size())
    , sorted_data(label_.size())
    , sorted_label(label_.size())
    , m(label_.size())
    , one_side(one_side_)
    , eps( 0.5 )
{
    // =====
    // sort data to handle easily
    // =====
    std::iota(this->sorted_index.begin(), this->sorted_index.end(), 0);
    std::sort(this->sorted_index.begin(),
              this->sorted_index.end(),
              [&data_] (size_t i, size_t j) {return data_[i] < data_[j];});
    for(int i=0; i<this->m; ++i) {
        this->sorted_data[i]  = data_[ this->sorted_index[i] ];
        this->sorted_label[i] = label_[ this->sorted_index[i] ];
    }

    // =====
    // keep threshoulds, and hypotheses
    // =====
    std::vector<double> tmp = this->sorted_data;
    tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
    double thr = tmp[0] - this->eps;
    this->threshoulds.push_back(thr);
    std::function<int(double)> f, g;
    f = std::bind(predict, thr, true, std::placeholders::_1);
    this->classifiers.push_back(f);
    if (!one_side) {
        g = std::bind(predict, thr, false, std::placeholders::_1);
        this->classifiers.push_back(g);
    }
    for(int i=0; i<tmp.size()-1; ++i) {

        double thr = (tmp[i] + tmp[i+1]) / 2.0;
        this->threshoulds.push_back(thr);
        f = std::bind(predict, thr, true, std::placeholders::_1);
        this->classifiers.push_back(f);
        if (!one_side) {
            g = std::bind(predict, thr, false, std::placeholders::_1);
            this->classifiers.push_back(g);
        }
    }
    thr = tmp[tmp.size()-1] + this->eps;
    this->threshoulds.push_back(thr);
    f = std::bind(predict, thr, true, std::placeholders::_1);
    this->classifiers.push_back(f);
    if (!one_side) {
        g = std::bind(predict, thr, false, std::placeholders::_1);
        this->classifiers.push_back(g);
    }

    // =====
    // calculate number of stumps
    // =====
    this->number_of_stump = this->threshoulds.size();
    if (!one_side) this->number_of_stump *= 2;
}

std::pair<double, std::function<int(double)> >
decisionstump_1d::
stump_one_side(std::vector<double> dist)
{
    auto f_ptr = this->classifiers.begin();
    std::function<int(double)> f = *f_ptr;
    std::function<int(double)> max_edge_classifier = *f_ptr;
    double max_edge = this->naive_calc_edge(dist, f);

    double edge = max_edge;
    size_t number_of_classifiers = this->classifiers.size();
    auto thr_ptr = this->threshoulds.begin()+1;
    int k = 0;
    f_ptr += 2;
    while(k < this->m && f_ptr < this->classifiers.end()) {
        f = *f_ptr;
        while (k < this->m && this->sorted_data[k] < *thr_ptr) {
            edge += 2 * dist[this->sorted_index[k]]
                      * this->sorted_label[k]
                      * f(this->sorted_data[k]);
            ++k;
        }
        if (edge > max_edge) {
            max_edge = edge;
            max_edge_classifier = f;
        }
        f_ptr += 2;
        ++thr_ptr;
    }
    return std::make_pair(max_edge, max_edge_classifier);
}

std::pair<double, std::function<int(double)> >
decisionstump_1d::
stump(std::vector<double> dist)
{
    if (this->one_side)
        return this->stump_one_side(dist);
    int    k = 0;
    auto   f_ptr  = this->classifiers.begin();
    auto   f = *f_ptr;
    double edge_1 = this->naive_calc_edge(dist, f);
    double edge_2 = -1 * edge_1;
    double max_edge            = (edge_1 > edge_2) ? edge_1 : edge_2;
    auto   max_edge_classifier = (edge_1 > edge_2) ? *f_ptr : *(f_ptr+1);
    auto   thr_ptr = this->threshoulds.begin();

    while(k < this->m && f_ptr < this->classifiers.end()) {
        f = *f_ptr;
        while(k < this->m && this->sorted_data[k] < *thr_ptr) {
            edge_1 += 2 * dist[this->sorted_index[k]]
                        * this->sorted_label[k]
                        * f(this->sorted_data[k]);
            ++k;
        }
        edge_2 = -1 * edge_1;
        if (std::max(edge_1, edge_2) > max_edge) {
            max_edge_classifier = (edge_1 > edge_2) ? *f_ptr : *(f_ptr+1);
            max_edge            = (edge_1 > edge_2) ? edge_1 : edge_2;
        }
        ++thr_ptr;
        f_ptr += 2;
    }

    return std::make_pair(max_edge, max_edge_classifier);
}


double
decisionstump_1d::
naive_calc_edge(std::vector<double> dist,
                std::function<int(double)> f)
{
    double edge = 0.0;
    for (int k=0; k<this->m; ++k)
        edge += dist[this->sorted_index[k]] *
                this->sorted_label[k] * f(this->sorted_data[k]);
    return edge;
}

std::vector<std::vector<double> >
decisionstump_1d::
all_edge_vector(const std::vector<std::vector<double> >&dat,
                const std::vector<int> &lab,
                const size_t &idx )
{
    double theta;
    if (this->ev_mat.size() > 0)
        return this->ev_mat;
    std::vector<double> v(this->m), w(this->m);
    for(std::function<int(double)> f : this->classifiers) {
        for(int i=0; i<this->m; ++i) {
            v[i] = lab[i] * f(dat[i][idx]);
        }
        this->ev_mat.push_back(v);
    }
    return this->ev_mat;
}

