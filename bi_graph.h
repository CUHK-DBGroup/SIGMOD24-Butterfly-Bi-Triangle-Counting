#ifndef BI_GRAPH
#define BI_GRAPH

#include <vector>
#include <string>
#include <queue>
#include <fstream>
#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <bitset>
#include "def.h"
#include "util.h"

template <typename vertex_t>
struct Edge {
    Edge() = default;
    Edge(vertex_t u, vertex_t v): u(u), v(v) {}

    bool operator < (const Edge &o) const {
        if (u == o.u) return v < o.v;
        return u < o.u;
    }
    
    friend std::ostream &operator <<(std::ostream &out, const Edge &e) {
        out << "(" << e.u << ", " << e.v << ")";
        return out; 
    }
    vertex_t u, v;
};
using vertex_t = int;
using edge_t = Edge<vertex_t>;

namespace std{
template<> struct hash<edge_t>{
    size_t operator()(const edge_t &p) const {
        auto res = (static_cast<size_t>(p.u) << (sizeof(int) * 8)) + static_cast<size_t>(p.v);
        return res;
    }
};
template<> struct equal_to<edge_t>{
    bool operator()(const edge_t &p1, const edge_t &p2) const {
        return p1.u == p2.u && p1.v == p2.v;
    }
};
}

class BiGraph {
public:
    BiGraph() = default;
    BiGraph(const string &alias, bool cleaned);
    ~BiGraph() = default;

    void read_graph(const string &alias, const string &file_path);
    void read_clean(const string &alias, const string &file_path);
    void write_cleaned(const string &output_path);
    void calc_graph();
    ll compute_n_wedges();
    void clear();
    void resize();


    vector<vertex_t> csr_vec;
    vector<vertex_t> csr_offset_vec;
    void ObtainCsrDegDec(int core_num);

    string _alias;
    vertex_t n_vertex{};
    ll n_edges{};
    vertex_t l_vertex_num{};
    vertex_t r_vertex_num{};
    vertex_t largest_index_in[2]{};
    ll n_wedge_in[2]{};
 


    vector <edge_t> edges;
    vector <ll> sum_deg_neighbors;
    vector <ll> sum_wedges;
    vector <vector <vertex_t> > adj;

    void vertex_range(Side side, vertex_t &l, vertex_t &r) const;
    vertex_t side_vertex_num(Side side) const;

    // neighbor intersection method
    vertex_t exact_neighbor_intersection(vertex_t a, vertex_t b) const;
    double   fast_neighbor_intersection(vertex_t a, vertex_t b, int iter) const; // adj must be ordered
    vertex_t ordered_neighbor_intersection(vertex_t a, vertex_t b) const;

    // graph sampling for butterfly counting
    ld bfc_wedge_sampling               (URD<double> &urd, std::mt19937 &gen);
    ld bfc_fast_wedge_sampling          (UID<ll> &uid, mt19937_64 &gen, int edge_sample_num);
    ld bfc_fast_edge_sampling           (UID<ll> &uid, mt19937_64 &gen, int vertex_sample_num);
    ld bfc_vertex_sampling              (UID<int> &uid, mt19937 &gen);
    ld bfc_weighted_pair_sampling       (UID<ll> &uid, mt19937_64 &gen);
    ld bfc_pair_sampling                (UID<vertex_t> &uid, mt19937 &gen);


    // graph sampling for bi-triangle counting
    ld btc_triple_sampling              (UID<vertex_t> &uid, mt19937 &gen);
    ld btc_weighted_triple_sampling     (UID<ll> &uid, mt19937_64 &gen);
    ld btc_vertex_sampling              (UID<vertex_t> &uid, mt19937 &gen);
    ld btc_edge_sampling                (UID<ll> &uid, mt19937_64 &gen);
    ld btc_wedge_sampling               (URD<double> &urd, mt19937& gen);
    ld btc_wedge_sampling2              (UID<ll> &uid, mt19937_64& gen);


    // local butterfly counting
    ld bfc_per_edge(vertex_t u, vertex_t v); // by neighbor_intersection
    ld bfc_per_edge_hashmap(vertex_t u, vertex_t v); // by hashmap
    ld bfc_per_edge_approx(vertex_t u, vertex_t v, int iter);
    ld bfc_per_vertex(vertex_t u);
    ld bfc_per_pair_fast(vertex_t u, vertex_t v);



    // local bi-triangle counting
    ld btc_per_triple_fast(vertex_t u, vertex_t v, vertex_t w);
    ld btc_per_vertex_gt(vertex_t u);
    ld btc_per_vertex(vertex_t u);
    ld btc_per_vertex2(vertex_t u);
    ld btc_per_edge(vertex_t u2, vertex_t l);
    ld btc_per_wedge(vertex_t u, vertex_t x, vertex_t y);



    // local 3/4-path counting
    ld path3_per_edge   (vertex_t u, vertex_t v);
    ld path4_per_wedge  (vertex_t u, vertex_t x, vertex_t y);
    ld path4_per_pair   (vertex_t u, vertex_t v);

    // 3/4-path sampling methods
    ld path3_edge_sampling  (UID<ll> &uid, mt19937_64 &gen);
    ld path4_wedge_sampling (URD<double> &urd, mt19937& gen);
    ld path4_pair_sampling  (UID<vertex_t> &uid1, UID<vertex_t> &uid2, mt19937 &gen, double ratio);
    ld path4_weighted_pair_sampling (UID<ll> &uid, mt19937_64 &gen);



    // weighted sampling methods
    // call this function first before call weighted_sample
    ld weighted_vertex_setup(Side side);
    ld weighted_pair_setup(Side side);
    ld weighted_pair_bias(Side side);
    vertex_t pair_sample_setup(Side side);
    ld inv_weighted_pair_setup(Side side);
    ld weighted_triple_setup(Side side);
    ld wedge_setup(Side side);

    vertex_t weighted_sample_v(mt19937_64 &gen, UID<ll> &uid) const;
    vertex_t alias_sample_v(mt19937 &gen, URD<double> &urd) const;
    vertex_t pair_sample_v(mt19937 &gen, UID<vertex_t> &uid) const;

private:
    void add_edge(vertex_t u, vertex_t v);
    void ComputeCore(int core_num, vector<vertex_t> &deg_vec);
    void BuildCsrVec(vector<vertex_t> &nid_oid_vec, vector<vertex_t> &oid_nid_vec);


    std::unordered_map <vertex_t, vertex_t> vertex2id[2];
    std::unordered_set <edge_t> edges_set;

    vector <vertex_t> index_map;
    vector <vertex_t> index_map2;

    HashMap<vertex_t, ld> hm1, hm2, hm3;
    std::unordered_map<vertex_t, ld> um1, um2;


    // for weighted sample method
    vector <vertex_t> weighted_sample_table;
    vector <vertex_t> J;
    vector <double> q;
    vector <vertex_t> vertex_table;
};
#endif // BI_GRAPH