#ifndef BFC_TRACKER_H
#define BFC_TRACKER_H

#include "bi_graph.h"
#include "def.h"

namespace tracker {

    void set_tracker(const string &tracker);

    void path4_exact_algorithm_tracker(BiGraph &graph, Side side, const string &gt_path);

    // exact bfc time tracker
    void bfc_exact_algorithm_tracker(BiGraph &graph, const string &gt_path);

    void bfc_fast_edge_sampling_tracker(BiGraph &graph, int rnd, int vertex_sample_num);
    

    void bfc_vertex_sampling_tracker(BiGraph &graph, int rnd, Side side);

    void bfc_wedge_sampling_tracker(BiGraph &graph, int rnd, Side side);

    // pair sampling
    void bfc_pair_sampling_tracker(BiGraph &graph, int rnd, Side side = L_SIDE);

    void bfc_weighted_pair_sampling_tracker(BiGraph &graph, int rnd, Side side = L_SIDE);

    // btc
    void btc_triple_sampling_tracker(BiGraph &graph, int rnd, Side side = L_SIDE);

    void btc_weighted_triple_sampling_tracker(BiGraph &graph, int rnd, Side side = L_SIDE);

    void btc_vertex_sampling_tracker(BiGraph &graph, int rnd, Side side = BOTH);

    void btc_edge_sampling_tracker(BiGraph &graph, int rnd);

    void btc_wedge_sampling_tracker(BiGraph &graph, int rnd, Side side = BOTH);

    // 3/4-path
    void path3_edge_sampling_tracker(BiGraph &graph, int rnd);

    void path4_wedge_sampling_tracker(BiGraph &graph, int rnd);

    void path4_pair_sampling_tracker(BiGraph &graph, int rnd);

    void path4_weighted_pair_sampling_tracker(BiGraph &graph, int rnd, Side side);

    void btc_exact_algorithm_tracker(BiGraph &graph, const string &btc_path);
}


#endif //BFC_TRACKER_H