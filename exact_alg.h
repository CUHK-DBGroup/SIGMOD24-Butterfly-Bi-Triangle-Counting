#ifndef BUTTERFLY_COUNTING_H
#define BUTTERFLY_COUNTING_H

#include "bi_graph.h"
#include "def.h"

namespace tracker {
    ll bfc_exact_pair(const BiGraph &graph);

    ld bfc_exact(const BiGraph &graph);

    ld btc_exact_swja(BiGraph &graph);

    ld calc_4path_opt(BiGraph &graph, Side side);
}


#endif // BUTTERFLY_COUNTING_H
