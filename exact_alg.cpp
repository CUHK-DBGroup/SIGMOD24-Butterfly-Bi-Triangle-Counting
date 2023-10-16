#include "exact_alg.h"
#include <map>

namespace tracker {
    vector<ll> hashmap_C;
    vector<vertex_t> aux;
    std::map<std::pair<int, int>, int> mp_cnt;
}


ll tracker::bfc_exact_pair(const BiGraph &graph) {
  int side = graph.n_wedge_in[0] > graph.n_wedge_in[1]; // ?
  ll res = 0;
  for (vertex_t u = (side == 0 ? 0 : graph.l_vertex_num); u < graph.largest_index_in[side]; u++) {
    auto ver_size = (vertex_t) graph.adj[u].size();
    for (vertex_t j = 0; j < ver_size; j++) {
      for (vertex_t k = 0; k < j; k++) {
        mp_cnt[std::make_pair(graph.adj[u][j], graph.adj[u][k])]++;
      }
    }
  }
  for (const auto &it: mp_cnt) {
    res = res + (ll) (it.second) * (it.second - 1) / 2;
  }
  mp_cnt.clear();
  return res;
}

ld tracker::bfc_exact(const BiGraph &graph) {
  hashmap_C.resize(graph.n_vertex);
  aux.resize(graph.n_vertex);

  const auto &adj = graph.adj;

  ll n_wedge_in[2] = {0, 0};
  vertex_t l = 0, r = 0;
  graph.vertex_range(L_SIDE, l, r);
  for (vertex_t u = l; u < r; ++u)
    n_wedge_in[0] += cn2(static_cast<ll>(adj[u].size()));

  graph.vertex_range(R_SIDE, l, r);
  for (vertex_t u = l; u < r; ++u)
    n_wedge_in[1] += cn2(static_cast<ll>(adj[u].size()));

  Side side = n_wedge_in[0] < n_wedge_in[1] ? L_SIDE : R_SIDE;
  graph.vertex_range(side, l, r);

  ld res = 0;
  for (vertex_t u = l; u < r; ++u) {
    vertex_t idx = 0;
    for (vertex_t u1hop: adj[u]) {
      for (vertex_t u2hop: adj[u1hop]) {
        if (u <= u2hop) break;
        if (!hashmap_C[u2hop]) aux[idx++] = u2hop;
        res += hashmap_C[u2hop];
        ++hashmap_C[u2hop];
      }
    }
    for (int j = 0; j < idx; ++j) {
      hashmap_C[aux[j]] = 0;
    }
  }

  return res;
}

//
ld tracker::btc_exact_swja(BiGraph &graph) {
  graph.ObtainCsrDegDec(0);
  // step 1 : compute the number of bi-triangles with invalid ones included.
  auto n_vertices = graph.n_vertex;
  vector<size_t> h2_nb_w_vec(n_vertices), h3_nb_w_vec(n_vertices);
  vector<int> h2_nb_id_vec(n_vertices), h3_nb_id_vec(n_vertices);
  int h2_nb_num = 0, h3_nb_num = 0;
  const auto &csr_vec = graph.csr_vec;
  const auto &csr_offset_vec = graph.csr_offset_vec;
  ld num_invalid_g1 = 0, num_invalid_g2 = 0, num_invalid_g3 = 0, num_invalid_g4 = 0; //T4.3 included
  ld num_bi_triangle = 0;

  int size_vec = 10;
  vector<ld> num_bi_triangle_vec(size_vec);
  ld g1_base = 0, g2_base = 0, g4_base = 0;
  vector<size_t> g2_base_vec(n_vertices);
  for (vertex_t vid = 0; vid < n_vertices; vid++) {
    g4_base = 0;
    h2_nb_num = 0;
    h3_nb_num = 0;
    int nb_stt_pos = csr_offset_vec[vid];
    int nb_end_pos = csr_offset_vec[vid + 1];
    // obtain hop-2 neighbor weight
    for (int nb_pos = nb_stt_pos; nb_pos < nb_end_pos; nb_pos++) {
      int nb_vid = csr_vec[nb_pos];
      if (nb_vid <= vid) {
        g4_base = nb_end_pos - nb_pos;
        break;
      }
      int h2_nb_stt_pos = csr_offset_vec[nb_vid];
      int h2_nb_end_pos = csr_offset_vec[nb_vid + 1];
      // obtain its hop-2 neighbors weight;
      g2_base = h2_nb_end_pos - h2_nb_stt_pos - 2;
      for (int h2_nb_pos = h2_nb_stt_pos; h2_nb_pos < h2_nb_end_pos; h2_nb_pos++) {
        int h2_nb_vid = csr_vec[h2_nb_pos];
        if (h2_nb_vid <= vid) {
          g1_base = h2_nb_pos - h2_nb_stt_pos;
          num_invalid_g1 += (g1_base) * (g1_base - 1) / 2;
          break;
        }
        //if (h2_nb_w_vec[h2_nb_vid]++ == 0) h2_nb_id_vec[h2_nb_num++] = h2_nb_vid;
        if (h2_nb_w_vec[h2_nb_vid] == 0) {
          h2_nb_id_vec[h2_nb_num] = h2_nb_vid;
          h2_nb_num++;
        }
        h2_nb_w_vec[h2_nb_vid]++;
        g2_base_vec[h2_nb_vid] += g2_base;
      }
    }
    // obtain its hop-3 neighbors weight;
    for (int h2_nb_pos = 0; h2_nb_pos < h2_nb_num; h2_nb_pos++) {
      int h2_nb_vid = h2_nb_id_vec[h2_nb_pos];
      ld h2_nb_w = h2_nb_w_vec[h2_nb_vid];
      h2_nb_w_vec[h2_nb_vid] = 0;
      num_invalid_g2 += g2_base_vec[h2_nb_vid] * (h2_nb_w - 1);
      g2_base_vec[h2_nb_vid] = 0;
      int h3_nb_stt_pos = csr_offset_vec[h2_nb_vid];
      int h3_nb_end_pos = csr_offset_vec[h2_nb_vid + 1];
      ld num_butterfly = h2_nb_w * (h2_nb_w - 1) / 2;
      num_invalid_g3 += num_butterfly * (h3_nb_end_pos - h3_nb_stt_pos);
      num_invalid_g4 += num_butterfly * g4_base;
      for (int h3_nb_pos = h3_nb_stt_pos; h3_nb_pos < h3_nb_end_pos; h3_nb_pos++) {
        int h3_nb_vid = csr_vec[h3_nb_pos];
        if (h3_nb_vid <= vid) {
          break;
        }
        if (h3_nb_w_vec[h3_nb_vid] == 0) h3_nb_id_vec[h3_nb_num++] = h3_nb_vid;
        h3_nb_w_vec[h3_nb_vid] += h2_nb_w;
      }
    }
    for (int h3_nb_pos = 0; h3_nb_pos < h3_nb_num; h3_nb_pos++) {
      int h3_nb_vid = h3_nb_id_vec[h3_nb_pos];
      ld h3_nb_w = h3_nb_w_vec[h3_nb_vid];
      h3_nb_w_vec[h3_nb_vid] = 0;
      num_bi_triangle_vec[vid % size_vec] += h3_nb_w * (h3_nb_w - 1) / 2;
    }
  }
  for (int i = 0; i < size_vec; i++) {
    num_bi_triangle += num_bi_triangle_vec[i];
  }
  num_bi_triangle -= num_invalid_g1;
  num_bi_triangle -= num_invalid_g2;
  num_bi_triangle -= num_invalid_g3;
  num_bi_triangle -= num_invalid_g4;
  return num_bi_triangle;
}

ld tracker::calc_4path_opt(BiGraph &graph, Side side) {
  ld res = 0;
  int l = 0, r = 0;
  graph.vertex_range(side, l, r);

  for (int i = l; i < r; ++i) {
    const auto &us = graph.adj[i];
    ld s = 0, s2 = 0;
    for (vertex_t v : us) {
      auto dv = graph.adj[v].size() - 1;
      s += dv;
      s2 += dv * dv;
    }
    res += s * s - s2;
  }
  return res / 2;
}