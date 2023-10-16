#include <sstream>
#include <ctime>
#include <algorithm>
#include <queue>
#include "bi_graph.h"
#include "util.h"


BiGraph::BiGraph(const string &alias, bool cleaned) {
  clear();
  std::cerr << " Processing the graph ... (please wait)" << std::endl;

  string file_path = get_dataset_path(alias, cleaned);
  if (cleaned) {
    read_clean(alias, file_path);
  } else {
    read_graph(alias, file_path);
  }
  std::cerr << " End Read graph!" << std::endl;
  std::cerr << " -n: " << n_vertex << ", -m: " << n_edges
            << ", -l: " << l_vertex_num << ", -r: " << r_vertex_num << std::endl;

  resize();
  calc_graph();
  std::cerr << " End calc the graph!" << std::endl;
  std::cerr << " n_wedge_in[0]: " << n_wedge_in[0] << ", n_wedges_in[1]: " << n_wedge_in[1] << "\n";
}

void BiGraph::read_graph(const string &alias, const string &file_path) {
  _alias = alias;

  std::ifstream fin(file_path);
  string buf;
  while (std::getline(fin, buf)) {
    const char *s = buf.c_str();
    const int len = (int) buf.length();
    vertex_t u, v;
    int ind = to_int(s, len, u);
    if (ind != -1) {
      ind = to_int(s + ind, len - ind, v);
      if (ind == -1) continue;
    } else continue;

    add_edge(u, v);
  }

  n_vertex = l_vertex_num + r_vertex_num;
  n_edges = static_cast<ll>(edges_set.size());
  largest_index_in[0] = l_vertex_num;
  largest_index_in[1] = n_vertex;

  adj.resize(n_vertex, vector<vertex_t>());

  for (const edge_t &edge: edges_set) {
    vertex_t u = vertex2id[0][edge.u];
    vertex_t v = vertex2id[1][edge.v] + l_vertex_num;
    adj[u].push_back(v);
    adj[v].push_back(u);
    edges.emplace_back(u, v);
  }

  fin.close();
}

void BiGraph::read_clean(const string &alias, const string &file_path) {
  _alias = alias;
  std::ifstream fin(file_path);

  string buf;
  std::getline(fin, buf);

  std::stringstream ss;
  ss << buf;
  while (ss >> buf) {
    if (buf == "-n") {
      ss >> n_vertex;
    } else if (buf == "-m") {
      ss >> n_edges;
    } else if (buf == "-l") {
      ss >> l_vertex_num;
    } else if (buf == "-r") {
      ss >> r_vertex_num;
    }
  }

  largest_index_in[0] = l_vertex_num;
  largest_index_in[1] = n_vertex;

  edges.resize(n_edges);
  adj.resize(n_vertex, vector<vertex_t>());

  for (ll i = 0; i < n_edges; ++i) {
    vertex_t u, v;
    fin >> u >> v;
    v += l_vertex_num;
    edges[i] = {u, v};
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  for (auto &us: adj) std::sort(us.begin(), us.end());
  fin.close();
}

void BiGraph::calc_graph() {
  n_wedge_in[0] = 0;
  for (vertex_t v = 0; v < largest_index_in[0]; v++) {
    n_wedge_in[0] += cn2(static_cast<ll> (adj[v].size()));
  }
  n_wedge_in[1] = 0;
  for (vertex_t v = largest_index_in[0]; v < largest_index_in[1]; v++) {
    n_wedge_in[1] += cn2(static_cast<ll>(adj[v].size()));
  }

  for (vertex_t v = 0; v < n_vertex; v++) {
    sort(adj[v].begin(), adj[v].end());

    sum_deg_neighbors[v] = 0;
    for (const vertex_t &neighbor: adj[v]) {
      sum_deg_neighbors[v] += static_cast<ll>(adj[neighbor].size());
    }
  }
}

void BiGraph::write_cleaned(const string &output_path) {
  std::cerr << " Write cleaned graph to: " << output_path << std::endl;

  std::ofstream fout(output_path);
  fout << " -n " << n_vertex
       << " -m " << n_edges
       << " -l " << l_vertex_num
       << " -r " << r_vertex_num << std::endl;

  std::sort(edges.begin(), edges.end());
  for (const auto &edge: edges) {
    fout << edge.u << " " << (edge.v - l_vertex_num) << "\n";
  }
  fout.close();
}

void BiGraph::add_edge(const vertex_t u, const vertex_t v) {
  if (vertex2id[0].find(u) == vertex2id[0].end()) {
    vertex2id[0][u] = l_vertex_num++;
  }
  if (vertex2id[1].find(v) == vertex2id[1].end()) {
    vertex2id[1][v] = r_vertex_num++;
  }
  edges_set.insert({u, v});
}

void BiGraph::vertex_range(Side side, vertex_t &l, vertex_t &r) const {
  if (side == Side::L_SIDE) l = 0, r = largest_index_in[0];
  else if (side == Side::R_SIDE) l = largest_index_in[0], r = largest_index_in[1];
  else if (side == Side::BOTH) l = 0, r = largest_index_in[1];
}

vertex_t BiGraph::side_vertex_num(Side side) const {
  if (side == L_SIDE) return l_vertex_num;
  if (side == R_SIDE) return r_vertex_num;
  return n_vertex;
}

void BiGraph::ObtainCsrDegDec(int core_num) {
  vector<vertex_t> deg_vec;
  ComputeCore(core_num, deg_vec);

  vector<vertex_t> nid_oid_vec = vector<int>(n_vertex, -1);
  vector<vertex_t> oid_nid_vec = vector<int>(n_vertex, -1);
  vector<vertex_t> deg_offset_vec = vector<int>(n_vertex);
  int max_deg = 0;
  for (int vid = 0; vid < n_vertex; vid++) {
    int deg = deg_vec[vid];
    if (deg < core_num) continue;
    max_deg = max_deg > deg ? max_deg : deg;
    deg_offset_vec[deg]++;
  }
  for (int deg = max_deg - 1; deg >= 2; deg--) deg_offset_vec[deg] += deg_offset_vec[deg + 1];
  for (int oid = 0; oid < n_vertex; oid++) {
    int deg = deg_vec[oid];
    if (deg < core_num) continue;
    nid_oid_vec[deg_offset_vec[deg + 1]++] = oid;
  }

  for (int nid = 0; nid < n_vertex; nid++) {
    int oid = nid_oid_vec[nid];
    if (oid == -1) break;
    oid_nid_vec[oid] = nid;
  }
  BuildCsrVec(nid_oid_vec, oid_nid_vec);
}

void BiGraph::ComputeCore(int core_num, vector<vertex_t> &deg_vec) {
  vector<bool> is_pushed = vector<bool>(n_vertex, false);
  std::queue<int> deleting_q;
  deg_vec.resize(n_vertex);
  for (int l_id = 0; l_id < largest_index_in[0]; l_id++) {
    deg_vec[l_id] = static_cast<vertex_t>(adj[l_id].size());
    if (deg_vec[l_id] < core_num) {
      deleting_q.push(l_id);
      is_pushed[l_id] = true;
    }
  }
  for (int r_id = largest_index_in[0]; r_id < largest_index_in[1]; r_id++) {
    deg_vec[r_id] = static_cast<vertex_t>(adj[r_id].size());
    if (deg_vec[r_id] < core_num) {
      deleting_q.push(r_id);
      is_pushed[r_id] = true;
    }
  }
  while (!deleting_q.empty()) {
    vertex_t vid = deleting_q.front();
    deleting_q.pop();
    for (vertex_t nb_id: adj[vid]) {
      if (is_pushed[nb_id]) continue;
      deg_vec[nb_id]--;
      if (deg_vec[nb_id] < core_num) {
        deleting_q.push(nb_id);
        is_pushed[nb_id] = true;
      }
    }
  }
}

void BiGraph::BuildCsrVec(vector<vertex_t> &nid_oid_vec, vector<vertex_t> &oid_nid_vec) {
  csr_vec.resize(2 * n_edges);
  csr_offset_vec.resize(n_vertex + 1);

  vertex_t new_l_v_num_ = 0;
  vertex_t new_r_v_num_ = 0;
  vertex_t new_vertex_num = 0;
  vertex_t new_edge_num = 0;
  vertex_t csr_offset_pos = 0;
  for (vertex_t nid = 0; nid < n_vertex; nid++) {
    vertex_t oid = nid_oid_vec[nid];
    if (oid == -1) break;
    if (oid < l_vertex_num) {
      new_l_v_num_++;
    } else {
      new_r_v_num_++;
    }
    new_vertex_num++;
    csr_offset_vec[nid] = csr_offset_pos;
    for (int nb_oid: adj[oid]) {
      int nb_nid = oid_nid_vec[nb_oid];
      if (nb_nid == -1) continue;
      new_edge_num++;
      csr_vec[csr_offset_pos++] = nb_nid;
    }
    sort(csr_vec.begin() + csr_offset_vec[nid],
         csr_vec.begin() + csr_offset_pos,
         [](auto a, auto b) -> bool { return a > b; }
    );
  }
  csr_offset_vec[new_vertex_num] = csr_offset_pos;
  csr_offset_vec.erase(csr_offset_vec.begin() + new_vertex_num + 1, csr_offset_vec.end());
  csr_vec.erase(csr_vec.begin() + new_edge_num, csr_vec.end());
  new_edge_num /= 2;
}


vertex_t BiGraph::exact_neighbor_intersection(vertex_t a, vertex_t b) const {
  if (adj[a].size() > adj[b].size()) std::swap(a, b);
  vertex_t res = 0;
  std::unordered_set<vertex_t> s;
  for (const vertex_t v: adj[a]) s.insert(v);
  for (const vertex_t v: adj[b]) {
    if (s.find(v) != s.end())
      ++res;
  }
  return res;
}

double BiGraph::fast_neighbor_intersection(vertex_t a, vertex_t b, int iter) const {
  if (adj[a].size() > adj[b].size()) std::swap(a, b);
  std::random_device rdw;
  std::mt19937 gene(rdw());
  std::uniform_int_distribution<vertex_t> dis(0, (int) adj[a].size() - 1);

  double num_wedges = 0;
  for (int i = 0; i < iter; i++) {
    int c = adj[a][dis(gene)];
    if (binary_search(adj[b].begin(), adj[b].end(), c)) {
      num_wedges += static_cast<double>(adj[a].size());
    }
  }
  num_wedges /= iter;
//	return num_wedges * (num_wedges - 1) / 2;
  return num_wedges; // note this!
}

vertex_t BiGraph::ordered_neighbor_intersection(vertex_t a, vertex_t b) const {
  size_t size_a = adj[a].size(), size_b = adj[b].size();
  size_t i = 0, j = 0;
  vertex_t cnt = 0;
  while (i < size_a && j < size_b) {
    if (adj[a][i] == adj[b][j]) {
      ++cnt;
      ++i;
      ++j;
    } else if (adj[a][i] < adj[b][j]) ++i;
    else ++j;
  }
  return cnt;
}

ll BiGraph::compute_n_wedges() {
  ll wedges = 0;
  for (int i = 0; i < n_vertex; i++) {
    wedges += ((ll) adj[i].size()) * ((ll) adj[i].size() - 1) >> 1;
    sum_wedges[i] = wedges;
  }
  return wedges;
}

void BiGraph::clear() {
  largest_index_in[0] = largest_index_in[1] = 0;
  n_vertex = 0;
  n_edges = 0;
  adj.clear();
  sum_wedges.clear();
  index_map.clear();
  index_map2.clear();
  hm1.clear();
  hm2.clear();
  hm3.clear();
}

void BiGraph::resize() {
  sum_deg_neighbors.resize(n_vertex);
  sum_wedges.resize(n_vertex);
  index_map.resize(n_vertex);
  index_map2.resize(n_vertex);
  hm1.resize(n_vertex);
  hm2.resize(n_vertex);
  hm3.resize(n_vertex);
}


/*========= bfc sampling method =========*/
ld BiGraph::bfc_wedge_sampling(URD<double> &urd, std::mt19937 &gen) {
  vertex_t u = alias_sample_v(gen, urd);
  while (adj[u].size() < 2) {
    u = alias_sample_v(gen, urd);
  }

  vertex_t v1 = 0, v2 = 1;
  if (adj[u].size() > 2) {
    UID<vertex_t> uid(0, static_cast<vertex_t>(adj[u].size() - 1));
    v1 = uid(gen), v2 = uid(gen);
    while (v1 == v2) {
      v1 = uid(gen), v2 = uid(gen);
    }
  }

  v1 = adj[u][v1];
  v2 = adj[u][v2];
  return exact_neighbor_intersection(v1, v2) - 1;
}

ld BiGraph::bfc_fast_wedge_sampling(UID<ll> &uid, mt19937_64 &gen,  int edge_sample_num) {
  ll ran = uid(gen);
  vertex_t lo = 0, hi = n_vertex - 1;
  while (lo < hi) {
    int mid = (lo + hi) >> 1;
    if (sum_wedges[mid] < ran) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  std::random_device rdev_wedge_low;
  std::mt19937 eng_wedge(rdev_wedge_low());
  std::uniform_int_distribution<vertex_t> uid1(0, int(adj[lo].size()) - 1);
  vertex_t l1 = uid1(eng_wedge);
  std::uniform_int_distribution<vertex_t> uid2(0, int(adj[lo].size()) - 2);
  vertex_t l2 = uid2(eng_wedge);

  if (l2 >= l1) l2++;
  l1 = adj[lo][l1];
  l2 = adj[lo][l2];
  return fast_neighbor_intersection(l1, l2, edge_sample_num) - 1; // -1 ?
}

ld BiGraph::bfc_fast_edge_sampling(UID<ll> &uid, mt19937_64 &gen, int vertex_sample_num) {
  auto random_edge = uid(gen);
  vertex_t u = edges[random_edge].u;
  vertex_t v = edges[random_edge].v;
  if (sum_deg_neighbors[u] > sum_deg_neighbors[v]) {
    std::swap(u, v);
  }
  if ((ll) adj[v].size() + sum_deg_neighbors[u] * 2 > vertex_sample_num) {
    return bfc_per_edge_approx(u, v, vertex_sample_num);
  } else {
    return bfc_per_edge_hashmap(u, v);
  }
}

ld BiGraph::bfc_vertex_sampling(UID<int> &uid, mt19937 &gen) {
  vertex_t random_vertex = uid(gen);
  return bfc_per_vertex(random_vertex);
}

ld BiGraph::bfc_weighted_pair_sampling(UID<ll> &uid, std::mt19937_64 &gen) {
  vertex_t u = weighted_sample_v(gen, uid);
  vertex_t v = weighted_sample_v(gen, uid);
  while (u == v) {
    u = weighted_sample_v(gen, uid);
    v = weighted_sample_v(gen, uid);
  }
  ll du = (ll) adj[u].size();
  ll dv = (ll) adj[v].size();
  return bfc_per_pair_fast(u, v) / du / dv;
}

ld BiGraph::bfc_pair_sampling(UID<vertex_t> &uid, std::mt19937 &gen) {
  vertex_t u = pair_sample_v(gen, uid);
  vertex_t v = pair_sample_v(gen, uid);
  while (u == v) {
    u = pair_sample_v(gen, uid);
    v = pair_sample_v(gen, uid);
  }
  return bfc_per_pair_fast(u, v);;
}



/*========= btc sampling method =========*/
ld BiGraph::btc_triple_sampling(UID<vertex_t> &uid, mt19937 &gen) {
  vertex_t u1 = uid(gen);
  vertex_t u2 = uid(gen);
  vertex_t u3 = uid(gen);
  while (u1 == u2 || u2 == u3 || u1 == u3) {
    u1 = uid(gen);
    u2 = uid(gen);
    u3 = uid(gen);
  }
  return btc_per_triple_fast(u1, u2, u3);
}

ld BiGraph::btc_weighted_triple_sampling(UID<ll> &uid, mt19937_64 &gen) {
  vertex_t u = weighted_sample_v(gen, uid);
  vertex_t v = weighted_sample_v(gen, uid);
  vertex_t w = weighted_sample_v(gen, uid);
  while (u == v || v == w || u == w) {
    u = weighted_sample_v(gen, uid);
    v = weighted_sample_v(gen, uid);
    w = weighted_sample_v(gen, uid);
  }
  return btc_per_triple_fast(u, v, w) / adj[u].size() / adj[v].size() / adj[w].size();
}

ld BiGraph::btc_vertex_sampling(UID<vertex_t> &uid, mt19937 &gen) {
  vertex_t u = uid(gen);
  return btc_per_vertex(u);
}

ld BiGraph::btc_edge_sampling(UID<ll> &uid, mt19937_64 &gen) {
  ll random_edge = uid(gen);
  vertex_t u = edges[random_edge].u;
  vertex_t v = edges[random_edge].v;
  return btc_per_edge(u, v);
}

ld BiGraph::btc_wedge_sampling(URD<double> &urd, mt19937 &gen) {
  vertex_t u = alias_sample_v(gen, urd);
  while (adj[u].size() < 2) {
    u = alias_sample_v(gen, urd);
  }
  vertex_t v1 = 0, v2 = 1;
  if (adj[u].size() > 2) {
    UID<vertex_t> uid(0, static_cast<vertex_t>(adj[u].size() - 1));
    v1 = uid(gen), v2 = uid(gen);
    while (v1 == v2) {
      v1 = uid(gen), v2 = uid(gen);
    }
  }

  v1 = adj[u][v1];
  v2 = adj[u][v2];
  return btc_per_wedge(u, v1, v2);
}

ld BiGraph::btc_wedge_sampling2(UID<ll> &uid, mt19937_64 &gen) {
  ll ran = uid(gen);
  vertex_t lo = 0, hi = n_vertex - 1;
  while (lo < hi) {
    int mid = (lo + hi) >> 1;
    if (sum_wedges[mid] < ran) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  std::random_device rdev_wedge_low;
  std::mt19937 eng_wedge(rdev_wedge_low());
  std::uniform_int_distribution<vertex_t> uid1(0, int(adj[lo].size()) - 1);
  vertex_t l1 = uid1(eng_wedge);
  std::uniform_int_distribution<vertex_t> uid2(0, int(adj[lo].size()) - 2);
  vertex_t l2 = uid2(eng_wedge);

  if (l2 >= l1) l2++;
  l1 = adj[lo][l1];
  l2 = adj[lo][l2];
  return btc_per_wedge(lo, l1, l2); // -1 ?
}


/*========= local butterfly counting =========*/
ld BiGraph::bfc_per_edge(vertex_t u, vertex_t v) {
  ld bfc_per_edge = 0;
  for (const vertex_t &w: adj[u])
    if (w != v) {
      bfc_per_edge += exact_neighbor_intersection(v, w) - 1;
    }
  return bfc_per_edge;
}

ld BiGraph::bfc_per_edge_hashmap(vertex_t u, vertex_t v) {
  for (const vertex_t &u1hop: adj[u])
    if (u1hop != v) {
      for (const vertex_t &u2hop: adj[u1hop])
        if (u2hop != u) {
          ++index_map[u2hop];
        }
    }
  ld bfc_per_edge = 0;
  for (const vertex_t &v1hop: adj[v]) {
    bfc_per_edge += index_map[v1hop];
  }

  // clear index_map
  for (const vertex_t &u1hop: adj[u]) {
    for (const vertex_t &u2hop: adj[u1hop]) {
      index_map[u2hop] = 0;
    }
  }

  return bfc_per_edge;
}

ld BiGraph::bfc_per_edge_approx(vertex_t u, vertex_t v, int iter) {
  static mt19937 gen(std::random_device{}());

  if (adj[u].size() <= 1 || adj[v].size() <= 1) {
    return 0;
  }
  UID<vertex_t> uid_u(0, (int) adj[u].size() - 1);
  UID<vertex_t> uid_v(0, (int) adj[v].size() - 1);
  ld res = 0;
  for (int i = 0; i < iter; ++i) {
    vertex_t x = adj[u][uid_u(gen)];
    vertex_t y = adj[v][uid_v(gen)];
    if (x != v && y != u && std::binary_search(adj[x].begin(), adj[x].end(), y)) {
      ++res;
    }
  }
  res *= (ld(adj[u].size()) * adj[v].size() / iter);
  return res;
}

ld BiGraph::bfc_per_vertex(vertex_t u) {
  ld res = 0;
  for (const vertex_t &u1hop: adj[u]) {
    for (const vertex_t &u2hop: adj[u1hop]) {
      if (u2hop != u) {
        res += index_map[u2hop];
        ++index_map[u2hop];
      }
    }
  }

  // clear index_map
  for (const vertex_t &u1hop: adj[u]) {
    for (const vertex_t &u2hop: adj[u1hop]) {
      index_map[u2hop] = 0;
    }
  }

  return res;
}

ld BiGraph::bfc_per_pair_fast(vertex_t u, vertex_t v) { // check
  for (const vertex_t &u1hop: adj[u]) {
    index_map[u1hop] = 1;
  }
  vertex_t co_neighbors = 0;
  for (const vertex_t &v1hop: adj[v]) {
    co_neighbors += index_map[v1hop];
  }

  // clear index_map
  for (const vertex_t &u1hop: adj[u]) {
    index_map[u1hop] = 0;
  }
  if (co_neighbors < 2) return 0;
  return static_cast<ld>(co_neighbors) * (co_neighbors - 1) / 2;
}


/*========= local bi-triangle counting =========*/
ld BiGraph::btc_per_triple_fast(vertex_t u, vertex_t v, vertex_t w) {
  ll c1 = 0; // |N(u) \cap N(v)|
  ll c2 = 0; // |N(v) \cap N(w)|
  ll c3 = 0; // |N(w) \cap N(u)|
  ll c4 = 0; // |N(u) \cap N(v) \cap N(w)|
  for (const vertex_t &u1hop: adj[u]) {
    index_map[u1hop] = u;
  }
  for (const vertex_t &v1hop: adj[v]) {
    // if (index_map[v1hop] == 1) ++c1;
    c1 += (index_map[v1hop] == u);
    index_map2[v1hop] = v;
  }
  for (const vertex_t &w1hop: adj[w]) {
    bool a = (index_map2[w1hop] == v);
    bool b = (index_map[w1hop] == u);
    if (a) ++c2;
    if (b) ++c3;
    if (a & b) ++c4;
  }


  return static_cast<ld>(c1) * c2 * c3 - static_cast<ld>(c1 + c2 + c3 - 2) * c4;
}

ld BiGraph::btc_per_vertex(vertex_t u) {
  ld res = 0, ex = 0;
  for (const vertex_t &u1: adj[u]) {
    index_map[u1] = 1; // to check if a vertex is a neighbor of u quickly
    for (const vertex_t &u2: adj[u1]) {
      if (u2 == u) continue;
      hm1.inc(u2, 1);
    }
  }
  // hm1[u2]: the number of wedges from u to u2
  hm1.foreach([&](vertex_t u2, auto c) {
      ex += cn2<ld>(c) * static_cast<ld>(adj[u2].size() - 2);
      for (const vertex_t &u3: adj[u2]) {
        hm2.inc(u3, c);
      }
  });
  // hm2[u3]: the number of 3-hop paths (super-wedges) from u to u3
  hm2.foreach([&](vertex_t u3, auto &c) {
      if (index_map[u3] == 1) {
        c -= static_cast<ld>(adj[u3].size() - 1);
      }
      res += cn2<ld>(c);
  });

  for (const vertex_t &v: adj[u]) {
    hm1.clear();
    // hm1[v2]: the number of wedges from v to v2 (not pass u)
    for (const vertex_t &v1: adj[v]) {
      if (v1 == u) continue;
      for (const vertex_t &v2: adj[v1]) {
        if (v2 == v) continue;
        hm1.inc(v2, 1);
      }
    }
    hm1.foreach([&](vertex_t, auto &c) { ex += cn2<ld>(c); });
  }

  // clear hm1, hm2 and index_map
  hm1.clear(), hm2.clear();
  for (const vertex_t &u1: adj[u]) index_map[u1] = 0;
  return res - ex;
}


ld BiGraph::btc_per_vertex_gt(vertex_t u) {
  ld res = 0;
  for (const vertex_t &u1: adj[u]) {
    for (const vertex_t &u2: adj[u1]) {
      if (u2 == u) continue;
      for (const vertex_t &u3: adj[u2]) {
        if (u3 == u1) continue;
        for (const vertex_t &u4: adj[u3]) {
          if (u4 == u2 || u4 == u) continue;
          for (const vertex_t &u5: adj[u4]) {
            if (u5 == u3 || u5 == u1) continue;
            for (const vertex_t &u6: adj[u5]) {
              if (u6 == u) ++res;
            }
          }
        }
      }
    }
  }

  return res;
}

ld BiGraph::btc_per_vertex2(vertex_t u) {
  ld res = 0, ex = 0;
  for (const vertex_t &u1: adj[u]) {
    index_map[u1] = u; // to check if a vertex is a neighbor of u quickly
    for (const vertex_t &u2: adj[u1]) {
      if (u2 == u) continue;
      um1[u2] += 1;
    }
  }
  for (const auto &[u2, c]: um1) {
    ex += cn2<ld>(c) * static_cast<ld>(adj[u2].size() - 2);
    for (const vertex_t &u3: adj[u2]) {
      um2[u3] += c;
    }
  }
  for (const auto &[u3, c]: um2) {
    if (index_map[u3] == u) {
      um2[u3] -= static_cast<ld>(adj[u3].size()) - 1;
    }
    res += cn2<ld>(c);
  }

  for (const auto &v: adj[u]) {
    std::unordered_map<vertex_t, ld> t;
    for (const vertex_t &v1: adj[v]) {
      if (v1 == u) continue;
      for (const vertex_t &v2: adj[v1]) {
        if (v2 == v) continue;
        t[v2] += 1;
      }
    }
    for (const auto &[_, c]: t) {
      ex += cn2<ld>(c);
    }
  }

  // clear hm1, hm2 and index_map
  um1.clear(), um2.clear();
  return res - ex;
}

ld BiGraph::btc_per_edge(vertex_t u, vertex_t l) {
  for (const vertex_t &u1: adj[u]) {
    index_map[u1] |= 1;
    if (u1 == l) continue;
    for (const vertex_t &u2: adj[u1]) {
      if (u2 == u) continue;
      hm1.inc(u2, 1);
    }
  }
  for (const vertex_t &l1: adj[l]) {
    index_map[l1] |= 2;
    if (l1 == u) continue;
    for (const vertex_t &l2: adj[l1]) {
      if (l2 == l) continue;
      hm2.inc(l2, 1);
    }
  }
  ld res = 0, ex = 0;
  hm2.foreach([&](vertex_t l2, auto c) {
      if (index_map[l2] & 1) {
        ex += c * (static_cast<ld>(adj[l2].size()) - 2);
      }
      for (const vertex_t &l3: adj[l2]) {
        hm3.inc(l3, c);
      }
  });
  hm1.for_keys([&](vertex_t u2) {
      if (!hm3.contains(u2)) return;
      if (index_map[u2] & 2)
        hm3[u2] -= static_cast<ld>(adj[u2].size() - 1);
      res += hm1[u2] * hm3[u2];
  });
  hm1.clear(), hm2.clear(), hm3.clear();
  for (const vertex_t &u1: adj[u]) index_map[u1] = 0;
  for (const vertex_t &l1: adj[l]) index_map[l1] = 0;
  return res - ex;
}


ld BiGraph::btc_per_wedge(vertex_t u, vertex_t x, vertex_t y) {
  if (adj[x].size() < adj[y].size()) std::swap(x, y);
  for (const vertex_t &x1: adj[x]) {
    if (x1 == u) continue;
    index_map[x1] = 1;
    for (const vertex_t &x2: adj[x1]) {
      if (x2 == x || x2 == y) continue;
      hm1.inc(x2, 1);
    }
  }
  ld res = 0, ex = 0;
  for (const vertex_t &y1: adj[y]) {
    if (y1 == u) continue;
    if (index_map[y1] == 1) ex += static_cast<ld>(adj[y1].size() - 2);
    for (const vertex_t &y2: adj[y1]) {
      if (y2 == y || y2 == x) continue;
      hm2.inc(y2, 1);
    }
  }
  hm2.foreach([&](auto v, auto c) {
      res += hm1[v] * c;
  });
  hm1.clear(), hm2.clear();
  for (const vertex_t &x1: adj[x]) index_map[x1] = 0;
  return res - ex;
}


// local 3-path counting
ld BiGraph::path3_per_edge(vertex_t u, vertex_t v) {
  auto du = adj[u].size();
  auto dv = adj[v].size();
  return (du - 1) * (dv - 1);
}

ld BiGraph::path4_per_wedge(vertex_t u, vertex_t x, vertex_t y) {
  UNUSED_VARIABLE(u);
  auto dx = adj[x].size();
  auto dy = adj[y].size();
  ld co = ordered_neighbor_intersection(x, y);
  return (dx - 1) * (dy - 1) - co;
}

ld BiGraph::path4_per_pair(vertex_t u, vertex_t v) {
  ld co = ordered_neighbor_intersection(u, v);
  ld du = adj[u].size();
  ld dv = adj[v].size();
  return co * ((du - 1) * (dv - 1) - co + 1);
}

ld BiGraph::path3_edge_sampling(UID<ll> &uid, mt19937_64 &gen) {
  ll random_edge = uid(gen);
  vertex_t u = edges[random_edge].u;
  vertex_t v = edges[random_edge].v;
  return path3_per_edge(u, v);
}

ld BiGraph::path4_wedge_sampling(URD<double> &urd, mt19937 &gen) {
  vertex_t u = alias_sample_v(gen, urd);
  while (adj[u].size() < 2) {
    u = alias_sample_v(gen, urd);
  }
  vertex_t v1 = 0, v2 = 1;
  if (adj[u].size() > 2) {
    UID<vertex_t> uid(0, static_cast<vertex_t>(adj[u].size() - 1));
    v1 = uid(gen), v2 = uid(gen);
    while (v1 == v2) {
      v1 = uid(gen), v2 = uid(gen);
    }
  }

  v1 = adj[u][v1];
  v2 = adj[u][v2];
  return path4_per_wedge(u, v1, v2);
}

ld BiGraph::path4_pair_sampling(UID<vertex_t> &uid1, UID<vertex_t> &uid2, mt19937 &gen, double ratio) {
  auto &uid = random_double() < ratio ? uid1 : uid2;
  vertex_t u = uid(gen);
  vertex_t v = uid(gen);
  while (u == v) {
    u = uid(gen);
    v = uid(gen);
  }
  return path4_per_pair(u, v);
}

ld BiGraph::path4_weighted_pair_sampling(UID<ll> &uid, mt19937_64 &gen) {
  vertex_t u = weighted_sample_v(gen, uid);
  vertex_t v = weighted_sample_v(gen, uid);
  while (u == v) {
    u = weighted_sample_v(gen, uid);
    v = weighted_sample_v(gen, uid);
  }
  auto du = adj[u].size();
  auto dv = adj[v].size();
  return path4_per_pair(u, v) / du / dv;
}


/*========= weighted sampling method =========*/
ld BiGraph::weighted_vertex_setup(Side side) {
  int t = (side == BOTH ? 2 : 1);
  weighted_sample_table.resize(n_edges * t);

  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  ld sum_deg = n_edges * t;
  ll cur = 0;
  for (vertex_t u = l; u < r; ++u) {
    ll deg = static_cast<ll>(adj[u].size());

    ll nxt = cur + deg;
    while (cur < nxt) {
      weighted_sample_table[cur] = u;
      ++cur;
    }
  }
  return sum_deg;
}


ld BiGraph::weighted_pair_setup(Side side) {
  int t = (side == BOTH ? 2 : 1);
  weighted_sample_table.resize(n_edges * t);

  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  ld bias = static_cast<ld>(n_edges) * n_edges * t;

  ll cur = 0;
  for (vertex_t u = l; u < r; ++u) {
    ll deg = static_cast<ll>(adj[u].size());
    bias -= deg * deg;

    ll nxt = cur + deg;
    while (cur < nxt) {
      weighted_sample_table[cur] = u;
      ++cur;
    }
  }
  bias = bias / 2;
  return bias;
}

ld BiGraph::weighted_pair_bias(Side side) {
  int t = (side == BOTH ? 2 : 1);

  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  ld bias = static_cast<ld>(n_edges) * n_edges * t;
  for (vertex_t u = l; u < r; ++u) {
    ll deg = static_cast<ll>(adj[u].size());
    bias -= deg * deg;
  }
  bias = bias / 2;
  return bias;
}

vertex_t BiGraph::pair_sample_setup(Side side) {
  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);
  vertex_t n = 0;

  for (vertex_t u = l; u < r; ++u) {
    if (!adj[u].empty()) {
      ++n;
    }
  }

  vertex_table.resize(n);
  int i = 0;
  for (vertex_t u = l; u < r; ++u) {
    if (!adj[u].empty()) {
      vertex_table[i++] = u;
    }
  }

  return n;
}

ld BiGraph::inv_weighted_pair_setup(Side side) {
  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  vector<double> prob(n_vertex);
  double s1 = 0, s2 = 0;

  for (vertex_t u = l; u < r; ++u) {
    if (adj[u].empty()) {
      prob[u] = 0;
      continue;
    }
    auto deg = static_cast<double>(adj[u].size());
    double inv_d = 1.0 / deg;
    prob[u] = inv_d;
    s1 += inv_d;
    s2 += inv_d * inv_d;
  }
  for (vertex_t u = l; u < r; ++u) {
    prob[u] /= s1;
  }

  alias_setup(prob, J, q);

  ld bias = static_cast<ld>(s1) * s1 - s2;
  bias /= 2;
  return bias;
}

ld BiGraph::weighted_triple_setup(Side side) {
  weighted_sample_table.resize(n_edges);

  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  ld bias = static_cast<ld>(n_edges) * n_edges * n_edges;
  ld s1 = 0, s2 = 0;
  ll cur = 0;
  for (vertex_t u = l; u < r; ++u) {
    ll deg = (ll) adj[u].size();

    bias -= deg * (s1 * s1 - s2);
    s1 += deg;
    s2 += deg * deg;

    ll nxt = cur + deg;
    while (cur < nxt) {
      weighted_sample_table[cur] = u;
      ++cur;
    }
  }
  bias /= 2;
  return bias;
}

ld BiGraph::wedge_setup(Side side) {
  vertex_t l = 0, r = 0;
  vertex_range(side, l, r);

  vector<double> prob(n_vertex);
  ld n_wedges = 0;
  for (vertex_t u = l; u < r; ++u) {
    double wedges = cn2(static_cast<double>(adj[u].size()));
    n_wedges += wedges;
    prob[u] = wedges;
  }
  for (vertex_t u = l; u < r; ++u) {
    prob[u] /= static_cast<double>(n_wedges);
  }
  alias_setup(prob, J, q);
  return n_wedges;
}

vertex_t BiGraph::weighted_sample_v(mt19937_64 &gen, UID<ll> &uid) const {
  return weighted_sample_table[uid(gen)];
}

vertex_t BiGraph::alias_sample_v(mt19937 &gen, URD<double> &urd) const {
  return alias_draw(gen, urd, J, q);
}

vertex_t BiGraph::pair_sample_v(mt19937 &gen, UID<vertex_t> &uid) const {
  return vertex_table[uid(gen)];
}