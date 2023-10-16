#include "bi_graph.h"
#include "util.h"
#include "tracker.h"
#include "gt_manager.h"
#include <cstring>

struct Params {
    Params() = default;

    Params(int argc, char *args[]) {
      int cnt = 0;
      while (cnt < argc) {
        if (strcmp(args[cnt], "-alias") == 0) {
          alias = args[++cnt];
        } else if (strcmp(args[cnt], "-alg") == 0) {
          alg = std::stoi(args[++cnt]);
        } else if (strcmp(args[cnt], "-rnd") == 0) {
          rnd = std::stoi(args[++cnt]);
        } else if (strcmp(args[cnt], "-cleaned") == 0) {
          cleaned = 1;
        } else if (strcmp(args[cnt], "-tracker") == 0) {
          tracker = args[++cnt];
        } else if (strcmp(args[cnt], "-side") == 0) {
          const string s = args[++cnt];
          if (s == "l") side = L_SIDE;
          else if (s == "r") side = R_SIDE;
          else if (s == "b") side = BOTH;
          else {
            std::cerr << " Error side: " << s << "!\n";
            exit(0);
          }
        } else if (strcmp(args[cnt], "-me") == 0) {
          me = args[++cnt];
        }
        cnt++;
      }
    }

    int alg{-1}, rnd{RND_NUM}, cleaned{0};
    Side side{L_SIDE};
    std::string alias{"dbpedia"}, tracker{"t-e"};
    std::string me{"bfc"};
};

void init_gt(const string &key, const string &alias) {
  GT::read_gt(get_gt_path(alias, "gt."));
  GT::show_gts(std::cout);

  ld exact_n_bf = GT::get_gt(key);
  if (exact_n_bf > 0) {
    set_exact_n_bf(exact_n_bf);
    std::cerr << " [" << alias << "-" << key << "] ground truth: " << exact_n_bf << std::endl;
  } else {
    std::cerr << " Not prepare ground truth for " << " [" << alias << "-" << key << "]!" << std::endl;
  }
}

void test() {}

void bfc(BiGraph &graph, int alg, int rnd, Side side) {
  init_gt("bfc", graph._alias);

  if (alg == 0) {
    auto gt_path = get_gt_path(graph._alias, "gt.");
    tracker::bfc_exact_algorithm_tracker(graph, gt_path);
  } else if (alg == 1) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::bfc_vertex_sampling_tracker(graph, i, side);
    }
  } else if (alg == 2) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::bfc_fast_edge_sampling_tracker(graph, i, N_FAST_EDGE_BFC_ITERATIONS);
    }
  } else if (alg == 3) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::bfc_wedge_sampling_tracker(graph, i, side);
    }
  } else if (alg == 4) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::bfc_pair_sampling_tracker(graph, i, side);
    }
  } else if (alg == 5) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::bfc_weighted_pair_sampling_tracker(graph, i, side);
    }
  }
}

void btc(BiGraph &graph, int alg, int rnd, Side side) {
  init_gt("btc", graph._alias);

  if (alg == 0) {
    auto gt_path = get_gt_path(graph._alias, "gt.");
    tracker::btc_exact_algorithm_tracker(graph, gt_path);
  } else if (alg == 1) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::btc_vertex_sampling_tracker(graph, i, side);
    }
  } else if (alg == 2) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::btc_edge_sampling_tracker(graph, i);
    }
  } else if (alg == 3) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::btc_wedge_sampling_tracker(graph, i, side);
    }
  } else if (alg == 4) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::btc_triple_sampling_tracker(graph, i, side);
    }
  } else if (alg == 5) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::btc_weighted_triple_sampling_tracker(graph, i, side);
    }
  }
}

void path3(BiGraph &graph, int alg, int rnd, Side side) {
  UNUSED_VARIABLE(alg);
  UNUSED_VARIABLE(side);

  init_gt("3path", graph._alias);
  for (int i = 1; i <= rnd; ++i) {
    tracker::path3_edge_sampling_tracker(graph, i);
  }
}

void path4(BiGraph &graph, int alg, int rnd, Side side) {
  init_gt("4path_" + to_string(side), graph._alias);

  if (alg == 0) {
    auto gt_path = get_gt_path(graph._alias, "gt.");
    tracker::path4_exact_algorithm_tracker(graph, side, gt_path);
  } else if (alg == 3) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::path4_wedge_sampling_tracker(graph, i);
    }
  } else if (alg == 4) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::path4_pair_sampling_tracker(graph, i);
    }
  } else if (alg == 5) {
    for (int i = 1; i <= rnd; ++i) {
      tracker::path4_weighted_pair_sampling_tracker(graph, i, side);
    }
  }

}


int main(int argc, char *args[]) {
  std::ios::sync_with_stdio(false);

  Params params(argc, args);
  int alg = params.alg;
  int rnd = params.rnd;
  int cleaned = params.cleaned;
  Side side = params.side;
  string alias = params.alias;
  auto me = params.me;

  tracker::set_tracker(params.tracker);

  BiGraph graph(alias, cleaned);
  if (me == "clean") {
    graph.write_cleaned(get_dataset_path(alias, true));
  } else if (me == "test") {
    test();
  } else if (me == "bfc") {
    bfc(graph, alg, rnd, side);
  } else if (me == "btc") {
    btc(graph, alg, rnd, side);
  } else if (me == "path3") {
    path3(graph, alg, rnd, side);
  } else if (me == "path4") {
    path4(graph, alg, rnd, side);
  } else {
    std::cerr << " error me: [" << me << "]!" << std::endl;
  }

  std::cerr << " Take a look at the output file ..." << std::endl;
  return 0;
}