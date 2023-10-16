#include "tracker.h"
#include <random>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include "exact_alg.h"
#include <sstream>
#include "util.h"
#include "gt_manager.h"

namespace tracker {
    string tk;

    void set_tracker(const string &tracker) {
      tk = tracker;
    }

    template<typename F1>
    void time_error_tracker(const F1 &sampler, ld bias, std::ofstream &rout) {
      double total_elapsed_time = 0;
      ld res_sum = 0;
      int sample_num = 1;
      int last_sn = 0;
      int inc_num = 0;
      while (total_elapsed_time <= LOCAL_TIME_LIMIT) {
        int cur_sample_num = sample_num - last_sn;
        clock_t op = clock();
        auto res = sampler();
        while (--cur_sample_num > 0) res += sampler();
        clock_t ed = clock();
        auto elapsed_time = get_elapsed_time(op, ed);

        res_sum += res;
        total_elapsed_time += elapsed_time;

        auto ans = res_sum / sample_num * bias;
        auto error = error_percent(ans);

        rout << std::fixed
             << std::left << std::setw(10) << std::setprecision(6) << total_elapsed_time << " "
             << std::left << std::setw(10) << std::setprecision(0) << sample_num << " "
             << std::left << std::setw(10) << std::setprecision(6) << error << " "
             << std::left << std::setw(30) << std::setprecision(6) << ans << std::endl;

        last_sn = sample_num;
        if (!inc_num) {
          if (elapsed_time > 0.0005) inc_num = sample_num;
          sample_num *= 2;
        } else sample_num += inc_num;
      }
    }

    template<typename F1>
    void variance_tracker(const F1 &sampler, ld bias, std::ofstream &rout) {
      double total_elapsed_time = 0;
      int sample_num = 1;
      int last_sn = 0;
      int inc_num = 0;
      ld s_var = 0;
      while (sample_num <= 10000 || total_elapsed_time <= LOCAL_TIME_LIMIT) {
        int cur_sample_num = sample_num - last_sn;

        clock_t op = clock();
        while (--cur_sample_num > 0) {
          auto res = sampler();
          auto std_dev = bias_with_gt(res * bias);
          s_var += std_dev * std_dev;
        }
        clock_t ed = clock();
        auto elapsed_time = get_elapsed_time(op, ed);
        total_elapsed_time += elapsed_time;

        rout << std::fixed
             << std::left << std::setw(10) << std::setprecision(6) << total_elapsed_time << " "
             << std::left << std::setw(10) << std::setprecision(0) << sample_num << " "
             << std::left << std::setw(30) << std::setprecision(6) << s_var / sample_num << std::endl;

        last_sn = sample_num;
        if (!inc_num) {
          if (elapsed_time > 0.05) inc_num = sample_num;
          sample_num *= 2;
        } else sample_num += inc_num;
      }
    }


    template<typename F1>
    void sampling_tracker(const F1 &sampler, ld bias, std::ofstream &rout) {
      if (tk == "t-e") {
        time_error_tracker(sampler, bias, rout);
      } else if (tk == "var") {
        variance_tracker(sampler, bias, rout);
      }
    }

    void path4_exact_algorithm_tracker(BiGraph &graph, Side side, const string &gt_path) {
      clock_t op = clock();
      ld exact_n_bf = calc_4path_opt(graph, side);
      clock_t ed = clock();
      double elapsed_time = get_elapsed_time(op, ed);

      set_exact_n_bf(exact_n_bf);

      std::cout << std::fixed << " exact_n_bf: " << exact_n_bf << std::endl;
      GT::set_gt("4path_" + to_string(side), exact_n_bf, elapsed_time);
      GT::save_gt(gt_path);
    }

    void bfc_exact_algorithm_tracker(BiGraph &graph, const string &bfc_path) {
      clock_t op = clock();
      ld exact_n_bf = bfc_exact(graph);
      clock_t ed = clock();
      double elapsed_time = get_elapsed_time(op, ed);

      set_exact_n_bf(exact_n_bf);

      std::cout << std::fixed << " exact_n_bf: " << exact_n_bf << std::endl;
      GT::set_gt("bfc", exact_n_bf, elapsed_time);
      GT::save_gt(bfc_path);
    }

    void bfc_vertex_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      auto result_path = get_result_path(graph._alias, "bfc_vertex", rnd, side);
      std::ofstream rout(result_path);

      vertex_t l = 0, r = 0;
      graph.vertex_range(side, l, r);

      mt19937 gen(std::random_device{}());
      UID<vertex_t> uid(l, r - 1);

      const auto &sampler = [&]() -> ld {
          return graph.bfc_vertex_sampling(uid, gen);
      };

      ld bias = static_cast<ld>(r - l) / 2.0 / (1 + (side == BOTH));
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void bfc_fast_edge_sampling_tracker(BiGraph &graph, int rnd, int vertex_sample_num) {
      auto result_path = get_result_path(graph._alias, "bfc_fast_edge", rnd);
      std::ofstream rout(result_path);

      mt19937_64 gen(std::random_device{}());
      auto n_edges = graph.n_edges;
      UID<ll> uid(0, n_edges - 1);

      const auto &sampler = [&]() -> ld {
          return graph.bfc_fast_edge_sampling(uid, gen, vertex_sample_num);
      };

      ld bias = static_cast<ld>(n_edges) / 4.0;
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void bfc_wedge_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      auto result_path = get_result_path(graph._alias, "bfc_wedge", rnd, side);
      std::ofstream rout(result_path);

      ld n_wedges = graph.wedge_setup(side);
      mt19937 gen(std::random_device{}());
      URD<double> urd(0, 1);

      const auto &sampler = [&]() -> ld {
          return graph.bfc_wedge_sampling(urd, gen);
      };

      ld bias = n_wedges / 2.0 / (1 + (side == BOTH));
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void bfc_pair_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      auto result_path = get_result_path(graph._alias, "bfc_pair2", rnd, side);
      std::ofstream rout(result_path);

      vertex_t n = graph.pair_sample_setup(side);
      mt19937 gen(std::random_device{}());
      UID<vertex_t> uid(0, n);

      const auto &sampler = [&]() -> ld {
          return graph.bfc_pair_sampling(uid, gen);
      };

      ld bias = cn2<ld>(n);
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void bfc_weighted_pair_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      auto result_path = get_result_path(graph._alias, "bfc_weighted_pair", rnd, side);
      std::ofstream rout(result_path);

      ld bias = graph.weighted_pair_setup(side);
      std::mt19937_64 gen(std::random_device{}());
      UID<ll> uid(0, graph.n_edges - 1);

      const auto &sampler = [&]() -> ld {
          return graph.bfc_weighted_pair_sampling(uid, gen);
      };

      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

// bi-triangle count time tracker
    void btc_vertex_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      string result_path = get_result_path(graph._alias, "btc_vertex", rnd, side);
      std::ofstream rout(result_path);

      vertex_t l = 0, r = 0;
      graph.vertex_range(side, l, r);
      ld bias = static_cast<ld>(r - l) / 3.0 / (1 + (side == BOTH));
      mt19937 gen(std::random_device{}());
      UID<vertex_t> uid(l, r - 1);

      const auto &sampler = [&]() -> ld {
          return graph.btc_vertex_sampling(uid, gen);
      };
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void btc_edge_sampling_tracker(BiGraph &graph, int rnd) {
      string result_path = get_result_path(graph._alias, "btc_edge", rnd);
      std::ofstream rout(result_path);

      ll n_edges = graph.n_edges;
      mt19937_64 gen(std::random_device{}());
      UID<ll> uid(0, n_edges - 1);

      const auto &sampler = [&]() -> ld {
          return graph.btc_edge_sampling(uid, gen);
      };

      ld bias = static_cast<ld>(n_edges) / 6.0;
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void btc_wedge_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      string result_path = get_result_path(graph._alias, "btc_wedge", rnd, side);
      std::ofstream rout(result_path);

      ld n_wedges = graph.wedge_setup(side);
      mt19937 gen(std::random_device{}());
      URD<double> urd(0, 1);

      const auto &sampler = [&]() -> ld {
          return graph.btc_wedge_sampling(urd, gen);
      };

      ld bias = n_wedges / 3.0 / (1 + (side == BOTH));
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void btc_triple_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      string result_path = get_result_path(graph._alias, "btc_triple", rnd, side);
      std::ofstream rout(result_path);

      vertex_t l = 0, r = 0;
      graph.vertex_range(side, l, r);
      std::mt19937 gen(std::random_device{}());
      UID<vertex_t> uid(l, r - 1);
      ld bias = static_cast<ld>(r - l) * (r - l - 1) * (r - l - 2) / 6.0;

      const auto &sampler = [&]() -> ld {
          return graph.btc_triple_sampling(uid, gen);
      };

      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void btc_weighted_triple_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      string result_path = get_result_path(graph._alias, "btc_weighted_triple", rnd, side);
      std::ofstream rout(result_path);

      ld bias = graph.weighted_triple_setup(side) / 2;
      std::mt19937_64 gen(std::random_device{}());
      UID<ll> uid(0, graph.n_edges - 1);

      const auto &sampler = [&]() -> ld {
          return graph.btc_weighted_triple_sampling(uid, gen);
      };

      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void path3_edge_sampling_tracker(BiGraph &graph, int rnd) {
      auto result_path = get_result_path(graph._alias, "path3_edge", rnd);
      std::ofstream rout(result_path);

      mt19937_64 gen(std::random_device{}());
      auto n_edges = graph.n_edges;
      UID<ll> uid(0, n_edges - 1);

      const auto &sampler = [&]() -> ld {
          return graph.path3_edge_sampling(uid, gen);
      };

      auto bias = static_cast<ld>(n_edges);
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void path4_wedge_sampling_tracker(BiGraph &graph, int rnd) {
      auto result_path = get_result_path(graph._alias, "path4_wedge", rnd);
      std::ofstream rout(result_path);

      ld n_wedges = graph.wedge_setup(BOTH);
      mt19937 gen(std::random_device{}());
      URD<double> urd(0, 1);

      const auto &sampler = [&]() -> ld {
          return graph.path4_wedge_sampling(urd, gen);
      };

      auto bias = static_cast<ld>(n_wedges);
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void path4_pair_sampling_tracker(BiGraph &graph, int rnd) {
      auto result_path = get_result_path(graph._alias, "path4_pair", rnd);
      std::ofstream rout(result_path);

      mt19937 gen(std::random_device{}());
      UID<vertex_t> uid1(0, graph.l_vertex_num - 1);
      UID<vertex_t> uid2(graph.l_vertex_num, graph.n_vertex - 1);
      auto n_pair_l = cn2<double>(graph.l_vertex_num);
      auto n_pair_r = cn2<double>(graph.r_vertex_num);
      auto ratio = n_pair_l / (n_pair_l + n_pair_r);

      const auto &sampler = [&]() -> ld {
          return graph.path4_pair_sampling(uid1, uid2, gen, ratio);
      };

      ld bias = n_pair_l + n_pair_r;
      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void path4_weighted_pair_sampling_tracker_b(BiGraph &graph, int rnd) {
      auto result_path = get_result_path(graph._alias, "path4_weighted_pair", rnd, BOTH);
      std::ofstream rout(result_path);

      graph.weighted_pair_setup(BOTH);
      auto bias_l = graph.weighted_pair_bias(L_SIDE);
      auto bias_r = graph.weighted_pair_bias(R_SIDE);

      auto n_edges = graph.n_edges;
      std::mt19937_64 gen(std::random_device{}());
      UID<ll> uid_l(0, n_edges - 1);
      UID<ll> uid_r(n_edges, n_edges * 2 - 1);
      const auto &sampler_l = [&]() -> ld {
          return graph.path4_weighted_pair_sampling(uid_l, gen);
      };
      const auto &sampler_r = [&]() -> ld {
          return graph.path4_weighted_pair_sampling(uid_r, gen);
      };

      double total_elapsed_time = 0;
      ld res_sl = 0, res_sr = 0;
      int sample_num = 1, last_sn = 0, inc_num = 0;
      while (total_elapsed_time <= LOCAL_TIME_LIMIT) {
        int cur_sample_num = sample_num - last_sn;
        clock_t op = clock();
        auto res_l = sampler_l();
        auto res_r = sampler_r();
        for (int i = 1; i < cur_sample_num; ++i) {
          res_l += sampler_l();
          res_r += sampler_r();
        }
        clock_t ed = clock();
        auto elapsed_time = get_elapsed_time(op, ed);

        res_sl += res_l;
        res_sr += res_r;
        total_elapsed_time += elapsed_time;

        auto ans = res_sl / sample_num * bias_l + res_sr / sample_num * bias_r;
        auto error = error_percent(ans);

        rout << std::fixed
             << std::left << std::setw(10) << std::setprecision(6) << total_elapsed_time << " "
             << std::left << std::setw(10) << std::setprecision(0) << sample_num << " "
             << std::left << std::setw(10) << std::setprecision(6) << error << " "
             << std::left << std::setw(30) << std::setprecision(6) << ans << std::endl;

        last_sn = sample_num;
        if (!inc_num) {
          if (elapsed_time > 0.0005) inc_num = sample_num;
          sample_num *= 2;
        } else sample_num += inc_num;
      }

      rout.close();
    }

    void path4_weighted_pair_sampling_tracker(BiGraph &graph, int rnd, Side side) {
      if (side == BOTH) {
        return path4_weighted_pair_sampling_tracker_b(graph, rnd);
      }
      auto result_path = get_result_path(graph._alias, "path4_weighted_pair", rnd, side);
      std::ofstream rout(result_path);

      auto bias = graph.weighted_pair_setup(side);
      auto n_edges = graph.n_edges;
      std::mt19937_64 gen(std::random_device{}());
      UID<ll> uid(0, n_edges - 1);
      const auto &sampler = [&]() -> ld {
          return graph.path4_weighted_pair_sampling(uid, gen);
      };

      sampling_tracker(sampler, bias, rout);

      rout.close();
    }

    void btc_exact_algorithm_tracker(BiGraph &graph, const string &btc_path) {
      clock_t op = clock();
      ld exact_n_bf = btc_exact_swja(graph);
      clock_t ed = clock();
      double elapsed_time = get_elapsed_time(op, ed);

      set_exact_n_bf(exact_n_bf);

      std::cout << std::fixed << " exact_n_bf: " << exact_n_bf << std::endl;
      GT::set_gt("btc", exact_n_bf, elapsed_time);
      GT::save_gt(btc_path);
    }

}