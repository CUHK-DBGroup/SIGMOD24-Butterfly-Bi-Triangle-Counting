#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <iostream>
#include <ctime>
#include <vector>
#include <random>
#include <queue>
#include "def.h"

template<typename key_t, typename value_t>
struct HashMap { // assume that _value[i] never be 0 after the first assignment for all i
    vector<key_t> _key;
    vector<value_t> _value;
    int _top{0};

    inline void resize(size_t new_size) {
      _key.resize(new_size);
      _value.resize(new_size);
    }

    inline void clear() { while (_top > 0) _value[_key[--_top]] = 0; }

    inline void set(key_t k, value_t v) {
      if (_value[k] == 0) _key[_top++] = k;
      _value[k] = v;
    }

    inline void inc(key_t k, value_t v) {
      if (_value[k] == 0)_key[_top++] = k;
      _value[k] += v;
    }

    inline value_t operator[](int i) const { return _value[i]; }

    inline value_t &operator[](int i) { return _value[i]; }

    inline bool contains(key_t k) const { return _value[k] != 0; }

    template<typename F>
    inline void for_keys(const F &f) const {
      for (int i = 0; i < _top; ++i) {
        f(_key[i]);
      }
    }

    template<typename F>
    inline void foreach(const F &f) {
      for (int i = 0; i < _top; ++i) {
        f(_key[i], _value[_key[i]]);
      }
    }
};

enum Side {
    NO_CARE, L_SIDE, R_SIDE, BOTH
};

std::string to_string(Side side);

std::string get_dataset_path(const string &alias, bool cleaned);

std::string get_result_path(const std::string &alias, const std::string &alg, Side side = NO_CARE);

std::string get_result_path(const std::string &alias, const std::string &alg, int rnd, Side side = NO_CARE);

std::string get_gt_path(const std::string &alias, const std::string &prefix);


int to_int(const char *s, int n, int &res);

// check if all the chars in the string are numbers
bool all_num(const std::string &s);

bool file_exist(const std::string &tmp_file);

ld read_gt(const std::string &gt_path);

void set_exact_n_bf(ld exact_n_bf);

ld error_percent(ld res, bool to_abs = true);

ld bias_with_gt(ld res, bool to_abs = true);

double get_elapsed_time(clock_t op, clock_t ed);

inline double random_double() {
  static URD<double> urd(0, 1);
  static mt19937 gen(std::random_device{}());
  return urd(gen);
}

template<typename T>
inline T cn2(T n) {
  return n * (n - 1) / 2;
}

template<typename T>
inline T cnk(T n, int k) {
  if (n < k) return 0;
  T res = 1;
  for (int i = 0; i < k; ++i) {
    res *= n - i;
  }
  for (int i = 1; i <= k; ++i) {
    res /= i;
  }
  return res;
}

template<typename T>
void alias_setup(const std::vector<double> &prob, std::vector<T> &J, std::vector<double> &q) {
  auto k = (T) prob.size();
  J.resize(k);
  q.resize(k);

  std::queue<T> smaller, larger;
  for (T i = 0; i < k; ++i) {
    q[i] = k * prob[i];
    if (q[i] < 1.0) smaller.push(i);
    else larger.push(i);
  }

  while (!smaller.empty() && !larger.empty()) {
    auto si = smaller.front();
    auto li = larger.front();
    smaller.pop(), larger.pop();
    J[si] = li;

    q[li] = q[li] + q[si] - 1.0;
    if (q[li] < 1.0) smaller.push(li);
    else larger.push(li);
  }
}

template<typename T>
T alias_draw(std::mt19937 &gen, std::uniform_real_distribution<double> &urd, const std::vector<T> &J,
             const std::vector<double> &q) {
  auto k = (size_t) (J.size() * urd(gen));
  return (urd(gen) < q[k] ? k : J[k]);
}

#endif // UTIL_H