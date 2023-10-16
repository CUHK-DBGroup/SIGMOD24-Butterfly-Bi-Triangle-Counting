#include "util.h"
#include <unistd.h>
#include <vector>
#include <queue>
#include <random>
#include <fstream>
#include <iostream>


std::string to_string(Side side) {
  static const string _s2s[4] = {"", "l", "r", "b"};
  return _s2s[side];
}


std::string get_dataset_path(const string &alias, bool cleaned) {
  const std::string tmp = cleaned ? "clean." : "out.";
  std::string res = DATASET_ROOT + alias + "/" + tmp + alias;
  std::cout << " Data path: " << res << std::endl;
  return res;
}

std::string get_result_path(const string &alias, const string &alg, Side side) {
  std::string suffix;
  if (side != NO_CARE) suffix += "_" + to_string(side);

  std::string res = RESULT_ROOT + alias + "/" + alg + "_"  + alias + suffix + ".out";
  std::cout << " Result path: " << res << std::endl;
  return res;
}


std::string get_result_path(const string &alias, const string &alg,  int rnd, Side side) {
    std::string suffix =  "";
    suffix += (rnd < 10 ? "0" : "") + std::to_string(rnd);
    if (side != NO_CARE) suffix += "_" + to_string(side);

    std::string res = RESULT_ROOT + alias + "/" + alg + "." + alias + suffix;
    std::cout << " Result path: " << res << std::endl;
    return res;
}

std::string get_gt_path(const std::string &alias, const std::string &prefix) {
  std::string res = GT_ROOT + alias + "/" + prefix + alias;
  return res;
}

int to_int(const char *s, int n, int &res) {
  if (n <= 0) return -1;
  res = 0;
  for (int i = 0; i < n; ++i) {
    char c = s[i];
    if (c == ' ') return i + 1;
    if (c > '9' || c < '0') return -1;
    res = res * 10 + c - '0';
  }
  return res;
}

bool all_num(const std::string &s) {
  for (char c: s) if (!(c >= '0' && c <= '9')) return false;
  return true;
}

bool file_exist(const std::string &tmp_file) {
  return access(tmp_file.c_str(), 0) == 0;
}

ld read_gt(const std::string &gt_path) {
  ld exact_n_bf = 0;
  if (file_exist(gt_path)) {
    std::ifstream infile(gt_path.c_str());
    infile >> exact_n_bf;
    infile.close();
    return exact_n_bf;
  } else {
    return -1;
  }
}

ld exact_n_bf;

void set_exact_n_bf(ld _exact_n_bf) {
  exact_n_bf = _exact_n_bf;
  std::cout << " set_exact_n_bf: " << exact_n_bf << std::endl;
}

ld error_percent(const ld res, bool to_abs) {
  if (exact_n_bf == 0) return 0;
  ld error = (res - exact_n_bf) / exact_n_bf * 100.0;
  if (to_abs && error < 0) error *= -1.0;
  return error;
}

ld bias_with_gt(const ld res, bool to_abs) {
  if (exact_n_bf == 0) return 0;
  ld error = (res - exact_n_bf);
  if (to_abs && error < 0) error *= -1.0;
  return error;
}

double get_elapsed_time(clock_t op, clock_t ed) {
  return double(ed - op) / CLOCKS_PER_SEC;
}
