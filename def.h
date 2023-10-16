#ifndef DEF
#define DEF

#include <utility>
#include <random>

#define UNUSED_VARIABLE(x) ((void)(x))
using std::vector;
using std::string;
using std::mt19937;
using std::mt19937_64;


// typedef
using ll = long long;
using ull = unsigned long long;
using ld = long double;
using us = unsigned short;

template<typename T>
using UID = std::uniform_int_distribution<T>;
template<typename T>
using URD = std::uniform_real_distribution<T>;

// constexpr def
constexpr double LOCAL_TIME_LIMIT = 40;
constexpr int LOCAL_ITERATIONS_NUM = 1;
constexpr int N_FAST_EDGE_BFC_ITERATIONS = 2100;

constexpr int SPARS_ITERATIONS_NUM = 1;
constexpr double SPARS_TIME_LIMIT = 16;
constexpr double SPARS_PROB_OP = 1e-6;
constexpr double SPARS_PROB_MAX = 0.00500;
constexpr int SPARS_PROB_OUT_BIAS = 1e6;

constexpr int RND_NUM = 50;
constexpr int RND_THRESHOLD = 2;
constexpr int TIME_SAMPLE_NUM = 1000;
constexpr double ERROR_THRESHOLD = 1;
constexpr double TIME_THERSHOLD = 100;


const std::string GT_ROOT = R"(GOUND_TRUTH_PATH/)";
const std::string DATASET_ROOT = R"(DATASET_PATH/)";
const std::string RESULT_ROOT = R"(RESULT_PATH/)";

const std::string BFC_PREFIX = "bfc_";
const std::string BTC_PREFIX = "swja_";


// const std::string DATASET_ROOT = R"(D:\Project\SIGMOD23_Bi-Triangle\kdd_weighted_butterfly\data)";

#endif // DEF