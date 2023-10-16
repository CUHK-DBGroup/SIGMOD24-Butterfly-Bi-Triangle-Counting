#include "gt_manager.h"
#include <iostream>
#include <sstream>
#include <map>
#include <fstream>
#include <iomanip>

namespace GT {
    struct Item {
        long double gt;
        double time;
    };
    std::map<std::string, Item> m;

    void show_gts(std::ostream &out) {
      for (const auto &[k, v] : m) {
        out << std::fixed
            << std::left << std::setw(10) << k << " "
            << std::left << std::setw(30) << std::setprecision(6) << v.gt << " "
            << std::setprecision(6) << v.time << "\n";
      }
    }

    void read_gt(const std::string &path) {
      std::ifstream fin(path);
      if (!fin.is_open()) {
        std::cerr << " error when open [" << path << "]" << std::endl;
        return;
      }
      std::stringstream ss;
      std::string buf;
      std::string name;
      long double gt;
      double time;
      while (!fin.eof()) {
        std::getline(fin, buf);
        ss.clear();
        ss << buf;
        if (!(ss >> name)) continue;
        if (!(ss >> gt)) continue;
        if (!(ss >> time)) continue;
        m[name] = {gt, time};
      }
      fin.close();
    }

    void set_gt(const std::string &key, long double gt, double time) {
      m[key] = {gt, time};
    }

    void save_gt(const std::string &path) {
      std::ofstream fout(path);
      if (!fout.is_open()) {
        std::cerr << " error when open [" << path << "]" << std::endl;
        return;
      }
      show_gts(fout);
      fout.close();
    }

    long double get_gt(const std::string &key) {
      if (!m.count(key)) return -1;
      return m[key].gt;
    }
}
