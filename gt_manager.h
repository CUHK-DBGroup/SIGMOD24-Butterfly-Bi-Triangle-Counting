#ifndef GT_MANAGER_H
#define GT_MANAGER_H

#include <string>

namespace GT {

    void show_gts(std::ostream &out);

    void read_gt(const std::string &path);

    void set_gt(const std::string &key, long double gt, double time = 0.0);

    void save_gt(const std::string &path);

    long double get_gt(const std::string &key);
}


#endif //GT_MANAGER_H
