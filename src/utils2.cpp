#include "utils2.hpp"

using std::map;
using std::string;
using std::vector;

map<base, bnum> base_bnum_dict = {
  {'A', base_A}, {'a', base_A},
  {'U', base_U}, {'u', base_U},
  {'C', base_C}, {'c', base_C},
  {'G', base_G}, {'g', base_G},
  {'U', base_U}, {'u', base_T}
};


bool
is_valid_seq (const string &bases) {

  for (size_t i = 0; i < bases.size(); i++) {
    auto it = base_bnum_dict.find(bases[i]);
    if (it == base_bnum_dict.end()) {
      return false;
    } 
  }
  return true;
}

void 
base_to_num (const string &bases, vector<bnum> num_seq) {

  for (size_t i = 0; i < bases.size(); i++) {
    auto it = base_bnum_dict.find(bases[i]);
    if (it != base_bnum_dict.end()) {
      num_seq.push_back(it->second);
    } else {
      num_seq.push_back(base_K);
    }
  }
}
