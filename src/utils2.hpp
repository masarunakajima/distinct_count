#ifndef UTILS
#define UTILS


#include <map>
#include <string>
#include <vector>




typedef char base;
typedef char bnum;

const bnum base_A = 0;
const bnum base_C = 1;
const bnum base_G = 2;
const bnum base_U = 3;
const bnum base_T = 3;
const bnum base_D = 6;
const bnum base_K = 9; // unknown base

extern std::map<base, bnum> base_bnum_dict;

bool
is_valid_seq (const std::string &bases);

void 
base_to_num (const std::string &bases, std::vector<bnum> num_seq);


#endif
