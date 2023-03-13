#ifndef CONSTANTS
#define CONSTANTS

#include "ss_names.hpp"


#include <string>
#include <utility>
#include <limits>
#include <math.h>

const unsigned max_unpaired_bases = 30;
const double loop_coef = 1.75;
// const unsigned alphabet_size = 4;
const double multi_const_a = 4.6;
const double multi_const_b = 0.4;
const double multi_const_c = 0.1;
const unsigned asym_max_leg_size = 4;
const int asym_max_leg_size_int = 4;
const double bimol = 161.8649;
// const double asymmetric_inc = 0.3;
// const double asymmetric_max = 3.0;


const char base_A = 0;
const char base_C = 1;
const char base_G = 2;
const char base_U = 3;
const char base_T = 3;
const char base_D = 6;
const char base_K = 9;

const bool GAIL = true;

const double bulge_slope = 7.0496;
const double bulge_interc = 478.17;

const double interior_slope = 7.2;
const double interior_interc = 546.31;

const unsigned undef_base = 9;
const std::pair<unsigned, unsigned> undef_base_pair = {9, 9};

const unsigned max_seq_size = 10000;

const double gas_const = 0.19872; // deca cal/mol/K
const double freezing_k = 273.15; // freezing point in kelvin.
const double ref_temp = 310.15; // reference temperature.

const std::string hloop_file = "hloop";
const std::string bloop_file = "bloop";
const std::string iloop_file = "iloop";
const std::string interior1x1_file = "interior1x1";
const std::string interior2x2_file = "interior2x2";
const std::string interior1x2_file = "interior1x2";
const std::string dangle3_file = "dangle3";
const std::string dangle5_file = "dangle5";
const std::string stack_file = "stack";
const std::string hmismatch_file = "hmismatch";
const std::string imismatch_file = "imismatch";
const std::string tetraloop_file = "tetraloop";
const std::string triloop_file = "triloop";
const std::string asym_file = "asym";
const std::string multicoef_file = "multicoef";
const std::string nonGC_file = "non_GC_penalty";
const std::string polyC_file = "polyC";

const std::string json_ext = ".json";
const std::string json_dG = "_dG.json";
const std::string json_dH = "_dH.json";

namespace hairpin_property
{
  extern unsigned size_3;
  extern unsigned size_4;
  extern unsigned max_unpaired_bases;
  extern unsigned min_loop_size;
  extern int min_loop_size_int;
}

namespace interior_property
{
  extern unsigned size_0;
  extern unsigned size_1;
  extern unsigned size_2;
  extern unsigned size_3;
  extern unsigned max_unpaired_bases;
  extern unsigned min_loop_size;
}

namespace loop_property
{
  extern unsigned size_0;
  extern unsigned size_1;
}

namespace multi_ptoperty
{
  extern unsigned size_0;
  extern unsigned size_1;
  extern unsigned size_8;
}


const unsigned filename_max = 1000;
const unsigned path_max = 1000;


#endif
