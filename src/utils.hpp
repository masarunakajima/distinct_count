#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#ifndef UTILS
#define UTILS

#include "ss_names.hpp"
#include <boost/multiprecision/cpp_int.hpp>

#include <unistd.h>
#include <thread>
#include <string>
#include <vector>
#include <stack>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <limits.h>
#include <map>
#include <sstream>
#include <cstdlib>
#include <algorithm>

//#include "nlohmann/json.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "constants.hpp"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct ss_node{
  bool paired;
  bool nick = false;

  unsigned start;
  unsigned end;
  unsigned u; // number of unpaired bases
  unsigned p; // number of base pairs
  std::vector<struct ss_node *> children;
  std::vector<unsigned> pchild_indices;

};


void
parse_pfunc_input_file(const std::string &input_file,
                       unsigned &num_strands,
                       std::vector<std::vector<char> > &strands,
                       std::vector <unsigned> &strand_ids);

void
parse_energy_input_file(const std::string &input_file,
                        unsigned int &num_strands,
                        std::vector<std::vector<char> > &strands,
                        std::vector <unsigned> &strand_ids,
                        std::string & ss);


static inline char base_to_num(const char b){
  if (b == 'A') return base_A;
  else if (b == 'C') return base_C;
  else if (b == 'G') return base_G;
  else if (b == 'T') return base_U;
  else if (b == 'U') return base_T;
  else if (b == 'D') return base_D;
  else{
    printf("Invalid base \"%c\".\n", b);
    return -1;
  }
  return -1;
}

static inline char num_to_base(const unsigned n){
  if (n == 0) return 'A';
  else if (n == 1) return 'C';
  else if (n == 2) return 'G';
  else if (n == 3) return 'U';
  else if (n == 6) return 'D';
  else if (n == 9) return 'X';
  else{
    printf("Invalid number \"%d\".\n", n);
    return ' ';
  }
  return ' ';
}

static inline
unsigned get_pos_pair_index(const unsigned pos1, const unsigned pos2) {
  return pos1 + max_seq_size*pos2;
}

inline std::string get_directory(const std::string &fpath){
  std::string dir = fpath;
  for (unsigned i = 0; i < fpath.size(); i++){
    if (dir.back() != '/') dir.pop_back();
    else {
      dir.pop_back();
      break;
    }

  }
  return dir;
}


bool is_dangle5(const unsigned i, const ss_node &parent);
bool is_dangle3(const unsigned i, const ss_node &parent);
bool is_dangle5_closing_pair(const ss_node &pnode);
bool is_dangle3_closing_pair(const ss_node &pnode);
void build_ss_tree(const std::string &ss, ss_node *node,
                   std::stack<ss_node*> &bp_stack,
                   unsigned &num_nicks);
void walk_ss_tree(const ss_node *node);


bool is_hairpin(const ss_node &node);
bool is_exterior(const ss_node &node);
bool is_interior(const ss_node &node);
bool is_multi(const ss_node &node);

bool is_polyC(const std::vector<char> &seq, const unsigned i, const unsigned j);
bool can_pair(const std::vector<char> &seq, const unsigned i,
              const unsigned j);
bool is_nonGC(const std::vector<char> &seq, const unsigned i, const unsigned j);


std::string GetCurrentWorkingDir();

std::string
convertToString(char* a, int size);


inline
bool
bases_can_pair(const char a, const char b) {
  return (a + b == 3) || (a + b == 5);
}


inline
bool
rna_bases_can_pair(const char a, const char b) {
  return (a + b == 3) || (a + b == 5);
}


inline
bool
bases_can_pair(const std::vector<char> &seq,
               const unsigned i, const unsigned j) {
  if (j <= i + hairpin_property::min_loop_size)
    return false;
  else
    return (seq[i] + seq[j] == 3) ||
      (seq[i] + seq[j] == 5);
}



inline
bool
bases_can_pair(const std::vector<char> &seq,
               const std::vector<std::vector<bool> > &eta,
               const unsigned i, const unsigned j) {

  if ((eta[i][j]) && (j <= i + hairpin_property::min_loop_size))
    return false;
  else
    return (seq[i] + seq[j] == 3) ||
      (seq[i] + seq[j] == 5);
}


void
get_divisors(const unsigned n, std::vector<unsigned> &sym);


void
get_orb_sizes(const std::vector<unsigned> &strand_ids, 
    unsigned &sym,
    std::vector<unsigned> &orb_sizes);


bool
is_valid(const std::vector<char> &seq,
         const std::vector<std::vector<bool> > &eta,
         const std::vector<std::pair<unsigned, unsigned> > & s,
         const std::string &s_str);

uint
get_sym_num (const std::vector<char> &seq,
             const std::vector<std::vector<bool> > &eta,
             const std::vector<std::pair<unsigned, unsigned> > & s,
             const std::string &s_str,
             const unsigned sym,
             const std::vector<unsigned> &orb_sizes);


void
get_pair_matrix(const std::vector<char> &seq,
                const std::vector<std::vector<bool> > &eta,
                std::vector<std::vector<unsigned> > & pair_matrix);


void
get_eta(const std::vector<char> & seq,
        const std::vector<unsigned> & nicks,
        std::vector<std::vector<bool> > & eta);

void
get_p(const std::vector<char> & seq,
    const std::vector<unsigned> & nicks,
    const std::vector<std::vector<bool> > & eta,
    std::vector<std::vector<bool> > & p,
    std::vector<std::vector<bool> > & p_red);
    

void
get_seq(const  std::vector<std::vector<char> > &strands,
        const std::vector<unsigned> &strand_ids,
        std::vector<char> &seq,
        std::vector<unsigned> &nicks);

void
get_seq(const  std::vector<std::vector<char> > &strands,
        const std::vector<unsigned> &strand_ids,
        std::vector<char> &seq,
        std::vector<unsigned> &nicks,
        std::vector<bool> &gam);


void
get_base_ss(const std::vector<unsigned> &nicks,
            const unsigned size,
    const std::vector<unsigned> &displacement,
            std::string &base_s);



void
ss2string (const std::vector<std::pair<unsigned, unsigned> > &pairs,
    const std::vector<unsigned> &displacement,
           std::string &s);





// // types used in secondary structure analysis
// namespace ss_names{
//   typedef unsigned long long ull;
//   typedef char base;
//   typedef std::pair<unsigned, unsigned> index_pair;
//   typedef std::pair<unsigned, unsigned> pos_pair;
//   typedef std::pair<base, base> base_pair;
//   typedef std::vector<std::vector<double> > matrix;
//   typedef std::vector<std::vector<ull> > s_matrix;
//   typedef std::vector<std::vector<std::vector<double> > > vec_3d;
//   typedef std::vector<std::vector<std::vector<size_t> > > s_vec_3d;
//   typedef std::vector<std::vector<std::reference_wrapper<matrix> > > matrix_block;
//   typedef std::vector<std::reference_wrapper<matrix> > matrix_vector;
//   typedef unsigned pos;
//   typedef std::unordered_map<unsigned, double> datamap;
//   typedef unsigned int uint;

// }





#endif
