
#ifndef SS_NAMES
#define SS_NAMES
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
#include <unordered_map>
#include <boost/multiprecision/cpp_int.hpp>




// types used in secondary structure analysis
namespace ss_names{
  typedef boost::multiprecision::uint1024_t ull;
  typedef char base;
  typedef std::pair<unsigned, unsigned> index_pair;
  typedef std::pair<unsigned, unsigned> pos_pair;
  typedef std::pair<base, base> base_pair;
  typedef std::vector<std::vector<double> > matrix;
  typedef std::vector<std::vector<ull> > s_matrix;
  typedef std::vector<std::vector<std::vector<double> > > vec_3d;
  typedef std::vector<std::vector<std::vector<size_t> > > s_vec_3d;
  typedef std::vector<std::vector<std::reference_wrapper<matrix> > > matrix_block;
  typedef std::vector<std::reference_wrapper<matrix> > matrix_vector;
  typedef unsigned pos;
  typedef std::unordered_map<unsigned, double> datamap;
  typedef unsigned int uint;

}





#endif
