#include "utils.hpp"

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;
using std::make_pair;
using std::runtime_error;
using std::to_string;
using std::map;
using std::tuple;
using std::make_tuple;
using std::get;
using std::stringstream;

//using json = nlohmann::json;

using namespace ss_names;

void
get_divisors(const uint n, vector<uint> &sym) {
  for (uint i = 1; i <= n; i++) {
    if (n % i == 0) {
      sym.push_back(i);
    }
  }
}


void
get_orb_sizes(const vector<uint> &strand_ids, 
    unsigned &sym,
    vector<uint> &orb_sizes) {
  uint num_strands = strand_ids.size();
  vector<uint> divisors;
  get_divisors(num_strands, divisors);

  for(uint d : divisors){
    bool same = true;
    for (uint i = 0; i < num_strands; i++) {
      if (strand_ids[i] != strand_ids[(i+d)%num_strands]) {
        same = false;
        break; 
      }
    }
    if (same) {
      sym = num_strands / d;
      break;
    }
  }

  vector<uint> sym_div;
  get_divisors(sym, orb_sizes);
  
}


void
parse_pfunc_input_file(const string &input_file, uint &num_strands,
                       vector<vector<base> > &strands,
                       vector <uint> &strand_ids) {
  std::ifstream ifile;
  ifile.open(input_file, std::ios::in);
  string line, seq_str, ss;
  getline(ifile, line);
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  num_strands = std::stoi(line);
  for (uint i = 0; i < num_strands; i++) {
    getline(ifile, line);
    // line = string(line.begin(), line.end()-1);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    // line = std::regex_replace(line, newlines_re, "");
    vector<base> seq;
    for (char c: line)
      seq.push_back(base_to_num(c));
    strands.push_back(seq);
  }
  getline(ifile, line);
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  // line = std::regex_replace(line, newlines_re, "");
  stringstream check1(line);
  string word;
  while(getline(check1, word, ' ')) {
    strand_ids.push_back(std::stoi(word));
  }
  // vector<unsigned> seq;
  // for (char c: seq_str) {
  //   seq.push_back(base_to_num(c));
  // }

}

void
parse_energy_input_file(const string &input_file, uint &num_strands,
                        vector<vector<base> > &strands,
                        vector <uint> &strand_ids,
                        string & ss) {
  std::ifstream ifile;
  ifile.open(input_file, std::ios::in);
  string line, seq_str;
  getline(ifile, line);
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  num_strands = std::stoi(line);
  for (uint i = 0; i < num_strands; i++) {
    getline(ifile, line);
    // line = string(line.begin(), line.end()-1);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    // line = std::regex_replace(line, newlines_re, "");
    vector<base> seq;
    for (char c: line)
      seq.push_back(base_to_num(c));
    strands.push_back(seq);
  }
  getline(ifile, line);
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  // line = std::regex_replace(line, newlines_re, "");
  stringstream check1(line);
  string word;
  while(getline(check1, word, ' ')) {
    strand_ids.push_back(std::stoi(word));
  }
  getline(ifile, line);
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  ss = line;
  // vector<unsigned> seq;
  // for (char c: seq_str) {
  //   seq.push_back(base_to_num(c));
  // }

}



bool
is_dangle5(const unsigned i, const ss_node &parent){
  return (i > 0 && !parent.children[i-1]->paired &&
          parent.children[i-1]->u > 1);
}
bool
is_dangle3(const unsigned i, const ss_node &parent){
  return (i < parent.children.size()-1 && !parent.children[i+1]->paired);
}
bool
is_dangle5_closing_pair(const ss_node &pnode){
  return pnode.children.size() > 0 && !(pnode.children.back()->paired) &&
    pnode.children.back()->u > 1;
}
bool
is_dangle3_closing_pair(const ss_node &pnode){
  return pnode.children.size() > 0 && !pnode.children[0]->paired;
}


void
build_ss_tree(const string &ss, ss_node *node, std::stack<ss_node*> &bp_stack,
              unsigned &num_nicks){
  unsigned i = node->start;
  unsigned left = i;
  // int right = left-1;

  while(i != node->end - 1){
    if (ss[i+1] == '.'){
      i++;
      node->u++;
    }
    else if (ss[i+1] == '+') {
      i++;
      node->nick = true;
      num_nicks += 1;
    }
    else if (ss[i+1] == '('){
      node->p++;
      if (left != i){
        ss_node *unode = new ss_node;
        unode->paired = 0;
        unode->start = left + 1;
        unode->end = i;
        unode->u = unode->end - unode->start + 1;
        node->children.push_back(unode);
      }
      unsigned end;
      unsigned balance = 1;

      for (unsigned y = i + 2; y < node->end; y++){
        if (ss[y] == '(') balance++;
        else if (ss[y] == ')') balance--;
        if (balance == 0){
          end = y;
          break;
        }
      }
      ss_node *pnode = new ss_node;
      pnode->paired = 1;
      pnode->start = i + 1;
      pnode->end = end;
      pnode->u = 0;
      pnode->p = 0;
      node->pchild_indices.push_back(node->children.size());
      node->children.push_back(pnode);
      // printf("pair start at %d\n", pnode->start);

      i = end;
      left = i;
    }
  }
  if (left != i){
    ss_node *unode = new ss_node;
    unode->paired = 0;
    unode->start = left + 1;
    unode->end = i;
    node->children.push_back(unode);
  }
  for (vector<ss_node* >::reverse_iterator it = node->children.rbegin();
       it != node->children.rend(); it++){
    if ((*it)->paired) bp_stack.push(*it);
    // bp_stack.push(*it);
  }
  for (vector<ss_node* >::reverse_iterator it = node->children.rbegin();
       it != node->children.rend(); it++){
    if ((*it)->paired) build_ss_tree(ss, *it, bp_stack, num_nicks);
  }

}

void walk_ss_tree(const ss_node *node){

  if (node->paired){
    cout << '(';
    for (unsigned i = 0; i < node->children.size(); i++){
      walk_ss_tree(node->children[i]);
    }
    cout << ')';
  }
  else{
    for (unsigned i = node->start; i <= node->end; i++){
      cout << '.';
    }
  }
}

bool is_hairpin(const ss_node &node){
  return node.paired && (node.p == 0) && (!node.nick);
}

bool is_exterior(const ss_node &node){
  return node.paired && (node.nick);
}

bool is_interior(const ss_node &node){
  return node.paired && (node.p == 1) && (!node.nick);
}

bool is_multi(const ss_node &node){
  return node.paired && (node.p > 1) && (!node.nick);
}

bool is_polyC(const vector<char> &seq, const unsigned i, const unsigned j) {
  for (unsigned k = i + 1; k < j; k++) {
    if (seq[k] != base_C) return false;
  }
  return true;
}




bool can_pair(const vector<char> &seq, const unsigned i,
              const unsigned j) {
  return (seq[i] + seq[j] == 3) ||
    (seq[i] + seq[j] == 5);
}

bool is_nonGC(const vector<char> &seq, const unsigned i, const unsigned j) {
  return (seq[i] != base_C) && (seq[j] != base_C);
}

std::string GetCurrentWorkingDir() {
  char buff[filename_max];
  GetCurrentDir( buff, filename_max );
  std::string current_working_dir(buff);
  return current_working_dir;
}


string
convertToString(char* a, int size)
{
  int i;
  string s = "";
  for (i = 0; i < size; i++) {
    s = s + a[i];
  }
  return s;
}


bool
is_valid(const vector<base> &seq,
         const vector<vector<bool> > &eta,
         const vector<pair<unsigned, unsigned> > & s,
         const string &s_str){


  // Check if basems can pair
  for (const pair<unsigned, unsigned> &p : s) {
    unsigned i = p.first;
    unsigned j = p.second;
    if ((eta[i][j]) && (j <= i + hairpin_property::min_loop_size)) {
      return false;
    }
    else if (!bases_can_pair(seq[i], seq[j])) {
      return false;
    }
  }

  // Check that the structure is connected.
  std::stack<unsigned> vstack ({0}) ;
  for (unsigned i = 0; i < s_str.size(); i++) {
    if (s_str[i] == '+') vstack.top() -= 1;
    else if (s_str[i] == '(') {
      vstack.push(1);
    }
    else if (s_str[i] == ')') {
      if (vstack.top() < 0) {
        return false;
      }
      else {
        vstack.pop();
      }
    }
  }
  if (vstack.top() < 0) return false;
  return true;

}


uint
get_sym_num (const vector<base> &seq,
             const vector<vector<bool> > &eta,
             const vector<pos_pair> & s,
             const string &s_str,
             const uint sym,
             const vector<uint> &orb_sizes) {
  if (s.size() == 0) {
    //cout << 1 << endl;
    //cout << s_str << endl;
    
    return 1;
  }
  vector<base> seq_c;
  for (auto b : seq) {
    if (b != base_K) seq_c.push_back(b);
  }
  vector<pos_pair> s_co;
  for (auto pair : s) {
    pos i = pair.first;
    pos j = pair.second;
    s_co.push_back({i, j});
  }

  for (unsigned i = 0;  i < orb_sizes.size(); i++) {
    uint c = sym/orb_sizes[i];
    // printf ("try sum number: %d\n", c);

    if (s_co.size() % c != 0) continue;
    uint inc = seq_c.size()/c;
    vector<pos_pair> s_c (s_co.begin(), s_co.end());
    // printf ("length: %d\n", s_co.size());

    for (pos_pair p = s_c[0]; s_c.size() > 0;) {

      bool valid = true;
      for (uint d = 0; d < c; d++) {
        unsigned pos1 = (p.first + d*inc) % seq_c.size();
        unsigned pos2 = (p.second + d*inc) % seq_c.size();
        unsigned min_pos = min(pos1, pos2);
        unsigned max_pos = max(pos1, pos2);
        pos_pair p_pair = {min_pos, max_pos};
        const vector<pos_pair>::const_iterator p_it =
          find(s_c.begin(), s_c.end(), p_pair);
        if (p_it == s_c.end()) {
          valid = false;
          break;
        }
        else {
          // printf ("found (%d, %d)\n", min_pos, max_pos);
          s_c.erase(p_it);
        }
      }
      if (!valid) {
        break;
      }
      else {
        if (s_c.size() > 0) {

          p = s_c[0];
        }
      }
    }
    if (s_c.size() == 0){
      //cout << c << endl;
      //cout << s_str << endl;
      return c;
    }

  }
  //cout << 1 << endl;
  //cout << s_str << endl;
  return 1;
}


void
get_pair_matrix(const vector<base> &seq,
                const vector<vector<bool> > &eta,
                vector<vector<uint> > & pair_matrix){
  uint size = seq.size();
  for (uint i = 0; i < size; i++) {
    vector<uint> v;
    for (uint j = i+1; j < size; j++) {
      if (eta[i][j]) {
        if (i + hairpin_property::min_loop_size < j)
          if (bases_can_pair(seq[i], seq[j]))
            v.push_back(j);
      }
      else {
        if (bases_can_pair(seq[i], seq[j]))
          v.push_back(j);
      }
    }
    pair_matrix.push_back(v);
  }

  //for (uint i = 0; i < pair_matrix.size(); i++) {
    //for (uint j = 0; j < pair_matrix[i].size(); j++) {
      //cout << pair_matrix[i][j] << " ";
    //}
    //cout << endl;
  //}

}


void
get_eta(const vector<base> & seq,
        const vector<pos> &nicks,
        vector<vector<bool> > & eta) {
  eta.resize(seq.size(), vector<bool>(seq.size()));
  for (unsigned i = 0; i < seq.size(); i++) {
    for (unsigned j = i; j < seq.size(); j++) {
      eta[i][j] = 0;
    }
  }
  pos start = 0;
  pos end = 0;
  for (unsigned k = 0; k < nicks.size(); k++) {
    end = nicks[k];
    for (pos i = start; i <= end; i++) {
      for (pos j = i; j <= end; j++) {
        eta[i][j] = 1;
      }
    }
    start = nicks[k]+1;
  }
  end = seq.size()-1;
  for (pos i = start; i <= end; i++) {
    for (pos j = i; j <= end; j++) {
      eta[i][j] = 1;
    }
  }

}

void
get_p(const vector<base> & seq,
    const vector<pos> &nicks,
    const vector<vector<bool> > & eta,
    vector<vector<bool> > & p,
    vector<vector<bool> > & p_red) {
  p.resize(seq.size(), vector<bool>(seq.size(), false));
  p_red.resize(seq.size(), vector<bool>(seq.size(), false));
  for (unsigned i = 0; i < seq.size(); i++) {
    for (unsigned j = i+1; j < seq.size(); j++) {
      if ((!eta[i][j]) || (j - i > 3)){
        p[i][j] = rna_bases_can_pair(seq[i], seq[j]);       
      }
      p_red[i][j] = rna_bases_can_pair(seq[i], seq[j]);       
    }
  }
}

void
get_seq(const vector<vector<base> > &strands,
        const vector<uint> &strand_ids,
        vector<base> &seq,
        vector<uint> & nicks){
  uint n_strands = strand_ids.size();
  for (uint i = 0; i < n_strands; i++) {
    for (base b : strands[strand_ids[i]-1]) {
      seq.push_back(b);
    }
    nicks.push_back((uint)seq.size()-1);
  }
  nicks.pop_back();

}


void
get_seq(const vector<vector<base> > &strands,
        const vector<uint> &strand_ids,
        vector<base> &seq,
        vector<uint> & nicks,
        vector<bool> & gam){
  uint n_strands = strand_ids.size();
  for (uint i = 0; i < n_strands; i++) {
    for (base b : strands[strand_ids[i]-1]) {
      seq.push_back(b);
    }
    nicks.push_back((uint)seq.size()-1);
  }
  nicks.pop_back();
  gam.resize(seq.size()+1, true);
  for (uint d : nicks){
    gam[d] = false; 
  }
}

void
get_base_ss(const vector<uint> &nicks,
            const uint size,
            const vector<uint> &displacement,
            string &base_s) {
  uint d = 0;
  for (uint i = 0; i < size; i++) {
    if (d != displacement[i]) {
      d = displacement[i];
      base_s.push_back('+');
    }  
    base_s.push_back('.');
  }
  //uint ind = 0;
  //for (uint n : nicks) {
    //for (;ind < n+1; ind++) {
      //base_s.push_back('.');
    //}
    //base_s.push_back('+');
    //ind = n+1;
  //}
  //for (;ind < size;ind++) {
    //base_s.push_back('.');
  //}
}


void
ss2string (const vector<pos_pair> &pairs,
    const vector<uint> &displacement,
           string &s) {
  for (auto p : pairs) {
    s[p.first+displacement[p.first]] = '(';
    s[p.second+displacement[p.second]] = ')';
  }
}
