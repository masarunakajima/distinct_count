#include "os.hpp"



using std::filesystem::path;
using std::vector;
using std::string;
using std::cerr;
using std::endl;
using std::stoi;

void
parse_input(const path &input_path,
    size_t &n_seqs,
    vector<vector<bnum> > seqs, 
    vector<size_t> &seq_ids) {

  std::ifstream ifile(input_path);
  if (!ifile.is_open()) {
    throw std::runtime_error("Error: could not open file: " + 
        input_path.string() + " in " + __FILE__ + " at line " + 
        std::to_string(__LINE__));
  }


  string line, seq_str, ss;
  size_t line_num = 0;
  while (std::getline(ifile, line)) {
    if (line_num == 0) {
      n_seqs = stoi(line);
    }
    else if (line_num <= n_seqs) {
      if (!line.empty() && line.back() == '\n') {
        line.pop_back();
      }
    }
    line_num++;
  }
  //getline(ifile, line);

  //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  //num_strands = std::stoi(line);
  //for (uint i = 0; i < num_strands; i++) {
    //getline(ifile, line);
    //// line = string(line.begin(), line.end()-1);
    //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    //// line = std::regex_replace(line, newlines_re, "");
    //vector<base> seq;
    //for (char c: line)
      //seq.push_back(base_to_num(c));
    //strands.push_back(seq);
  //}
  //getline(ifile, line);
  //line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  //// line = std::regex_replace(line, newlines_re, "");
  //stringstream check1(line);
  //string word;
  //while(getline(check1, word, ' ')) {
    //strand_ids.push_back(std::stoi(word));
  //}


}
