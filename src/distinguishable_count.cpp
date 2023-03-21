#include "distinguishable_count.hpp"

using std::string;
using std::vector;
using std::stringstream;
using std::map;
using std::cerr;
using std::runtime_error;
using std::cout;
using std::endl;
using std::filesystem::path;

//void
//parse_input(input_path, n_uniq_path, seqs, seq_ids);



//void
//parse_input_file(const string &input_file, size_t &num_strands,
                       //vector<vector<base> > &strands,
                       //vector<size_t> &strand_ids) {
  //std::ifstream ifile;
  //ifile.open(input_file, std::ios::in);
  //string line, seq_str, ss;
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

//}





int main(int argc, const char **argv) {

  try {
    // Get input file
    if (argc == 1) {
      throw std::runtime_error("Missing input file in " __FILE__ " at line " +
          std::to_string(__LINE__));
    }
    path input_path = argv[1];
    
    // Parse input file
    size_t n_uniq_seqs;
    vector<vector<bnum> > seqs;
    vector<size_t> seq_ids;
    parse_input(input_path, n_uniq_seqs, seqs, seq_ids);
    cout << "made it!" << endl;
    






    

  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
