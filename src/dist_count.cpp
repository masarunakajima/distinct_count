//#include "read_constants_nupack.hpp"
#include "constants.hpp"
#include "utils.hpp"
//#include "parameters.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::runtime_error;
using namespace ss_names;


struct
count_matrices{
  s_matrix S;
  s_matrix Se;
  s_matrix Sb;

  s_matrix R;
  s_matrix Re;
  s_matrix Rb;

  void
  initialize(size_t m_size){
    S = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));
    Se = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));
    Sb = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));

    R = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));
    Re = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));
    Rb = s_matrix(m_size+1, vector<ull>(m_size + 1, 0));

    for (pos i = 0; i <= m_size; i++) S[i][i] = 1;
    for (pos i = 1; i <= m_size; i++) S[i][i-1] = 1;
  }
};





void
update_Se(const vector<bool> &gam,
    const vector<pos> &nicks,
    const pos i, const pos j,
    count_matrices &cm){
   
  // Do nothinkg if nicks are at the edge 
  if ((gam[i] + gam[j-1] == 0) && (i < j-1)) {
  }
  // Only account for single unit at i or j-1
  else if (gam[i] + gam[j-1] == 1) {
    cm.Se[i][j] = cm.S[i+1][j-1]; 
  }
  // General case
  else {
    for (pos d : nicks) {
      if ((i <= d) && (d < j)) {
        cm.Se[i][j] += cm.S[i+1][d]*cm.S[d+1][j-1]; 
      }
    }
  }
}

void
update_Sb(const vector<bool> &gam,
    const vector<vector<bool> > &p,
    const pos i, const pos j,
    count_matrices &cm){
  cm.Sb[i][j] = p[i][j]*(gam[i]*gam[j-1]*cm.S[i+1][j-1] 
    + cm.Se[i][j]);
}

void
update_S(const vector<bool> &gam,
    const pos i, const pos j,
    count_matrices &cm){
  cm.S[i][j] = gam[i]*cm.S[i+1][j];
  for (pos k = i+1; k <= j; k++){
    cm.S[i][j] += gam[k]*cm.Sb[i][k]*cm.S[k+1][j];
  }
}



void
update_Re(const vector<bool> &gam,
    const vector<pos> &nicks,
    const pos i, const pos j,
    count_matrices &cm){
  // do nothing if there are nicks on both ends
  if ((gam[i] + gam[j-1] == 0) ) {

  }
  // if there is nick on only one end
  else if (gam[i] + gam[j-1] == 1) {
    cm.Re[i][j] = cm.R[i+1][j-1]; 
  }
  // general case
  else {
    for (pos d : nicks) {
      if ((i <= d) && (d < j)) {
        cm.Re[i][j] += cm.R[i+1][d]*cm.S[d+1][j-1]
          + cm.S[i+1][d]*cm.R[d+1][j-1]; 
      }
    }
  }
}


void
update_Rb(const vector<bool> &gam,
    const vector<vector<bool> > &p_red,
    const pos i, const pos j,
    count_matrices &cm){
  cm.Rb[i][j] = p_red[i][j]*
    (gam[i]*gam[j-1]*(cm.S[i+1][j-1]+cm.R[i+1][j-1]) 
    + cm.Re[i][j]);
}

void
update_R(const vector<bool> &gam,
    const pos i, const pos j,
    count_matrices &cm){
  cm.R[i][j] = gam[i]*cm.R[i+1][j] + cm.Rb[i][j];
  for (pos k = i+1; k < j; k++){
    cm.R[i][j] += gam[k]*cm.Sb[i][k]*cm.R[k+1][j];
    cm.R[i][j] += gam[k]*cm.Rb[i][k]*cm.S[k+1][j];
  }
}


int
main(int argc, const char **argv) {
  try {

    string outfile;
    string param_name = "rna1995";
    string param_dir;
    // run mode flags
    bool VERBOSE = false;
    double T = 37.0;

    const string description =
      "Calculate free energy of secondary structure for a particular\
      sequence. The input file must contain the sequence in the first line\
      and the secondary structure in the dot-parenthesis notation in the\
      seconda line.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           description, "<sequence-structure-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("params", 'p', "parameter name (default: rna1995)",
                      false, param_name);
    opt_parse.add_opt("paramdir", 'd', "parameter directory",
                      false, param_dir);
    opt_parse.add_opt("temp", 't', "temeprature in celcius (default: 37)",
                      false, T);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string input_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    uint num_strands;
    vector<vector<base> > strands;
    vector<uint> strand_ids;
    parse_pfunc_input_file(input_file, num_strands, strands, strand_ids);
    uint sym;
    vector<uint> orb_sizes;
    get_orb_sizes(strand_ids, sym, orb_sizes);
    vector<base> tot_seq;
    vector<pos> nicks;
    vector<bool> gam;
    get_seq(strands, strand_ids, tot_seq, nicks, gam);
    unsigned size = tot_seq.size();

    vector<vector<bool> > eta(size,  vector<bool>(size,0));
    get_eta(tot_seq, nicks, eta);
    vector<vector<bool> > p(size,  vector<bool>(size,0));
    vector<vector<bool> > p_red(size,  vector<bool>(size,0));
    get_p(tot_seq, nicks, eta, p, p_red);

    size_t seq_size = size;
    count_matrices cm;

    cm.initialize(seq_size);
    for (pos j = 1; j < (pos)seq_size ; j++) {
      for (pos i = j - 1; i != (unsigned)-1 ; i--) {
        update_Se(gam, nicks, i, j, cm);
        update_Sb(gam, p, i, j, cm);
        update_S(gam, i, j, cm);
        update_Re(gam, nicks, i, j, cm);
        update_Rb(gam, p_red, i, j, cm);
        update_R(gam, i, j, cm);

        //printf("s[%d][%d] =  ", i,j);
        //std::cout << cm.S[i][j] << std::endl;
        // printf("sb[%d][%d] = ", i,j);
        // std::cout << cm.Sb[i][j] << std::endl;
        // printf("se[%d][%d] = ", i,j);
        // std::cout << cm.Se[i][j] << std::endl;
        //printf("re[%d][%d] = ", i,j);
        //std::cout << cm.Re[i][j] << std::endl;
        //printf("rb[%d][%d] = ", i,j);
        //std::cout << cm.Rb[i][j] << std::endl;
        //printf("r[%d][%d] =  ", i,j);
        //std::cout << cm.R[i][j] << std::endl;
        // printf("rb[%d][%d] = ", i,j);
        // std::cout << cm.Rb[i][j] << std::endl;
        // printf("re[%d][%d] = ", i,j);
        // std::cout << cm.Re[i][j] << std::endl;

      }
    }
    ull count = cm.S[0][seq_size-1];

    vector<ull> sym_structures;
    uint n_sym = size/sym; 
    for (uint k = 0; k < orb_sizes.size()-1; k++){
      uint b = orb_sizes[k];
      sym_structures.push_back(cm.R[0][b*n_sym-1]);
    }
    sym_structures.push_back(cm.S[0][seq_size-1]);

    vector<ull> sym_cor_structures;
    for (uint k = 0; k < orb_sizes.size(); k++){
      uint b = orb_sizes[k];
      ull count = sym_structures[k];
      for (uint j = 0; j < k; j++){
        uint a = orb_sizes[j];
        if (b % a == 0) {
          count -= sym_cor_structures[j];
        }
      }
      sym_cor_structures.push_back(count);
    }

    ull dist_count = 0;
    for (uint k = 0; k < orb_sizes.size(); k++){
      uint b = orb_sizes[k];
      dist_count += sym_cor_structures[k]/b;
    }



    // output
    if (outfile == "") {
      
      printf("Number of structures\n");
      std::cout << "Total\t: " << count << std::endl;
      // printf("Total\t: %lld\n",  count);

      //for (uint k = 0; k < sym; k++) {
        //std::cout << "b = " << orb_sizes[k] << "\t: " <<
          //sym_structures[k]*max_sym/sym[i] << std::endl;
      //}
      printf("Number of distinguishable structures\n");
      std::cout << "Total\t: " << dist_count << std::endl;

      //for (uint i = 0; i < sym.size(); i++) {
        //std::cout << "sym = " << sym[i] << "\t: " <<
          //sym_structures[i] << std::endl;

      //}

    }
    else {
      std::fstream fs;
      fs.open (outfile, std::fstream::out);
      fs << "Number of structures\n";
      fs << "Total\t: " << count << std::endl;
      // printf("Total\t: %lld\n",  count);

      //for (uint i = 0; i < sym.size(); i++) {
        //fs << "sym = " << sym[i] << "\t: " <<
          //sym_structures[i]*max_sym/sym[i] << std::endl;
      //}
      fs <<"Number of distinguishable structures\n";
      fs << "Total\t: " << dist_count << std::endl;

      //for (uint i = 0; i < sym.size(); i++) {
        //fs << "sym = " << sym[i] << "\t: " <<
          //sym_structures[i] << std::endl;

      //}
      fs.close();

      // FILE *fp;
      // fp = fopen(outfile.c_str(),"w");
      // for (uint i = 0; i < sym.size(); i++) {
      //   fprintf(fp, "sym = %d: %ld\n", sym[i], sym_structures[i]);
      // }
      // fclose(fp);
    }
  }
  catch (runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}



//TODO: Make function to check if the input is valid
