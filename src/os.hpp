#ifndef OS
#define OS


#include "utils2.hpp"
#include "fstream"
#include "filesystem"

#include <vector>
#include <string>
#include <iostream>


void
parse_input(const std::filesystem::path &input_path,
    size_t &n_seqs,
    std::vector<std::vector<bnum> > seqs, 
    std::vector<size_t> &seq_ids);


#endif
