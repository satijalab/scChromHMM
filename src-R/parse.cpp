#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP get_state(std::string fpath, 
               std::string chr_name,
               std::string cell_names_file) {
    std::ifstream file (fpath, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
      std::cout << "ERROR opening file";
    }

    std::vector<uint32_t> mat_size (2);
    file.read(reinterpret_cast<char*>(mat_size.data()), sizeof(uint32_t) * 2);

    size_t num_rows { mat_size[0] };
    size_t num_cols { mat_size[1] };

    std::vector<uint32_t> indptr (num_cols);
    file.read(reinterpret_cast<char*>(indptr.data()), sizeof(uint32_t) * num_cols);

    num_cols -= 1;
    size_t nnz = indptr[num_cols];
    std::vector<uint8_t> values (nnz);
    std::vector<uint32_t> indices (nnz);
    file.read(reinterpret_cast<char*>(values.data()), sizeof(uint8_t) * nnz);
    file.read(reinterpret_cast<char*>(indices.data()), sizeof(uint32_t) * nnz);
    
    uint8_t check;
    if (file >> check) {
      std::cout << "END of File not reached" << check << std::endl;
    } else {
      file.close();
    }
    
    std::vector<std::string> cell_names {num_cols}; 
    {
      std::ifstream file (cell_names_file, std::ios::in | std::ios::binary);
      if (!file.is_open()) {
        std::cout << "ERROR";
      }
      for (size_t i=0; i<num_cols; i++) {
        file >> cell_names[i];
      }
      file.close();
    } // filling cell names
    
    std::vector<std::string> region_names {num_rows}; 
    {
      for (size_t i=0; i<num_rows; i++) {
        region_names[i] = chr_name + '-' + std::to_string(i*200) + '-' + std::to_string((i+1)*200);
      }
    } // filling region names

    Rcpp::S4 mat("dgCMatrix");
    Rcpp::IntegerVector rcpp_indptr(indptr.begin(), indptr.end());
    Rcpp::IntegerVector rcpp_indices(indices.begin(), indices.end());
    Rcpp::NumericVector rcpp_values(values.begin(), values.end());

    //std::cout << "Mat Size: " << num_rows << " x " << num_cols << "; with nnz: " << nnz;
    mat.slot("Dim") = Rcpp::IntegerVector::create(num_rows, num_cols);

    // Setting p
    mat.slot("p") = rcpp_indptr;

    // Setting 'x'.
    mat.slot("x") = rcpp_values;

    // Setting 'i'.
    mat.slot("i") = rcpp_indices;
    
    // setting dimnames
    mat.slot("Dimnames") = Rcpp::List::create(region_names, cell_names);

    return SEXP(mat);
}
