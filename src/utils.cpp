#include <Rcpp.h>
#include <vector>
#include <string>

//' @export
// [[Rcpp::export]]
std::vector<int> pairwise_hamming(const std::vector<std::string> &strs) {
  std::vector<int> distances;
  for (int i = 0; i < strs.size(); ++i) {
    for (int j = i + 1; j < strs.size(); ++j) {
      const std::string &str1 = strs[i];
      const std::string &str2 = strs[j];
      if (str1.length() != str2.length())
        Rcpp::stop("Strings must have the same length");

      int ed = 0;
      for (int c_id = 0; c_id < str1.length(); ++c_id) {
        if (str1[c_id] != str2[c_id]) {
          ed++;
        }
      }
      distances.push_back(ed);
    }
  }

  return distances;
}
