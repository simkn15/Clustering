#include <Rcpp.h>
#include "transclust/TransClust.hpp"

using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

List cppCluster(TransClust& tc){
   clustering res = tc.cluster();

   return List::create(
      Named("labels") = res.id2object,
      Named("costs") = res.cost,
      Named("thresholds") = res.threshold,
      Named("clusters") = res.clusters
   );
}

// [[Rcpp::export]]
List cppTransClustFile(
      std::string filename,
      std::string file_type,
      List params
   ) {
   TransClust tc(
         filename,
         TCC::TransClustParams()
            .set_file_type(file_type)
            .set_sim_fallback(as<float>(params["sim_fallback"]))
            .set_threshold(as<float>(params["threshold"]))
            .set_p(as<float>(params["p"]))
            .set_f_att(as<float>(params["f_att"]))
            .set_f_rep(as<float>(params["f_rep"]))
            .set_R(as<unsigned>(params["R"]))
            .set_dim(as<unsigned>(params["dim"]))
            .set_start_t(as<float>(params["start_t"]))
            .set_d_init(as<float>(params["d_init"]))
            .set_d_maximal(as<float>(params["d_maximal"]))
            .set_s_init(as<float>(params["s_init"]))
            .set_f_s(as<float>(params["f_s"]))
            .set_fpt_time_limit(as<float>(params["fpt_time_limit"]))
            .set_fpt_max_cost(as<float>(params["fpt_max_cost"]))
            .set_fpt_step_size(as<float>(params["fpt_step_size"]))
            .set_disable_force(as<bool>(params["disable_force"]))
            .set_disable_fpt(as<bool>(params["disable_fpt"]))
            .set_seed(as<unsigned>(params["seed"]))
   );
   return cppCluster(tc);
}

// [[Rcpp::export]]
List cppTransClustDist(
      NumericVector sim_matrix_1d,
      unsigned num_o,
      List params) {
   std::vector<float> sm(sim_matrix_1d.begin(),sim_matrix_1d.end());
   TransClust tc(
         sm,
         num_o,
         TCC::TransClustParams()
            .set_sim_fallback(as<float>(params["sim_fallback"]))
            .set_threshold(as<float>(params["threshold"]))
            .set_p(as<float>(params["p"]))
            .set_f_att(as<float>(params["f_att"]))
            .set_f_rep(as<float>(params["f_rep"]))
            .set_R(as<unsigned>(params["R"]))
            .set_dim(as<unsigned>(params["dim"]))
            .set_start_t(as<float>(params["start_t"]))
            .set_d_init(as<float>(params["d_init"]))
            .set_d_maximal(as<float>(params["d_maximal"]))
            .set_s_init(as<float>(params["s_init"]))
            .set_f_s(as<float>(params["f_s"]))
            .set_fpt_time_limit(as<float>(params["fpt_time_limit"]))
            .set_fpt_max_cost(as<float>(params["fpt_max_cost"]))
            .set_fpt_step_size(as<float>(params["fpt_step_size"]))
            .set_disable_force(as<bool>(params["disable_force"]))
            .set_disable_fpt(as<bool>(params["disable_fpt"]))
            .set_seed(as<unsigned>(params["seed"]))
   );
   return cppCluster(tc);
}
