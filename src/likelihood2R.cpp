
#include <ggdmcHeaders/common_type_casting.h>
#include <ggdmcHeaders/likelihood.h>
#include <ggdmcHeaders/likelihood_type_casting.h>

//' Compute Likelihood for Behavioural Models
//'
//' These functions compute likelihoods for behavioural models, with
//' \code{compute_subject_likelihood} handling a single subject and
//' \code{compute_likelihood} handling multiple subjects.
//'
//' @description Computes the likelihood for multiple subjects by
//' aggregating results from individual subject computations.
//'
//' @param dmis A List of S4 data model instances (one per subject)
//' @param dmi One S4 data model instance (for one subject)
//' @param parameter_r A List (one per subject) of or one NumericVectors
//'                    containing model parameters
//' @param debug Logical flag for debug mode (default = FALSE)
//'
//' @return A List object containing aggregated results across subjects,
//' with:
//' \itemize{
//'   \item \code{likelihood} - The total likelihood value over multiple
//'   subject
//'   \item Additional debugging information when \code{debug = TRUE}
//' }
//'
//' @details These function expose the internal mechanism of the design-based
//' likelihood computation.
//' @examples
//' \dontrun{
//' model <- ggdmcModel::BuildModel(
//'     p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1",
//'                 t0 = "1"),
//'     match_map = list(M = list(s1 = "r1", s2 = "r2")),
//'     factors = list(S = c("s1", "s2")),
//'     constants = c(sd_v = 1, st0 = 0),
//'     accumulators = c("r1", "r2"),
//'     type = "lba" )
//'
//' dmis <- ggdmcModel::BuildDMI(hdat, model)
//' nsubject <- length(unique(hdat$s))
//' parameters <- list()
//' for (i in seq_len(nsubject)) {
//'     new_p_vector <- p_vector[model@pnames]
//'     parameters[[i]] <- new_p_vector
//' }
//'
//' result <- compute_subject_likelihood(sub_dmis[[1]], parameters[[1]], F)
//' print(result)
//'
//' result <- compute_likelihood(pop_dmis, parameters, F)
//' print(result)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List compute_likelihood(const Rcpp::List &dmis,
                              const Rcpp::List &parameter_r,
                              bool debug = false) {
  auto l_ptr = new_likelihoods(dmis);
  unsigned int n_subject = l_ptr.size();

  Rcpp::List out(n_subject);

  for (size_t subject_idx = 0; subject_idx < n_subject; ++subject_idx) {
    if (debug) {
      l_ptr[subject_idx]->print_is_empty_cell("Empty cells: ");
      l_ptr[subject_idx]->print_rt();
    }

    auto parameters = Rcpp::as<std::vector<double>>(parameter_r[subject_idx]);

    size_t n_cell = l_ptr[subject_idx]->m_model->m_n_cell;
    Rcpp::List cell_out(n_cell);

    for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx) {

      l_ptr[subject_idx]->likelihood(parameters, debug);
      // l_ptr[subject_idx]->lba_likelihood(parameters, debug);
      cell_out[cell_idx] = l_ptr[subject_idx]->m_density[cell_idx];
    }
    out[subject_idx] = cell_out;
  }

  return out;
}
//' @rdname compute_likelihood
//' @export
// [[Rcpp::export]]
Rcpp::List compute_subject_likelihood(const Rcpp::S4 &dmi,
                                      const Rcpp::NumericVector &parameter_r,
                                      bool debug = false) {

  Rcpp::S4 model_r = dmi.slot("model");
  Rcpp::List data_r = dmi.slot("data");
  std::string model_str = model_r.slot("type");

  auto l_ptr = new_likelihood(dmi);
  if (debug) {
    l_ptr->print_is_empty_cell("Empty cells: ");
    l_ptr->print_rt();
  }

  auto parameters = Rcpp::as<std::vector<double>>(parameter_r);
  size_t n_cell = l_ptr->m_model->m_n_cell;

  // This step computes likelihoods and create an output, m_density.
  l_ptr->likelihood(parameters, debug);

  // if (model_str == "lba") {
  //   l_ptr->lba_likelihood(parameters, debug);
  // } else if (model_str == "fastdm") {
  //   l_ptr->ddm_likelihood(parameters, debug);
  // } else {
  //   throw std::runtime_error("Undefined model type");
  // }

  Rcpp::List cell_out(n_cell);

  for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx) {
    cell_out[cell_idx] = l_ptr->m_density[cell_idx];
  }

  return cell_out;
}
