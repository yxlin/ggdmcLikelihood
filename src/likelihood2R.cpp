
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
//' @param dmis A list of S4 data model instances (one per subject)
//' @param dmi One S4 data model instance (for one subject)
//' @param parameter_r A list (one per subject) of or one numeric vector
//'                    containing model parameters
//' @param debug Logical flag for debug mode (default = FALSE)
//'
//' @return
//' with:
//' \itemize{
//'   \item \code{compute_likelihood} returns a list. Each element is
//'          the likelihood for a subject. The element in the inner list
//'          is the likelihood for a condition.
//'   \item \code{compute_subject_likelihood} returns also a list. Each
//'          element is the likelihood for a condition.
//' }
//'
//' @details These functions provide access to the internal mechanism of 
//' the design-based likelihood computation. They primarily intended to 
//' initialise new 'samples' or to verify that the likelihood evaluations, 
//' when associated with a particular design, are computed accurately.
//'
//' @examples
//' # Example dataset
//' hdat <- data.frame(
//'   RT = round(runif(15, min = 0.4, max = 1.2), 7),
//'   R  = sample(c("r1", "r2", "r3"), size = 15, replace = TRUE),
//'   s  = rep(1:3, each = 5),
//'   S  = rep(c("s1", "s2", "s3"), each = 5),
//'   stringsAsFactors = FALSE
//' )
//' dat <- hdat[hdat$s==1, ]
//'
//' p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)
//' nsubject <- length(unique(hdat$s))
//'
//' if(requireNamespace("ggdmcModel", quietly = TRUE)) {
//'     BuildModel <- getFromNamespace("BuildModel", "ggdmcModel")
//'     BuildDMI   <- getFromNamespace("BuildDMI", "ggdmcModel")
//' 
//'     model <- BuildModel(
//'         p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
//'         match_map = list(M = list(s1 = "r1", s2 = "r2")),
//'         factors = list(S = c("s1", "s2")),
//'         constants = c(st0 = 0, sd_v = 1),
//'         accumulators = c("r1", "r2"),
//'         type = "lba")
//'     pop_dmis <- BuildDMI(hdat, model)
//'     sub_dmis <- BuildDMI(dat, model)
//'
//'     parameters <- list()
//'     for (i in seq_len(nsubject)) {
//'         new_p_vector <- p_vector[model@pnames]
//'         parameters[[i]] <- new_p_vector
//'     }
//'
//'     result1 <- compute_subject_likelihood(sub_dmis[[1]], parameters[[1]], FALSE)
//'     result2 <- compute_likelihood(pop_dmis, parameters, FALSE)
//' }
//'
//' print(result1)
//' print(result2)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List compute_likelihood(const Rcpp::List &dmis,
                              const Rcpp::List &parameter_r, bool debug = false)
{
    auto l_ptr = new_likelihoods(dmis);
    unsigned int n_subject = l_ptr.size();

    Rcpp::List out(n_subject);

    for (size_t subject_idx = 0; subject_idx < n_subject; ++subject_idx)
    {
        if (debug)
        {
            l_ptr[subject_idx]->print_is_empty_cell("Empty cells: ");
            l_ptr[subject_idx]->print_rt();
        }

        auto parameters =
            Rcpp::as<std::vector<double>>(parameter_r[subject_idx]);

        size_t n_cell = l_ptr[subject_idx]->m_model->m_n_cell;
        Rcpp::List cell_out(n_cell);

        for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx)
        {
            l_ptr[subject_idx]->likelihood(parameters, debug);
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
                                      bool debug = false)
{
    Rcpp::S4 model_r = dmi.slot("model");
    Rcpp::List data_r = dmi.slot("data");
    std::string model_str = model_r.slot("type");

    auto l_ptr = new_likelihood(dmi);
    if (debug)
    {
        l_ptr->print_is_empty_cell("Empty cells: ");
        l_ptr->print_rt();
    }

    auto parameters = Rcpp::as<std::vector<double>>(parameter_r);
    size_t n_cell = l_ptr->m_model->m_n_cell;

    // This step computes likelihoods and creates an output, named 'm_density'.
    l_ptr->likelihood(parameters, debug);
    Rcpp::List cell_out(n_cell);

    for (size_t cell_idx = 0; cell_idx < n_cell; ++cell_idx)
    {
        cell_out[cell_idx] = l_ptr->m_density[cell_idx];
    }

    return cell_out;
}
