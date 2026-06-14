

#' ADAPT-ME: Longitudinal Microbiome Differential Abundance Analysis
#'
#' @description
#' Runs the full ADAPT-ME pipeline: BUM-guided reference selection with
#' within-subject stability scoring, followed by mixed-effects Tobit
#' count-ratio tests with a subject random intercept.
#'
#' @details
#' ADAPT-ME takes in a longitudinal metagenomics count table as a phyloseq object.
#' The phyloseq object must have metadata containing at least one condition variable
#' (`cond.var`) and one subject identifier (`subject.var`) for repeated-measures data.
#' The condition variable can be numeric (continuous) or character (categorical).
#' ADAPT-ME does not support multigroup comparison. If there are multiple groups,
#' specify the baseline through `base.cond`; ADAPT-ME then carries out DAA between
#' `base.cond` and all others. Additional covariates can be adjusted for via `adj.var`.
#'
#' Zero counts are treated as left-censored observations. Mixed-effects Tobit models
#' are fitted for log count ratios between each taxon and the summed reference counts,
#' with a subject-level random intercept integrated out via 15-node Gauss-Hermite
#' quadrature. This removes between-subject baseline variation from the residual,
#' substantially improving power in paired and longitudinal designs.
#'
#' Reference taxa are selected iteratively using a BUM p-value model. The reference
#' score combines median fold-change distance with within-subject standard deviation,
#' so that temporally unstable taxa are not mistaken for stable references.
#'
#' For cross-sectional data (no repeated measures), use `adapt_crosssectional()`
#' instead, which fits fixed-effects Tobit models and does not require `subject.var`.
#'
#' Rare taxa and low-depth samples can be filtered with `prev.filter` and
#' `depth.filter`. DA taxa are called at BH-adjusted p-value < `alpha` (default 0.05).
#'
#' The returned value is a `DAresult` S4 object. Use `summary()` and `plot()` to
#' explore the output.
#' 
#' 
#' @param input_data a phyloseq object
#' @param cond.var the variable representing the conditions to compare, a character string
#' @param base.cond the condition chosen as baseline. This is only used when the condition is categorical.
#' @param adj.var the names of the variables to be adjusted, a vector of character strings
#' @param subject.var the metadata variable identifying repeated measurements from the same subject
#' @param censor the value to censor at for zero counts, default 1
#' @param prev.filter taxa whose prevalences are smaller than the cutoff will be excluded from analysis, default 0.05
#' @param depth.filter a sample would be discarded if its library size is smaller than the threshold
#' @param alpha the cutoff of the adjusted p values
#' @importFrom stats optim median p.adjust qchisq
#' @returns a `DAresult` type object contains the input and the output. Use summary and plot to explore the output
#' @export
#' 
#' @examples
#' \dontrun{
#' longitudinal_results <- adapt(input_data=longitudinal_phyloseq,
#'        cond.var="Timepoint", base.cond="baseline",
#'        adj.var="Site", subject.var="SubjectID")
#' }
adapt <- function(input_data, cond.var, base.cond = NULL, adj.var=NULL, subject.var, censor=1,
                  prev.filter=0.05, depth.filter=1000, alpha=0.05){
  
  # preprocess input phyloseq object
  preprocessed_output <- preprocess(input_data, cond.var, base.cond, adj.var, subject.var, prev.filter, depth.filter)
  
  # check if the arguments are valid
  stopifnot("The zero counts need to be censored at a positive number!" = censor > 0)
  stopifnot("The cutoff for adjusted p-values needs to be between 0 and 0.5!" = alpha > 0 & alpha < 0.5)
  
  count_table <- preprocessed_output$count_table
  complete_design_matrix <- preprocessed_output$design_matrix
  DAAname <- preprocessed_output$DAAname
  
  taxa_names <- colnames(count_table)

  
  reftaxa <- taxa_names # initially all the taxa are reference taxa(relative abundance)
  cat("Selecting Reference Set... ")
  while(1){
    relabd_result <- count_ratio(count_table = count_table, design_matrix = complete_design_matrix,
                                 subject = preprocessed_output$subject, censor = censor,
                                 reftaxa = reftaxa, test_all=FALSE)
    
    estimated_effect <- relabd_result$log10foldchange
    pvals <- relabd_result$pval
    names(pvals) <- relabd_result$Taxa
    names(estimated_effect) <- relabd_result$Taxa
    # check distribution of p values
    bumfit <- BUM_fit(pvals)
    loglik <- BUM_llk(bumfit$estim_params, pvals[!is.na(pvals)])
    if (2*loglik > qchisq(0.95, 1)){ # need to continue shrinking reference taxa set
      distance2med <- abs(estimated_effect - median(estimated_effect, na.rm=TRUE))
      stability <- relabd_result$within_subject_sd
      names(stability) <- relabd_result$Taxa
      distance_scale <- median(distance2med, na.rm=TRUE)
      stability_scale <- median(stability, na.rm=TRUE)
      if (!is.finite(distance_scale) || distance_scale <= 0) distance_scale <- 1
      if (!is.finite(stability_scale) || stability_scale <= 0) stability_scale <- 1
      scaled_distance <- distance2med / distance_scale
      scaled_stability <- stability / stability_scale
      reference_score <- scaled_distance + scaled_stability
      reference_score[!is.finite(reference_score)] <- Inf
      sorted_score <- sort(reference_score)
      ordered_taxanames <- names(sorted_score)
      reftaxa <- ordered_taxanames[seq_len(length(ordered_taxanames)/2)]
    } else{
      break
    }
  }
  cat(sprintf("%d taxa selected as reference\n", length(reftaxa)))
  all_CR_results <- count_ratio(count_table=count_table, design_matrix=complete_design_matrix,
                                subject = preprocessed_output$subject, censor = censor,
                                reftaxa=reftaxa, test_all=TRUE)

  all_pvals <- all_CR_results$pval
  names(all_pvals) <- all_CR_results$Taxa
  all_adjusted_pvals <- p.adjust(all_pvals, method="BH")
  all_CR_results$adjusted_pval <- all_adjusted_pvals

  significant_pvals <- all_adjusted_pvals[all_adjusted_pvals < alpha & !is.na(all_adjusted_pvals)]
  if (length(significant_pvals) > 0){
    DiffTaxa <- names(significant_pvals)
  } else{
    DiffTaxa <- c()
  }
  cat(sprintf("%d differentially abundant taxa detected\n", length(DiffTaxa)))
  if(is.null(DiffTaxa)) DiffTaxa <- ""
  output <- new("DAresult",
                DAAname=DAAname,
                reference=reftaxa, 
                signal=DiffTaxa,
                details=all_CR_results,
                input=input_data)

  invisible(output)
}
