


#' @importFrom phyloseq filter_taxa prune_samples otu_table sample_data taxa_are_rows sample_sums
#' @importFrom stats model.matrix
preprocess <- function(input_data, cond.var, base.cond, adj.var, subject.var, prev.filter, depth.filter){
  
  # check if input data type is phyloseq
  stopifnot("Input data isn't a phyloseq object!" = is(input_data, 'phyloseq'))
  # check the prevalence threshold and the sequencing depth threshold
  stopifnot("The prevalence filter (prev.filter) has to be a positive number between 0 and 1!" = prev.filter >=0 & prev.filter < 1)
  stopifnot("The sequencing depth filter (depth.filter) has to be a nonnegative number!" = depth.filter >= 0)
  # filter phyloseq object based on taxa prevalence and sequencing depth
  subset_data <- filter_taxa(input_data, function(x) mean(x>0) > prev.filter, TRUE)
  subset_data <- prune_samples(sample_sums(subset_data) > depth.filter, subset_data)
  
  # check if the variables exist in the metadata
  metadata <- data.frame(sample_data(subset_data))
  allcols <- colnames(metadata)
  stopifnot("The main variable name for conditions need to be a string." = 
              is(cond.var, "character") && length(cond.var) == 1)
  stopifnot("The subject variable name needs to be a string." =
              is(subject.var, "character") && length(subject.var) == 1)

  stopifnot("The variables for adjustments should be either NULL or a vector of character strings." = 
              is(adj.var, "character") | is(adj.var, "NULL"))

  if (subject.var == cond.var || subject.var %in% adj.var) {
    stop("The subject variable cannot also be the condition or an adjustment variable!")
  }

  selected_cols <- unique(c(cond.var, adj.var, subject.var))
  if (!all(selected_cols %in% allcols)){
    unavailable_cols <- selected_cols[!selected_cols %in% allcols]
    stop(sprintf("Some columns are not available in the metadata! (%s)", 
                 paste(unavailable_cols, collapse=",")))
  }
  subset_metadata <- metadata[, selected_cols, drop=FALSE]
  
  # dichotomize categorical variables if selected, set up design matrix
  if (any(is.na(subset_metadata))){
    stop("No missing data allowed in the metadata!")
  }
  main_variable <- subset_metadata[, cond.var]
  if (length(unique(main_variable)) == 1){
    stop("All samples share the same condition!")
  }
  if (!is.numeric(main_variable)){
    if (is.null(base.cond)){
      base.cond <- unique(main_variable)[1]
    } else if (!base.cond %in% main_variable) {
      base.cond <- unique(main_variable)[1]
    }
    cat(sprintf("Choose '%s' as the baseline condition\n", base.cond))
    if (length(unique(main_variable)) == 2){
      others <- setdiff(main_variable, base.cond)
    } else{
      others <- "others"
    }
    main_variable <- as.integer(main_variable != base.cond)
    cond.var <- sprintf("%s (%s VS %s)", cond.var, others, base.cond)
  }
  adjustments <- NULL
  if (!is.null(adj.var)){
    adjustments <- subset_metadata[, adj.var, drop=FALSE]
    adjustments<- model.matrix(~., data=adjustments)
    adjustments<- adjustments[, -1] # remove intercept
  }
  complete_design_matrix <- cbind(1, main_variable, adjustments)
  colnames(complete_design_matrix)[1:2] <- c("(Intercept)", "condition")
  subject <- factor(subset_metadata[, subject.var])
  if (nlevels(subject) < 2) {
    stop("The subject variable needs at least two unique subjects!")
  }
  if (!any(tabulate(subject) > 1)) {
    stop("At least one subject must have repeated observations for longitudinal analysis!")
  }
  
  # parse the count matrix and edit the count_ratio function
  count_table <- otu_table(subset_data)
  if (taxa_are_rows(subset_data)){
    count_table <- t(count_table)
  }
  count_table <- as(count_table, "matrix")
  
  if (any(is.na(count_table))){
    stop("No missing data allowed in the count table!")
  }
  
  cat(sprintf("%d taxa and %d samples being analyzed...\n", 
              ncol(count_table), nrow(count_table)))
  
  output <- list(count_table=count_table, design_matrix=complete_design_matrix,
                 DAAname = cond.var, subject=subject, clean_metadata = subset_metadata)
  
  return(output)
  
}
