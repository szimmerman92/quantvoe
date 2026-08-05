#' Run initial association 
#'
#' Run initial association for a single feature
#' @param j dependent feature name
#' @param independent_variables A tibble containing the information for your independent variables (e.g. bacteria relative abundance, age). The columns should correspond to different variables, and the rows should correspond to different units,  (e.g. individual1, individual2, etc). If passing multiple datasets, pass a named list of the same length and in the same order as the dependent_variables parameter. If running from the command line, pass the path (or comma separated paths) to .rds files containing the data, one per dataset.
#' @param dependent_variables A tibble containing the information for your dependent variables (e.g. bacteria relative abundance, age). The columns should correspond to different variables, and the rows should correspond to different units, such as individuals (e.g. individual1, individual2, etc). If passing multiple datasets, pass a named list of the same length and in the same order as the independent_variables parameter. If running from the command line, pass the path (or comma separated paths) to .rds files containing the data, one per dataset.
#' @param primary_variable The column name from the independent_variables tibble containing the key variable you want to associate with disease in your first round of modeling (prior to vibration). For example, if you are interested fundamentally identifying how well age can predict height, you would make this value a string referring to whatever column in said dataframe refers to "age."
#' @param constant_adjusters A character vector (or just one string) corresponding to column names in your dataset to include in every vibration. (default = NULL)
#' @param model_type Specifies regression type -- "glm", "survey", or "negative_binomial". Survey regression will require additional parameters (at least weight, nest, strata, and ids). Any model family (e.g. gaussian()), or any other parameter can be passed as the family argument to this function.
#' @param proportion_cutoff Float between 0 and 1. Filter out dependent features that are this proportion of zeros or more (default = 1, so no filtering done.)
#' @param family GLM family (default = gaussian()). For help see help(glm) or help(family).
#' @param ids Name of column in dataframe specifying cluster ids from largest level to smallest level. Only relevant for survey data. (Default = NULL).
#' @param strata Name of column in dataframe with strata. Relevant for survey data. (Default = NULL).
#' @param weights Name of column containing sampling weights.
#' @param nest If TRUE, relabel cluster ids to enforce nesting within strata.
#' @param num_knots Numeric vector specifying the values of the primary variable of where the knots should be. If length of num_knots is greater than 0 generate the B-spline basis matrix for a natural cubic spline on primary_variable
#' @param spline_type if num_knots is 0 then this is NULL. Otherwise define the spline function as a string. Options are c("bs","pspline","ns")
#' @param quantile_bounds if your primary variable is a factor corresponding to quantiles set this, otherwise default is NULL.
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @keywords regression, initial association
regression <- function(j,independent_variables,dependent_variables,primary_variable,constant_adjusters,model_type,proportion_cutoff,family,ids,strata,weights,nest,num_knots,spline_type,quantile_bounds){
  # print("In regression()")
  # print(weights)
  feature_name = colnames(dependent_variables)[j+1]
  if(primary_variable %in% colnames(dependent_variables %>% dplyr::select(.data$sampleID,c(feature_name)))){
    return('Regression not run due to presence of primary variable in dependent variable dataframe')
  }
  regression_df = suppressMessages(dplyr::left_join(dependent_variables %>% dplyr::select(.data$sampleID,c(feature_name)),independent_variables)) #%>% dplyr::mutate_if(is.factor, as.character)) %>% dplyr::mutate_if(is.character, as.factor))
  regression_df = regression_df %>% dplyr::select(-.data$sampleID)
  #final check to confirm ready for regressions
  
  if(length(num_knots) > 1 | num_knots[1] > 0) {

    primary_variable_treated = paste("splines::", spline_type, "(",primary_variable,",knots=",ifelse(length(num_knots) >0, paste("c(", paste(num_knots,collapse=","), ")", sep = ""), num_knots),")",sep="") # add this code
    
  } else {
    primary_variable_treated = primary_variable
  }
  #run regression
  if(!is.null(constant_adjusters)){
    primary_variable_formodel = paste(paste(constant_adjusters,sep='+',collapse='+'),'+',primary_variable_treated)
  } 
  else if(is.null(constant_adjusters)) {
    primary_variable_formodel = primary_variable_treated
  }
  
  # if dependent variable is of type survival.
  if(class(dependent_variables[[feature_name]])[1]=="Surv" & model_type != "survey") { 
    # we will do cox regression
    myformula = stats::as.formula(paste(feature_name, "~ ",primary_variable_formodel))
    # print(myformula)
    weights = regression_df %>% dplyr::select(weights) %>% unlist %>% unname
    #regression_model = coxph(formula=myformula,weights=weights,data=regression_df)
    regression_model = tryCatch({
      coxph(formula = myformula, weights = weights, data = regression_df)
    }, error = function(e) {
      NULL
    })

    # return(regression_model)
    if(is.null(regression_model)) {
      return(NULL)
    }
    # print(summary(regression_model))
    #return(broom::tidy(coxph(formula=myformula,weights=weights,data=regression_df)) %>% dplyr::mutate(feature=feature_name))
    #return(tryCatch(broom::tidy(coxph(formula=myformula,weights=weights,data=regression_df)) %>% dplyr::mutate(feature=feature_name),
                   #warning = function(w) w,
                    #error = function(e) e
    #))
  }
  
  
  if(model_type=='negative_binomial'){
    #regression_model = MASS::glm.nb(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)

    regression_model = tryCatch({
      MASS::glm.nb(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)
    }, error = function(e) {
      NULL
    })
    
    # return(regression_model)
    if(is.null(regression_model)) {
      return(NULL)
    }
    
    #return(tryCatch(broom::tidy(MASS::glm.nb(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)) %>% dplyr::mutate(feature=feature_name),
    #                warning = function(w) w, 
    #                error = function(e) e
    #))
  }
  if(model_type=='survey'){
    options(survey.lonely.psu="adjust")
    dsn=survey::svydesign(weights = regression_df %>% 
                            dplyr::pull(weights) 
                          , ids = regression_df %>%
                            dplyr::pull(ids)
                          , nest = as.logical(nest)
                          , strata = regression_df %>% 
                            dplyr::pull(strata)
                          , data = regression_df)
    # View(regression_df)
    if(class(dependent_variables[[feature_name]])[1]=="Surv") {
      # print("I'm in survival and survey")
      myformula = stats::as.formula(paste(feature_name, "~ ",primary_variable_formodel))
      
      # print(nrow(regression_df))
      # myformula = stats::as.formula("cancer_death_status_surv ~ splines::bs(LBXCOT, df = 6)")
      # print(myformula)
   
      # X <- splines::bs(dsn$variables$LBXCOT, df = 6)
      # 
      # # Interior knots (the ones bs() chose internally)
      # print(attr(X, "knots"))
      # 
      # # Boundary knots (min/max used, or whatever you specified via Boundary.knots)
      # print(attr(X, "Boundary.knots"))
      # 
      # # Degree (default cubic = 3)
      # print(attr(X, "degree"))
      
      # X <- splines::bs(dsn$variables$LBXCOT, knots = c(0.011, 0.026, 0.07, 3.23400000000001))
      # View(apply(X, 2, sd))          # look for a column with ~0 variance -> that's your #7
      # print(range(dsn$variables$LBXCOT))  # compare to the outer knot, 3.234
      # print(sum(dsn$variables$LBXCOT > 3.234))  # if this is 0 or tiny, that's the smoking gun
      
      
      # df <- dsn$variables
      # print(colnames(df))
      
      # df_counts <- df %>%
      #   filter(LBXCOT > 3.23400000000001) %>%
      #   count(cancer_death_status_surv)
      # print(df_counts)
       
      # View(regression_df)
      
      # X <- model.matrix(myformula, data = dsn$variables)
      # View(X)
      # View(cor(X))

      # print(qr(X)$rank)
      # print(ncol(X))
      # hist(X[,ncol(X)])
      # qr_X <- qr(X)
      # # Columns qr_X$pivot gives you the order; the last (ncol(X) - rank) 
      # # pivoted columns are the "extra" ones not contributing new rank
      # rank <- qr_X$rank
      # redundant_cols <- qr_X$pivot[(rank + 1):ncol(X)]
      # print(colnames(X)[redundant_cols])
      
      
      
      # To view the matrix with all variables in the regression formula along with the dummy columns for the categorical variables
      # View(model.matrix(myformula, regression_df))
 
      # View(dsn$variables)
      #regression_model = survey::svycoxph(formula=myformula,design=dsn)
      
      regression_model <- tryCatch({
        survey::svycoxph(formula=myformula,design=dsn)
       }, error = function(e) {
         NULL
       })
      # View(regression_model)
      # return(regression_model)
      if(is.null(regression_model)) {
        return(NULL)
      }
      
      # print(regression_model)
      #return(tryCatch(broom::tidy(survey::svycoxph(formula=myformula,design=dsn)) %>% dplyr::mutate(feature=feature_name),
      #                warning = function(w) w,
      #                error = function(e) e
      #))
    } else {
      #regression_model = survey::svyglm(family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),design=dsn)
      regression_model = tryCatch({
        survey::svyglm(family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),design=dsn)
      }, error = function(e) {
        NULL
      })
      
      # return(regression_model)
      if(is.null(regression_model)) {
        return(NULL)
      }
      
      #return(tryCatch(broom::tidy(survey::svyglm(family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),design=dsn)) %>% dplyr::mutate(feature=feature_name),
      #                warning = function(w) w,
      #                error = function(e) e
      #)) 
    }
  }
  if(model_type=='glm'){
    #regression_model = stats::glm(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)
    regression_model = tryCatch({
      stats::glm(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)
    }, error = function(e) {
      NULL
    })
    
    # return(regression_model)
    if(is.null(regression_model)) {
      return(NULL)
    }
    
    #return(tryCatch(broom::tidy(stats::glm(weights=regression_df %>% dplyr::select(weights) %>% unlist %>% unname,family=family,formula=stats::as.formula(stringr::str_c("I(`", feature_name,"`) ~ ",primary_variable_formodel)),data = regression_df)) %>% dplyr::mutate(feature=feature_name),
    #                warning = function(w) w,
    #                error = function(e) e
    #))
  }
  # print(primary_variable)
  
  # View(regression_model)
  
  # print("Right before termplot")
  
  if(model_type == "survey" & length(num_knots) > 0)
  {
    # print("Running manual_termplot")
    
    df_predicted_risk_primary_variable <- termplot(regression_model
             , se = TRUE
             , data   = dsn$variables 
             , plot = FALSE)[[primary_variable]] 
    
    # df_predicted_risk_primary_variable <- manual_termplot(
    #   model = regression_model,
    #   data   = dsn$variables,
    #   var    = primary_variable,
    #   knots  = num_knots,
    #   boundary_knots = range(dsn$variables %>%
    #                            pull(primary_variable), na.rm = TRUE),
    #   use_observed_x = TRUE,
    #   center_method  = "mean_basis"
    # )[[primary_variable]]
      
  } else {
    
    df_predicted_risk_primary_variable <- termplot(regression_model
                                                   , se = TRUE
                                                   , plot = FALSE)[[primary_variable]]
    
  }
  # print(names(termplot(regression_model, se = TRUE, plot = FALSE)))
  # View(tidy(regression_model))
  # View(glance(regression_model))
  # View(df_predicted_risk_primary_variable)
  
  data_type_primary_variable <- class(df_predicted_risk_primary_variable$x)
  
  # print(quantile_bounds)
  if(is.null(quantile_bounds) == TRUE)
  {
    # print("is.null(quantile_bounds)")
    min_primary_variable = min(df_predicted_risk_primary_variable$x, na.rm = TRUE)
    # print(min_primary_variable)
    center = with(df_predicted_risk_primary_variable, y[df_predicted_risk_primary_variable$x == min_primary_variable])
    # View(summary(regression_model))
    degree_of_freedom = summary(regression_model)$logtest["df"] %>% as.numeric(.)
    df_predicted_risk_primary_variable <- df_predicted_risk_primary_variable %>%
      mutate(y_centered = y - center) %>%
      mutate(ci_lower_y = y_centered - 1.96*se) %>%
      mutate(ci_upper_y = y_centered + 1.96*se) %>%
      mutate(statistic = (y_centered/se)) %>%
      mutate(p.value = 2*pt(q = -abs(statistic), df = degree_of_freedom)) %>% 
      mutate(primary_variable_name=primary_variable) %>% 
      mutate(feature=feature_name)
    
  } else {
    
    print("compute associations to quantile models")
    # print(quantile_bounds)
    lower_bound_quantiles <- quantile_bounds[1:length(quantile_bounds) - 1]
    # print(lower_bound_quantiles)
    upper_bound_quantiles <- quantile_bounds[2:length(quantile_bounds)]
    
    min_primary_variable <- "_1"
    center = with(df_predicted_risk_primary_variable, y[df_predicted_risk_primary_variable$x == min_primary_variable])
    # View(summary(regression_model))
    degree_of_freedom = summary(regression_model)$logtest["df"] %>% as.numeric(.)
    # View(df_predicted_risk_primary_variable)
    
    df_predicted_risk_primary_variable <- df_predicted_risk_primary_variable %>%
      mutate(lower_bound_quantiles = lower_bound_quantiles) %>%
      mutate(upper_bound_quantiles = upper_bound_quantiles) %>%
      mutate(y_centered = y - center) %>%
      mutate(ci_lower_y = y_centered - 1.96*se) %>%
      mutate(ci_upper_y = y_centered + 1.96*se) %>%
      mutate(statistic = (y_centered/se)) %>%
      mutate(p.value = 2*pt(q = -abs(statistic), df = degree_of_freedom)) %>% 
      mutate(primary_variable_name=primary_variable) %>% 
      mutate(feature=feature_name)

  }
  
 
  # View(df_predicted_risk_primary_variable)
  
  if(class(dependent_variables[[feature_name]])[1]=="Surv")
  {
    df_predicted_risk_primary_variable <- df_predicted_risk_primary_variable %>%
      mutate(effect_size = exp(y_centered)
             , effect_size_ci_lower = exp(ci_lower_y)
             , effect_size_ci_upper = exp(ci_upper_y))
  } else {
    df_predicted_risk_primary_variable <- df_predicted_risk_primary_variable %>%
      mutate(effect_size = y_centered
             , effect_size_ci_lower = ci_lower_y
             , effect_size_ci_upper = ci_upper_y)
  }
  regression_model_df_stats = broom::tidy(regression_model) %>% dplyr::mutate(feature=feature_name) %>% dplyr::mutate(primary_variable_name=primary_variable)
  return(list(regression_model_df_stats,df_predicted_risk_primary_variable))
}

#' Run all associations for a given dataset
#'
#' Function to run all associations for dependent and indepenent features.
#' @param x merged independent an dependent data for a given dataset
#' @param primary_variable The column name from the independent_variables tibble containing the key variable you want to associate with disease in your first round of modeling (prior to vibration). For example, if you are interested fundamentally identifying how well age can predict height, you would make this value a string referring to whatever column in said dataframe refers to "age."
#' @param constant_adjusters A character vector (or just one string) of column names corresponding to column names in your dataset to include in every vibration. (default = NULL)
#' @param vibrate TRUE/FALSE -- run vibrations (default=TRUE)
#' @param model_type Specifies regression type -- "glm", "survey", or "negative_binomial". Survey regression will require additional parameters (at least weight, nest, strata, and ids). Any model family (e.g. gaussian()), or any other parameter can be passed as the family argument to this function.
#' @param proportion_cutoff Float between 0 and 1. Filter out dependent features that are this proportion of zeros or more (default = 1, so no filtering will be done.)
#' @param family GLM family (default = gaussian()). For help see help(glm) or help(family).
#' @param ids Name of column in dataframe specifying cluster ids from largest level to smallest level. Only relevant for survey data. (Default = NULL).
#' @param strata Name of column in dataframe with strata. Relevant for survey data. (Default = NULL).
#' @param weights Name of column containing sampling weights.
#' @param nest If TRUE, relabel cluster ids to enforce nesting within strata.
#' @param num_knots If num_knots is greater than 0 generate the B-spline basis matrix for a natural cubic spline on primary_variable
#' @param spline_type if num_knots is 0 then this is NULL. Otherwise define the spline function as a string. Options are c("bs","pspline","ns")
#' @param quantile_bounds if your primary variable is a factor corresponding to quantiles set this, otherwise default is NULL.
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @keywords regression, initial association
run_associations <- function(x,primary_variable,constant_adjusters,model_type,proportion_cutoff,vibrate,family,ids,strata,weights,nest,num_knots,spline_type,quantile_bounds){
  dependent_variables <- dplyr::as_tibble(x[[1]])
  classes = dependent_variables %>% summarise_all(class)
  classes = as.character(classes)
  colnames(dependent_variables)[[1]]='sampleID'
  # get the column names that are of survival class
  surv_vars = colnames(dependent_variables)[classes=="Surv"]
  toremove = which(colSums(dependent_variables %>% dplyr::select(-.data$sampleID,-all_of(surv_vars)) == 0,na.rm=TRUE)/nrow(dependent_variables)>proportion_cutoff)
  print(paste("Removing",length(toremove),"features that are at least",proportion_cutoff*100,"percent zero values."))
  dependent_variables=dependent_variables %>% dplyr::select(-(toremove+1))
  if(ncol(dependent_variables)==1){
    print('After filtering your data, you had nothing left. Try changing your filtering threshold for zero-value data and running again.')
    #quit()
    return(NULL)
  }
  independent_variables <- dplyr::as_tibble(x[[2]])
  # print(independent_variables)
  
  colnames(independent_variables)[[1]]='sampleID'
  print(paste('Computing',as.character(ncol(dependent_variables)-1),'associations for dataset',as.character(unname(unlist(x[[3]])))))
  colnames(dependent_variables)[1]='sampleID'
  colnames(independent_variables)[1]='sampleID'
  tokeep = independent_variables %>% dplyr::select_if(~ length(unique(.)) > 1) %>% colnames
  todrop = setdiff(colnames(independent_variables),tokeep)
  
  # print("Right here")
  # print(tokeep)
  if(length(todrop)>1){
    print('Note -- The following variables lack multiple levels and will be dropped should you run a vibration:')
    print(todrop)
    if(primary_variable %in% todrop){
      print('One of the variables being dropped is your variable of interethis will result in the pipeline failing. Please adjust your independent variables and try again.')
      return(NULL)
      #quit()
    }
  }
  independent_variables=independent_variables %>% dplyr::select(-tidyselect::all_of(todrop))
  constant_adjusters = setdiff(constant_adjusters,todrop)
  if(ncol(independent_variables)==2){
    vibrate=FALSE
  }
  overlap = intersect(colnames(independent_variables %>% dplyr::select(-.data$sampleID)),colnames(dependent_variables %>% dplyr::select(-.data$sampleID)))
  
  if(length(overlap)>0){
    print('The following variables are in both the dependent and independent datasets. This may cause some some regressions to fail, though the pipeline will still run to completion.')
    print(overlap)
  }
  out = purrr::map(seq_along(dependent_variables %>% dplyr::select(-.data$sampleID)), function(j) regression(j,independent_variables,dependent_variables,primary_variable,constant_adjusters,model_type,proportion_cutoff,family,ids,strata,weights,nest,num_knots,spline_type,quantile_bounds))
  out_success = out[unlist(purrr::map(out,function(x) class(x)=="list"))]
  if(length(out_success)!=length(out)){
    print(paste('Dropping',length(out)-length(out_success),'features with regressions that failed to converge.'))
  }
  if(length(out_success)==0){
    print(paste("All of your regression output failed. Printing error messages to screen."))
    Sys.sleep(3)
    print(out)
    return(NULL)
    #quit()
  }
  out_success_model_out = purrr::map(out_success, function(x) x[[1]])
  out_success_model_out_termplot = purrr::map(out_success, function(x) x[[2]])
  out_success_model_out_termplot = dplyr::bind_rows(out_success_model_out_termplot)
  
  out_success_model_out = out_success_model_out %>% dplyr::bind_rows() %>% dplyr::filter(grepl(primary_variable,.data$term)==TRUE) %>% group_by(term) %>% dplyr::mutate( bonferroni = stats::p.adjust(.data$p.value, method = "bonferroni"), BH = stats::p.adjust(.data$p.value, method = "BH"), BY = stats::p.adjust(.data$p.value, method = "BY"))
  # out_success = out_success %>% dplyr::mutate(dataset_id=x[[3]])
  return(list('output' = out_success_model_out,'vibrate' = vibrate,'termplot'=out_success_model_out_termplot))
}

#' Deploy associations across datasets
#'
#' Top-level function to run all associations for all datasets.
#' @param bound_data merged independent an dependent data for all datasets
#' @param primary_variable The column name from the independent_variables tibble containing the key variable you want to associate with disease in your first round of modeling (prior to vibration). For example, if you are interested fundamentally identifying how well age can predict height, you would make this value a string referring to whatever column in said dataframe refers to "age."
#' @param constant_adjusters A character vector (or just one string) of column names corresponding to column names in your dataset to include in every vibration. (default = NULL)
#' @param vibrate TRUE/FALSE -- run vibrations (default=TRUE)
#' @param model_type Specifies regression type -- "glm", "survey", or "negative_binomial". Survey regression will require additional parameters (at leaset weight, nest, strata, and ids). Any model family (e.g. gaussian()), or any other parameter can be passed as an additional argument to this function.
#' @param proportion_cutoff Float between 0 and 1. Filter out dependent features that are this proportion of zeros or more (default = 1, so no filtering done.)
#' @param family GLM family (default = gaussian()). For help see help(glm) or help(family).
#' @param ids Name of column in dataframe specifying cluster ids from largest level to smallest level. Only relevant for survey data. (Default = NULL).
#' @param strata Name of column in dataframe with strata. Relevant for survey data. (Default = NULL).
#' @param weights Name of column containing sampling weights.
#' @param nest If TRUE, relabel cluster ids to enforce nesting within strata.
#' @param num_knots If num_knots is greater than 0 generate the B-spline basis matrix for a natural cubic spline on primary_variable
#' @param spline_type if num_knots is 0 then this is NULL. Otherwise define the spline function as a string. Options are c("bs","pspline","ns")
#' @param quantile_bounds if your primary variable is a factor corresponding to quantiles set this, otherwise default is NULL.
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @keywords regression, initial association
#' @export
compute_initial_associations <- function(bound_data,primary_variable, constant_adjusters = NULL,model_type = 'glm', proportion_cutoff = 1,vibrate = TRUE,family = gaussian(),ids = NULL,strata =NULL,weights =NULL,nest = NULL,num_knots=0,spline_type=NULL,quantile_bounds){
  # print("In compute_initial_associations")
  # print(weights)
  
  output = apply(bound_data, 1, function(x) run_associations(x,primary_variable,constant_adjusters,model_type,proportion_cutoff,vibrate,family,ids,strata,weights,nest,num_knots,spline_type,quantile_bounds))
  # View(output)
  output_regs = purrr::map(output, function(x) x[[1]])
  output_vib = unlist(unname(unique(purrr::map(output, function(x) x[[2]]))))
  output_term = purrr::map(output, function(x) x[[3]])
  if(FALSE %in% output_vib & vibrate!=FALSE){
    output_vib=FALSE
    print('For at least one dataset, we dropped all the variables that you could possible vibrate over due to lacking multiple levels. Vibrate parameter being set to FALSE.')
  }
  output_regs = dplyr::bind_rows(output_regs)
  output_term = dplyr::bind_rows(output_term)
  return(list('output'=output_regs,'vibrate'=output_vib,'termplot'=output_term))
}
