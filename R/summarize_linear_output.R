
#' Unnest vibration data
#'
#' Unnest regression output for vibrations.
#' @param vib_df VoE dataframe with columns for each adjuster (output of get_adjuster_expanded_vibrations).
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
#' @keywords voe analysis
filter_unnest_feature_vib <- function(vib_df) {
  return(vib_df %>% dplyr::slice(which(purrr::map_lgl(vib_df$feature_fit, ~class(.)[[1]] == "tbl_df"))) %>% tidyr::unnest(.data$feature_fit))
}

#' Expand vibration data
#'
#' Add column to vibration output for each adjuster (indicating presence or absence in a model).
#' @param voe_df Raw vibration output, the first entry in the output of compute_vibrations.
#' @param adjusters A dataframe, each column corresponding to the adjusters used in each dataset for vibrations. This is the second entry in the compute_vibrations output. 
#' @param constant_adjusters A character vector (or just one string) corresponding to column names in your dataset to include in every vibration. (default = NULL)
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
#' @keywords voe analysis
get_adjuster_expanded_vibrations <- function(voe_df, adjusters, constant_adjusters) {
  copy_voe_df <- rlang::duplicate(voe_df, shallow = FALSE)
  adjusters= unique(unlist(unname(purrr::map(adjusters, function(x) unlist(x)))))
  for (variable in adjusters) {
    copy_voe_df  = copy_voe_df %>% dplyr::mutate(newcol = purrr::map_int(copy_voe_df$vars, ~(variable %in% .)))
    colnames(copy_voe_df)[length(colnames(copy_voe_df))] <- variable
  }
  copy_voe_df = copy_voe_df %>% dplyr::select(-tidyselect::all_of(constant_adjusters))
  return(copy_voe_df)
}

#' Find confounders
#'
#' Model confounding from vibration analysis.
#' @param voe_list_for_reg A dataframe of expanded VoE output (output of filter_unnest_feature_vib)
#' @keywords voe analysis
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @export
find_confounders_linear <- function(voe_list_for_reg){
  trylinear=FALSE
  ptype=unique(voe_list_for_reg$term)
  voe_adjust_for_reg_ptype <- voe_list_for_reg %>% dplyr::select(-.data$dataset_id,-.data$independent_feature) %>% dplyr::select_if(~ length(unique(.)) > 1) %>% dplyr::select(-c(.data$full_fits,.data$std.error,.data$statistic,.data$termplot_fit))
  voe_adjust_for_reg_ptype$estimate = abs(voe_adjust_for_reg_ptype$estimate)
  if('dependent_feature' %in% colnames(voe_adjust_for_reg_ptype)){ 
    if(!(1 %in% unique(unlist(unname(table(voe_adjust_for_reg_ptype$dependent_feature)))))){
      tryCatch({
        if("term"%in% colnames(voe_adjust_for_reg_ptype)) {
          voe_adjust_for_reg_ptype_split = split(voe_adjust_for_reg_ptype,voe_adjust_for_reg_ptype$term)
          fit_estimate = purrr::map(voe_adjust_for_reg_ptype_split, ~ lme4::lmer(data=.x %>% dplyr::select_if(~ length(unique(.)) > 1),stats::as.formula(estimate ~ . +(1|dependent_feature) - dependent_feature - estimate - p.value),control = lme4::lmerControl(optimizer = "bobyqa")))
          fullmod=fit_estimate
          fit_estimate_forplot = purrr::map(fit_estimate, ~ broom.mixed::tidy(.x) %>% dplyr::mutate(sdmin=(.data$estimate - .data$std.error),sdmax=(.data$estimate + .data$std.error)))
        } else {
          fit_estimate=lme4::lmer(data=voe_adjust_for_reg_ptype,stats::as.formula(estimate ~ . +(1|dependent_feature) - dependent_feature - estimate - p.value),control = lme4::lmerControl(optimizer = "bobyqa"))
          fullmod=fit_estimate
          fit_estimate_forplot=broom.mixed::tidy(fit_estimate) %>% dplyr::mutate(sdmin=(.data$estimate-.data$std.error),sdmax=(.data$estimate+.data$std.error))
        }
      },
      error = function(e){
        print('Note: Mixed effect modeling to identify sources of confounding failed. Running a simple linear model instead. If you want to try this analysis yourself, you can access the raw data for this yourself in the output and follow the methodological layout in the docs.')
        trylinear <<- TRUE
      })
    }
    if(trylinear==TRUE){
      tryCatch({
        if("term"%in% colnames(voe_adjust_for_reg_ptype)) {
          voe_adjust_for_reg_ptype_split = split(voe_adjust_for_reg_ptype,voe_adjust_for_reg_ptype$term)
          fit_estimate = purrr::map(voe_adjust_for_reg_ptype_split, ~ stats::lm(data=.x %>% dplyr::select_if(~ length(unique(.)) > 1),stats::as.formula(estimate ~ . - estimate - dependent_feature - p.value)))
          fullmod=fit_estimate
          fit_estimate_forplot = purrr::map(fit_estimate, ~ broom::tidy(.x) %>% dplyr::mutate(sdmin=(.data$estimate - .data$std.error),sdmax=(.data$estimate + .data$std.error)))
        } else {
          fit_estimate=stats::lm(data=voe_adjust_for_reg_ptype,stats::as.formula(estimate ~ . - estimate - dependent_feature - p.value))
          fullmod=fit_estimate
          fit_estimate_forplot=broom::tidy(fit_estimate) %>% dplyr::mutate(sdmin=(.data$estimate - .data$std.error),sdmax=(.data$estimate + .data$std.error))
        }
      },
      error = function(e){
        fit_estimate_forplot = 'Confounder analysis failed.'
        fullmod='Confounder analysis failed.'
        print('Confounder analysis failed. We recommend looking at the raw vibration output to see what the issue may be.')
      })
    }
  } else{
    tryCatch({
      print('Note: Using regular linear model for confounder analysis instead of a mixed effect one. See the GitHub README for more details.')
      if("term"%in% colnames(voe_adjust_for_reg_ptype)) {
        voe_adjust_for_reg_ptype_split = split(voe_adjust_for_reg_ptype
                                               , voe_adjust_for_reg_ptype$term)
        
        fit_estimate = lapply(voe_adjust_for_reg_ptype_split
                              , function(term_df) {
                                term_df = term_df[,-1]
                                if(length(unique(term_df))>1 
                                   & sum(!is.na(term_df$estimate)) > 0) 
                                {
                                  fit_res_temp = lm(estimate ~ . - estimate - p.value,data = term_df) 
                                }
                                })
        fit_estimate = fit_estimate[!sapply(fit_estimate, is.null)]
        # fit_estimate = purrr::map(voe_adjust_for_reg_ptype_split
        #                           , ~ stats::lm(data=.x %>% 
        #                                           dplyr::select_if(~ (length(unique(.)) > 1 & sum(!is.na(.$estimate)) > 0))
        #                                         ,stats::as.formula(estimate ~ . - estimate - p.value)))
   
        fullmod=fit_estimate
        fit_estimate_forplot = purrr::map(fit_estimate, ~ broom::tidy(.x) %>% dplyr::mutate(sdmin=(.data$estimate - .data$std.error),sdmax=(.data$estimate + .data$std.error)))
      } else {
        fit_estimate=stats::lm(data=voe_adjust_for_reg_ptype,stats::as.formula(estimate ~ . - estimate - p.value))
        fullmod=fit_estimate
        fit_estimate_forplot=broom::tidy(fit_estimate) %>% dplyr::mutate(sdmin=(.data$estimate - .data$std.error),sdmax=(.data$estimate + .data$std.error))
      }
    },
    error = function(e){
      fit_estimate_forplot = 'Confounder analysis failed.'
      fullmod='Confounder analysis failed.'
      print('Confounder analysis failed. We recommend looking at the raw vibration output to see what the issue may be.')
    })
  }
  return(list(summarized_output = fit_estimate_forplot, full_model = fullmod))
}

#' summarize_vibration_data_by_feature_termplot
#'
#' Summarize output of vibrations for each dependent feature of interest for non linear voe.
#' @param df A dataframe of expanded VoE output (output of filter_unnest_feature_vib)
#' @param center_for_effect_size determines how to calculate janus effect. 0 if values are linear, 1 if values are hazard ratio or risk ratios or fold differences 

#' @keywords voe analysis
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
summarize_vibration_data_by_feature_termplot <- function(df,center_for_effect_size){
  # View(df)
  summarized_voe_data <- df %>% 
    group_by(independent_feature,dependent_feature,x) %>%
    summarise(estimate_quantile_1 = quantile(effect_size, probs = 0.01, na.rm = TRUE)
              , estimate_quantile_50 = quantile(effect_size, probs = 0.5, na.rm = TRUE)
              , estimate_quantile_99 = quantile(effect_size, probs = 0.99, na.rm = TRUE)
              , estimate_diff_99_1 = estimate_quantile_99 - estimate_quantile_1
              , num_models = sum(!is.na(effect_size))
              , janus_effect = sum(effect_size>center_for_effect_size,na.rm=TRUE)/num_models
              , pval_quantile_1 = quantile(p.value, probs = 0.01, na.rm = TRUE)
              , pval_quantile_50 = quantile(p.value, probs = 0.5, na.rm = TRUE)
              , pval_quantile_99 = quantile(p.value, probs = 0.99, na.rm = TRUE)
              , pvalue_diff_99_1 = pval_quantile_99 - pval_quantile_1
              , estimate_mean_over_vibration = mean(effect_size, na.rm = TRUE)
              , estimate_ci_lower_over_vibration = mean(effect_size_ci_lower, na.rm = TRUE)
              , estimate_ci_upper_over_vibration = mean(effect_size_ci_upper, na.rm = TRUE)
              , estimate_se_over_vibration = (mean(ci_upper_y, na.rm = TRUE) - mean(ci_lower_y, na.rm = TRUE))/(2*1.96)
              , t_statistic_estimate_over_vibration = mean(y_centered, na.rm = TRUE)/estimate_se_over_vibration
              , crosses_the_center = ifelse(center_for_effect_size > estimate_ci_lower_over_vibration 
                                            & center_for_effect_size < estimate_ci_upper_over_vibration
                                            , "yes"
                                            , "no")
              # , p.value_over_vibration = 2*pt(q = (t_statistic_estimate_over_vibration)
              #                                 , df = num_models - 1)
              , p.value_over_vibration =  exp(-0.717*(abs(t_statistic_estimate_over_vibration)) - 0.416*abs(t_statistic_estimate_over_vibration)^2)
              # , exact_results = fisher.test(x=factor(as.numeric(effect_size>center_for_effect_size),levels=c(0,1))
              #                               , y=factor(as.numeric(p.value < 0.05),levels=c(0,1)))$p.value
              ) %>%
    ungroup(.)
  
  if("lower_bound_quantiles" %in% colnames(df)) {
    summarized_voe_data = merge(summarized_voe_data
                                , df %>% 
                                  select(x,lower_bound_quantiles,upper_bound_quantiles) %>%
                                  unique(.)
                                , by="x")
  }
  
    # contigency_table = df %>% group_by(independent_feature,dependent_feature,x) %>%
    #   summarise(exact_results = fisher.test(x=factor(effect_size>center_for_effect_size,levels=c(0,1)), y=factor(p.value < 0.05,levels=c(0,1)))$p.value)
    # 
    # contigency_table = df %>% group_by(independent_feature,dependent_feature,x) %>%
    #   summarise(above_center_sig = sum(effect_size>center_for_effect_size & p.value < 0.05,na.rm=TRUE),
    #             above_center_not_sig = sum(effect_size>center_for_effect_size & p.value >= 0.05,na.rm=TRUE),
    #             below_center_sig = sum(effect_size<center_for_effect_size & p.value < 0.05,na.rm=TRUE),
    #             below_center_not_sig = sum(effect_size<center_for_effect_size & p.value >= 0.05,na.rm=TRUE))
    # 
    # contigency_table$percent_sig_above_center = (contigency_table$above_center_sig/(contigency_table$above_center_sig+contigency_table$above_center_not_sig)) * 100
    # contigency_table$percent_sig_below_center = (contigency_table$below_center_sig/(contigency_table$below_center_sig+contigency_table$below_center_not_sig)) * 100
    # contigency_table$percent_sig_above_center[is.na(contigency_table$percent_sig_above_center)] = 0
    # contigency_table$percent_sig_below_center[is.na(contigency_table$percent_sig_below_center)] = 0
    # 
    # exact_res = apply(contigency_table, 1, function(counts) {
    #   temp_df = data.frame(above_center=c(counts["above_center_sig"],counts["above_center_not_sig"]),below_center=c(counts["below_center_sig"],counts["below_center_not_sig"]))
    #   temp_df$above_center = as.numeric(temp_df$above_center)
    #   temp_df$below_center = as.numeric(temp_df$below_center)
    #   rownames(temp_df) = c("sig","not sig")
    #   fishers_res = fisher.test(temp_df,alternative="two.sided",conf.int=TRUE)
    #   pvalue = fishers_res$p.value
    #   effect_size = fishers_res$estimate
    #   confidence_interval_lower = fishers_res$conf.int[1]
    #   confidence_interval_upper = fishers_res$conf.int[2]
    #   return(c(pvalue,effect_size,confidence_interval_lower,confidence_interval_upper))
    # })
    # exact_res = t(exact_res)
    # colnames(exact_res) = c("pvalue","effect_size","confidence_interval_lower","confidence_interval_upper")
    # exact_res = as.data.frame(exact_res)
    # # do some basic checks because sometimes p-values are incorrect if we lack variability
    # exact_res$pvalue[rowSums(contigency_table[,c("above_center_sig","below_center_sig")])==0] = 1
    # exact_res$pvalue[(contigency_table$percent_sig_above_center == 100 & contigency_table$percent_sig_below_center == 0) |
    #                    (contigency_table$percent_sig_below_center == 100 & contigency_table$percent_sig_above_center == 0)] = 0
    # contigency_table = cbind(contigency_table,exact_res)
  return(summarized_voe_data)
}
  
#' summarize_vibration_data_by_feature
#'
#' Summarize output of vibrations for each dependent feature of interest.
#' @param df A dataframe of expanded VoE output (output of filter_unnest_feature_vib)
#' @param center_for_effect_size determines how to calculate janus effect. 0 if values are linear, 1 if values are hazard ratio or risk ratios or fold differences 
#' @keywords voe analysis
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
summarize_vibration_data_by_feature <- function(df,center_for_effect_size){
  # p <- c(0.01,.5,.99)
  # p_names <- purrr::map_chr(p, ~paste0('estimate_quantile_',.x*100, "%"))
  # p_funs <- purrr::map(p, ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% purrr::set_names(nm = p_names)
  # model_counts = df %>% group_by(dependent_feature,term) %>% dplyr::count(.data$dependent_feature) %>% dplyr::rename(number_of_models=.data$n) %>% ungroup()
  # 
  # janus_effect = df %>% dplyr::group_by(.data$dependent_feature,.data$term) %>% 
  # dplyr::summarise(janus_effect = sum(.data$estimate > center_for_effect_size, na.rm = TRUE)/sum(is.finite(.data$estimate), na.rm = TRUE)) %>% dplyr::ungroup()
  # 
  # 
  # df_estimates = suppressMessages(df %>% 
  #                                 dplyr::group_by(.data$dependent_feature,.data$term) %>%
  #                                 dplyr::summarize_at(dplyr::vars(.data$estimate), tibble::lst(!!!p_funs)) %>%
  #                                 dplyr::mutate(estimate_diff_99_1 = .data$`estimate_quantile_99%`-.data$`estimate_quantile_1%`))%>% ungroup()
  # 
  # df_estimates = merge(df_estimates,janus_effect)
  # 
  # 
  # p_names <- purrr::map_chr(p, ~paste0('pval_quantile_',.x*100, "%"))
  # p_funs <- purrr::map(p, ~purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>% purrr::set_names(nm = p_names)
  # df_pval = df %>% dplyr::group_by(.data$dependent_feature,.data$term) %>% dplyr::summarize_at(dplyr::vars(.data$p.value), tibble::lst(!!!p_funs)) %>% dplyr::mutate(pvalue_diff_99_1 = .data$`pval_quantile_99%`-.data$`pval_quantile_1%`) %>% ungroup()
  # summarized_voe_data=dplyr::bind_cols(model_counts, df_estimates %>% dplyr::select(-c(dependent_feature,term)),df_pval %>% dplyr::select(-c(dependent_feature,term)))
  print(colnames(df))
  summarized_voe_data <- df %>% 
    group_by(dependent_feature
             , term) %>%
    summarise(estimate_quantile_1 = quantile(estimate, probs = 0.01)
              , estimate_quantile_50 = quantile(effect_size, probs = 0.5)
              , estimate_quantile_99 = quantile(effect_size, probs = 0.99)
              , estimate_diff_99_1 = estimate_quantile_99 - estimate_quantile_1
              , num_models = sum(!is.na(effect_size))
              , janus_effect = sum(effect_size>center_for_effect_size,na.rm=TRUE)/num_models
              , pval_quantile_1 = quantile(p.value, probs = 0.01)
              , pval_quantile_50 = quantile(p.value, probs = 0.5)
              , pval_quantile_99 = quantile(p.value, probs = 0.99)
              , pvalue_diff_99_1 = pval_quantile_99 - pval_quantile_1
              , estimate_mean_over_vibration = mean(effect_size)
              , estimate_ci_lower_over_vibration = mean(effect_size_ci_lower)
              , estimate_ci_upper_over_vibration = mean(effect_size_ci_upper)
              , estimate_se_over_vibration = (mean(ci_upper_y) - mean(ci_lower_y))/(2*1.96)
              , t_statistic_estimate_over_vibration = mean(y_centered)/estimate_se_over_vibration
              , crosses_the_center = ifelse(center_for_effect_size > estimate_ci_lower_over_vibration 
                                            & center_for_effect_size < estimate_ci_upper_over_vibration
                                            , "yes"
                                            , "no")
              # , p.value_over_vibration = 2*pt(q = (t_statistic_estimate_over_vibration)
              #                                 , df = num_models - 1)
              , p.value_over_vibration =  exp(-0.717*(abs(t_statistic_estimate_over_vibration)) - 0.416*abs(t_statistic_estimate_over_vibration)^2)
              # , exact_results = fisher.test(x=factor(as.numeric(effect_size>center_for_effect_size),levels=c(0,1))
              #                               , y=factor(as.numeric(p.value < 0.05),levels=c(0,1)))$p.value
    ) %>%
    ungroup(.)
  
  return(summarized_voe_data)
}


#' Analyze VoE data
#'
#' Post-process vibration output. Only works with data immediately out of the compute_vibrations function.
#' @param vibration_output Output list from the compute vibrations function.
#' @param confounder_analysis TRUE/FALSE -- run confounder analysis (default = TRUE).
#' @param constant_adjusters A character vector (or just one string) corresponding to column names in your dataset to include in every vibration. (default = NULL)
#' @param num_knots If num_knots is greater than 0 generate the B-spline basis matrix for a natural cubic spline on primary_variable
#' @param center_for_effect_size determines how to calculate janus effect. 0 if values are linear, 1 if values are hazard ratio or risk ratios or fold differences 
#' @keywords voe analysis
#' @importFrom rlang .data
#' @importFrom dplyr "%>%"
#' @export
analyze_voe_data <- function(vibration_output,confounder_analysis,constant_adjusters,num_knots,center_for_effect_size){
  voe_annotated =get_adjuster_expanded_vibrations(vibration_output[[1]], vibration_output[[2]],constant_adjusters)
  # View(voe_annotated)
  voe_unnested_annotated = filter_unnest_feature_vib(voe_annotated) %>% dplyr::select(-.data$vars)
  # View(voe_unnested_annotated)
  # View(voe_annotated$termplot_fit)
  voe_termplot_unnested = voe_annotated %>% 
    dplyr::slice(which(purrr::map_lgl(voe_annotated$termplot_fit, ~class(.)[[1]] == "data.frame"))) %>% 
    tidyr::unnest(.data$termplot_fit) %>% 
    dplyr::select(-c(vars, full_fits, feature_fit))
  # View(voe_termplot_unnested)
  # print("analyze voe data")
  
  if(length(num_knots) > 1 | num_knots[1] > 0) {
    summarized = summarize_vibration_data_by_feature_termplot(df=voe_termplot_unnested,center_for_effect_size=center_for_effect_size)
  } else {
    # print("Here!")
    #summarized = summarize_vibration_data_by_feature(voe_unnested_annotated,center_for_effect_size)
    summarized = summarize_vibration_data_by_feature_termplot(df=voe_termplot_unnested,center_for_effect_size=center_for_effect_size)
  }
  c_analysis='No confounder analysis completed.'
  # if(confounder_analysis==TRUE & num_knots[1] == 0 & length(num_knots) == 1){ # NOTE! WE NEED TO DO THIS FOR NON LINEAR BUT NOT NOW
  #   if(nrow(voe_unnested_annotated)>=10){
  #     c_analysis = find_confounders_linear(voe_unnested_annotated)
  #   }
  #   else{
  #     print('Skipping confounder analysis, as not enough vibrations (under 10) completed to make it worthwhile.') 
  #   }
  # } 
  return(list('summarized_vibration_output'= summarized,'confounder_analysis'=c_analysis,'data'=voe_unnested_annotated,'termplot_data'=voe_termplot_unnested))
}

