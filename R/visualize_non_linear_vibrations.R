#' Visualization of Vibration of Effects
#'
#' Angel Hair Plot
#' Visualize the associations between outcome and main predictor across all model specifications. 
#' Each model specification is a strand of angel hair
#' @param list_voe_results list of the association statistics for all the specifications
#' 
angel_hair_plot <- function(list_voe_results)
{
  df_summarized_voe <- list_voe_results$vibration_output$summarized_vibration_output
  
  df_all_voe_models <- list_voe_results$vibration_output$termplot_data
  
  df_combination <- df_summarized_voe %>%
    select(dependent_feature
           , independent_feature) %>%
    unique(.) %>%
    mutate(combination = paste(dependent_feature 
                               , independent_feature
                               , sep = " - "))
  
  all_combinations <- df_combination %>%
    pull(combination)
  
  num_combinations <- length(all_combinations)
  
  list_combination <- list()
  
  for(i in seq(num_combinations))
  {
    combination_i <- all_combinations[i]
    print(combination_i)
    subset_combination_i <- df_combination %>%
      filter(combination == combination_i)
    
    dependent_feature_i <- subset_combination_i %>%
      pull(dependent_feature)
    
    independent_feature_i <- subset_combination_i %>%
      pull(independent_feature)
    
    subset_voe_models_i <- df_all_voe_models %>%
      filter(dependent_feature == dependent_feature_i 
             & independent_feature == independent_feature_i)
    # View(subset_voe_models_i)
    
    binary_var_pres_absent = subset_voe_models_i[,list_voe_results$vibration_variables[,1]]
    binary_var_characters = lapply(colnames(binary_var_pres_absent), function(mycolname) {
      vector_vibrate_variable <- binary_var_pres_absent %>% 
        pull(all_of(mycolname))
     
     string_vibration_variable <- ifelse(vector_vibrate_variable  == 0
                                         , ""
                                         , mycolname)
     })
    binary_var_characters = bind_cols(binary_var_characters)
    new_column_names = paste("character",colnames(binary_var_pres_absent),sep="_")
    colnames(binary_var_characters) = new_column_names
    subset_voe_models_i = cbind(subset_voe_models_i,binary_var_characters)
    
    formula_covariates <- apply(subset_voe_models_i[,new_column_names],1,function(x) paste(x[x!=""],collapse=" + "))
    
    subset_voe_models_i <- subset_voe_models_i %>%
      mutate(formula_covariates = formula_covariates) %>%
      mutate(regression_formula = paste(dependent_feature
                                        , " ~ "
                                        , independent_feature
                                        , ifelse(formula_covariates == ""
                                                 , ""
                                                 , " + ")
                                        , formula_covariates
                                        , sep = "")) 
    percent_confounders = rowSums(subset_voe_models_i[,list_voe_results$vibration_variables[,1]])/length(list_voe_results$vibration_variables[,1])
    subset_voe_models_i$confounder_category = "vibration model"
    subset_voe_models_i$confounder_category[percent_confounders==0] = "univariate model"
    subset_voe_models_i$confounder_category[percent_confounders==1] = "full model"
    # View(subset_voe_models_i)
    
   
    
    if("lower_bound_quantiles" %in% colnames(df_all_voe_models))
    {
      
      subset_voe_models_i <- subset_voe_models_i %>%
        mutate(label = ifelse(confounder_category == "full model"
                              , gsub("_", "", x)
                              , "")) %>%
        mutate(mean_bound = sqrt(lower_bound_quantiles*upper_bound_quantiles))
      
      angel_hair_plot <- ggplot(data = subset_voe_models_i) +
        geom_rect(mapping = aes(xmin = lower_bound_quantiles
                                , xmax = upper_bound_quantiles
                                , ymin = effect_size_ci_lower
                                , ymax = effect_size_ci_upper
                                , color = confounder_category
                                , group = regression_formula
                                , fill = p.value)) +
        geom_text(mapping = aes(x = mean_bound
                                , y = effect_size
                                , label = label)
                  , color = "black"
                  , size = 5) +
        scale_fill_gradientn(colors = c("gray10", "gray50", "gray85", "gray95")
                             , limits = c(0, 1)
                             , breaks = c(0.01, 0.05, 0.1, 1)
                             , values = c(0, 0.01
                                          , 0.011, 0.05
                                          , 0.051, 0.1
                                          , 0.101, 1)) +
        guides( fill = guide_colorbar(title = "p-values"
                                       , barwidth = 30
                                       , label.theme = element_text(size = 7)
                                       , nbin = 200
                                       , order = 2)) +
        scale_x_log10() +
        scale_y_log10() +
        geom_hline(yintercept = 1) +
        theme(legend.position = "top"
              , legend.direction = "horizontal"
              , legend.box = "vertical")
      
    } else {
      angel_hair_plot <- ggplot(data = subset_voe_models_i) +
        geom_line(mapping = aes(x = x
                                , y = effect_size
                                , color = confounder_category
                                , group = regression_formula)) +
        scale_x_log10() +
        geom_hline(yintercept = 1)
    }
    
    list_combination[[combination_i]] <- angel_hair_plot
    
  }
  return(list_combination)
}

#' Fettuccine, linguine, or bucatini plot
#' Visualize the summary of the associations across all model specifications
#' One curve to represent the mean association across all model specifications
#' Confidence interval represent the 95% CI across all model specifications
fettuccine_linguine_bucatini_plot <- function(list_voe_results)
{
  df_summarized_voe <- list_voe_results$vibration_output$summarized_vibration_output
  
  df_combination <- df_summarized_voe %>%
    select(dependent_feature
           , independent_feature) %>%
    unique(.) %>%
    mutate(combination = paste(dependent_feature 
                               , independent_feature
                               , sep = " - "))
  
  all_combinations <- df_combination %>%
    pull(combination)
  
  num_combinations <- length(all_combinations)
  
  list_combination <- list()
  
  for(i in seq(num_combinations))
  {
    combination_i <- all_combinations[i]
    print(combination_i)
    subset_combination_i <- df_combination %>%
      filter(combination == combination_i)
    
    dependent_feature_i <- subset_combination_i %>%
      pull(dependent_feature)
    
    independent_feature_i <- subset_combination_i %>%
      pull(independent_feature)
    
    subset_summarized_voe_i <- df_summarized_voe %>%
      filter(dependent_feature == dependent_feature_i 
             & independent_feature == independent_feature_i)
    
    subset_summarized_voe_i$sig_stars = ""
    subset_summarized_voe_i$sig_stars[subset_summarized_voe_i$p.value_over_vibration < 0.05] = "*"
    # View(subset_summarized_voe_i)
    
    subset_summarized_voe_i <- subset_summarized_voe_i %>%
      mutate(mean_bound = sqrt(lower_bound_quantiles*upper_bound_quantiles)) %>%
      mutate(label = gsub("_", "", x))
      
    if("lower_bound_quantiles" %in% colnames(subset_summarized_voe_i))
    {
      fettuccine_linguine_bucatini_plot <- ggplot(data = subset_summarized_voe_i) +
        geom_rect(mapping = aes(xmin = lower_bound_quantiles
                                , xmax = upper_bound_quantiles
                                , ymin = estimate_ci_lower_over_vibration
                                , ymax = estimate_ci_upper_over_vibration
                                , fill = p.value_over_vibration)
                  , color = "black") +
        geom_text(mapping = aes(x = mean_bound
                                , y = estimate_mean_over_vibration
                                , label = label)
                  , color = "black"
                  , size = 5) +
        scale_fill_gradientn(colors = c("gray10", "gray50", "gray85", "gray95")
                             , limits = c(0, 1)
                             , breaks = c(0.01, 0.05, 0.1, 1)
                             , values = c(0, 0.01
                                          , 0.011, 0.05
                                          , 0.051, 0.1
                                          , 0.101, 1)) +
        guides( fill = guide_colorbar(title = "p-values"
                                      , barwidth = 30
                                      , label.theme = element_text(size = 7)
                                      , nbin = 200
                                      , order = 2)) +
        scale_x_log10() +
        scale_y_log10() +
        geom_hline(yintercept = 1) +
        theme(legend.position = "top"
              , legend.direction = "horizontal"
              , legend.box = "vertical") 
      
    } else {
      fettuccine_linguine_bucatini_plot <- ggplot(data = subset_summarized_voe_i ) +
        geom_hline(yintercept = 1) +
        geom_ribbon(mapping = aes(x = x
                                  , ymin = estimate_ci_lower_over_vibration
                                  , ymax = estimate_ci_upper_over_vibration
        )
        , fill = "gray70"
        , alpha = 0.5) +
        geom_line(mapping = aes(x = x
                                , y = estimate_mean_over_vibration
                                , group = "mean")) +
        scale_x_log10() +
        scale_y_log10() +
        geom_text(aes(x=x,y=estimate_ci_upper_over_vibration,label=sig_stars)) 
    }
    

    
    list_combination[[combination_i]] <- fettuccine_linguine_bucatini_plot
    
  }
  return(list_combination)
}


#' Dot plot
#' Visualize summary of the association over deciles of concentrations for all main predictors
#' A dot represnt the mean associatio for values within a decile and a main predictor

make_dot_plot <- function(list_voe_results,method) {
  df_summarized_voe <- list_voe_results$vibration_output$summarized_vibration_output
  
  df_combination <- df_summarized_voe %>%
    select(dependent_feature
           , independent_feature) %>%
    unique(.) %>%
    mutate(combination = paste(dependent_feature 
                               , independent_feature
                               , sep = " - "))
  
  all_combinations <- df_combination %>%
    pull(combination)
  
  num_combinations <- length(all_combinations)
  
  for(i in seq(num_combinations))
  {
    combination_i <- all_combinations[i]
    print(combination_i)
    subset_combination_i <- df_combination %>%
      filter(combination == combination_i)
    
    dependent_feature_i <- subset_combination_i %>%
      pull(dependent_feature)
    
    independent_feature_i <- subset_combination_i %>%
      pull(independent_feature)
    
    subset_summarized_voe_i <- df_summarized_voe %>%
      filter(dependent_feature == dependent_feature_i 
             & independent_feature == independent_feature_i)
    
    subset_summarized_voe_i = subset_summarized_voe_i %>% 
      group_by(independent_feature,dependent_feature) %>%
      mutate(quantile_bins = cut(x,
                                 breaks=quantile(x,seq(0,1,by=0.1)))) %>% 
      mutate(quantile_bins = as.numeric(quantile_bins))

    # first bin of cut is exclusive so lets put values that are equal to minimum number into first bin
    subset_summarized_voe_i$quantile_bins[is.na(subset_summarized_voe_i$quantile_bins)] = 1

    subset_summarized_voe_i$quantile_groups = paste(subset_summarized_voe_i$independent_feature,subset_summarized_voe_i$dependent_feature,subset_summarized_voe_i$quantile_bins,sep="split_here")
    subset_summarized_voe_i_list = split(subset_summarized_voe_i,subset_summarized_voe_i$quantile_groups)
    if(method == "meta") {
      
      meta_analysis_out = lapply(subset_summarized_voe_i_list, function(tempdf) {
        quantile_temp = unique(tempdf$quantile_bins)
        independent_feature_temp = unique(tempdf$independent_feature)
        dependent_feature_temp = unique(tempdf$dependent_feature)
        meta_out = metagen(TE=estimate_mean_over_vibration,
                         seTE=estimate_se_over_vibration,
                         studlab=quantile_bins,
                         data=tempdf,
                         common = FALSE, random = TRUE,
                         method.tau = 'REML',
                         method.random.ci = FALSE,
                         prediction = TRUE,
                         sm = "HR")
        meta_pval = meta_out$pval.random
        meta_estimate = meta_out$TE.random
        return(c(meta_pval,meta_estimate,quantile_temp,independent_feature_temp,dependent_feature_temp))
      })
      meta_analysis_out = do.call("rbind",meta_analysis_out)
      colnames(meta_analysis_out) = c("pval","estimate","quantile","independent_feature","dependent_feature")
      meta_analysis_out = as.data.frame(meta_analysis_out)
      meta_analysis_out$pval = as.numeric(meta_analysis_out$pval)
      meta_analysis_out$neg_log10pval = -log10(meta_analysis_out$pval)
      meta_analysis_out$estimate = as.numeric(meta_analysis_out$estimate)
      meta_analysis_out$quantile = factor(meta_analysis_out$quantile,levels=1:10)
    # make plot
      ggplot(meta_analysis_out,aes(x=independent_feature,y=quantile)) + geom_point(aes(size=neg_log10pval,color=estimate))
    } else if(method="mean") {
      quantile_summarized_data = lapply(subset_summarized_voe_i_list, function(tempdf) {
        quantile_temp = unique(tempdf$quantile_bins)
        independent_feature_temp = unique(tempdf$independent_feature)
        dependent_feature_temp = unique(tempdf$dependent_feature)
        estimate_temp = mean(tempdf$estimate_mean_over_vibration)
        pvalue_temp = mean(tempdf$p.value_over_vibration)
        return(c(pvalue_temp,estimate_temp,quantile_temp,independent_feature_temp,dependent_feature_temp))
      })
      quantile_summarized_data = do.call("rbind",quantile_summarized_data)
      colnames(quantile_summarized_data) = c("pval","estimate","quantile","independent_feature","dependent_feature")
      quantile_summarized_data = as.data.frame(quantile_summarized_data)
      quantile_summarized_data$pval = as.numeric(quantile_summarized_data$pval)
      quantile_summarized_data$neg_log10pval = -log10(quantile_summarized_data$pval)
      quantile_summarized_data$estimate = as.numeric(quantile_summarized_data$estimate)
      quantile_summarized_data$quantile = factor(quantile_summarized_data$quantile,levels=1:10)
      ggplot(quantile_summarized_data,aes(x=independent_feature,y=quantile)) + geom_point(aes(size=neg_log10pval,color=estimate))
      else if(method="i dont know yet") {
      quantile_summarized_data = lapply(subset_summarized_voe_i_list, function(tempdf) {
        mean_of_vibration_mean = mean(tempdf$estimate_mean_over_vibration)
        ci_lower_of_ci_lower = mean(tempdf$estimate_ci_lower_over_vibration)
        ci_upper_of_ci_upper = mean(tempdf$estimate_ci_upper_over_vibration)
        estimate_se_of_vibration_se = (mean(ci_upper_of_ci_upper) - mean(ci_lower_of_ci_lower))/(2*1.96)
        t_statistic_estimate_of_vibrations = mean(tempdf$y_centered)/estimate_se_of_vibration_se
        p.value_over_vibration2 =  exp(-0.717*(abs(t_statistic_estimate_of_vibrations)) - 0.416*abs(t_statistic_estimate_of_vibrations)^2)
        return(c(p.value_over_vibration2,mean_of_vibration_mean,quantile_temp,independent_feature_temp,dependent_feature_temp))
      }
  }
}
