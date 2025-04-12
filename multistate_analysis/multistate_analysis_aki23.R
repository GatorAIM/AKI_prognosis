################################################################################
#####       Multi-State Analysis of AKI Recovery (Stage 2&3 at Onset)     ######
################################################################################

####  Step 1: Load and Preprocess Data 
# Load libraries
library('arrow')
library('dplyr')
library('tidyr')
library('survival')
library('stringr')
library('ggplot2')

# Load data
base_path = './'
aki_subgrp = 'aki23'

df <- read_parquet(paste0(base_path,"data_msm_", aki_subgrp,".parquet"))
df <- subset(df, select = -`__index_level_0__`)
colnames(df) <- c("ENCOUNTERID", 
                  "FROM_STATUS",  
                  "TO_STATUS", 
                  'TRANS',
                  'tstart',
                  'tstop',
                  'EVENT_STATUS',
                  "SCR_BASELINE", "SCR_FD", "SCR_MEAN",
                  "SYSTOLIC_FD", "SYSTOLIC_MEAN", "DIASTOLIC_FD",  "DIASTOLIC_MEAN", 
                  "AGE",  "ONSET_SINCE_ADMIT",
                  'ALBUMIN', "ALP", 'AST',  "BILIRUBIN", "BNP", "CALCIUM", 
                  "CO2", "HEMATOCRIT", "POTASSIUM", "RBC", 
                  "site"
                  )

# Mask site names
site_mapping <- c(
  "... MASKED_FOR_ANONYMITY": "... MASKED_FOR_ANONYMITY"
)
site_names <- c("Site1", "Site2", "Site3", "Site4")  
df$site <- recode(df$site, !!!site_mapping)
df$site <- factor(df$site)
df$ENCOUNTERID <- factor(df$ENCOUNTERID)
site_data <- split(df, df$site)


# Trim outliers and scale features
trim_outliers <- function(x) {
  upper_bound <- quantile(x, 0.99, na.rm = TRUE)
  lower_bound <- quantile(x, 0.01, na.rm = TRUE)
  x[x > upper_bound & !is.na(x)] <- upper_bound
  x[x < lower_bound & !is.na(x)] <- lower_bound
  return(x)
}


model_datasets <- list()

for (site_name in names(site_data)){
  data_transformed <- site_data[[site_name]]
  nonscaled_vars <- c('ENCOUNTERID',
                      'EVENT_STATUS',
                      'FROM_STATUS',
                      'TO_STATUS',
                      'TRANS',
                      'tstart',
                      'tstop',
                      'site'
                      ) 
  continuous_vars <- setdiff(colnames(data_transformed), nonscaled_vars)  
  data_transformed[, continuous_vars] <- lapply(data_transformed[, continuous_vars], trim_outliers)
   
  data_scaled = data_transformed 
  continuous_vars <- setdiff(colnames(data_transformed), nonscaled_vars) 
  data_scaled[, continuous_vars] <- scale(data_transformed[, continuous_vars])

  data_scaled$TRANS = as.factor(data_scaled$TRANS)
  model_datasets[[site_name]] <- data_scaled
}


####  Step 2: Fit  Multi-State Cox PH Models 

covs = c("SCR_BASELINE", "SCR_FD", "SCR_MEAN", 
         "SYSTOLIC_FD", "SYSTOLIC_MEAN",  "DIASTOLIC_FD", "DIASTOLIC_MEAN", 
         "AGE", "ONSET_SINCE_ADMIT",
         'ALBUMIN', "ALP", 'AST',  "BILIRUBIN", "BNP", "CALCIUM", 
         "CO2", "HEMATOCRIT", "POTASSIUM", "RBC")


mscox_fit <- list()
surv_fit <- list()
for (site_name in  site_names){
  if (site_name =='Site1'){
    cox_formula = as.formula(paste('Surv(tstart, tstop, EVENT_STATUS) ~ ', 
                                   paste(paste(covs[!covs %in% c("CO2")], 
                                               collapse = ":TRANS + "), 
                                         ':TRANS'), '+strata(TRANS)'))
  } else if (site_name =='Site2'){
    cox_formula = as.formula(paste('Surv(tstart, tstop, EVENT_STATUS) ~ ', 
                                   paste(paste(covs[!covs %in% c("BNP")], 
                                               collapse = ":TRANS + "), 
                                         ':TRANS'), '+strata(TRANS)'))
  } else if (site_name =='Site3'){
    cox_formula = as.formula(paste('Surv(tstart, tstop, EVENT_STATUS) ~ ', 
                                   paste(paste(covs, 
                                               collapse = ":TRANS + "), 
                                         ':TRANS'), '+strata(TRANS)'))
  } else{
    cox_formula = as.formula(paste('Surv(tstart, tstop, EVENT_STATUS) ~ ', 
                                   paste(paste(covs, 
                                               collapse = ":TRANS + "), 
                                         ':TRANS'), '+strata(TRANS)'))
  }
    
  mscox_site <- coxph(formula = cox_formula, 
                      data = model_datasets[[site_name]], 
                      id=ENCOUNTERID,
                      control = coxph.control(
                        iter.max = 100,        
                        outer.max = 100 
                      )
  )
  
  mscox_fit[[site_name]] <- mscox_site
}



#### Step 3:  Plot Log Hazard Ratios   
get_legend <- function(a_gplot, return_all = FALSE) {
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  
  if (length(leg) > 1 && return_all) {
    return(tmp$grobs[leg])
  } else {
    return(tmp$grobs[[leg[1]]])  
  }
}

coef_data <- data.frame()
for (site in site_names) {
  cox_model <- mscox_fit[[site]]
  coefs <- summary(cox_model)$coefficients
  transition_regex_1 <- "^(.*?):TRANS([0-9]+)$"
  transition_regex_2 <- "^TRANS([0-9]+):(.*)$"
  parsed_1 <- str_match(rownames(coefs)[1:12], transition_regex_1)
  parsed_2 <- str_match(rownames(coefs)[13:length(rownames(coefs))], transition_regex_2)
  parsed <- rbind(parsed_1[,c(2,3)], parsed_2[,c(3,2)])
  coefs_df <- data.frame(
    Variable = parsed[, 1],  
    Transition = parsed[, 2], 
    Coefficient = coefs[, "coef"],
    LowerCI = coefs[, "coef"] - 1.96 * coefs[, "se(coef)"],
    UpperCI = coefs[, "coef"] + 1.96 * coefs[, "se(coef)"],
    Site = site,
    stringsAsFactors = FALSE
  )
  coef_data <- rbind(coef_data, coefs_df)
}

coef_data <- na.omit(coef_data)
var_common <- c('SCR_BASELINE', 'SCR_FD', 
                'SCR_MEAN', 'SYSTOLIC_FD',
                'SYSTOLIC_MEAN', 'ALBUMIN')
var_name_list <- c('SCR_BASELINE' = 'SCr (baseline)', 
                   'SCR_FD' = 'SCr (slope)',
                   'SCR_MEAN' = 'SCr (mrv)',
                   'SYSTOLIC_FD'= 'SysBP (slope)',
                   'SYSTOLIC_MEAN'= 'SysBP (mrv)',
                   'ALBUMIN' = 'Albumin',
                   'DIASTOLIC_FD'= 'DiaBP (slope)',  
                   'DIASTOLIC_MEAN'= 'DiaBP (mrv)',
                   'AGE'= 'Age',
                   'ONSET_SINCE_ADMIT'= 'Days Since Adm.',
                   'ALP' = 'ALP',
                   'AST' = 'AST',
                   'BILIRUBIN'= 'Bilirubin',
                   'BNP' = 'BNP',
                   'CALCIUM' = 'Calcium',
                   'CO2' =  'CO2',
                   "HEMATOCRIT" = 'Hematocrit',
                   "POTASSIUM" = 'Potassium',
                   'RBC' = 'RBC'
)
trans_name_list <- c(
  '1' = 'AKI0 to Discharge',
  '5' = 'AKI1 to Discharge',
  '9' = 'AKI2&3 to Discharge',
  '6' = 'AKI1 to AKI0',
  '10' = 'AKI2&3 to AKI0', 
  '2' = 'AKI0 to AKI1',
  '11' = 'AKI2&3 to AKI1',
  '3' = 'AKI0 to AKI2&3',
  '7' = 'AKI1 to AKI2&3',
  '4' = 'AKI0 to Death',
  '8' = 'AKI1 to Death',
  '12' = 'AKI2&3 to Death'
)

all_vars <- unique(coef_data$Variable)
all_transitions <- unique(coef_data$Transition)
expanded_data <- expand.grid(Variable = all_vars, 
                             Transition = all_transitions, 
                             Site = unique(coef_data$Site))
full_coef_data <- merge(expanded_data, coef_data, 
                        by = c("Variable", "Transition", "Site"), 
                        all.x = TRUE)
full_coef_data$Variable <- factor(full_coef_data$Variable, levels = all_vars)
full_coef_data$Transition <- factor(full_coef_data$Transition, levels = names(trans_name_list))

custom_labeller <- function(var_name_list, var_common) {
  function(variable) {
    sapply(variable, function(x) {
      if (x %in% var_common) {
        bquote(bold(.(var_name_list[[x]])))  
      } else {
        bquote(plain(.(var_name_list[[x]])))  
      }
    })
  }
}

plot_forest <- function(data, variables, var_name_list, trans_name_list, var_common) {
  data_subset <- data[(data$Variable %in% variables) & (data$Transition %in% names(trans_name_list)), ]
  
  ggplot(data_subset, aes(x = Coefficient, y = Transition, color = Site)) +
    geom_point() +
    geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.6) +
    geom_hline(yintercept = 3.5, linetype = "dashed", color = "grey", size = 0.5) +
    geom_hline(yintercept = 5.5, linetype = "dashed", color = "grey", size = 0.5) +
    geom_hline(yintercept = 7.5, linetype = "dashed", color = "grey", size = 0.5) +
    geom_hline(yintercept = 9.5, linetype = "dashed", color = "grey", size = 0.5) +
    facet_wrap(
      ~ Variable,
      scales = "free_x",
      nrow = 1,
      labeller = as_labeller(custom_labeller(var_name_list, var_common), label_parsed)
    ) +
    scale_y_discrete(labels = trans_name_list) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10),  # Apply consistent size here
      legend.position = "bottom",
      panel.grid.major = element_blank()
    ) +
    labs(x = "Coefficient (log hazard ratio)", y = NULL, color = "Site")
}

set_variable_order <- function(data, variable_order) {
  data$Variable <- factor(data$Variable, levels = variable_order)
  return(data)
}

variables_group_1 <- names(var_name_list)[1:5]
variables_group_2 <- names(var_name_list)[6:10]
variables_group_3 <- names(var_name_list)[11:15]
variables_group_4 <- names(var_name_list)[16:19]

full_coef_data_group_1 <- set_variable_order(full_coef_data, variables_group_1)
full_coef_data_group_2 <- set_variable_order(full_coef_data, variables_group_2)
full_coef_data_group_3 <- set_variable_order(full_coef_data, variables_group_3)
full_coef_data_group_4 <- set_variable_order(full_coef_data, variables_group_4)

plot1 <- plot_forest(full_coef_data_group_1, variables_group_1, var_name_list, trans_name_list, var_common)
plot2 <- plot_forest(full_coef_data_group_2, variables_group_2, var_name_list, trans_name_list, var_common)
plot3 <- plot_forest(full_coef_data_group_3, variables_group_3, var_name_list, trans_name_list, var_common)
plot4 <- plot_forest(full_coef_data_group_4, variables_group_4, var_name_list, trans_name_list, var_common)

mscox_plot_final <-grid.arrange(
  arrangeGrob(
    plot1 + theme(legend.position = "none"),
    plot2 + theme(legend.position = "none"),
    plot3 + theme(legend.position = "none"),
    plot4 + theme(legend.position = "none"),
    ncol = 1
  ),
  get_legend(plot1),
  ncol = 1, heights = c(10, 1)
)

# Save the log hazard ratio plot
png(paste0(base_path,"/result/plot_mscox_aki23.png"), 
    width = 12, 
    height = 10, 
    units = "in", 
    res = 300)
grid.draw(mscox_plot_final)  
dev.off()




#### Step 4: Estimate the Instantaneous Hazard Functions

get_boot_insthaz <- function(surv_fit, boot_iter, site){
  time_range <- 0:7
  cumhaz <- surv_fit$cumhaz
  time <- surv_fit$time
  strata_counts <- surv_fit$strata  
  strata_list <- split(cumhaz, rep(names(strata_counts), strata_counts))
  time_list <- split(time, rep(names(strata_counts), strata_counts))
  padded_cumhaz_list <- list()
  padded_time_list <- list()
  instantaneous_hazards <- list()
  
  for (strata_name in names(strata_list)) {
    original_time <- time_list[[strata_name]]
    original_cumhaz <- strata_list[[strata_name]]
    padded_cumhaz <- numeric(length(time_range))
    padded_time <- time_range
    padded_cumhaz[1] <- 0
    padded_cumhaz[match(original_time, time_range)] <- original_cumhaz
    missing_times <- setdiff(time_range, original_time)
    padded_cumhaz[match(missing_times, time_range)] <- 0
    padded_cumhaz_list[[strata_name]] <- padded_cumhaz
    padded_time_list[[strata_name]] <- padded_time
    instantaneous_hazards[[strata_name]] <- diff(padded_cumhaz) / diff(padded_time)
  }
  
  results <- data.frame(
    Strata = rep(names(instantaneous_hazards), each = length(time_range) - 1),
    Time = rep(time_range[-1], times = length(instantaneous_hazards)),
    Instantaneous_Hazard = unlist(instantaneous_hazards)
  )
  
  results$Strata <- as.integer(gsub("strata\\(TRANS\\)=([0-9]+)", "\\1", results$Strata))
  results$boot_iter <- boot_iter
  results$site <- site
  return(results)
}

# Bootstrapping the instantaneous hazard functions
n_boot <- 99  
inst_hazard <- data.frame()
for (site_name in site_names) {
  data_orig <- model_datasets[[site_name]]
  enc_list <- unique(data_orig$ENCOUNTERID)
  n_enc <- length(enc_list)
  for (boot_iter in 0:n_boot) {
    if (boot_iter == 0) {
      surv_fit_site <- survfit(Surv(tstart, tstop, EVENT_STATUS) ~ strata(TRANS), 
                               data = data_orig)
    } else {
      set.seed(123 + boot_iter)
      enc_boot_list <- sample(enc_list, n_enc, replace = TRUE)
      data_boot <- merge(data.frame(ENCOUNTERID = enc_boot_list), 
                         data_orig, 
                         by = "ENCOUNTERID", 
                         all.x = TRUE)
      surv_fit_site <- survfit(Surv(tstart, tstop, EVENT_STATUS) ~ strata(TRANS), 
                               data = data_boot)
    }
    
    inst_hazard_boot <- get_boot_insthaz(surv_fit_site, boot_iter, site_name)
    inst_hazard <- rbind(inst_hazard, inst_hazard_boot)
  }
}



####  Step 5:  Plot the baseline hazard functions  
inst_hazard_summary <- inst_hazard %>%
  mutate(Strata = factor(Strata, levels = names(trans_name_list), labels = trans_name_list)) %>%
  group_by(Strata, Time, site) %>%
  summarize(
    mean_hazard = mean(Instantaneous_Hazard, na.rm = TRUE),
    lower_ci = quantile(Instantaneous_Hazard, 0.025, na.rm = TRUE),
    upper_ci = quantile(Instantaneous_Hazard, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

strata_site_combinations <- inst_hazard_summary %>%
  distinct(Strata, site)

padding_rows <- strata_site_combinations %>%
  mutate(
    Time = 0,
    mean_hazard = 0,
    lower_ci = 0,
    upper_ci = 0
  )

inst_hazard_summary_padded <- inst_hazard_summary %>%
  bind_rows(padding_rows) %>%
  arrange(Strata, site, Time)

create_subplot <- function(transitions, inst_hazard_data, subtitle, nrow) {
  filtered_data <- inst_hazard_data %>% filter(Strata %in% transitions)
  ggplot(filtered_data, aes(x = Time, y = mean_hazard, color = site, linetype = Strata)) +
    geom_line(size = 0.8, alpha = 0.7) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = site), 
                color = NA, alpha = 0.1, show.legend = FALSE) +
    scale_linetype_discrete() +
    scale_x_continuous(breaks = 0:7, limits = c(0, 7)) + 
    labs(title = paste("Transitions", subtitle),
         x = "Days", y = "Instantaneous Hazard Rate") +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      legend.position = c(0.01, 1.03), 
      legend.justification = c(0, 1)
    ) +
    guides(
      linetype = guide_legend(nrow = nrow),
      color = "none"  
    )
}

transitions_list <- list(
  subplot1 = c("AKI0 to Discharge", "AKI1 to Discharge", "AKI2&3 to Discharge"),
  subplot2 = c("AKI1 to AKI0", "AKI2&3 to AKI0"),
  subplot3 = c("AKI0 to AKI1", "AKI2&3 to AKI1"),
  subplot4 = c("AKI0 to AKI2&3", "AKI1 to AKI2&3"),
  subplot5 = c("AKI0 to Death", "AKI1 to Death", "AKI2&3 to Death")
)

plot1 <- create_subplot(transitions_list$subplot1, inst_hazard_summary_padded, 'to Discharge', 2)
plot2 <- create_subplot(transitions_list$subplot2, inst_hazard_summary_padded, 'to No AKI', 1)
plot3 <- create_subplot(transitions_list$subplot3, inst_hazard_summary_padded, 'to AKI1', 1)
plot4 <- create_subplot(transitions_list$subplot4, inst_hazard_summary_padded, 'to AKI2&3', 1)
plot5 <- create_subplot(transitions_list$subplot5, inst_hazard_summary_padded, 'to Death', 2)

legend_site <- get_legend(
  ggplot(inst_hazard_summary, aes(x = Time, y = mean_hazard, color = site)) +
    geom_line() +
    labs(color = "Site") +
    theme(legend.position = "bottom")
)

basehaz_plot_final <- grid.arrange(
  plot1, plot2, plot3, plot4, plot5, legend_site,
  layout_matrix = rbind(
    c(1, 2, 3),
    c(4, 5, 6)
  ),
  heights = c(8, 8),
  nrow = 2,
  ncol = 3
)

# Save the baseline hazard plot
png(paste0(base_path,"/result/plot_basehaz_aki23.png"), 
    width = 11, 
    height = 8, 
    units = "in", 
    res = 300)
grid.draw(basehaz_plot_final)  
dev.off()

