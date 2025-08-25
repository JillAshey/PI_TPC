## Run as background job in R

library(dplyr)
library(purrr)
library(minpack.lm)  # for nlsLM
library(tidyverse)
library(nls.multstart)
library(broom)
library(LoLinR)
library(readr)
library(lubridate)
library(fuzzyjoin)
library(stringr)
library(purrr)
library(future)
library(furrr)
library(minpack.lm)

# 1️⃣ Load the metadata ----

run_meta <- read_csv("data/pi_curves/pi_curve_metadata_run.csv") %>%
  mutate(
    Start = as_datetime(Date) + hours(hour(Start.time)) + minutes(minute(Start.time)),
    Stop  = as_datetime(Date) + hours(hour(Stop.time)) + minutes(minute(Stop.time))
  )

sample_meta <- read_csv("data/pi_curves/pi_curve_metadata_sample.csv")

# Make sure Chamber well # is numeric
sample_meta <- sample_meta %>%
  mutate(`Chamber` = as.numeric(`Chamber`))


# 2️⃣ List all raw data files ----

files <- list.files(path = "data/pi_curves/raw_data/Oct2024", pattern = "*.csv", full.names = TRUE)


# 3️⃣ Read, process & join each file ----

all_data <- map_dfr(files, function(f) {
  # Extract channel number from filename, e.g. Ch1 → 1
  channel <- str_extract(f, "Ch\\d+") %>% str_remove("Ch") %>% as.numeric()
  
  # Read the file
  raw <- read_csv(f) %>%
    mutate(
      # Parse datetime (adjust format if needed)
      DateTime = dmy_hms(Date),
      Channel = channel
    )
  
  # Join to run metadata by fuzzy time match
  raw_joined <- fuzzy_left_join(
    raw,
    run_meta,
    by = c("DateTime" = "Start", "DateTime" = "Stop"),
    match_fun = list(`>=`, `<=`)
  )
  
  # Join to sample metadata by Date, Run #, and Channel
  raw_joined <- raw_joined %>%
    left_join(
      sample_meta,
      by = c(
        "Date.y" = "Date",
        "Run" = "Run",
        "Channel" = "Chamber"
      )
    )
  
  return(raw_joined)
})


# 4️⃣ Inspect and save ----

all_data %>% count(is.na(`Colony_ID`))
all_data %>% count(is.na(Run))

write_csv(all_data, "output/all_joined_data_Oct2024.csv")

# Make sure you have all_data loaded already

# 1️⃣ Make list of all Colony_ID / Run combinations (including blanks)
colony_runs <- all_data %>%
  filter(!is.na(Colony_ID)) %>%
  distinct(Colony_ID, Run)

# 2️⃣ Function: blank correction + slope extraction + save/return interval summary
fit_colony_intervals <- function(colony_id, run_num) {
  
  # Sample data for this colony/run
  sample_data <- all_data %>%
    filter(Colony_ID == colony_id, Run == run_num)
  
  # Skip if no data
  if (nrow(sample_data) == 0) return(NULL)
  
  # Blank data for same run
  blank_data <- all_data %>%
    filter(Colony_ID == "blank", Run == run_num)
  
  # Blank correction if possible
  if (nrow(blank_data) > 0) {
    blank_means <- blank_data %>%
      group_by(Light_Value) %>%
      summarise(blank_O2 = mean(Oxygen, na.rm = TRUE), .groups = "drop")
    
    sample_data <- sample_data %>%
      left_join(blank_means, by = "Light_Value") %>%
      mutate(Oxygen = Oxygen - blank_O2)
  }
  
  # Split into intervals by light step
  intervals <- sample_data %>%
    group_by(Light_Value) %>%
    group_split()
  
  # Local regression
  fit_reg <- function(df) {
    rankLocReg(
      xall = as.numeric(df$`Delta T [min]`),
      yall = df$Oxygen,
      alpha = 0.2,
      method = "pc",
      verbose = FALSE
    )
  }
  
  fits <- map(intervals, fit_reg)
  
  # Slopes per interval
  micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
  
  # Interval summary table
  interval_df <- map2_dfr(
    intervals,
    seq_along(intervals),
    ~ tibble(
      Colony_ID   = colony_id,
      Run         = run_num,
      Interval    = .y,
      Light_Value = unique(.x$Light_Value),
      N_points    = nrow(.x),
      O2_start    = head(.x$Oxygen, 1),
      O2_end      = tail(.x$Oxygen, 1),
      DeltaT_min  = min(.x$`Delta T [min]`),
      DeltaT_max  = max(.x$`Delta T [min]`),
      O2_min      = min(.x$Oxygen, na.rm = TRUE),
      O2_max      = max(.x$Oxygen, na.rm = TRUE),
      Slope       = micromol.L.s[.y]
    )
  )
  
  # Save file for this colony/run
  fn <- paste0("output/intervals_", colony_id, "_run", run_num, ".csv")
  write_csv(interval_df, fn)
  
  return(interval_df)
}

# 3️⃣ Run for all colony/run combos in parallel
plan(multisession)

all_interval_summaries <- future_pmap_dfr(
  list(colony_runs$Colony_ID, colony_runs$Run),
  fit_colony_intervals
)

# Save combined CSV
write_csv(all_interval_summaries, "output/all_interval_summaries.csv")



# # # Define function to fit and extract parameters for one Colony_ID and Run
# # fit_colony <- function(colony_id, run_num = 1) {
# #   # Filter data for this colony and run, optionally exclude Light_Value == 0 if desired
# #   sample_data <- all_data %>%
# #     filter(Colony_ID == colony_id, Run == run_num)
# #   
# #   # Split by Light_Value intervals
# #   intervals <- sample_data %>%
# #     group_by(Light_Value) %>%
# #     group_split()
# #   
# #   # Fit local linear regression for each interval
# #   fit_reg <- function(df) {
# #     rankLocReg(
# #       xall = as.numeric(df$`Delta T [min]`),
# #       yall = df$Oxygen,
# #       alpha = 0.2,
# #       method = "pc",
# #       verbose = FALSE
# #     )
# #   }
# #   
# #   fits <- map(intervals, fit_reg)
# #   
# #   # Extract slopes
# #   micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
# #   light_levels <- map_dbl(intervals, ~ unique(.x$Light_Value))
# #   
# #   # Prepare data for nonlinear model fitting
# #   pi_data <- tibble(I = light_levels, P = micromol.L.s)
# #   
# #   # Fit Platt model with Rd parameter (catch errors)
# #   fit_result <- tryCatch({
# #     nlsLM(
# #       P ~ Pmax * tanh(alpha * I / Pmax) + Rd,
# #       data = pi_data,
# #       start = list(Pmax = max(pi_data$P), alpha = 0.01, Rd = -0.1),
# #       control = nls.lm.control(maxiter=100)
# #     )
# #   }, error = function(e) NULL)
# #   
# #   if(is.null(fit_result)) {
# #     return(tibble(
# #       Colony_ID = colony_id,
# #       Pmax = NA_real_,
# #       alpha = NA_real_,
# #       Rd = NA_real_,
# #       Ik = NA_real_,
# #       Pgross = NA_real_
# #     ))
# #   }
# #   
# #   pars <- coef(fit_result)
# #   Pmax <- pars["Pmax"]
# #   alpha <- pars["alpha"]
# #   Rd <- pars["Rd"]
# #   Ik <- Pmax / alpha
# #   Pgross <- Pmax + Rd
# #   
# #   tibble(
# #     Colony_ID = colony_id,
# #     Pmax = Pmax,
# #     alpha = alpha,
# #     Rd = Rd,
# #     Ik = Ik,
# #     Pgross = Pgross
# #   )
# # }
# # 
# # # Run for all unique Colony_ID values (remove blanks and NAs)
# # colony_ids <- unique(all_data$Colony_ID)
# # colony_ids <- colony_ids[!is.na(colony_ids) & colony_ids != "blank"]
# # 
# # # Apply the function over all colonies
# # results <- map_dfr(colony_ids, fit_colony)
# # 
# # # Save results to file
# # write.csv(results, "output/photosynthesis_fit_results_Oct2024.csv", row.names = FALSE)
# 
# # Create all combinations of Colony_ID and Run — INCLUDING blanks
# # We'll remove NAs but keep "blank"
# colony_runs <- all_data %>%
#   filter(!is.na(Colony_ID)) %>%
#   distinct(Colony_ID, Run)
# 
# # Function to fit model for one Colony_ID & Run, subtracting blank baseline
# fit_colony_intervals <- function(colony_id, run_num) {
#   
#   # Get sample data
#   sample_data <- all_data %>%
#     filter(Colony_ID == colony_id, Run == run_num)
#   
#   # Get blank data for same run
#   blank_data <- all_data %>%
#     filter(Colony_ID == "blank", Run == run_num)
#   
#   # Skip if no sample data
#   if (nrow(sample_data) == 0) return(NULL)
#   
#   # Blank correction if blank exists
#   if (nrow(blank_data) > 0) {
#     blank_means <- blank_data %>%
#       group_by(Light_Value) %>%
#       summarise(blank_O2 = mean(Oxygen, na.rm = TRUE), 
#                 .groups = "drop")
#     
#     sample_data <- sample_data %>%
#       left_join(blank_means, by = "Light_Value") %>%
#       mutate(Oxygen = Oxygen - blank_O2)
#   }
#   
#   # Split into intervals by light step
#   intervals <- sample_data %>%
#     group_by(Light_Value) %>%
#     group_split()
#   
#   # Local regression for each interval
#   fit_reg <- function(df) {
#     rankLocReg(
#       xall = as.numeric(df$`Delta T [min]`),
#       yall = df$Oxygen,
#       alpha = 0.2,
#       method = "pc",
#       verbose = FALSE
#     )
#   }
#   
#   fits <- map(intervals, fit_reg)
#   
#   # Slopes for each interval
#   micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
#   
#   # Build interval summary table
#   interval_df <- map2_dfr(
#     intervals,
#     seq_along(intervals),
#     ~ tibble(
#       Colony_ID   = colony_id,
#       Run         = run_num,
#       Interval    = .y,
#       Light_Value = unique(.x$Light_Value),
#       N_points    = nrow(.x),
#       O2_start    = head(.x$Oxygen, 1),
#       O2_end      = tail(.x$Oxygen, 1),
#       DeltaT_min  = min(.x$`Delta T [min]`),
#       DeltaT_max  = max(.x$`Delta T [min]`),
#       O2_min      = min(.x$Oxygen, na.rm = TRUE),
#       O2_max      = max(.x$Oxygen, na.rm = TRUE),
#       Slope       = micromol.L.s[.y]
#     )
#   )
#   
#   # Save to CSV for this colony/run
#   interval_file <- paste0("output/intervals_", colony_id, "_run", run_num, ".csv")
#   write_csv(interval_df, interval_file)
#   
#   return(interval_df)
# }
# 
# 
# fit_colony <- function(colony_id, run_num) {
#   
#   # Get sample data
#   sample_data <- all_data %>%
#     filter(Colony_ID == colony_id, Run == run_num)
#   
#   # Get blank data for same run
#   blank_data <- all_data %>%
#     filter(Colony_ID == "blank", Run == run_num)
#   
#   # If no sample data, return NA row
#   if (nrow(sample_data) == 0) {
#     return(tibble(
#       Colony_ID = colony_id, Run = run_num,
#       Pmax = NA_real_, alpha = NA_real_, Rd = NA_real_,
#       Ik = NA_real_, Pgross = NA_real_
#     ))
#   }
#   
#   # If blank exists, subtract its O2 from sample O2 (matching by Light_Value)
#   if (nrow(blank_data) > 0) {
#     # average blank oxygen per light level
#     blank_means <- blank_data %>%
#       group_by(Light_Value) %>%
#       summarise(blank_O2 = mean(Oxygen, na.rm = TRUE), .groups = "drop")
#     
#     # join and adjust sample O2
#     sample_data <- sample_data %>%
#       left_join(blank_means, by = "Light_Value") %>%
#       mutate(Oxygen = Oxygen - blank_O2)
#   }
#   
#   # get intervals of data 
#     intervals <- sample_data %>%
#       group_by(Light_Value) %>%
#       group_split()
#     
#     fit_reg <- function(df) {
#       rankLocReg(
#         xall = as.numeric(df$`Delta T [min]`),
#         yall = df$Oxygen,
#         alpha = 0.2, method = "pc", verbose = FALSE
#       )
#     }
#     
#     fits <- map(intervals, fit_reg)
#     micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
#     light_levels <- map_dbl(intervals, ~ unique(.x$Light_Value))
#     
#     # --- 1️⃣ CREATE & SAVE INTERVAL-LEVEL SUMMARY BEFORE FITTING MODEL
#     
#     interval_df <- map2_dfr(
#       intervals,
#       seq_along(intervals),
#       ~ tibble(
#         Colony_ID = colony_id,
#         Run       = run_num,
#         Interval  = .y,
#         Light_Value = unique(.x$Light_Value),
#         N_points    = nrow(.x),
#         O2_start    = head(.x$Oxygen, 1),
#         O2_end      = tail(.x$Oxygen, 1),
#         DeltaT_min  = min(.x$`Delta T [min]`),
#         DeltaT_max  = max(.x$`Delta T [min]`),
#         O2_min      = min(.x$Oxygen),
#         O2_max      = max(.x$Oxygen),
#         Slope       = micromol.L.s[.y]
#       )
#     )
#     
#     # SAVE to CSV (appends if file exists, or use write_csv for single sample/run)
#     interval_file <- paste0("output/intervals_", colony_id, "_run", run_num, ".csv")
#     write_csv(interval_df, interval_file)
#     
#     # --- 2️⃣ THEN FIT MODEL AS BEFORE
#     
#     pi_data <- tibble(I = light_levels, P = micromol.L.s)
#     
#     fit_result <- tryCatch({
#       nlsLM(
#         P ~ Pmax * tanh(alpha * I / Pmax) + Rd,
#         data = pi_data,
#         start = list(Pmax = max(pi_data$P), alpha = 0.01, Rd = -0.1),
#         control = nls.lm.control(maxiter = 100)
#       )
#     }, error = function(e) NULL)
#     
#     if (is.null(fit_result)) {
#       return(tibble(
#         Colony_ID = colony_id, Run = run_num,
#         Pmax = NA_real_, alpha = NA_real_, Rd = NA_real_,
#         Ik = NA_real_, Pgross = NA_real_
#       ))
#     }
#     
#     pars <- coef(fit_result)
#     Pmax   <- pars["Pmax"]
#     alpha  <- pars["alpha"]
#     Rd     <- pars["Rd"]
#     Ik     <- Pmax / alpha
#     Pgross <- Pmax + Rd
#     
#     tibble(
#       Colony_ID = colony_id, Run = run_num,
#       Pmax = Pmax, alpha = alpha, Rd = Rd,
#       Ik = Ik, Pgross = Pgross
#     )
#   }
#   
# 
# #   
# #   # Split by light step
# #   intervals <- sample_data %>%
# #     group_by(Light_Value) %>%
# #     group_split()
# #   
# #   # Local linear regression for each light step
# #   fit_reg <- function(df) {
# #     rankLocReg(
# #       xall = as.numeric(df$`Delta T [min]`),
# #       yall = df$Oxygen,
# #       alpha = 0.2, method = "pc", verbose = FALSE
# #     )
# #   }
# #   
# #   fits <- map(intervals, fit_reg)
# #   
# #   # Extract slopes for each step
# #   micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
# #   light_levels <- map_dbl(intervals, ~ unique(.x$Light_Value))
# #   pi_data <- tibble(I = light_levels, P = micromol.L.s)
# #   
# #   # Fit Platt model with Rd
# #   fit_result <- tryCatch({
# #     nlsLM(
# #       P ~ Pmax * tanh(alpha * I / Pmax) + Rd,
# #       data = pi_data,
# #       start = list(Pmax = max(pi_data$P), alpha = 0.01, Rd = -0.1),
# #       control = nls.lm.control(maxiter = 100)
# #     )
# #   }, error = function(e) NULL)
# #   
# #   # If fit fails
# #   if (is.null(fit_result)) {
# #     return(tibble(
# #       Colony_ID = colony_id, Run = run_num,
# #       Pmax = NA_real_, alpha = NA_real_, Rd = NA_real_,
# #       Ik = NA_real_, Pgross = NA_real_
# #     ))
# #   }
# #   
# #   # If fit succeeds, extract parameters
# #   pars <- coef(fit_result)
# #   Pmax   <- pars["Pmax"]
# #   alpha  <- pars["alpha"]
# #   Rd     <- pars["Rd"]
# #   Ik     <- Pmax / alpha
# #   Pgross <- Pmax + Rd
# #   
# #   tibble(
# #     Colony_ID = colony_id, Run = run_num,
# #     Pmax = Pmax, alpha = alpha, Rd = Rd,
# #     Ik = Ik, Pgross = Pgross
# #   )
# # }

# 🚀 Run in parallel — still includes blanks in the loop so blank correction works
# plan(multisession)
# 
# results <- future_pmap_dfr(
#   list(colony_runs$Colony_ID, colony_runs$Run),
#   fit_colony
# )
# 
# # Save
# write.csv(results, "output/photosynthesis_fit_results_all_runs_blank_corrected_Oct2024.csv", row.names = FALSE)
# 
