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

getwd()

# Load the metadata ----
run_meta <- read_csv("data/pi_curves/pi_curve_metadata_run.csv") %>%
  mutate(
    Start = as_datetime(Date) + hours(hour(Start.time)) + minutes(minute(Start.time)),
    Stop  = as_datetime(Date) + hours(hour(Stop.time)) + minutes(minute(Stop.time)))

sample_meta <- read_csv("data/pi_curves/pi_curve_metadata_sample.csv")

# Make sure Chamber well # is numeric
sample_meta <- sample_meta %>%
  mutate(`Chamber` = as.numeric(`Chamber`))

# List all raw data files ----

files <- list.files(path = "data/pi_curves/raw_data/Jan2024/", pattern = "*.csv", full.names = TRUE)

# List all raw data files ----

all_data <- map_dfr(files, function(f) {
  # Extract channel number from filename, e.g. Ch1 → 1
  channel <- str_extract(f, "channel\\d+") %>% str_remove("channel") %>% as.numeric()
  
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

# Inspect and save ----
all_data %>% count(is.na(`Colony_ID`))
all_data %>% count(is.na(Run))

write_csv(all_data, "output/pi_curves/Jan2024/all_joined_data_Jan2024.csv")

# Make list of all Colony_ID / Run combinations (including blanks)
colony_runs <- all_data %>%
  filter(!is.na(Colony_ID)) %>%
  distinct(Colony_ID, Run)

# Calculate interval slopes and QC for all blanks (ignore run) - at this timepoint, only two out of five runs had blanks. So I am going to average the blanks and then apply the average to all runs 
blank_data <- all_data %>%
  filter(Colony_ID == "blank")

blank_intervals <- blank_data %>%
  group_by(Run, Light_Value) %>%
  group_split()

# Fit regressions and extract slopes for each blank interval
fit_reg <- function(df) {
  rankLocReg(
    xall = as.numeric(df$`Delta T [min]`),
    yall = df$Oxygen,
    alpha = 0.2,
    method = "pc",
    verbose = FALSE
  )
}

blank_fits <- map(blank_intervals, fit_reg)
blank_slopes <- map_dbl(blank_fits, ~ pluck(., "allRegs", "b1", 1))
blank_summary <- map2_dfr(
  blank_intervals,
  blank_slopes,
  ~ tibble(
    Light_Value = unique(.x$Light_Value),
    Blank_Slope = .y
  )
)

# Average blank slopes per light value, across all blanks/runs
avg_blank_slopes <- blank_summary %>%
  group_by(Light_Value) %>%
  summarise(Avg_Blank_Slope = mean(Blank_Slope, na.rm = TRUE), .groups = "drop")






# Function: blank correction + slope extraction + save/return interval summary
fit_colony_intervals_avgblank <- function(colony_id, run_num) {
  # Sample data for this colony/run
  sample_data <- all_data %>%
    filter(Colony_ID == colony_id, Run == run_num)
  if (nrow(sample_data) == 0) return(NULL)
  
  # Split into intervals by light step
  intervals <- sample_data %>%
    group_by(Light_Value) %>%
    group_split()
  
  # Fit local regression for each light step
  fits <- map(intervals, fit_reg)
  micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
  light_levels <- map_dbl(intervals, ~ unique(.x$Light_Value))
  
  # Build interval summary table
  interval_df <- tibble(
    Colony_ID   = colony_id,
    Run         = run_num,
    Light_Value = light_levels,
    Raw_Slope   = micromol.L.s
  ) %>%
    left_join(avg_blank_slopes, by = "Light_Value") %>%
    mutate(Blank_Corrected_Slope = Raw_Slope - Avg_Blank_Slope)
  
  # Save per colony/run (optional)
  fn <- paste0("output/pi_curves/Jan2024/intervals_", colony_id, "_run", run_num, "_Jan2024.csv")
  write_csv(interval_df, fn)
  
  return(interval_df)
}

# Process all samples and save combined output
plan(multisession)

all_interval_summaries <- future_pmap_dfr(
  list(colony_runs$Colony_ID, colony_runs$Run),
  fit_colony_intervals_avgblank
)

write_csv(all_interval_summaries, "output/pi_curves/Jan2024/all_interval_summaries_Jan2024.csv")


