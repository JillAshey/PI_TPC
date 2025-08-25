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

files <- list.files(path = "data/pi_curves/raw_data/May2024/", pattern = "*.csv", full.names = TRUE)

# Read, process & join each file ----
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

write_csv(all_data, "output/pi_curves/May2024/all_joined_data_May2024.csv")

# Make list of all Colony_ID / Run combinations (including blanks)
colony_runs <- all_data %>%
  filter(!is.na(Colony_ID)) %>%
  distinct(Colony_ID, Run)

# Function: blank correction + slope extraction + save/return interval summary
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
  fn <- paste0("output/pi_curves/May2024/intervals_", colony_id, "_run", run_num, "_May2024.csv")
  write_csv(interval_df, fn)
  
  return(interval_df)
}

# Run for all colony/run combos in parallel
plan(multisession)

all_interval_summaries <- future_pmap_dfr(
  list(colony_runs$Colony_ID, colony_runs$Run),
  fit_colony_intervals
)

# Save combined CSV
write_csv(all_interval_summaries, "output/pi_curves/May2024/all_interval_summaries_May2024.csv")
