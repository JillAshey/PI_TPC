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

files <- list.files(path = "data/pi_curves/raw_data/Sept2023/", pattern = "*.csv", full.names = TRUE)

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

write_csv(all_data, "output/pi_curves/Sept2023/all_joined_data_Sept2023.csv")

# ggplot(blank_data, aes(x = `Delta T [min]`, y = Oxygen)) +
#   geom_line() +
#   #facet_wrap(~DateTime) +
#   labs(
#     x = "Time",
#     y = paste0("Oxygen (", unique(all_data$`Oxygen Unit`), ")"),
#     color = "Colony"
#   ) +
#   theme_minimal()  

# Make list of all Colony_ID / Run combinations (including blanks)
colony_runs <- all_data %>%
  filter(!is.na(Colony_ID)) %>%
  distinct(Colony_ID, Run)

# Calculate interval slopes and QC for all blanks (ignore run) - at this timepoint, only two out of five runs had blanks. So I am going to average the blanks and then apply the average to all runs 
blank_data <- all_data %>%
  filter(Colony_ID == "blank")

blank_intervals <- blank_data %>%
  group_by(Run, Light_Value) %>%
  filter(!is.na(Colony_ID)) %>%
  group_split()

# Fit regressions and extract slopes for each blank interval
fit_reg <- function(df) {
  rankLocReg(
    xall = as.numeric(df$`Delta T [min]`),
    yall = df$Oxygen,
    alpha = 0.1,
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

# Main function with detailed interval summary and avg blank correction
fit_colony_intervals_avgblank <- function(colony_id, run_num) {
  # Subset sample data
  sample_data <- all_data %>%
    filter(Colony_ID == colony_id, Run == run_num)
  if (nrow(sample_data) == 0) return(NULL)
  
  # Split data by light step
  intervals <- sample_data %>%
    group_by(Light_Value) %>%
    group_split()
  
  fits <- map(intervals, fit_reg)
  micromol.L.s <- map_dbl(fits, ~ pluck(., "allRegs", "b1", 1))
  
  # Build detailed interval summary
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
  ) %>%
    left_join(avg_blank_slopes, by = "Light_Value") %>%
    #mutate(Blank_Corrected_Slope = Slope - Avg_Blank_Slope)
  
  fn <- paste0("output/pi_curves/Sept2023/intervals_", colony_id, "_run", run_num, "_Sept2023.csv")
  write_csv(interval_df, fn)
  
  interval_df
}

# Run for all samples and save big CSV
plan(multisession)

all_interval_summaries <- future_pmap_dfr(
  list(colony_runs$Colony_ID, colony_runs$Run),
  fit_colony_intervals_avgblank
)

write_csv(all_interval_summaries, "output/pi_curves/Sept2023/all_interval_summaries_Sept2023.csv")





###### troubleshooting since the blanks wont run 
print(length(blank_intervals))
print(lapply(blank_intervals, dim))

blank_data <- all_data %>% filter(Colony_ID == "blank")
blank_intervals <- blank_data %>% group_by(Run, Light_Value) %>% group_split()
result <- lapply(blank_intervals, fit_reg)  # Try sequentially, NOT parallel -- started running at 1:18pm

## the result df did not get produced -- I ended it with no output and no errors after running for 20 mins 
test2 <- fit_reg(blank_intervals[[2]])
print(summary(blank_intervals[[2]]$Oxygen))
print(summary(blank_intervals[[2]]$`Delta T [min]`))

test5 <- fit_reg(blank_intervals[[5]])
print(summary(blank_intervals[[5]]$Oxygen))
print(summary(blank_intervals[[5]]$`Delta T [min]`))

test1 <- fit_reg(blank_intervals[[1]]) # does not run 
print(summary(blank_intervals[[1]]$Oxygen))
print(summary(blank_intervals[[1]]$`Delta T [min]`))

test8 <- fit_reg(blank_intervals[[9]])
print(summary(blank_intervals[[9]]$Oxygen))
print(summary(blank_intervals[[9]]$`Delta T [min]`))

## It may be an issue of too many data points in some of the intervals of the blanks. maybe thin the data... There is clearly something weird going on with the blanks. 

