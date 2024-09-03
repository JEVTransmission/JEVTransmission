library(dplyr)

# Function to calculate key metrics for a single simulation
calculate_metrics <- function(simulation_data) {
  metrics <- simulation_data %>%
    summarise(
      peak_infected_mosquitoes = max(Im1 + Im2 + Im3),
      peak_infected_birds = max(Ib1 + Ib2 + Ib3),
      peak_severe_human_cases = max(Sev1 + Sev2 + Sev3),
      total_infected_mosquitoes = sum(Im1 + Im2 + Im3),
      total_infected_birds = sum(Ib1 + Ib2 + Ib3),
      total_severe_human_cases = sum(Sev1 + Sev2 + Sev3),
      time_to_peak_infected_mosquitoes = which.max(Im1 + Im2 + Im3),
      time_to_peak_infected_birds = which.max(Ib1 + Ib2 + Ib3),
      time_to_peak_severe_human_cases = which.max(Sev1 + Sev2 + Sev3)
    )
  return(metrics)
}

# Calculate metrics for each simulation in all three regions
metrics_P1_1 <- bind_rows(lapply(results_P1_1, calculate_metrics))
metrics_P2_1 <- bind_rows(lapply(results_P2_1, calculate_metrics))
metrics_P3_1 <- bind_rows(lapply(results_P3_1, calculate_metrics))

# Aggregate metrics across all simulations for each region
aggregate_metrics <- function(metrics) {
  summary_stats <- metrics %>%
    summarise(
      mean_peak_infected_mosquitoes = mean(peak_infected_mosquitoes),
      sd_peak_infected_mosquitoes = sd(peak_infected_mosquitoes),
      mean_peak_infected_birds = mean(peak_infected_birds),
      sd_peak_infected_birds = sd(peak_infected_birds),
      mean_peak_severe_human_cases = mean(peak_severe_human_cases),
      sd_peak_severe_human_cases = sd(peak_severe_human_cases),
      mean_time_to_peak_infected_mosquitoes = mean(time_to_peak_infected_mosquitoes),
      sd_time_to_peak_infected_mosquitoes = sd(time_to_peak_infected_mosquitoes),
      mean_time_to_peak_infected_birds = mean(time_to_peak_infected_birds),
      sd_time_to_peak_infected_birds = sd(time_to_peak_infected_birds),
      mean_time_to_peak_severe_human_cases = mean(time_to_peak_severe_human_cases),
      sd_time_to_peak_severe_human_cases = sd(time_to_peak_severe_human_cases)
    )
  return(summary_stats)
}

summary_P1_1 <- aggregate_metrics(metrics_P1_1)
summary_P2_1 <- aggregate_metrics(metrics_P2_1)
summary_P3_1 <- aggregate_metrics(metrics_P3_1)

# Combine all metrics for ANOVA or t-tests
all_metrics <- bind_rows(
  metrics_P1_1 %>% mutate(region = "Region 1"),
  metrics_P2_1 %>% mutate(region = "Region 2"),
  metrics_P3_1 %>% mutate(region = "Region 3")
)

# Perform ANOVA to compare regions for each metric
anova_peak_infected_mosquitoes <- aov(peak_infected_mosquitoes ~ region, data = all_metrics)
anova_peak_infected_birds <- aov(peak_infected_birds ~ region, data = all_metrics)
anova_peak_severe_human_cases <- aov(peak_severe_human_cases ~ region, data = all_metrics)

# Perform t-tests if necessary (pairwise comparisons)
# For example, comparing peak infected mosquitoes between regions 1 and 2
t_test_peak_infected_mosquitoes_R1_R2 <- t.test(peak_infected_mosquitoes ~ region, 
                                                data = all_metrics %>% filter(region %in% c("Region 1", "Region 2")))

# peak amount ####
# Boxplot for peak infected mosquitoes
ggplot(all_metrics, aes(x = region, y = peak_infected_mosquitoes)) +
  geom_boxplot() +
  labs(title = "Peak Infected Mosquitoes by Initial Infection Region", y = "Peak Infected Mosquitoes", x = "Region")

# Boxplot for peak infected birds
ggplot(all_metrics, aes(x = region, y = peak_infected_birds)) +
  geom_boxplot() +
  labs(title = "Peak Infected Birds by Initial Infection Region", y = "Peak Infected Birds", x = "Region")

# Boxplot for peak severe human cases
ggplot(all_metrics, aes(x = region, y = peak_severe_human_cases)) +
  geom_boxplot() +
  labs(title = "Peak Severe Human Cases by Initial Infection Region", y = "Peak Severe Human Cases", x = "Region")
# peak time ####
# Boxplot for time to peak infected mosquitoes
ggplot(all_metrics, aes(x = region, y = time_to_peak_infected_mosquitoes)) +
  geom_boxplot() +
  labs(title = "Time to Peak Infected Mosquitoes by Initial Infection Region", y = "Time to Peak (Days)", x = "Region")

# Boxplot for time to peak infected birds
ggplot(all_metrics, aes(x = region, y = time_to_peak_infected_birds)) +
  geom_boxplot() +
  labs(title = "Time to Peak Infected Birds by Initial Infection Region", y = "Time to Peak (Days)", x = "Region")

# Boxplot for time to peak severe human cases
ggplot(all_metrics, aes(x = region, y = time_to_peak_severe_human_cases)) +
  geom_boxplot() +
  labs(title = "Time to Peak Severe Human Cases by Initial Infection Region", y = "Time to Peak (Days)", x = "Region")

