# Combine all simulation results
all_results_P1_1 <- bind_rows(results_P1_1)

# Combine parameters into a dataframe
param_df <- bind_rows(params_list_P1_1)

# cumulative incidence human ####
# Sensitivity analysis: Correlate parameters with key outputs (e.g., Cumulative Incidence in Humans)
output_of_interest <- all_results_P1_1 %>% 
  group_by(simulation) %>%
  summarise(CumulativeIncidence = max(CInc_human1 + CInc_human2 + CInc_human3))

# Correlation analysis
cor_results <- cor(param_df, output_of_interest$CumulativeIncidence)

# Output the correlation results
print(cor_results)

# Visualising the sensitivity analysis
cor_df <- as.data.frame(cor_results) %>%
  rownames_to_column(var = "Parameter") %>%
  gather(key = "Output", value = "Correlation", -Parameter)

ggplot(cor_df, aes(x = Parameter, y = Correlation, fill = Output)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Sensitivity Analysis: Correlation of Parameters with Cumulative Incidence",
       x = "Parameter", y = "Correlation with Cumulative Incidence") +
  theme_minimal()

# infected mosquito ####
output_of_interest <- all_results_P1_1 %>% 
  group_by(simulation) %>%
  summarise(InfectedMosquito = max(Im1 + Im2 + Im3))

# Example: Correlation analysis
cor_results <- cor(param_df, output_of_interest$InfectedMosquito)

# Output the correlation results
print(cor_results)

# Visualising the sensitivity analysis
cor_df <- as.data.frame(cor_results) %>%
  rownames_to_column(var = "Parameter") %>%
  gather(key = "Output", value = "Correlation", -Parameter)

ggplot(cor_df, aes(x = Parameter, y = Correlation, fill = Output)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Sensitivity Analysis: Correlation of Parameters with Infected Mosquitoes",
       x = "Parameter", y = "Correlation with Infected Mosquitoes") +
  theme_minimal()

# infected birds ####
output_of_interest <- all_results_P1_1 %>% 
  group_by(simulation) %>%
  summarise(InfectedBirds = max(Ib1 + Ib2 + Ib3))

# Correlation analysis
cor_results <- cor(param_df, output_of_interest$InfectedBirds)

# Output the correlation results
print(cor_results)

# Visualising the sensitivity analysis
cor_df <- as.data.frame(cor_results) %>%
  rownames_to_column(var = "Parameter") %>%
  gather(key = "Output", value = "Correlation", -Parameter)

ggplot(cor_df, aes(x = Parameter, y = Correlation, fill = Output)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Sensitivity Analysis: Correlation of Parameters with Infected Birds",
       x = "Parameter", y = "Correlation with Infected Birds") +
  theme_minimal()

# asymptomatic and mild human case ####
output_of_interest <- all_results_P1_1 %>% 
  group_by(simulation) %>%
  summarise(AsymptomaticMild = max(A1 + A2 + A3))

# Correlation analysis
cor_results <- cor(param_df, output_of_interest$AsymptomaticMild)

# Output the correlation results
print(cor_results)

# Visualising the sensitivity analysis
cor_df <- as.data.frame(cor_results) %>%
  rownames_to_column(var = "Parameter") %>%
  gather(key = "Output", value = "Correlation", -Parameter)

ggplot(cor_df, aes(x = Parameter, y = Correlation, fill = Output)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Sensitivity Analysis: Correlation of Parameters with Asymptomatic \n and Mild Human Cases",
       x = "Parameter", y = "Correlation with Asymptomatic and Mild Human Cases") +
  theme_minimal()

# severe human case ####
output_of_interest <- all_results_P1_1 %>% 
  group_by(simulation) %>%
  summarise(Severe = max(Sev1 + Sev2 + Sev3))

# Correlation analysis
cor_results <- cor(param_df, output_of_interest$Severe)

# Output the correlation results
print(cor_results)

# Visualising the sensitivity analysis
cor_df <- as.data.frame(cor_results) %>%
  rownames_to_column(var = "Parameter") %>%
  gather(key = "Output", value = "Correlation", -Parameter)

ggplot(cor_df, aes(x = Parameter, y = Correlation, fill = Output)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Sensitivity Analysis: Correlation of Parameters with Severe Human Cases",
       x = "Parameter", y = "Correlation with Severe Human Cases") +
  theme_minimal()
