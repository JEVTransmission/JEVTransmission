# Load packages ####
library(pacman)
library(parallel)
library(gridExtra)
library(scales)
p_load(deSolve, tidyverse, ggplot2)
install.packages(c("sensitivity", "lhs"))
library(sensitivity)
library(lhs)
#load("JEV.RData")


# install.packages("htmltools")
# remotes::install_github("rstudio/htmltools")

# definition ####
# Sm: Number of female Culex.tritaeniorhynchus mosquitoes who are susceptible to JEV infection
# Em: Number of female Culex.tritaeniorhynchus mosquitoes who are infected with JEV and not yet infectious to birds
# Im: Number of female Culex.tritaeniorhynchus mosquitoes who are infected with JEV and infectious to birds
# Sb: Number of birds who are susceptible to JEV infection
# Eb: Number of birds who are infected with JEV and not yet infectious to mosquitoes
# Ib: Number of birds who are infected with JEV and infectious to mosquitoes
# Rb: Number of birds who are no longer infectious and currently immune to JEV
# Sh: Number of humans who are susceptible to JEV infection
# Eh: Number of humans who are infected with JEV but still in latency period
# A: Number of humans who are infected with JEV and asymptomatic or have only mild symptoms
# Sev: Number of humans who are infected with JEV and have severe symptoms
# Rh: Number of humans who are no longer infectious and currently immune to JEV

# Define basic dynamic bird-Vector Model ####
seir_mbh_move <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    
    # population and lambda
    ## patch 1
    ### population
    M1 = Sm1+Em1+Im1          # mosquito population
    B1 = Sb1+Eb1+Ib1+Rb1      # bird population
    P1 = Sh1+Eh1+A1+Sev1+Rh1  # human population
    ### lambda
    lambda_m1 = a*c*Ib1/B1
    lambda_b1 = a*(M1/B1)*b*(Im1/M1)
    lambda_h1 = a_h*(M1/P1)*b_h*(Im1/M1)
    
    ## patch 2
    ### population
    M2 = Sm2+Em2+Im2          # mosquito population
    B2 = Sb2+Eb2+Ib2+Rb2      # bird population
    P2 = Sh2+Eh2+A2+Sev2+Rh2  # human population
    ### lambda
    lambda_m2 = a*c*Ib2/B2
    lambda_b2 = a*(M2/B2)*b*(Im2/M2)
    lambda_h2 = a_h*(M2/P2)*b_h*(Im2/M2)
    
    ## patch 3
    ### population
    M3 = Sm3+Em3+Im3          # mosquito population
    B3 = Sb3+Eb3+Ib3+Rb3      # bird population
    P3 = Sh3+Eh3+A3+Sev3+Rh3  # human population
    ### lambda
    lambda_m3 = a*c*(Ib3/B3)
    lambda_b3 = a*(M3/B3)*b*(Im3/M3)
    lambda_h3 = a_h*(M3/P3)*b_h*(Im3/M3)
    
    
    # equations
    ## patch 1
    ### mosquito equations
    dSm1 = mu_m*M1 - lambda_m1*Sm1 - mu_m*Sm1 
    dEm1 = lambda_m1*Sm1 - (gamma_m + mu_m)*Em1
    dIm1 = gamma_m*Em1 - mu_m*Im1
    
    ### bird equations
    dSb1 = max(mu_b * B1 - lambda_b1 * Sb1 + rho_b * Rb1 - mu_b * Sb1 
               - m12 * Sb1 + m21 * Sb2 
               - m13 * Sb1 + m31 * Sb3, 0)
    dEb1 = max(lambda_b1 * Sb1 - (gamma_b + mu_b) * Eb1 
               - m12 * Eb1 + m21 * Eb2 
               - m13 * Eb1 + m31 * Eb3, 0)
    dIb1 = max(gamma_b * Eb1 - (r_b + mu_b) * Ib1 
               - m12 * Ib1 + m21 * Ib2 
               - m13 * Ib1 + m31 * Ib3, 0)
    dRb1 = max(r_b * Ib1 - (rho_b + mu_b) * Rb1 
               - m12 * Rb1 + m21 * Rb2 
               - m13 * Rb1 + m31 * Rb3, 0)
    
    dCInc_bird1 = lambda_b1*Sb1
    
    ### human equations
    dSh1 = mu_h*P1 - lambda_h1*Sh1 + rho_h*Rh1 - mu_h*Sh1 
    dEh1 = lambda_h1*Sh1 - (gamma_h + mu_h)*Eh1 
    dA1  = (1-p)*gamma_h*Eh1 - (r_a + mu_h)*A1 
    dSev1= p*gamma_h*Eh1 - (r_sev + mu_h + mu_sev)*Sev1 
    dRh1 = r_a*A1 + r_sev*Sev1 - (rho_h + mu_h)*Rh1 
    dCInc_human1 = lambda_h1*Sh1
    
    ## patch 2
    ### mosquito equations
    dSm2 = mu_m*M2 - lambda_m2*Sm2 - mu_m*Sm2  
    dEm2 = lambda_m2*Sm2 - (gamma_m + mu_m)*Em2  
    dIm2 = gamma_m*Em2 - mu_m*Im2  
    ### bird equations
    dSb2 = max(mu_b * B2 - lambda_b2 * Sb2 + rho_b * Rb2 - mu_b * Sb2 
               - m21 * Sb2 + m12 * Sb1 
               - m23 * Sb2 + m32 * Sb3,0)
    dEb2 = max(lambda_b2 * Sb2 - (gamma_b + mu_b) * Eb2 
               - m21 * Eb2 + m12 * Eb1 
               - m23 * Eb2 + m32 * Eb3,0)
    dIb2 = max(gamma_b * Eb2 - (r_b + mu_b) * Ib2 
               - m21 * Ib2 + m12 * Ib1 
               - m23 * Ib2 + m32 * Ib3,0)
    dRb2 = max(r_b * Ib2 - (rho_b + mu_b) * Rb2 
               - m21 * Rb2 + m12 * Rb1 
               - m23 * Rb2 + m32 * Rb3,0)
    dCInc_bird2 = lambda_b2 * Sb2
    
    ### human equations
    dSh2 = mu_h*P2 - lambda_h2*Sh2 + rho_h*Rh2 - mu_h*Sh2  
    dEh2 = lambda_h2*Sh2 - (gamma_h + mu_h)*Eh2 
    dA2  = (1-p)*gamma_h*Eh2 - (r_a + mu_h)*A2 
    dSev2= p*gamma_h*Eh2 - (r_sev + mu_h + mu_sev)*Sev2 
    dRh2 = r_a*A2 + r_sev*Sev2 - (rho_h + mu_h)*Rh2 
    dCInc_human2 = lambda_h2*Sh2
    
    ## patch 3
    ### mosquito equations
    dSm3 = mu_m*M3 - lambda_m3*Sm3 - mu_m*Sm3 
    dEm3 = lambda_m3*Sm3 - (gamma_m + mu_m)*Em3 
    dIm3 = gamma_m*Em3 - mu_m*Im3  
    ### bird equations
    dSb3 = max(mu_b * B3 - lambda_b3 * Sb3 + rho_b * Rb3 - mu_b * Sb3 
               - m31 * Sb3 + m13 * Sb1 
               - m32 * Sb3 + m23 * Sb2,0)
    dEb3 = max(lambda_b3 * Sb3 - (gamma_b + mu_b) * Eb3 
               - m31 * Eb3 + m13 * Eb1 
               - m32 * Eb3 + m23 * Eb2,0)
    dIb3 = max(gamma_b * Eb3 - (r_b + mu_b) * Ib3 
               - m31 * Ib3 + m13 * Ib1 
               - m32 * Ib3 + m23 * Ib2,0)
    dRb3 = max(r_b * Ib3 - (rho_b + mu_b) * Rb3 
               - m31 * Rb3 + m13 * Rb1 
               - m32 * Rb3 + m23 * Rb2,0)
    dCInc_bird3 = lambda_b3 * Sb3    
    
    ### human equations
    dSh3 = mu_h*P3 - lambda_h3*Sh3 + rho_h*Rh3 - mu_h*Sh3  
    dEh3 = lambda_h3*Sh3 - (gamma_h + mu_h)*Eh3 
    dA3  = (1-p)*gamma_h*Eh3 - (r_a + mu_h)*A3 
    dSev3= p*gamma_h*Eh3 - (r_sev + mu_h + mu_sev)*Sev3 
    dRh3 = r_a*A3 + r_sev*Sev3 - (rho_h + mu_h)*Rh3 
    dCInc_human3 = lambda_h3*Sh3
    
    output <- c(dSm1, dEm1, dIm1,
                dSb1, dEb1, dIb1, dRb1, dCInc_bird1,
                dSh1, dEh1, dA1, dSev1, dRh1, dCInc_human1,
                dSm2, dEm2, dIm2,
                dSb2, dEb2, dIb2, dRb2, dCInc_bird2,
                dSh2, dEh2, dA2, dSev2, dRh2, dCInc_human2,
                dSm3, dEm3, dIm3,
                dSb3, dEb3, dIb3, dRb3, dCInc_bird3,
                dSh3, dEh3, dA3, dSev3, dRh3, dCInc_human3
    )
    list(output)
  })
  
}

# define parameter sampling function ####
sample_parameters <- function(){
  parms <- c(
        a = runif(1, min = 0.2, max = 0.5),  # Bird feeding rate
        b = runif(1, min = 0.3, max = 0.6),    # Transmission efficiency mosquito to bird
        c = runif(1, min = 0.3, max = 0.6),    # Transmission efficiency bird to mosquito
        
        a_h = runif(1, min = 0.2, max = 0.5),  # Human feeding rate
        b_h = runif(1, min = 0.2, max = 0.5),    # Transmission efficiency mosquito to human
        
        # Mosquito parameters
        mu_m = runif(1, min = 1/60, max = 1/10),  # Mosquito mortality rate
        gamma_m = runif(1, min = 1/10, max = 1/3),  # Mosquito incubation period
        
        # Bird parameters
        mu_b = runif(1, min = 1/(60*365), max = 1/(20*365)),  # Bird mortality rate
        gamma_b = runif(1, min = 1/21, max = 1),  # Bird incubation period
        r_b = runif(1, min = 1/30, max = 1/5),    # Bird recovery rate
        rho_b = runif(1, min = 1/(10*365), max = 1/(1*365)),  # Bird immunity loss rate
        
        # Human parameters
        mu_h = runif(1, min = 1/(86.3*365), max = 1/(82*365)),  # Human mortality rate
        p = runif(1, min = 0.1, max = 0.4),    # Probability of severe symptoms in humans
        gamma_h = runif(1, min = 1/26, max = 1/2),  # Human incubation period
        r_a = runif(1, min = 1/14, max = 1/5),    # Human recovery rate for asymptomatic/mild cases
        r_sev = runif(1, min = 1/(1*365), max = 1/(0.5*365)),  # Recovery rate for severe cases
        mu_sev = runif(1, min = 1/4, max = 2/5),  # Mortality rate for severe cases
        rho_h = runif(1, min = 1/(20*365), max = 1/(10*365)),  # Human immunity loss rate
        
        # Movement rates using gravity model
        m12 = runif(1, min = movement_rate_ranges[1, "Min"], max = movement_rate_ranges[1, "Max"]),
        m21 = runif(1, min = movement_rate_ranges[2, "Min"], max = movement_rate_ranges[2, "Max"]),
        m13 = runif(1, min = movement_rate_ranges[3, "Min"], max = movement_rate_ranges[3, "Max"]),
        m31 = runif(1, min = movement_rate_ranges[4, "Min"], max = movement_rate_ranges[4, "Max"]),
        m23 = runif(1, min = movement_rate_ranges[5, "Min"], max = movement_rate_ranges[5, "Max"]),
        m32 = runif(1, min = movement_rate_ranges[6, "Min"], max = movement_rate_ranges[6, "Max"])

      )
      return(parms)
    }
    


# Input definitions ####
#Initial value
# one first infection in each patch
start_P1_1<-c(Sm1=13826,   Em1=0, Im1=0, # mosquito population 13826
         Sb1=25180-1, Eb1=0, Ib1=1, Rb1=0,  CInc_bird1=0, # bird population 25180
         Sh1=110,     Eh1=0, A1=0 ,  Sev1=0, Rh1=0, CInc_human1=0,  # human population 110
         
         Sm2=349173,  Em2=0, Im2=0, # mosquito population 349173
         Sb2=18481,   Eb2=0, Ib2=0, Rb2=0,  CInc_bird2=0, # bird population 18481
         Sh2=10420,   Eh2=0, A2=0,  Sev2=0, Rh2=0, CInc_human2=0,  # human population 110420
         
         Sm3=3486916, Em3=0, Im3=0, # mosquito population 3486916
         Sb3=84561,   Eb3=0, Ib3=0, Rb3=0,  CInc_bird3=0, # bird population 84561
         Sh3=605420,  Eh3=0, A3=0,  Sev3=0, Rh3=0, CInc_human3=0  # human population 605420
)

start_P2_1<-c(Sm1=13826,   Em1=0, Im1=0, # mosquito population 13826
              Sb1=25180, Eb1=0, Ib1=0, Rb1=0,  CInc_bird1=0, # bird population 25180
              Sh1=110,     Eh1=0, A1=0 ,  Sev1=0, Rh1=0, CInc_human1=0,  # human population 110
              
              Sm2=349173,  Em2=0, Im2=0, # mosquito population 349173
              Sb2=18481-1,   Eb2=0, Ib2=1, Rb2=0,  CInc_bird2=0, # bird population 18481
              Sh2=10420,   Eh2=0, A2=0,  Sev2=0, Rh2=0, CInc_human2=0,  # human population 110420
              
              Sm3=3486916, Em3=0, Im3=0, # mosquito population 3486916
              Sb3=84561,   Eb3=0, Ib3=0, Rb3=0,  CInc_bird3=0, # bird population 84561
              Sh3=605420,  Eh3=0, A3=0,  Sev3=0, Rh3=0, CInc_human3=0  # human population 605420
)

start_P3_1<-c(Sm1=13826,   Em1=0, Im1=0, # mosquito population 13826
              Sb1=25180, Eb1=0, Ib1=0, Rb1=0,  CInc_bird1=0, # bird population 25180
              Sh1=110,     Eh1=0, A1=0 ,  Sev1=0, Rh1=0, CInc_human1=0,  # human population 110
              
              Sm2=349173,  Em2=0, Im2=0, # mosquito population 349173
              Sb2=18481,   Eb2=0, Ib2=0, Rb2=0,  CInc_bird2=0, # bird population 18481
              Sh2=10420,   Eh2=0, A2=0,  Sev2=0, Rh2=0, CInc_human2=0,  # human population 110420
              
              Sm3=3486916, Em3=0, Im3=0, # mosquito population 3486916
              Sb3=84561-1,   Eb3=0, Ib3=1, Rb3=0,  CInc_bird3=0, # bird population 84561
              Sh3=605420,  Eh3=0, A3=0,  Sev3=0, Rh3=0, CInc_human3=0  # human population 605420
)

# Vector of timesteps
times <- seq(0, 365*10, 1)

# run multiple simulations ####
# run simulations 
n_simulations <- 1000 # number of simulations

results_P1_1 <- list()
params_list_P1_1 <- list()
for (i in 1:n_simulations) {
  set.seed(4665 + i)
  parms <- sample_parameters()
  sim <- ode(y = start_P1_1, times = times, func = seir_mbh_move, parms = parms,
             method = "bdf", rtol = 1e-3, atol = 1e-4)
  results_P1_1[[i]] <- as_tibble(as.data.frame(sim)) %>%
    mutate(simulation = i)
  params_list_P1_1[[i]] <- parms
}
####

results_P2_1 <- list()
for (i in 1:n_simulations) {
  set.seed(4665 + i)
  parms <- sample_parameters()
  sim <- ode(y = start_P2_1, times = times, func = seir_mbh_move, parms = parms,
             method = "bdf", rtol = 1e-3, atol = 1e-4)
  results_P2_1[[i]] <- as_tibble(as.data.frame(sim)) %>%
    mutate(simulation = i)
}

results_P3_1 <- list()
for (i in 1:n_simulations) {
  set.seed(4665 + i)
  parms <- sample_parameters()
  sim <- ode(y = start_P3_1, times = times, func = seir_mbh_move, parms = parms,
             method = "bdf", rtol = 1e-3, atol = 1e-4)
  results_P3_1[[i]] <- as_tibble(as.data.frame(sim)) %>%
    mutate(simulation = i)
}

# combine all simulation results
all_results_P1_1 <- bind_rows(results_P1_1)
all_results_P2_1 <- bind_rows(results_P2_1)
all_results_P3_1 <- bind_rows(results_P3_1)

#save(list = ls(), file = "JEV.RData")

# plotting ####
# Function to calculate the mean line
calculate_mean <- function(data, variable) {
  data %>%
    group_by(time) %>%
    summarise(mean_value = mean(!!sym(variable), na.rm = TRUE))
}
# Generate the plots
## infected mosquitoes Ib
### together
ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(aes(y = Im1/1000, group = simulation, colour = "Patch 1"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Im2/1000, group = simulation, colour = "Patch 2"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Im3/1000, group = simulation, colour = "Patch 3"), alpha = 0.4, size = 0.5) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im1"), aes(y = mean_value/1000), colour = "darkred", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im2"), aes(y = mean_value/1000), colour = "darkgreen", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im3"), aes(y = mean_value/1000), colour = "darkblue", size = 1, linetype = "solid", show.legend = FALSE) +
  labs(title = "Infected Mosquitoes",
       x = "Time (days)",
       y = "Number of Infected Mosquitoes \n (in thousands)",
       colour = "Patch") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12),  
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

### separate
#### Patch 1
p1_Im1 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im1"), aes(y = mean_value/1000, colour = "Patch 1 Mean"), size = 0.6) +
  labs(title = "Infected Mosquitoes in Patch 1",
       x = "Time (days)",
       y = "Number of Infected Mosquitoes \n (in thousands)"
       ) +
  scale_colour_manual(values = c("Patch 1 Mean" = "darkred")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 2
p2_Im2 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im2"), aes(y = mean_value/1000, colour = "Patch 2 Mean"), size = 0.6) +
  labs(title = "Infected Mosquitoes in Patch 2",
       x = "Time (days)",
       y = "Number of Infected Mosquitoes \n (in thousands)") +
  scale_colour_manual(values = c("Patch 2 Mean" = "darkgreen")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 3
p3_Im3 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Im3"), aes(y = mean_value/1000, colour = "Patch 3 Mean"), size = 0.6) +
  labs(title = "Infected Mosquitoes in Patch 3",
       x = "Time (days)",
       y = "Number of Infected Mosquitoes \n (in thousands)") +
  scale_colour_manual(values = c("Patch 3 Mean" = "darkblue")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Display the three plots together
grid.arrange(p1_Im1, p2_Im2, p3_Im3, ncol = 1)

## infected birds Ib
### together
ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(aes(y = Ib1/1000, group = simulation, colour = "Patch 1"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Ib2/1000, group = simulation, colour = "Patch 2"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Ib3/1000, group = simulation, colour = "Patch 3"), alpha = 0.4, size = 0.5) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib1"), aes(y = mean_value/1000), colour = "darkred", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib2"), aes(y = mean_value/1000), colour = "darkgreen", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib3"), aes(y = mean_value/1000), colour = "darkblue", size = 1, linetype = "solid", show.legend = FALSE) +
  labs(title = "Infected Birds",
       x = "Time (days)",
       y = "Number of Infected Birds \n (in thousands)",
       colour = "Patch") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12),  
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )
### separate
#### Patch 1
p1_Ib1 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib1"), aes(y = mean_value/1000, colour = "Patch 1 Mean"), size = 0.6) +
  labs(title = "Infected Birds in Patch 1",
       x = "Time (days)",
       y = "Number of Infected Birds \n (in thousands)"
  ) +
  scale_colour_manual(values = c("Patch 1 Mean" = "darkred")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 2
p2_Ib2 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib2"), aes(y = mean_value/1000, colour = "Patch 2 Mean"), size = 0.6) +
  labs(title = "Infected Birds in Patch 2",
       x = "Time (days)",
       y = "Number of Infected Birds \n (in thousands)") +
  scale_colour_manual(values = c("Patch 2 Mean" = "darkgreen")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 3
p3_Ib3 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Ib3"), aes(y = mean_value/1000, colour = "Patch 3 Mean"), size = 0.6) +
  labs(title = "Infected Birds in Patch 3",
       x = "Time (days)",
       y = "Number of Infected Birds \n (in thousands)") +
  scale_colour_manual(values = c("Patch 3 Mean" = "darkblue")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Display the three plots together
grid.arrange(p1_Ib1, p2_Ib2, p3_Ib3, ncol = 1)

## A compartment human
### together
ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(aes(y = A1/1000, group = simulation, colour = "Patch 1"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = A2/1000, group = simulation, colour = "Patch 2"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = A3/1000, group = simulation, colour = "Patch 3"), alpha = 0.4, size = 0.5) +
  geom_line(data = calculate_mean(all_results_P2_1, "A1"), aes(y = mean_value/1000), colour = "darkred", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "A2"), aes(y = mean_value/1000), colour = "darkgreen", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "A3"), aes(y = mean_value/1000), colour = "darkblue", size = 1, linetype = "solid", show.legend = FALSE) +
  labs(title = "Asymptomatic or Mild Symptom Human Case",
       x = "Time (days)",
       y = "Asymptomatic or Mild Symptom Human Case \n (in thousands)",
       colour = "Patch") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12),  
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )
### separate
#### Patch 1
p1_A1 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "A1"), aes(y = mean_value/1000, colour = "Patch 1 Mean"), size = 0.6) +
  labs(title = "Asymptomatic or Mild Symptom Human Case in Patch 1",
       x = "Time (days)",
       y = "Number of Asymptomatic or \n Mild Symptom Human Case \n (in thousands)"
  ) +
  scale_colour_manual(values = c("Patch 1 Mean" = "darkred")) +
  theme_minimal() +
  #ylim(0, 0.0005) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 2
p2_A2 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "A2"), aes(y = mean_value/1000, colour = "Patch 2 Mean"), size = 0.6) +
  labs(title = "Asymptomatic or Mild Symptom Human Case in Patch 2",
       x = "Time (days)",
       y = "Number of Asymptomatic or \n Mild Symptom Human Case \n (in thousands)") +
  scale_colour_manual(values = c("Patch 2 Mean" = "darkgreen")) +
  theme_minimal() +
  #ylim(0, 0.05) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 3
p3_A3 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "A3"), aes(y = mean_value/1000, colour = "Patch 3 Mean"), size = 0.6) +
  labs(title = "Asymptomatic or Mild Symptom Human Case in Patch 3",
       x = "Time (days)",
       y = "Number of Asymptomatic or \n Mild Symptom Human Case \n (in thousands)") +
  scale_colour_manual(values = c("Patch 3 Mean" = "darkblue")) +
  theme_minimal() +
  #ylim(0, 2) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Display the three plots together
grid.arrange(p1_A1, p2_A2, p3_A3, ncol = 1)

## Sev compartment human
### together
ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(aes(y = Sev1/1000, group = simulation, colour = "Patch 1"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Sev2/1000, group = simulation, colour = "Patch 2"), alpha = 0.4, size = 0.5) +
  geom_line(aes(y = Sev3/1000, group = simulation, colour = "Patch 3"), alpha = 0.4, size = 0.5) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev1"), aes(y = mean_value/1000), colour = "darkred", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev2"), aes(y = mean_value/1000), colour = "darkgreen", size = 1, linetype = "solid", show.legend = FALSE) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev3"), aes(y = mean_value/1000), colour = "darkblue", size = 1, linetype = "solid", show.legend = FALSE) +
  labs(title = "Severe Human Case",
       x = "Time (days)",
       y = "Severe Human Case \n (in thousands)",
       colour = "Patch") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12),  
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )
### separate
#### Patch 1
p1_Sev1 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev1"), aes(y = mean_value/1000, colour = "Patch 1 Mean"), size = 0.6) +
  labs(title = "Severe Human Case in Patch 1",
       x = "Time (days)",
       y = "Number of Severe Human Case \n (in thousands)"
  ) +
  scale_colour_manual(values = c("Patch 1 Mean" = "darkred")) +
  theme_minimal() +
  ylim(0, 0.00005) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 2
p2_Sev2 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev2"), aes(y = mean_value/1000, colour = "Patch 2 Mean"), size = 0.6) +
  labs(title = "Severe Human Case in Patch 2",
       x = "Time (days)",
       y = "Number of Severe Human Case \n (in thousands)") +
  scale_colour_manual(values = c("Patch 2 Mean" = "darkgreen")) +
  theme_minimal() +
  ylim(0, 0.005) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Patch 3
p3_Sev3 <- ggplot(all_results_P2_1, aes(x = time)) +
  geom_line(data = calculate_mean(all_results_P2_1, "Sev3"), aes(y = mean_value/1000, colour = "Patch 3 Mean"), size = 0.6) +
  labs(title = "Severe Human Case in Patch 3",
       x = "Time (days)",
       y = "Number of Severe Human Case \n (in thousands)") +
  scale_colour_manual(values = c("Patch 3 Mean" = "darkblue")) +
  theme_minimal() +
  ylim(0, 0.3) +
  theme(
    plot.title = element_text(size = 12), 
    axis.title.x = element_text(size = 9), 
    axis.title.y = element_text(size = 9), 
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 6)
  )

#### Display the three plots together
grid.arrange(p1_Sev1, p2_Sev2, p3_Sev3, ncol = 1)
