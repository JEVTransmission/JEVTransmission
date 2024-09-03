install.packages("terra")
install.packages("sf")

library(terra)
library(sf)
library(ggplot2)

# load regions shapefile
regions <- vect("Singapore_3Regions_fixed.gpkg.shp")

# vegetation import
vegetation_raster <- rast("Vegetation_Map_10m_DRichards.tif")
print(vegetation_raster)
plot(vegetation_raster)
# Save the plot to a PNG file
dev.copy(png, file = "vegetation_raster_plot.png")  # Copies the current plot to a PNG file
dev.off()
mean_vegetation_by_region <- extract(vegetation_raster, regions, fun=mean, na.rm=TRUE)
print(mean_vegetation_by_region)

# human density import
mean_PopulationDensity_by_region <- read.csv("PopulationDensity.csv")
print(mean_PopulationDensity_by_region)

# region area import
region_area_by_region <- read.csv("RegionArea.csv")
print(region_area_by_region)

# normalize
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Normalize the factors
mean_vegetation_norm <- normalize(mean_vegetation_by_region$Vegetation_Map_10m_DRichards)
mean_human_density_norm <- normalize(mean_PopulationDensity_by_region$PopulationDensity)
region_area_norm <- normalize(region_area_by_region$Area)

# Combine the normalized data into a data frame
data_normalized <- data.frame(
  ID = mean_vegetation_by_region$ID,
  Vegetation = mean_vegetation_norm,
  HumanDensity = mean_human_density_norm,
  Area = region_area_norm
)

print(data_normalized)

# gravity model ####
# Assuming data_normalized already contains normalized data
# Bird populations in each region
B1 <- 25180   # Example population for Region 1
B2 <- 18481   # Example population for Region 2
B3 <- 84561   # Example population for Region 3

# Distances between regions
D12 <- 12
D13 <- 18
D23 <- 15

# Gravity model parameters
K <- 1/100000000  # Proportionality constant
b <- 2  # Distance exponent

# Number of simulations
n_simulations <- 1000

# Store results
movement_rates <- matrix(NA, nrow = n_simulations, ncol = 6)
colnames(movement_rates) <- c("m12", "m21", "m13", "m31", "m23", "m32")

for (i in 1:n_simulations) {
  # Assign random weights within the specified ranges
  w1 <- runif(1, min = 0, max = 1)
  w2 <- runif(1, min = -1, max = 0)
  w3 <- runif(1, min = 0, max = 1)
  
  # Calculate A_i for each region with a minimum threshold to avoid negative values
  data_normalized$A_i <- pmax(0, w1 * data_normalized$Vegetation +
                                w2 * data_normalized$HumanDensity +
                                w3 * data_normalized$Area)
  
  A12 <- data_normalized$A_i[2]
  A21 <- data_normalized$A_i[1]
  A13 <- data_normalized$A_i[3]
  A31 <- data_normalized$A_i[1]
  A23 <- data_normalized$A_i[3]
  A32 <- data_normalized$A_i[2]
  
  # # Calculate movement rates using gravity model with adjusted K and capping
  # movement_rates[i, "m12"] <- pmin(K * A12 * (B1 * B2) / (D12^b), 0.05)
  # movement_rates[i, "m21"] <- pmin(K * A21 * (B1 * B2) / (D12^b), 0.05)
  # movement_rates[i, "m13"] <- pmin(K * A13 * (B1 * B3) / (D13^b), 0.05)
  # movement_rates[i, "m31"] <- pmin(K * A31 * (B1 * B3) / (D13^b), 0.05)
  # movement_rates[i, "m23"] <- pmin(K * A23 * (B2 * B3) / (D23^b), 0.05)
  # movement_rates[i, "m32"] <- pmin(K * A32 * (B2 * B3) / (D23^b), 0.05)
  
  # Calculate movement rates using gravity model with adjusted K and capping
  movement_rates[i, "m12"] <- K * A12 * (B1 * B2) / (D12^b)
  movement_rates[i, "m21"] <- K * A21 * (B1 * B2) / (D12^b)
  movement_rates[i, "m13"] <- K * A13 * (B1 * B3) / (D13^b)
  movement_rates[i, "m31"] <- K * A31 * (B1 * B3) / (D13^b)
  movement_rates[i, "m23"] <- K * A23 * (B2 * B3) / (D23^b)
  movement_rates[i, "m32"] <- K * A32 * (B2 * B3) / (D23^b)
}

# Calculate the mean movement rates
mean_movement_rates <- colMeans(movement_rates)

# Print the mean movement rates
print(mean_movement_rates)

# Calculate the mean and standard deviation for each movement rate
movement_rate_stats <- data.frame(
  Mean = apply(movement_rates, 2, mean),
  SD = apply(movement_rates, 2, sd)
)

# Print the calculated means and standard deviations
print(movement_rate_stats)

# Define a reasonable range using mean ± 1 SD
movement_rate_ranges <- data.frame(
  Min = movement_rate_stats$Mean - movement_rate_stats$SD,
  Max = movement_rate_stats$Mean + movement_rate_stats$SD
)

# Ensure the minimum is not negative
movement_rate_ranges$Min <- pmax(0, movement_rate_ranges$Min)

print(movement_rate_ranges)
#save(data_normalized, file = "JEV.RData")


