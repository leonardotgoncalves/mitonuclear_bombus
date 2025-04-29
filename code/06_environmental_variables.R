# Load the required package
library(raster)

# Load the CSV with coordinates
coords_df <- read.csv("coordinates.csv")

# Load WorldClim bioclimatic variables (resolution 2.5 arc-min)
bio_stack <- getData('worldclim', var='bio', res=2.5)

# Subset the desired bio variables
selected_bios <- bio_stack[[c(1, 4, 5, 6, 7)]]  # bio1, bio4, bio5, bio6, bio7

# Extract bio values for the given coordinates
bio_values <- extract(selected_bios, coords_df)

# Combine the coordinates with the extracted values
results <- cbind(coords_df, bio_values)

# Write the results to a new CSV
write.csv(results, "bio_variables.csv", row.names = FALSE)

# Optional: View first few rows
head(results)
