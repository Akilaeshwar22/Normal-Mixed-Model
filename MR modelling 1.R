# MR modelling - 1 (Normal mixed model)

# Load necessary libraries
library(tidyverse) # For data manipulation and plotting
library(mixtools)  # For mixture modeling
library(lme4)      # For mixed-effects models
library(readxl)

# --- Step 1: Data Loading and Cleaning ---

# Read datasets
data_2013 <- read_excel("D:\\CMC-WTRL\\MR\\MR modelling\\Mixture modelling.xlsx", sheet = "2013")
data_2022 <- read_excel("D:\\CMC-WTRL\\MR\\MR modelling\\Mixture modelling.xlsx", sheet = "2022")


# Combine the two datasets into a single data frame
my_data <- bind_rows(data_2013, data_2022)

# Ensure no negative antibody titres
my_data <- my_data %>%
  mutate(antibodytitre = pmax(0.1, antibodytitre))

# --- Step 2: Define Serostatus Based on Manufacturer's Cutoffs ---

# Create a new serostatus variable with three categories (Negative, Equivocal, Positive)
my_data <- my_data %>%
  mutate(serostatus = case_when(
    disease == "Measles" & antibodytitre < 9 ~ "Negative",
    disease == "Measles" & antibodytitre >= 9 & antibodytitre <= 11 ~ "Equivocal",
    disease == "Measles" & antibodytitre > 11 ~ "Positive",
    disease == "Rubella" & antibodytitre < 8 ~ "Negative",
    disease == "Rubella" & antibodytitre >= 8 & antibodytitre <= 11 ~ "Equivocal",
    disease == "Rubella" & antibodytitre > 11 ~ "Positive",
    TRUE ~ NA_character_ # Catch any other cases
  ))

my_data %>%
  filter(year == 2013, is.na(serostatus)) %>%
  select(disease, antibodytitre)

my_data %>%
  filter(year == 2013, is.na(cluster)) %>%
  select(disease, antibodytitre, cluster)

my_data %>%
  filter(year == 2013, disease == "Rubella", is.na(antibodytitre))

# --- Step 3: Descriptive and Exploratory Analysis ---

# Summary statistics including the new serostatus proportions
cat("--- Summary Statistics by Serostatus ---\n")
print(my_data %>%
        group_by(year, disease, serostatus) %>%
        summarise(n = n()))

# Create a folder to store the plots if it doesn't already exist
folder_path <- "D:\\CMC-WTRL\\MR\\MR modelling\\Plots"
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
}

# Plot the distribution of antibody titres
distribution_abtitre <- my_data %>%
  ggplot(aes(x = antibodytitre, fill = as.factor(year))) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  facet_wrap(~ disease, scales = "free") +
  labs(
    title = "Distribution of antibody titres by disease",
    x = "Antibody titre (Log10 scale)",
    y = "Density",
    fill = "Year"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(face = "bold", colour = "black", size = 12),
    axis.title = element_text(face = "bold", colour = "black", size = 14),
    plot.title = element_text(face = "bold", colour = "black", size = 16),
    legend.title = element_text(face = "bold", colour = "black", size = 12),
    legend.text = element_text(face = "bold", colour = "black", size = 10),
    # Add styling for facet text
    strip.text = element_text(face = "bold", colour = "black", size = 12)
  ); distribution_abtitre

# Save the age group plot as a high-quality TIFF file
ggsave(
  filename = file.path(folder_path, "Antibody titres.tiff"),
  plot = distribution_abtitre,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300
)

# --- Step 4: Mixture Modeling and Re Calculation ---

# Define R0 values for each disease
R0_values <- c("Measles" = 12, "Rubella" = 5)

# Create a data frame to store the results
results_df <- tibble(
  year = numeric(),
  disease = character(),
  seronegative_prop = numeric(),
  Re = numeric()
)

# Loop through each disease and year to perform the analysis
for (d in c("Measles", "Rubella")) {
  for (y in c(2013, 2022)) {
    data_subset <- my_data %>%
      filter(disease == d, year == y)
    
    if (nrow(data_subset) > 5) {
      set.seed(123)
      mix_model <- normalmixEM(log(data_subset$antibodytitre), k = 2, fast = TRUE)
      seronegative_prop <- min(mix_model$lambda)
      immune_prop <- 1 - seronegative_prop
      Re <- R0_values[[d]] * (1 - immune_prop)
      
      results_df <- results_df %>%
        add_row(
          year = y,
          disease = d,
          seronegative_prop = seronegative_prop,
          Re = Re
        )
    }
  }
}

cat("\n--- Mixture Model and Re Results ---\n")
print(results_df)


# --- Step 5: Geographic Cluster Analysis with Multilevel Model ---

# For the multilevel model, we must exclude the equivocal cases
# as the model requires a binary outcome (Positive/Negative)
my_data_binary <- my_data %>%
  filter(serostatus != "Equivocal") %>%
  mutate(seronegative_binary = ifelse(serostatus == "Negative", 1, 0))

# Loop through each disease and fit a multilevel model
for (d in c("Measles", "Rubella")) {
  cat(paste0("\n--- Fitting multilevel model for ", toupper(d), " using manufacturer's cutoffs ---\n"))
  data_subset <- my_data_binary %>% filter(disease == d)
  
  if (nrow(data_subset) > 5) {
    model_with_clusters <- glmer(
      seronegative_binary ~ factor(year) + age + gender + (1 | cluster),
      data = data_subset,
      family = binomial
    )
    print(summary(model_with_clusters))
  }
}


# --- Step 6: Data Visualization ---

# --- Step 3: Corrected Data Visualization ---
# Plot serostatus proportions by age group (stacked bar chart)
my_data_binned <- my_data %>%
  mutate(age_group = cut(age,
                         breaks = seq(0, 80, by = 5),
                         right = FALSE,
                         labels = paste0(seq(0, 75, by = 5), "-", seq(4, 79, by = 5))))

serostatus_by_age <- my_data_binned %>%
  filter(!is.na(serostatus)) %>%
  group_by(year, disease, age_group, serostatus) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(year, disease, age_group) %>%
  mutate(prop = n / sum(n))

# Assign the plot to a variable named `plot_age`
plot_age <- ggplot(serostatus_by_age, aes(x = age_group, y = prop, fill = serostatus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(disease ~ year, scales = "free_x") +
  labs(
    title = "Serostatus proportions by age group, year, and disease",
    x = "Age group (years)",
    y = "Proportion",
    fill = "Serostatus"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Negative" = "tomato", "Equivocal" = "gold", "Positive" = "seagreen")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(face = "bold", colour = "black", size = 12),
    axis.title = element_text(face = "bold", colour = "black", size = 14),
    plot.title = element_text(face = "bold", colour = "black", size = 16),
    legend.title = element_text(face = "bold", colour = "black", size = 12),
    legend.text = element_text(face = "bold", colour = "black", size = 10),
    # Add styling for facet text
    strip.text = element_text(face = "bold", colour = "black", size = 12)
  );plot_age 

# Plot serostatus proportions by cluster
serostatus_by_cluster <- my_data %>%
  filter(!is.na(serostatus)) %>%
  group_by(year, disease, cluster, serostatus) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(year, disease, cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = paste0("C", 1:12)))

# Assign the plot to a variable named `plot_cluster`
plot_cluster <- ggplot(serostatus_by_cluster, aes(x = cluster, y = prop, fill = serostatus)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(disease ~ year, scales = "free_x") +
  labs(
    title = "Serostatus proportions by geographic cluster",
    x = "Cluster",
    y = "Proportion",
    fill = "Serostatus"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Negative" = "tomato", "Equivocal" = "gold", "Positive" = "seagreen")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(face = "bold", colour = "black", size = 12),
    axis.title = element_text(face = "bold", colour = "black", size = 14),
    plot.title = element_text(face = "bold", colour = "black", size = 16),
    legend.title = element_text(face = "bold", colour = "black", size = 12),
    legend.text = element_text(face = "bold", colour = "black", size = 10),
    # Add styling for facet text
    strip.text = element_text(face = "bold", colour = "black", size = 12)
  );plot_cluster


# Save the age group plot using the correct plot variable
ggsave(
  filename = file.path(folder_path, "proportions by age.tiff"),
  plot = plot_age,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300
)

# Save the geographic cluster plot using the correct plot variable
ggsave(
  filename = file.path(folder_path, "proportions by cluster.tiff"),
  plot = plot_cluster,
  device = "tiff",
  width = 12,
  height = 8,
  dpi = 300
)



# --- Simple Mixture Model Plots ---

# Make a folder for mixture plots
plot_folder <- "D:/CMC-WTRL/MR/MR modelling/Plots"
if (!dir.exists(plot_folder)) dir.create(plot_folder)

for (d in c("Measles", "Rubella")) {
  for (y in c(2013, 2022)) {
    data_subset <- my_data %>%
      filter(disease == d, year == y, !is.na(antibodytitre))
    
    if (nrow(data_subset) > 5) {
      set.seed(123)
      mix_model <- normalmixEM(log(data_subset$antibodytitre), k = 2, fast = TRUE)
      
      # Base histogram of log titres
      hist(log(data_subset$antibodytitre),
           breaks = 30, freq = FALSE, col = "lightgray",
           main = paste(d, y, "Mixture Model"),
           xlab = "Log antibody titre")
      
      # Add the two fitted normal curves
      curve(mix_model$lambda[1] * dnorm(x, mix_model$mu[1], mix_model$sigma[1]),
            add = TRUE, col = "red", lwd = 2)
      curve(mix_model$lambda[2] * dnorm(x, mix_model$mu[2], mix_model$sigma[2]),
            add = TRUE, col = "blue", lwd = 2)
      
      # Add the combined mixture density
      curve(mix_model$lambda[1] * dnorm(x, mix_model$mu[1], mix_model$sigma[1]) +
              mix_model$lambda[2] * dnorm(x, mix_model$mu[2], mix_model$sigma[2]),
            add = TRUE, col = "darkgreen", lwd = 2, lty = 2)
      
      # Save each plot
      fname <- file.path(plot_folder, paste0(d, "_", y, "_mixture.png"))
      dev.copy(png, fname, width = 800, height = 600)
      dev.off()
    }
  }
}

ggplot(my_data, aes(sample = log(antibodytitre))) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~ disease + year, scales = "free") +
  labs(title = "Q-Q plots of log antibody titres")

library(rstatix)

my_data %>%
  mutate(log_titre = log(antibodytitre)) %>%   # create a new column
  group_by(disease, year) %>%
  shapiro_test(log_titre)                      # run test on that column

library(moments)

my_data %>%
  group_by(disease, year) %>%
  summarise(
    skewness = skewness(log(antibodytitre), na.rm = TRUE),
    kurtosis = kurtosis(log(antibodytitre), na.rm = TRUE)
  )

ggplot(my_data, aes(x = log(antibodytitre))) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  facet_wrap(~ disease + year, scales = "free") +
  labs(title = "Histogram of log(antibody titres) with density",
       x = "Log(antibody titre)", y = "Density")

# Load necessary libraries
# You've already loaded tidyverse, so this is just a reminder
library(tidyverse)

# For Measles, ensure "Measles" is capitalized
measles_data <- my_data %>%
  filter(disease == "Measles") %>%
  mutate(log_titre = log(antibodytitre))

# Now, your plotting code will work as intended
ggplot(measles_data, aes(x = log_titre)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightblue", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(measles_data$log_titre), sd = sd(measles_data$log_titre)), color = "red", size = 1.2) +
  facet_wrap(~year) +
  labs(title = "Distribution of Log-Transformed Measles Antibody Titres",
       x = "Log(Antibody Titre)",
       y = "Density") +
  theme_minimal()

# For Rubella, the same principle applies
rubella_data <- my_data %>%
  filter(disease == "Rubella") %>%
  mutate(log_titre = log(antibodytitre))

ggplot(rubella_data, aes(x = log_titre)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.2, fill = "lightgreen", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(rubella_data$log_titre), sd = sd(rubella_data$log_titre)), color = "red", size = 1.2) +
  facet_wrap(~year) +
  labs(title = "Distribution of Log-Transformed Rubella Antibody Titres",
       x = "Log(Antibody Titre)",
       y = "Density") +
  theme_minimal()
