library(dplyr)
library(ggplot2)
setwd("Z:/home/users/grenova/K-tag/IHC/Spleen ASC architecture/Quantification/Distance to nearest PC")

# Read files and combine into one dataframe
samples <- list(
  Sample1 = "Sample1_WT.txt",
  Sample2 = "Sample2_WT.txt",
  Sample3 = "Sample3_Kika.txt",
  Sample4 = "Sample4_WT.txt",
  Sample5= "Sample5_Kika2.txt",
  Sample6 = "EHP309_Sample1_WT.txt",
  Sample7 = "EHP309_Sample2_Kika.txt"
)

data <- lapply(names(samples), function(x) {
  df <- read.delim(samples[[x]], header = TRUE, sep = "\t")
  data.frame(
    Sample = x,
    Distance = df$`Nearest.CD138.distance..µm.`
  )
}) %>%
  bind_rows() %>%
  filter(!is.na(Distance))


# Assign groups
data <- data %>%
  mutate(Group = case_when(
    Sample %in% c("Sample1", "Sample2", "Sample4", "Sample6") ~ "WT",
    Sample %in% c("Sample3", "Sample5", "Sample7") ~ "Kika"
  ))


# Create distance bins
bin_width <- 2

binned <- data %>%
  filter(Distance <= 50) %>%
  mutate(Bin = cut(
    Distance,
    breaks = seq(0, 50, by = bin_width),
    include.lowest = TRUE
  )) %>%
  group_by(Sample, Group, Bin) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  group_by(Sample) %>%
  mutate(Relative_frequency = Count / sum(Count))


# Calculate mean and SD between samples
summary_data <- binned %>%
  group_by(Group, Bin) %>%
  summarise(
    Mean = mean(Relative_frequency),
    SD = sd(Relative_frequency),
    .groups = "drop"
  )


# Extract midpoint of bins
summary_data$Distance <- 
  as.numeric(sub("\\((.+),(.+)\\]", "\\1", summary_data$Bin)) + bin_width/2


# Plot mean ± SD
ggplot(summary_data, aes(x = Distance, y = Mean, color = Group, fill = Group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    alpha = 0.25,
    color = NA
  ) +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 0.4)) +
  theme_classic() +
  labs(
    x = "Distance to nearest CD138 cell (µm)",
    y = "Relative frequency",
    title = "Distance to CD138: Mean ± SD"
  )

ggsave("CD138_distance.png", plot = last_plot(), width = 10, height = 9, dpi = 300)
ggsave("CD138_distance.pdf", plot = last_plot(), width = 10, height = 9, dpi = 300)
