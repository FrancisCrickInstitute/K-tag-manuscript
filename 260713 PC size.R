library(readr)
library(dplyr)
library(ggplot2)
setwd("Z:/home/users/grenova/K-tag/IHC/Spleen ASC architecture/Quantification/PC size")


# Read the QuPath .txt file
Sample5 <- read_tsv("EHP596_Sample5_Kika.txt")
Sample3 <- read_tsv("EHP596_Sample3_Kika.txt")
Sample2 <- read_tsv("EHP309_Sample2_Kika.txt")

all_samples <- bind_rows(
  Sample5 %>% mutate(Sample = "Sample5"),
  Sample3 %>% mutate(Sample = "Sample3"),
  Sample2 %>% mutate(Sample = "Sample2")
)

colnames(all_samples)

unique(Sample5$Parent)

# Keep only Edge and Centre objects, and create a grouping variable
Sample5 <- Sample5 %>%
  filter(grepl("^(Edge|Centre)", Parent)) %>%
  mutate(Group = ifelse(grepl("^Edge", Parent), "Edge", "Centre"))
Sample3 <- Sample3 %>%
  filter(grepl("^(Edge|Centre)", Parent)) %>%
  mutate(Group = ifelse(grepl("^Edge", Parent), "Edge", "Centre"))
Sample2 <- Sample2 %>%
  filter(grepl("^(Edge|Centre)", Parent)) %>%
  mutate(Group = ifelse(grepl("^Edge", Parent), "Edge", "Centre"))
all_samples <- all_samples %>%
  filter(grepl("^(Edge|Centre)", Parent)) %>%
  mutate(Group = ifelse(grepl("^Edge", Parent), "Edge", "Centre"))

colnames(Sample5)
unique(Sample3$Parent)
unique(Sample3$Parent)

# Plot histogram
Plot_S5 <- ggplot(Sample5, aes(x = `Nucleus: Area`, fill = Group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(
    x = "Nucleus: Area",
    y = "Density",
    title = "Distribution of Nucleus Area"
  ) +
  theme_classic()

Plot_S3 <- ggplot(Sample3, aes(x = `Nucleus: Area`, fill = Group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(
    x = "Nucleus: Area",
    y = "Density",
    title = "Distribution of Nucleus Area"
  ) +
  theme_classic()

Plot_S2 <- ggplot(Sample2, aes(x = `Nucleus: Area`, fill = Group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(
    x = "Nucleus: Area",
    y = "Density",
    title = "Distribution of Nucleus Area"
  ) +
  theme_classic()

Plot_all <- ggplot(all_samples, aes(x = `Nucleus: Area`, fill = Group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.5,
    position = "identity"
  ) +
  coord_cartesian(xlim = c(0, 200)) +
  labs(
    x = "Nucleus: Area",
    y = "Density",
    title = "Distribution of Nucleus Area"
  ) +
  theme_classic()

ggsave("PC_size_histogram_S2.pdf", Plot_S2)
ggsave("PC_size_histogram_S3.pdf", Plot_S3)
ggsave("PC_size_histogram_S5.pdf", Plot_S5)
ggsave("PC_size_histogram_all.pdf", Plot_all)


## KI67

ki67_summary_grouped <- Sample2 %>%
  group_by(Group) %>%
  summarise(
    Total_cells = n(),
    Ki67_positive = sum(`Nucleus: Ki67 mean` > 100, na.rm = TRUE),
    Percent_positive = 100 * Ki67_positive / Total_cells
  )

ki67_summary_grouped <- Sample3 %>%
  group_by(Group) %>%
  summarise(
    Total_cells = n(),
    Ki67_positive = sum(`Nucleus: Ki67 mean` > 100, na.rm = TRUE),
    Percent_positive = 100 * Ki67_positive / Total_cells
  )

ki67_summary_grouped <- all_samples %>%
  group_by(Group) %>%
  summarise(
    Total_cells = n(),
    Ki67_positive = sum(`Nucleus: Ki67 mean` > 100, na.rm = TRUE),
    Percent_positive = 100 * Ki67_positive / Total_cells
  )

ki67_summary_grouped

colnames()

