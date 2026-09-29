library(dplyr)
library(ggplot2)

threshold <- 100

setwd("Z:/home/users/grenova/K-tag/IHC/Spleen ASC architecture/Quantification/Cluster composition")

Sample2_kika <- read.delim("EHP309_Sample2_kika2.txt", header = TRUE, sep = "\t")
Sample3_kika <- read.delim("EHP596_Sample3_kika.txt", header = TRUE, sep = "\t")
Sample5_kika <- read.delim("EHP596_Sample5_kika.txt", header = TRUE, sep = "\t")

colnames(Sample2_kika)
Sample2_kika <- Sample2_kika %>%
  rename(
    Nucleus..Ig.mean = Nucleus..Ig.s.mean
  )

all_samples <- bind_rows(
  Sample2_kika %>% mutate(Sample = "Sample2"),
  Sample3_kika %>% mutate(Sample = "Sample3"),
  Sample5_kika %>% mutate(Sample = "Sample5")
)

colnames(all_samples)

summary(all_samples$Parent)
summary(all_samples$Sample)


composition <- all_samples %>%
  mutate(
    B220_pos = Cell..B220.mean > threshold,
    CD138_pos = Cell..CD138.mean > threshold,
    Ki67_pos = Nucleus..Ki67.mean > threshold,
    Ig_pos = Nucleus..Ig.mean > threshold
  ) %>%
  
  # Remove cells negative for both B220 and CD138
  filter(B220_pos | CD138_pos)


composition <- composition %>%
  mutate(
    Population = paste0(
      ifelse(B220_pos, "B220+", "B220-"), " ",
      ifelse(CD138_pos, "CD138+", "CD138-"), " ",
      ifelse(Ig_pos, "Ig+", "Ig-"), " ",
      ifelse(Ki67_pos, "Ki67+", "Ki67-")
    )
  )

# Count populations per annotation
composition <- composition %>%
count(Sample, Parent, Population)


# Convert to % within each annotation
composition <- composition %>%
group_by(Sample, Parent) %>%
  mutate(
    Percent = 100 * n / sum(n)
  ) %>%
  ungroup()



#Create a pie-chart function

#This will calculate the median composition across annotations within each sample:
  
make_pie <- function(comp, sample_name) {
  
  pie_data <- comp %>%
    filter(Sample == sample_name) %>%
    group_by(Population) %>%
    summarise(
      Median = median(Percent),
      .groups = "drop"
    ) %>%
    arrange(desc(Median)) %>%
    mutate(
      Population = factor(Population, levels = Population)
    )
  
  ggplot(pie_data, aes(x = "", y = Median, fill = Population)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = my_colors) +
    theme_void() +
    labs(title = sample_name)
}


pie2 <- make_pie(composition, "Sample2")
pie3 <- make_pie(composition, "Sample3")
pie5 <- make_pie(composition, "Sample5")

ggsave("Sample2_piechart.pdf", pie2)
ggsave("Sample3_piechart.pdf", pie3)
ggsave("Sample5_piechart.pdf", pie5)

# Calculate median composition across all annotations
overall_composition <- composition %>%
  group_by(Population) %>%
  summarise(
    Median = median(Percent),
    SD = sd(Percent),
    .groups = "drop"
  ) %>%
  arrange(desc(Median)) %>%
  mutate(
    Population = factor(Population, levels = Population)
  )


# Generate pie chart
common_pie <- ggplot(overall_composition,
       aes(x = "", y = Median, fill = Population)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = my_colors) +
  theme_void() +
  labs(
    title = "Median cellular composition across all samples"
  )

unique(composition$Population)
my_colors <- c(
  "B220+ CD138+ Ig- Ki67+" = "#E7B4BC",
  "B220+ CD138- Ig+ Ki67+" = "#E694A6",
  "B220+ CD138- Ig+ Ki67-" = "#7791EA",
  "B220+ CD138- Ig- Ki67+" = "#979FC6",
  "B220+ CD138- Ig- Ki67-" = "#D8E9E6",
  "B220- CD138+ Ig+ Ki67+" = "#CACDDC",
  "B220- CD138+ Ig+ Ki67-" = "#F2A843", 
  "B220- CD138+ Ig- Ki67-" = "#EAC4AF", 
  "B220+ CD138+ Ig+ Ki67-" = "#5BD285", 
  "B220- CD138+ Ig- Ki67+" = "#979797",
  "B220+ CD138+ Ig+ Ki67+" = "#9A3E67", 
  "B220+ CD138+ Ig- Ki67-" = "#AFE4DE"
)

My12Colours <- c('IgM PCs'="#E7B4BC",'ABCs'="#D8E9E6", 'Spleen IgM PCs'="#CACDDC",
                 'IgA PCs'="#F2A843", 'IgM ABC/PCs'="#979FC6", '11'="#EAC4AF", '10'="#5BD285",
                 'CS PCs'="#E694A6", '9'="#979797", 'IgD ABCs'="#7791EA", 'ABC/PB'= "#9A3E67",
                 'IgD ABC/PCs'="#AFE4DE" )

ggsave("piechart_all.pdf", common_pie)
