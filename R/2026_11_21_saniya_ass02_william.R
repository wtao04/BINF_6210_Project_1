############################################################
# BINF 6210 Assignment 2 - Revised Script
#
# Original Author: William Tao
#
# Revisions and Additional Analysis: Saniya Chadha
#
# Contributions:
# - William Tao: Developed the original script, performed
#   the initial TE analyses, statistical testing, and
#   visualizations.
#
# - Saniya Chadha: Reviewed and revised the script by
#   improving figure labels and readability, adding
#   automated figure export with ggsave(), creating a
#   reusable find_outliers() function, adding code comments
#   and formatting improvements, and incorporating a
#   chi-square test to compare TE class composition between
#   species.
############################################################

# 01 Load libraries + datasets ----

library(tidyverse)
library(viridis)
library(ggtext)
library(ggplot2)

t_spiralis <- read_csv("../data/GCA_001447595.1.csv")
t_spiralis
t_britovi <- read_csv("../data/GCA_001447585.1.csv")
t_britovi

#excluding non-TE matching classes, and adding row numbers to create unique identifier for each row, enabling exclusion of outliers

ls_non_te <- c("Simple_repeat", "Low_complexity", "Satellite", "snRNA", "rRNA", "ARTEFACT", "Unknown", "simple_repeat")

t_spiralis <- t_spiralis |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_spiralis) #1037

t_britovi <- t_britovi |> mutate(row = 1:n()) |> filter(matching_class %notin% ls_non_te)
dim(t_britovi) #1177


# 02 Identify, Plot, and Remove Outliers (using interquartile range method) ----

#EDIT: code is harder to read since its mostly one long line of statements with no comments.

#Suggested comment: "Calculate the first and third quartile of percent divergence"
Q1_spiralis <- quantile(t_spiralis$perc_div, .25)
Q3_spiralis <- quantile(t_spiralis$perc_div, .75)
#Suggested comment: "Calculate the interquartile range"
IQR_spiralis <- IQR(t_spiralis$perc_div)
#Suggested comment: "Find the observations out side of the 1.5 range of the IQR threshold 
outliers_spiralis <- subset(t_spiralis, t_spiralis$perc_div<(Q1_spiralis - 1.5*IQR_spiralis) | t_spiralis$perc_div>(Q3_spiralis + 1.5*IQR_spiralis))
#EDIT: can write nrow(outliers_spiralis)
nrow(outliers_spiralis)
dim(outliers_spiralis) #7

#EDIT: Since the code is repeated twice for both species maybe we can create a function to be used for both species
find_outliers <- function(df) {
Q1 <- quantile(df$perc_div, 0.25)
Q3 <- quantile(df$perc_div, 0.75)
IQR_val <- IQR(df$perc_div)

subset(
  df,
  perc_div < (Q1 - 1.5 * IQR_val) |
    perc_div > (Q3 + 1.5 * IQR_val)
)
}
# Calculate statistics for T. spiralis
Q1_spiralis <- quantile(t_spiralis$perc_div, 0.25)
Q3_spiralis <- quantile(t_spiralis$perc_div, 0.75)
IQR_spiralis <- IQR(t_spiralis$perc_div)
outliers_spiralis <- find_outliers(t_spiralis)

# Calculate statistics for T. britovi
Q1_britovi <- quantile(t_britovi$perc_div, 0.25)
Q3_britovi <- quantile(t_britovi$perc_div, 0.75)
IQR_britovi <- IQR(t_britovi$perc_div)
outliers_britovi <- find_outliers(t_britovi)
# Remove outliers
t_spiralis_out_removed <- t_spiralis |>
  filter(row %notin% outliers_spiralis$row)

t_britovi_out_removed <- t_britovi |>
  filter(row %notin% outliers_britovi$row)



#EDIT: Consider reformatting long ggplot statements across multiple lines and grouping related arguments together additional spacing and indentation would improve readability and make future modifications or debugging easier.
t_spiralis |> 
  ggplot(aes(matching_class,y = perc_div)) + 
  geom_boxplot() + 
  
  geom_hline(aes(
    yintercept = Q1_spiralis - 1.5*IQR_spiralis, 
    color = "Q1")) + 
  
  geom_hline(aes(
    yintercept = Q3_spiralis + 1.5*IQR_spiralis,  
    colour = "Q3")) + 
  
  scale_color_manual(name = "Outlier Thresholds", 
                     values = c("red", "blue")) + 
  guides(x = guide_axis(angle = 45)) + 
  
  labs(title = "Outlier Detection by TE Matching Class - *T. spiralis*", 
       x = "TE Matching Class", 
       y = "Percent Divergence") + 
  theme(plot.title = element_markdown())


#Suggested comment: "Calculate the first and third quartile of percent divergence"
Q1_britovi <- quantile(t_britovi$perc_div, .25)
Q3_britovi <- quantile(t_britovi$perc_div, .75)
#Suggested comment: "Calculate the interquartile range"
IQR_britovi <- IQR(t_britovi$perc_div)
#Suggested comment: "Find the observations out side of the 1.5 range of the IQR threshold 
outliers_britovi <- subset(t_britovi, t_britovi$perc_div<(Q1_britovi - 1.5*IQR_britovi) | t_britovi$perc_div>(Q3_britovi + 1.5*IQR_britovi)) 
#EDIT: can write nrow(outliers_spiralis)
dim(outliers_britovi) #1

t_britovi |> 
  ggplot(aes(matching_class,y = perc_div)) + 
  geom_boxplot() + 
  
  geom_hline(aes(
    yintercept = Q1_britovi - 1.5*IQR_britovi, 
    color = "Q1")) + 
  
  geom_hline(aes(
    yintercept = Q3_britovi + 1.5*IQR_britovi,  
    colour = "Q3" )) + 
  
  scale_color_manual(name = "Outlier Threshold",
                     values = c("red", "blue"))+ 
  guides(x = guide_axis(angle = 45)) + 
  
  labs(title = "Outlier Detection by TE Matching Class - *T. britovi*", 
       x = "TE Matching Class", 
       y = "Percent Divergence") + 
  theme(plot.title = element_markdown())

#Remove outliers from dataset
t_spiralis_out_removed <- t_spiralis |> filter(row %notin% outliers_spiralis$row)

dim(t_spiralis_out_removed) #1030

t_britovi_out_removed <- t_britovi |> filter(row %notin% outliers_britovi$row)
dim(t_britovi_out_removed) #1176

#EDIT: Summary of code after identifiying outliers, table created
summary_stats <- data.frame(
  Species = c("T. spiralis", "T. britovi"),
  Outliers = c(nrow(outliers_spiralis),
               nrow(outliers_britovi))
) |> view()

# 03 Perform Wilcoxon tests: percent divergence between species, within each TA mathcing_class ----
#EDIT: another additional test may be a chi-square test that may complement the Wilcoxon test by comparing the frequency distribution of TE classes between species 


te_combined <- bind_rows(
  t_spiralis |> mutate(species = "spiralis"),
  t_britovi  |> mutate(species = "britovi")
) |>
  filter(
    !matching_class %in% ls_non_te,
    !is.na(matching_class),
    is.finite(perc_div)
  )

wilcox_results <- te_combined |>
  group_by(matching_class) |>
  group_modify(~ {
    x <- .x$perc_div[.x$species == "spiralis"]
    y <- .x$perc_div[.x$species == "britovi"]
    
    # Skip classes with fewer than two observations in either species
    # or where every value is identical
    can_test <- length(x) >= 2 &&
      length(y) >= 2 &&
      n_distinct(c(x, y)) > 1
    
    tibble(
      n_spiralis = length(x),
      n_britovi = length(y),
      median_spiralis = if (length(x) > 0) median(x) else NA_real_,
      median_britovi = if (length(y) > 0) median(y) else NA_real_,
      p_value = if (can_test) {
        wilcox.test(
          x, y,
          alternative = "two.sided",
          paired = FALSE,
          exact = FALSE
        )$p.value
      } else NA_real_
    )
  }) |>
  ungroup() |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) |>
  arrange(p_adjusted)

print(wilcox_results, n = Inf)

#EDIT: a chi-square test to compare TE class composition between species. 
class_counts = te_combined |> 
  count(matching_class, species) |> 
  pivot_wider(
    names_from = species, 
    values_from = n, 
    values_fill = 0
  )
count_matrix = class_counts |> 
  select(spiralis, britovi) |> 
  as.matrix()
chi_result = chisq.test(count_matrix)
chi_result

# 04 Plot average percent divergence within matching_class between t_spiralis & t_britovi using boxplots ----

# find distribution of %div of classes in britovi
shared_classes <- intersect(t_spiralis$matching_class, t_britovi$matching_class)

te_shared <- bind_rows(
  t_britovi  |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. britovi"),
  t_spiralis |> 
    select(matching_class, perc_div) |> 
    mutate(species = "T. spiralis")) |>
  filter(matching_class %in% shared_classes)


# EDIT: adding a figure axis label that is well detailed and title would help make the plots more readable
# EDIT: The combined boxplot directly compares both species within a single visualization, the individual species plots may be redundant and could be removed
te_shared |> filter(species == "T. britovi") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot() +
  labs(
    title = "% Divergence of Shared TE Classes in T. britovi", 
    x = "% Divergence",
    y = "TE Class"
  )
#EDIT: adding ggsave to save it locally (Karl likes)
ggsave("T_britovi_boxplot.png", width = 8, height = 6)

te_shared |> filter(species == "T. spiralis") |>
  ggplot(aes(x = perc_div, y = matching_class)) + geom_boxplot() +
  labs(
    title = "% Divergence of Shared TE Classes in T. spiralis", 
    x = "% Divergence",
    y = "TE Class"
  )
#EDIT: adding ggsave to save it locally (Karl likes)
ggsave("T_spiralis_boxplot.png", width = 8, height = 6)

# combine datasets to make combined graph
ggplot(data = te_shared, mapping = aes(x = perc_div, y = matching_class, fill = species)) +
  geom_boxplot() +
  labs(
    title = "% Divergence in Shared TE Classes of T. britovi and T. spiralis",
    x = "% Divergence",
    y = "TE Class")
#EDIT: adding ggsave to save it locally (Karl likes)
ggsave("T_spiralis_britovi_boxplot.png", width = 8, height = 6)
