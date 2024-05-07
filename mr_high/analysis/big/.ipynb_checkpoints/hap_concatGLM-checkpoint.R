library("leaps")
library("ggplot2")
library('dplyr')
library("kableExtra")
library("knitr")

# cleaning data 
data <- big_scores
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

# concatGLM analysis
data$concat_seq_score <- 1 - data$concat_seq_score
concat_data <- data[data$concat_seq_score != 2, ]
concatglm <- glm(concat_seq_score~
                   no_of_taxa+
                   exponential_growth_rate+
                   individuals_per_taxa+
                   effective_population_size+
                   species_birth_rate+
                   tree_wide_sub_rate+
                   no_of_sites+
                   ado+
                   amplification_error+
                   sequencing_errors, 
                 data = concat_data)
summary(concatglm)




# making the concatGLM table
glm_summary <- summary(concatglm)

parameter_labels <- data.frame(
  Parameter = c("(Intercept)", "no_of_taxa", "exponential_growth_rate",
                "individuals_per_taxa", "effective_population_size",
                "species_birth_rate", "tree_wide_sub_rate", "no_of_sites",
                "ado", "amplification_error", "sequencing_errors"),
  Label = c("Intercept","Number of Tumours", "Exponential Growth Rate",
            "Cells per Tumour", "Effective Population Size",
            "Metastatic Rate", "Mutation Rate", "Number of Sites",
            "Allelic Dropout", "Amplification Error", "Sequencing Errors")
)

rownames(glm_summary$coefficients) <- parameter_labels$Label[
  match(rownames(glm_summary$coefficients), parameter_labels$Parameter)
]

glm_summary_with_stars <- data.frame(
  glm_summary$coefficients,
  stars = ifelse(glm_summary$coefficients[, "Pr(>|t|)"] < 0.001, "***",
                 ifelse(glm_summary$coefficients[, "Pr(>|t|)"] < 0.01, "**",
                        ifelse(glm_summary$coefficients[, "Pr(>|t|)"] < 0.05, "*",
                               "")))
)

colnames(glm_summary_with_stars) <- gsub("Std..Error", "Standard Error", colnames(glm_summary_with_stars))
colnames(glm_summary_with_stars) <- gsub("t.value", "t-value", colnames(glm_summary_with_stars))
colnames(glm_summary_with_stars) <- gsub("Pr...t..", "p-value", colnames(glm_summary_with_stars))
colnames(glm_summary_with_stars)[colnames(glm_summary_with_stars) == "stars"] <- ""

kable_table <- kable(glm_summary_with_stars, format = "html", align = "c", caption = "GLM on Phylogenetic Accuracy of Concatenated H-CS Trees (MR-High)") %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    font_size = 20, 
    full_width = FALSE, 
    position = "center", 
    fixed_thead = TRUE 
  ) %>%
  column_spec(1:ncol(glm_summary_with_stars), extra_css = "font-family: 'Noto Serif Sinhala', serif; color: black;")

temp_html <- tempfile(fileext = ".html")
save_kable(kable_table, file = temp_html)

save_path_png <- "12_high/results/big_scores/hap_concatGLM.png"
webshot(temp_html, file = save_path_png, zoom = 2)
unlink(temp_html)


# R^2 analysis

# Rename the parameters in the concat_data dataset
colnames(concat_data) <- c("Taxa", "CellsPerTumour", "EffPopSize", "MetastaticRate",
                           "MutationRate", "NumberOfSites", "ExpGrowthRate", "ADO", "AmpError", "SeqError",
                           "concat_seq_score","ASTRAL_score", "hap_FastRFS_score", "vcf_merge_seq_score","vcf_FastRFS_score")

# Perform the regression analysis with the renamed parameters
concat_r2 <- regsubsets(concat_seq_score ~ Taxa + CellsPerTumour +
                          EffPopSize + MetastaticRate + MutationRate + NumberOfSites +
                          ExpGrowthRate + ADO + AmpError + SeqError,
                        data = concat_data)

reg.summary=summary(concat_r2)
plot(concat_r2, scale = "r2", main = "Subset Selection - Phylogenetic Accuracy of Concatenated H-CS Trees (MR-High)")

