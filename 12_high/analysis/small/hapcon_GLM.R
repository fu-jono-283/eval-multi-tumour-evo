library("leaps")
library("webshot")
library("ggplot2")
library('dplyr')
library("kableExtra")
library("knitr")

# cleaning data 
data <- small_scores
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

# hapconGLM analysis
data$hapcon_seq_score <- 1 - data$hapcon_seq_score
hapcon_data <- data[data$hapcon_seq_score != 2, ]
hapconglm <- glm(hapcon_seq_score~
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
                 data = hapcon_data)
summary(hapconglm)


# making the hapconGLM table
glm_summary <- summary(hapconglm)

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

kable_table <- kable(glm_summary_with_stars, format = "html", align = "c", caption = "GLM on Phylogenetic Accuracy of Reconstructed H-CS Trees (MR-High)", escape = FALSE) %>%
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
save_path_png <- "12_high/results/small_scores/hapconGLM.png"
webshot(temp_html, file = save_path_png, zoom = 2)
unlink(temp_html)

# R^2 analysis

# Rename the parameters in the hapcon_data dataset
colnames(hapcon_data) <- c("Taxa", "CellsPerTumour", "EffectivePopSize", "MetastaticRate",
                           "MutationRate", "NumberOfSites", "ExpGrowthRate", "ADO", "AmpError", "SeqError",
                           "Replicate", "scSEQ_score", "hapcon_tree_scores", "vcf_tree_scores")

hapcon_r2 <- regsubsets(hapcon_tree_scores ~ Taxa + CellsPerTumour +
                          EffectivePopSize + MetastaticRate + MutationRate + NumberOfSites +
                          ExpGrowthRate + ADO + AmpError + SeqError,
                        data = hapcon_data)

reg.summary=summary(hapcon_r2)
plot(hapcon_r2, scale = "r2", main = "Subset Selection - Phylogenetic Accuracy of Reconstructed H-CS Trees (MR-High)")




