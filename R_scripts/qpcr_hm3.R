library(tidyverse)
library(pheatmap)
library(multcompView) # Multiple comparison
library(ggthemes)
library(forcats) # Rename levels
data <- read.csv("./tidy/qpcr_hmc3.csv")

data_summary <- data %>%
        group_by(treatment) %>%
        summarize(IL1B = log2(mean(IL1B, na.rm = TRUE)),
                  IL6 = log2(mean(IL6, na.rm = TRUE)),
                  TNFA = log2(mean(TNFA, na.rm = TRUE)),
                  APOE = log2(mean(APOE, na.rm = TRUE)),
                  CD33 = log2(mean(CD33, na.rm = TRUE)),
                  ST6GALNAC3 = log2(mean(ST6GALNAC3, na.rm = TRUE)),
                  TREM2 = log2(mean(TREM2, na.rm = TRUE))) %>%
        filter(!(treatment == "control"))


################### HEATMAP ###############################################
data_summary_df <- as.data.frame(data_summary)

rownames(data_summary_df) <- data_summary_df$treatment

rownames(data_summary_df) <- c("AβO", "AβO/Chol", "AβO/Chol/LPS", "AβO/Fruc",
                               "AβO/Fruc/Chol/LPS","AβO/LPS","Chol","Chol/LPS",
                               "Fruc","Fruc/Chol","Fruc/LPS","LPS")

data_summary_df$treatment <- NULL

pheatmap(data_summary_df, main = "Gene expression (log2FC): 2^–∆∆Ct",  display_numbers = TRUE,
         border_color = "black",
         number_color = "black")


########################### Stats ##########################
data$treatment <- as.factor(data$treatment)

data$treatment <- factor(data$treatment, levels = c("control", "abo", "abo_cho","abo_cho_lps","abo_fru","abo_fru_cho_lps",
                            "abo_lps","cho","cho_lps","fru","fru_cho","fru_lps","lps"))


data$treatment <- do.call(
        fct_recode,
        c(list(data$treatment), setNames(c("control", "abo", "abo_cho","abo_cho_lps","abo_fru","abo_fru_cho_lps",
                                           "abo_lps","cho","cho_lps","fru","fru_cho","fru_lps","lps"), c("Control","AβO", "AβO/Chol", "AβO/Chol/LPS", "AβO/Fruc",
                                                                                                         "AβO/Fruc/Chol/LPS","AβO/LPS","Chol","Chol/LPS",
                                                                                                         "Fruc","Fruc/Chol","Fruc/LPS","LPS")))
)

### IL1B
anova <- aov(IL1B ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(IL1B), sd = sd(IL1B)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean+sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "IL1B") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,12) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))


### IL6
anova <- aov(IL6 ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(IL6, na.rm = TRUE), sd = sd(IL6, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "IL6") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,55) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))


### TNFA
anova <- aov(TNFA ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(TNFA, na.rm = TRUE), sd = sd(TNFA, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "TNFA") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,20) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))

### APOE
anova <- aov(APOE ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(APOE, na.rm = TRUE), sd = sd(APOE, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "APOE") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,15) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))

### CD33
anova <- aov(CD33 ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(CD33, na.rm = TRUE), sd = sd(CD33, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "CD33") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))

### ST6GALNAC3
anova <- aov(ST6GALNAC3 ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(ST6GALNAC3, na.rm = TRUE), sd = sd(ST6GALNAC3, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "ST6GALNAC3") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))


### TREM2
anova <- aov(TREM2 ~ treatment, data = data)

# Tukey's test
tukey <- TukeyHSD(anova)

# compact letter display
cld <- multcompLetters4(anova, tukey)

dt <- group_by(data, treatment) %>%
        summarise(mean=mean(TREM2, na.rm = TRUE), sd = sd(TREM2, na.rm = TRUE)) %>%
        arrange(desc(mean))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$treatment)
dt$cld <- cld$Letters

ggplot(dt, aes(treatment, mean)) + 
        geom_bar(stat = "identity", aes(fill = mean), show.legend = FALSE) +
        geom_errorbar(aes(ymin = mean-sd, ymax=mean + sd), width = 0.2) +
        labs(x = "Treatment", y = "Fold Change", title = "TREM2") +
        geom_text(aes(label = cld, y = mean + sd), vjust = -0.5) +
        ylim(0,5) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, color = "black"))
