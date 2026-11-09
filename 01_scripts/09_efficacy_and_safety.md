# 9. Treatment efficacy and safety (Fig. 2)

Unlike every other module in this repository, this figure is built
directly from **clinical trial metadata** (egg counts, cure status, and
adverse-event records), not from a metagenomic read-based or MAG-based
pipeline — there is no "generation" step here, only analysis and
plotting. It is included as its own module because it doesn't fit under
any of the other numbered stages, and was previously missing from this
repository entirely.

Requires `efficacy_table.csv` (per-arm/country egg-count and cure-rate
summary, built from the parent trial's randomisation/outcome data) and the
`physeq_pre` phyloseq object's `sample_data` (built in module 2, which
carries the adverse-event columns `AE_any`/`AE_GID` alongside the
microbiome metadata).

## Step 0 — shared helper and colour palette

```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(scales)
library(patchwork)
library(ggsci)
library(gghalves)

arm_colours <- c("ALB" = "#E07A5F", "FDCx1" = "#3D5A80", "FDCx3" = "#293241")

add_binomial_ci <- function(dat, events, total) {
  dat %>%
    rowwise() %>%
    mutate(ci = list(binom.test(x = .data[[events]], n = .data[[total]])$conf.int),
           lower = ci[[1]], upper = ci[[2]]) %>%
    ungroup() %>%
    select(-ci)
}
```

## Panel A — cure rate by microscopy and qPCR

```r
efficacy_all <- read_csv('efficacy_table.csv')

# Microscopy-based cure rate: hookworm (Ethiopia) and T. trichiura (Mozambique)
efficacy_plot_data <- efficacy_all %>%
  filter((country == "Mozambique" & species == "TR") | (country == "Ethiopia" & species == "HKW")) %>%
  filter(treatment_arm %in% c("ALB", "FDCx1", "FDCx3")) %>%
  mutate(cured_n = n_pre - n_post) %>%
  add_binomial_ci(events = "cured_n", total = "n_pre") %>%
  mutate(
    cure_rate = cured_n / n_pre,
    treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3")),
    facet_label = case_when(
      country == "Ethiopia" ~ "Hookworm-Eth\nmicroscopy",
      country == "Mozambique" & species == "TR" ~ "T. trichiura-Moz\nmicroscopy"
    )
  )

# qPCR-based cure rate: T. trichiura, Mozambique only
qpcr_plot_data <- metadata_wide %>%
  filter(country == "Mozambique", TR_pre == 1) %>%
  filter(treatment_arm %in% c("ALB", "FDCx1", "FDCx3")) %>%
  mutate(cured = TR_Ct_Mean_post >= 50) %>%
  group_by(treatment_arm) %>%
  summarise(cured_n = sum(cured, na.rm = TRUE), n_pre = sum(!is.na(TR_Ct_Mean_post)), .groups = "drop") %>%
  add_binomial_ci(events = "cured_n", total = "n_pre") %>%
  mutate(cure_rate = cured_n / n_pre, treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3")),
         facet_label = "T. trichiura-Moz\nqPCR")

efficacy_plot_data_all <- bind_rows(
  efficacy_plot_data %>% select(treatment_arm, cured_n, n_pre, lower, upper, cure_rate, facet_label),
  qpcr_plot_data %>% select(treatment_arm, cured_n, n_pre, lower, upper, cure_rate, facet_label)
) %>%
  mutate(facet_label = factor(facet_label, levels = c("Hookworm-Eth\nmicroscopy", "T. trichiura-Moz\nmicroscopy", "T. trichiura-Moz\nqPCR")))

p_efficacy <- ggplot(efficacy_plot_data_all, aes(x = cure_rate, y = treatment_arm, colour = treatment_arm)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(cured_n, "/", n_pre)), nudge_y = 0.25, colour = "black", size = 3) +
  facet_wrap(~facet_label, nrow = 1) +
  scale_colour_manual(values = arm_colours) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = "Cure rate (95% CI)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(size = 9))
```

## Panel B — paired pre/post egg counts, ALB arm (raincloud plot)

```r
# tri_pre_post and hook_pre_post are built in module 2 (taxonomy), Part B,
# from the parent trial's egg-count data
raincloud_data <- bind_rows(
  hook_pre_post %>% filter(arm == "ALB") %>% mutate(group = "Hookworm - Ethiopia"),
  tri_pre_post %>% filter(arm == "ALB") %>% mutate(group = "T. trichiura - Mozambique")
) %>%
  mutate(
    log_epg = log(epg + 1),
    point = factor(point, levels = c("Pre-treatment", "Post-treatment")),
    x_num = as.numeric(point),
    x_jit = x_num + ifelse(point == "Pre-treatment", -0.12, 0.12)
  )

p_erri_raincloud <- ggplot(raincloud_data, aes(x = x_num, y = log_epg)) +
  geom_half_violin(data = filter(raincloud_data, point == "Pre-treatment"), aes(x = x_num),
                    side = "l", fill = "#E07A5F", colour = NA, alpha = 0.6, trim = FALSE) +
  geom_half_violin(data = filter(raincloud_data, point == "Post-treatment"), aes(x = x_num),
                    side = "r", fill = "#3D5A80", colour = NA, alpha = 0.6, trim = FALSE) +
  geom_line(aes(x = x_jit, group = id), colour = "grey60", alpha = 0.3, linewidth = 0.3) +
  geom_point(aes(x = x_jit, colour = point), size = 1.3, alpha = 0.7) +
  scale_colour_manual(values = c("Pre-treatment" = "#E07A5F", "Post-treatment" = "#3D5A80")) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Pre-treatment", "Post-treatment"), limits = c(0.4, 2.6)) +
  facet_wrap(~group, nrow = 1) +
  labs(x = NULL, y = "Log(EPG + 1)") +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_rect(fill = "grey92", colour = NA))
```

## Panel C — adverse-event forest plot (any AE and GID, by country and arm)

```r
safety_plot_data <- data.frame(sample_data(physeq_pre)) %>%
  rownames_to_column("Sample_ID") %>%
  distinct(Participant_ID, country, treatment_arm, AE_any, AE_GID) %>%
  pivot_longer(cols = c(AE_any, AE_GID), names_to = "outcome", values_to = "status") %>%
  mutate(
    event = status == "Yes",
    outcome = recode(outcome, AE_any = "Any adverse event", AE_GID = "Gastrointestinal adverse event")
  ) %>%
  group_by(country, treatment_arm, outcome) %>%
  summarise(events = sum(event, na.rm = TRUE), total = sum(!is.na(event)), .groups = "drop") %>%
  add_binomial_ci(events = "events", total = "total") %>%
  mutate(proportion = events / total, treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3")))

p_safety <- ggplot(safety_plot_data, aes(x = proportion, y = treatment_arm, colour = treatment_arm)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.15, linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0(events, "/", total)), nudge_y = 0.25, colour = "black", size = 3) +
  facet_grid(outcome ~ country) +
  scale_colour_manual(values = arm_colours) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(-0.05, 0.605), breaks = seq(0, 1, 0.2)) +
  labs(x = "Participants with an event (95% CI)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none")
```

## Panel D — specific adverse-event type, by country and treatment arm

```r
meta_ae <- as.data.frame(sample_data(physeq_pre)) %>%
  rownames_to_column("Sample_ID") %>%
  distinct(Participant_ID, country, treatment_arm, across(starts_with("AE_")))

# Specific AE term columns only (drop the aggregate any/GID/NSD/Other summary columns)
ae_term_cols <- meta_ae %>%
  select(starts_with("AE_")) %>%
  select(-any_of(c("AE_any", "AE_GID", "AE_NSD", "AE_OTHER"))) %>%
  colnames()

ae_type_data <- meta_ae %>%
  select(Participant_ID, country, treatment_arm, all_of(ae_term_cols)) %>%
  pivot_longer(cols = all_of(ae_term_cols), names_to = "AE", values_to = "status") %>%
  mutate(AE = str_remove(AE, "^AE_"), AE = str_replace_all(AE, "_", " "), event = status == "Yes" | status == 1) %>%
  filter(event) %>%
  count(country, treatment_arm, AE, name = "n") %>%
  mutate(treatment_arm = factor(treatment_arm, levels = c("ALB", "FDCx1", "FDCx3")))

ae_order <- ae_type_data %>%
  group_by(AE) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(AE)

ae_type_data <- ae_type_data %>% mutate(AE = factor(AE, levels = ae_order))

p_ae_type <- ggplot(ae_type_data, aes(x = AE, y = n, fill = treatment_arm)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~country, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = arm_colours, name = NULL) +
  labs(x = NULL, y = "Participants with event") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
```

## Assemble Figure 2

```r
figure_2 <- (p_efficacy / p_erri_raincloud) | (p_safety / p_ae_type)
figure_2

ggsave("Figures/Figure2.pdf", figure_2, width = 14, height = 10)
ggsave("Figures/Figure2.png", figure_2, width = 14, height = 10, dpi = 300)
```
