# ============================================================
# Demo: Relative abundance estimation for Bathynomus giganteus
# ============================================================
# Purpose:
#   - Demonstrate CPUE-based relative abundance estimation
#   - Integrate catch counts with trap effort
#   - Fit a negative binomial GLM with effort offset to estimate environmental preference
#
# Inputs:
#   - ../data/BTW_expedition_catch2022.csv
#   - ../data/BTW_expedition_effort2022.csv
#
# Outputs:
#   - ../data/cpue_summary_bahamas.csv
#   - ../figures/relative_abundance_bahamas.png
# ============================================================

# ------------------------
# 0. Libraries
# ------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(MASS)
library(broom)

# ------------------------
# 1. Load data
# ------------------------
catch_counts <- read.csv(
  "../data/BTW_expedition_catch2022.csv",
  stringsAsFactors = FALSE
)

effort_df <- read.csv(
  "../data/BTW_expedition_effort2022.csv",
  stringsAsFactors = FALSE
)

# ------------------------
# 2. Parse dates and compute trap-hours
# ------------------------
effort_df <- effort_df %>%
  mutate(
    date = mdy(date),
    start_dt = as.POSIXct(
      paste(date, deploy_start),
      format = "%Y-%m-%d %H:%M",
      tz = "UTC"
    ),
    end_dt = as.POSIXct(
      paste(date, deploy_end),
      format = "%Y-%m-%d %H:%M",
      tz = "UTC"
    ),
    duration_hrs = as.numeric(difftime(end_dt, start_dt, units = "hours")),
    trap_hours   = traps * duration_hrs
  )

# ------------------------
# 3. Join catch and effort, compute CPUE
# ------------------------
cpue_df <- catch_counts %>%
  left_join(
    effort_df %>% dplyr::select(set, traps, trap_hours),
    by = "set"
  ) %>%
  mutate(
    cpue_per_trap_hour = n_caught / trap_hours,
    cpue_per_trap      = n_caught / traps
  )

# ------------------------
# 4. Assign region and covariates
# ------------------------
cpue2 <- cpue_df %>%
  mutate(
    region = case_when(
      set %in% c("BTW_DW1","BTW_DW2","BTW_DW3","BTW_DW4") ~ "Nassau",
      set %in% c("BTW_DW5","BTW_DW6","BTW_DW7","BTW_DW8") ~ "Exuma",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(
    effort_df %>% dplyr::select(set, duration_hrs, depth, bottom_temp),
    by = "set"
  )

# ------------------------
# 5. Summarise CPUE and covariates by region
# ------------------------
summary_tbl <- cpue2 %>%
  group_by(region) %>%
  summarise(
    cpue_traphour_mean = mean(cpue_per_trap_hour),
    cpue_traphour_sd   = sd(cpue_per_trap_hour),
    cpue_trap_mean     = mean(cpue_per_trap),
    cpue_trap_sd       = sd(cpue_per_trap),
    depth_mean         = mean(depth, na.rm = TRUE),
    depth_sd           = sd(depth, na.rm = TRUE),
    temp_mean          = mean(bottom_temp, na.rm = TRUE),
    temp_sd            = sd(bottom_temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -region,
    names_to  = c("metric", "stat"),
    names_sep = "_(?=(mean|sd)$)"
  ) %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    metric = recode(metric,
                    cpue_traphour = "Catch per Trap Hour",
                    cpue_trap     = "Catch per Trap",
                    depth         = "Depth (m)",
                    temp          = "Temperature (°C)"
    ),
    ci_lower = pmax(mean - sd, 0),
    ci_upper = mean + sd
  )

# Save summary table
write.csv(
  summary_tbl,
  "../data/cpue_summary_bahamas.csv",
  row.names = FALSE
)

# ------------------------
# 6. Fit negative binomial GLM
# ------------------------
cpue2 <- cpue2 %>%
  mutate(
    depth_s = as.numeric(scale(depth)),
    temp_s  = as.numeric(scale(bottom_temp))
  )

m_nb <- glm.nb(
  n_caught ~ depth_s + temp_s + offset(log(trap_hours)),
  data = cpue2
)
summary(m_nb)

glm_tidy <- broom::tidy(
  m_nb,
  conf.int = TRUE,
  exponentiate = TRUE
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term,
                  depth_s = "Depth",
                  temp_s  = "Temperature"
    ),
    p_txt = format.pval(p.value, digits = 2, eps = 1e-3),
    irr_txt = sprintf(
      "%s: IRR = %.2f (95%% CI %.2f–%.2f), p = %s",
      term, estimate, conf.low, conf.high, p_txt
    )
  )

glm_subtitle <- paste(glm_tidy$irr_txt, collapse = "\n")

# ------------------------
# 7. Plot relative abundance
# ------------------------
# dir.create("../figures", showWarnings = FALSE)

p <- ggplot(summary_tbl, aes(x = region, y = mean, fill = region)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2
  ) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("Nassau" = "orange", "Exuma" = "turquoise")
  ) +
  labs(
    title = expression(paste("Relative Abundance of ", italic("B. giganteus"))),
    subtitle = glm_subtitle,
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),
    axis.text   = element_text(color = "black"),
    strip.text  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", size = 16)
  )

ggsave(
  "../figures/relative_abundance_bahamas.jpg",
  p,
  width = 9,
  height = 6,
  dpi = 300
)

message("Demo complete: relative abundance estimated and plotted.")
