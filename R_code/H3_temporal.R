# 05_H3_temporal.R
#
# H3: Species will show temporal activity shifts at roadside locations,
# with larger-bodied and herbivore species showing the greatest displacement.
#
# Approach: Overlap coefficient (Dhat1) with bootstrap CI from the overlap
# package. Body mass modelled as continuous log10(kg).



# Time in radians ---------------------------------------------------------

if (!"time_rad" %in% colnames(animals_site)) {
  animals_site <- animals_site |>
    mutate(time_rad = (hour(eventStart) * 3600 +
                         minute(eventStart) * 60 +
                         second(eventStart)) / 86400 * 2 * pi)
}


# Per-species overlap coefficients ----------------------------------------
# Minimum 5 detections per treatment to attempt estimation;
# 5–9 flagged as low-confidence in the forest plot.

testable_species <- animals_site |>
  count(scientificName, treatment) |>
  pivot_wider(names_from = treatment, values_from = n, values_fill = 0L) |>
  filter(Roadside >= 20, Interior >= 20) |>
  pull(scientificName)

cat(sprintf("Species with ≥20 detections per treatment: %d\n", length(testable_species)))

overlap_results <- map(testable_species, compute_overlap) |>
  list_rbind() |>
  mutate(
    n_road = map_int(scientificName, \(s)
                     sum(animals_site$scientificName == s & animals_site$treatment == "Roadside")),
    n_int = map_int(scientificName, \(s)
                    sum(animals_site$scientificName == s & animals_site$treatment == "Interior"))
  ) |>
  left_join(guild_lookup |> select(scientificName, diet, log10_mass_kg, activity),
            by = "scientificName")


# Overlap categories (Ridout & Linkie 2009 thresholds)
categorize_overlap <- function(delta) {
  case_when(
    is.na(delta)  ~ NA_character_,
    delta > 0.85  ~ "High (>0.85)",
    delta >= 0.60 ~ "Moderate (0.60–0.85)",
    TRUE          ~ "Low (<0.60)"
  )
}

overlap_results <- overlap_results |>
  mutate(
    sample_adequate  = n_road >= 10 & n_int >= 10,
    overlap_category = factor(
      if_else(!is.na(delta), categorize_overlap(delta), NA_character_),
      levels = c("High (>0.85)", "Moderate (0.60–0.85)", "Low (<0.60)")
    )
  )

overlap_results |>
  filter(!is.na(delta)) |>
  count(overlap_category, .drop = FALSE) |>
  mutate(pct = round(100 * n / sum(n), 1)) |>
  print()


# Guild and body mass summaries -------------------------------------------

guild_stats <- overlap_results |>
  filter(!is.na(delta), !is.na(diet)) |>
  group_by(diet) |>
  summarise(
    n_species    = n(),
    mean_delta   = mean(delta),
    median_delta = median(delta),
    sd_delta     = sd(delta),
    .groups = "drop"
  ) |>
  arrange(mean_delta)

print(guild_stats)

cor_mass <- cor.test(overlap_results$log10_mass_kg, overlap_results$delta,
                     method = "spearman", exact = FALSE)

cat(sprintf("Body mass vs overlap: ρ = %.3f, p = %.4f\n",
            cor_mass$estimate, cor_mass$p.value))


# Plots -------------------------------------------------------------------

OVERLAP_COLOURS <- c(
  "Low (<0.60)"          = "#D7263D",
  "Moderate (0.60–0.85)" = "#F4B400",
  "High (>0.85)"         = "#00ab41"
)

# Forest plot: species-level overlap with 95% CI
# Shape encodes estimator: Dhat4 for species with >= 50 detections per
# treatment (filled circle), Dhat1 otherwise (open circle).
p_sp_overlap_forest <- overlap_results |>
  filter(!is.na(delta)) |>
  mutate(
    scientificName = fct_reorder(scientificName, delta),
    alpha_val = if_else(sample_adequate, 1.0, 0.55)
  ) |>
  ggplot(aes(delta, scientificName, colour = overlap_category)) +
  geom_vline(xintercept = c(0.60, 0.85), linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high, alpha = I(alpha_val)),
                 height = 0.3, linewidth = 0.8) +
  geom_point(aes(size = log10_mass_kg, shape = estimator, alpha = I(alpha_val)),
             stroke = 1) +
  scale_colour_manual(values = OVERLAP_COLOURS, name = "Overlap category") +
  scale_size_body(range = c(2, 8)) +
  scale_shape_manual(
    values = c("Dhat4" = 19L, "Dhat1" = 21L),
    labels = c("Dhat4" = "Dhat4 (n ≥ 50)", "Dhat1" = "Dhat1 (n < 50)"),
    name = "Estimator"
  ) +
  scale_x_continuous(limits = c(0, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.60, 0.75, 0.85, 1.0)) +
  labs(
    title    = "H3: Species-Level Temporal Overlap",
    subtitle = "Left = substantial temporal shift | Right = similar timing\nPoint size = body mass (kg); shape = estimator used",
    x = "Overlap coefficient", y = NULL
  ) +
  theme_bwindi(legend_pos = "right")

print(p_sp_overlap_forest)
bwindi_save(p_sp_overlap_forest, "H3_species_overlap_forest", width = "double")


# Overlap by body mass tertile — violin + labelled points
# A scatter of delta ~ mass works poorly with n < 15; grouping into tertiles
# shows the distribution shape while keeping individual species identifiable.
if (!exists("make_mass_tertile")) {
  make_mass_tertile <- function(x) {
    brks <- unique(quantile(x, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
    if (length(brks) >= 4) {
      cut(x, breaks = brks,
          labels = c("Small-bodied", "Mid-bodied", "Large-bodied"),
          include.lowest = TRUE)
    } else {
      factor(dplyr::ntile(x, 3), labels = c("Small-bodied", "Mid-bodied", "Large-bodied"))
    }
  }
}

overlap_mass_plot <- overlap_results |>
  filter(!is.na(delta), !is.na(log10_mass_kg)) |>
  mutate(
    mass_tertile = make_mass_tertile(log10_mass_kg),
    mass_tertile = factor(mass_tertile, levels = c("Small-bodied", "Mid-bodied", "Large-bodied"))
  )

p_mass_overlap <- overlap_mass_plot |>
  ggplot(aes(mass_tertile, delta)) +
  geom_hline(yintercept = c(0.60, 0.85), linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_violin(aes(fill = mass_tertile), alpha = 0.20, colour = NA,
              width = 0.8, trim = TRUE, scale = "width") +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high, colour = diet),
    width = 0.08, linewidth = 0.7, alpha = 0.7,
    position = position_jitter(width = 0.08, seed = 42)
  ) +
  geom_point(
    aes(size = log10_mass_kg, colour = diet, shape = sample_adequate),
    alpha = 0.9, stroke = 0.8,
    position = position_jitter(width = 0.08, seed = 42)
  ) +
  geom_text_repel(
    aes(label = scientificName, colour = diet),
    size = 2.6, fontface = "italic", max.overlaps = 15,
    seed = 42, min.segment.length = 0.2
  ) +
  scale_fill_manual(
    values = c("Small-bodied" = "#008080", "Mid-bodied" = "#D7263D",
               "Large-bodied" = "#F4B400"),
    guide = "none"
  ) +
  scale_colour_diet() +
  scale_size_body(range = c(2, 7)) +
  scale_shape_manual(
    values = c("TRUE" = 19L, "FALSE" = 21L),
    labels = c("TRUE" = "≥10 detections", "FALSE" = "5–9 detections"),
    name = "Sample size"
  ) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = c(0, 0.25, 0.5, 0.60, 0.75, 0.85, 1.0)) +
  labs(
    title    = "H3: Temporal Overlap by Body Mass Class",
    subtitle = sprintf(
      "Dashed lines at Δ̂₁ = 0.60 and 0.85 mark category boundaries\nSpearman ρ = %.3f, p = %.4f",
      cor_mass$estimate, cor_mass$p.value
    ),
    x = "Body mass class", y = "Overlap coefficient (Δ̂₁)", colour = "Diet guild"
  ) +
  theme_bwindi(legend_pos = "right")

print(p_mass_overlap)
bwindi_save(p_mass_overlap, "H3_bodymass_overlap", width = "double")

# Mean overlap by guild
p_guild_mean <- guild_stats |>
  mutate(diet = fct_reorder(diet, mean_delta)) |>
  ggplot(aes(mean_delta, diet, fill = diet)) +
  geom_vline(xintercept = c(0.60, 0.85), linetype = "dashed",
             colour = "grey50", linewidth = 0.7) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_text(aes(label = sprintf("n=%d  Δ̂₁=%.2f", n_species, mean_delta)),
            hjust = -0.05, size = 3.5, fontface = "bold") +
  scale_fill_diet(guide = "none") +
  scale_x_continuous(limits = c(0, 1.1),
                     breaks = c(0, 0.25, 0.5, 0.60, 0.75, 0.85, 1.0)) +
  labs(
    title    = "H3: Mean Temporal Overlap by Diet Guild",
    subtitle = "H3 prediction: Herbivores < Carnivores",
    x = "Mean overlap coefficient (Δ̂₁)", y = "Diet guild"
  ) +
  theme_bwindi()

print(p_guild_mean)
bwindi_save(p_guild_mean, "H3_guild_mean_overlap", width = "single")


# KDE overlap curves by diet guild × treatment ----------------------------
# Evaluates the same von Mises kernel used by overlapEst() on a fine grid,
# so the shaded intersection area is a direct visual of the Dhat1 statistic.
# Each panel shows the two treatment curves and their overlap region.

# Grid of 512 equally-spaced points on [0, 2π], matching overlap::densityFit
TIME_GRID <- seq(0, 2 * pi, length.out = 512)

# Convert radians back to clock hours for axis labels
rad_to_hour <- function(x) x / (2 * pi) * 24

# Build KDE curves for one group (identified by a logical index into animals_site)
kde_curves <- function(data, facet_var) {
  data |>
    filter(!is.na({{ facet_var }}), !is.na(time_rad)) |>
    group_by({{ facet_var }}, treatment) |>
    group_map(\(d, g) {
      bw   <- getBandWidth(d$time_rad)
      dens <- densityFit(d$time_rad, TIME_GRID, bw)
      # g holds the group keys: g[[1]] = facet level, g[[2]] = treatment
      tibble(
        facet     = as.character(g[[1]]),
        treatment = as.character(g[[2]]),
        t_rad     = TIME_GRID,
        hour      = rad_to_hour(TIME_GRID),
        density   = dens / (2 * pi)
      )
    }) |>
    list_rbind()
}

# Returns the pointwise minimum between Roadside and Interior KDE curves.
# Kept separate from the long-format KDE data so geom_ribbon has a clean
# x + ymin + ymax without a treatment column causing ambiguity.
overlap_ribbon <- function(kde_df) {
  kde_df |>
    pivot_wider(names_from = treatment, values_from = density) |>
    mutate(overlap_density = pmin(Roadside, Interior, na.rm = TRUE)) |>
    select(-Roadside, -Interior)
}

# Guild-level KDE curves
guild_kde <- kde_curves(
  animals_site |>
    filter(!is.na(diet)) |>
    mutate(diet = factor(diet, levels = c("Herbivore",  "Omnivore"))),
  diet
) |>
  rename(diet = facet) |>
  mutate(diet = factor(diet, levels = c("Herbivore", 
                                        "Omnivore")))

guild_ribbon <- overlap_ribbon(guild_kde)

# Pull Dhat1 per guild for facet annotation
guild_delta_labels <- overlap_results |>
  filter(!is.na(delta), !is.na(diet)) |>
  group_by(diet) |>
  summarise(mean_delta = mean(delta), .groups = "drop")

p_diel_guild <- guild_kde |>
  ggplot(aes(hour, density, colour = treatment, fill = treatment)) +
  night_bands +
  geom_ribbon(
    data = guild_ribbon,
    aes(x = hour, ymin = 0, ymax = overlap_density),
    inherit.aes = FALSE,
    fill = "grey60", alpha = 0.35, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_text(
    data = guild_delta_labels,
    aes(x = 22, y = Inf, label = sprintf("\u0394\u02c61 = %.2f", mean_delta)),
    inherit.aes = FALSE, vjust = 1.4, hjust = 1, size = 3, colour = "grey30"
  ) +
  scale_colour_treatment() +
  scale_fill_treatment() +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                     labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  facet_wrap(~diet, nrow = 1) +
  labs(
    title    = "H3: Temporal Activity Overlap by Diet Guild",
    subtitle = "Von Mises KDE curves — grey fill = overlap region (Dhat1); shaded bands = night",
    x = "Hour of day", y = "Kernel density", colour = "Treatment", fill = "Treatment"
  ) +
  theme_bwindi(legend_pos = "top")

print(p_diel_guild)
bwindi_save(p_diel_guild, "H3_diel_guild", width = "double")


# KDE overlap curves by body mass tertile x treatment --------------------

mass_kde_data <- animals_site |>
  left_join(guild_lookup |> select(scientificName, log10_mass_kg), by = "scientificName") |>
  filter(!is.na(log10_mass_kg), !is.na(time_rad)) |>
  mutate(mass_tertile = cut(
    log10_mass_kg,
    breaks = quantile(log10_mass_kg, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
    labels = c("Small-bodied", "Mid-bodied", "Large-bodied"),
    include.lowest = TRUE
  ))

mass_kde <- kde_curves(mass_kde_data, mass_tertile) |>
  rename(mass_tertile = facet) |>
  mutate(mass_tertile = factor(mass_tertile,
                               levels = c("Small-bodied", "Mid-bodied", "Large-bodied")))

mass_ribbon <- overlap_ribbon(mass_kde)

# Mean Dhat1 per mass tertile for annotation
mass_delta_labels <- overlap_results |>
  filter(!is.na(delta), !is.na(log10_mass_kg)) |>
  mutate(mass_tertile = make_mass_tertile(log10_mass_kg)) |>
  group_by(mass_tertile) |>
  summarise(mean_delta = mean(delta), .groups = "drop") |>
  mutate(mass_tertile = factor(mass_tertile,
                               levels = c("Small-bodied", "Mid-bodied", "Large-bodied")))

p_diel_mass <- mass_kde |>
  ggplot(aes(hour, density, colour = treatment, fill = treatment)) +
  night_bands +
  geom_ribbon(
    data = mass_ribbon,
    aes(x = hour, ymin = 0, ymax = overlap_density),
    inherit.aes = FALSE,
    fill = "grey60", alpha = 0.35, colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_text(
    data = mass_delta_labels,
    aes(x = 22, y = Inf, label = sprintf("\u0394\u02c61 = %.2f", mean_delta)),
    inherit.aes = FALSE, vjust = 1.4, hjust = 1, size = 3, colour = "grey30"
  ) +
  scale_colour_treatment() +
  scale_fill_treatment() +
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24),
                     labels = c("00:00", "06:00", "12:00", "18:00", "24:00")) +
  facet_wrap(~mass_tertile, nrow = 1) +
  labs(
    title    = "H3: Temporal Activity Overlap by Body Mass Tertile",
    subtitle = "Von Mises KDE curves — grey fill = overlap region (Dhat1); tertiles for display only",
    x = "Hour of day", y = "Kernel density", colour = "Treatment", fill = "Treatment"
  ) +
  theme_bwindi(legend_pos = "top")

print(p_diel_mass)
bwindi_save(p_diel_mass, "H3_diel_body_mass", width = "double")


# Summary table and interpretation ----------------------------------------

h3_summary <- overlap_results |>
  filter(!is.na(delta)) |>
  select(scientificName, n_road, n_int, delta, ci_low, ci_high,
         overlap_category, diet, log10_mass_kg) |>
  arrange(delta) |>
  mutate(across(c(delta, ci_low, ci_high), round, 3),
         body_mass_kg = round(10^log10_mass_kg, 1))

print(h3_summary, n = Inf)

overall_mean <- mean(overlap_results$delta, na.rm = TRUE)
cat(sprintf("\nOverall community mean Δ̂₁: %.3f\n", overall_mean))

herb_mean <- guild_stats$mean_delta[guild_stats$diet == "Herbivore"]
carn_mean <- guild_stats$mean_delta[guild_stats$diet == "Carnivore"]
if (length(herb_mean) && length(carn_mean)) {
  cat(sprintf("Herbivore: %.3f | Carnivore: %.3f | diff: %.3f — %s\n",
              herb_mean, carn_mean, herb_mean - carn_mean,
              if (herb_mean < carn_mean) "SUPPORTS H3" else "DOES NOT SUPPORT H3"))
}

cat(sprintf("Body mass: ρ = %.3f, p = %.4f — %s\n",
            cor_mass$estimate, cor_mass$p.value,
            if (cor_mass$estimate < 0 && cor_mass$p.value < 0.05)
              "SUPPORTS H3" else "DOES NOT SUPPORT H3"))
