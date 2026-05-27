# 03_H1_richness.R
#
# H1: Species and guild richness will be significantly lower at roadside
# locations, with herbivore guilds showing the greatest reduction while
# carnivore guilds remain stable or increase.
#
# Approach: Poisson GLMMs (pair_id random intercept), iNEXT rarefaction
# to confirm sampling adequacy, and a continuous body mass × treatment
# interaction to test size-sensitivity.



# iNEXT rarefaction -------------------------------------------------------
# Checks that species accumulation curves approach asymptote at the observed
# effort before interpreting raw richness differences between treatments.

inext_input <- list(
  Roadside = build_incidence_freq(animals_site, "Roadside"),
  Interior = build_incidence_freq(animals_site, "Interior")
)

inext_out <- iNEXT(inext_input, q = c(0, 1, 2),
                   datatype = "incidence_freq", endpoint = 36, nboot = 1000)

p_inext <- ggiNEXT(inext_out, type = 1, facet.var = "Order.q") +
  scale_colour_treatment() +
  scale_fill_treatment() +
  labs(
    title    = "H1: Species Diversity Rarefaction & Extrapolation",
    subtitle = "Incidence-based Hill numbers — shaded bands = 95% CI",
    x = "Number of camera stations", y = "Species diversity (Hill number)",
    colour = NULL, fill = NULL
  ) +
  theme_bwindi()

print(p_inext)
bwindi_save(p_inext, "H1_inext_rarefaction", width = "double")


# Overall richness GLMM ---------------------------------------------------

richness_cam <- det_site |>
  select(cam_label, pair_id, treatment, n_species) |>
  mutate(treatment = relevel(factor(treatment), ref = "Interior"))

m_rich <- glmmTMB(n_species ~ treatment + (1 | pair_id),
                  data = richness_cam, family = poisson)

print(summary(m_rich))

emm_rich      <- emmeans(m_rich, ~treatment, type = "response")
contrast_rich <- pairs(emm_rich, reverse = TRUE)
print(emm_rich)
cat("\nContrast (Roadside / Interior):\n")
print(contrast_rich)

emm1_df <- as_tibble(emm_rich)
pct_reduction <- with(emm1_df,
                      (rate[treatment == "Interior"] - rate[treatment == "Roadside"]) /
                        rate[treatment == "Interior"] * 100
)
cat(sprintf("Roadside reduction: %.1f%%\n", pct_reduction))

sim_rich <- simulateResiduals(m_rich, n = 1000, plot = FALSE)
print(testDispersion(sim_rich))
print(testZeroInflation(sim_rich))


# Diet guild GLMMs --------------------------------------------------------

guild_richness_cam <- animals_site |>
  group_by(cam_label, pair_id, treatment, diet) |>
  summarise(n_sp = n_distinct(scientificName), .groups = "drop") |>
  complete(nesting(cam_label, pair_id, treatment), diet, fill = list(n_sp = 0)) |>
  mutate(treatment = relevel(factor(treatment), ref = "Interior"))

# Per-guild models; collect marginal means and contrasts into a summary table
guild_glmm_results <- guild_richness_cam |>
  group_by(diet) |>
  group_map(\(d, g) {
    m   <- glmmTMB(n_sp ~ treatment + (1 | pair_id), data = d, family = poisson)
    emm <- emmeans(m, ~treatment, type = "response")
    cont <- pairs(emm, reverse = TRUE)
    
    emm_df    <- as_tibble(emm)
    int_mean  <- emm_df$rate[emm_df$treatment == "Interior"]
    road_mean <- emm_df$rate[emm_df$treatment == "Roadside"]
    eff_pct   <- (int_mean - road_mean) / int_mean * 100
    
    cont_df <- as_tibble(cont)
    print(summary(m))
    
    tibble(
      diet      = g$diet,
      mean_int  = round(int_mean, 2),
      mean_road = round(road_mean, 2),
      ratio     = round(cont_df$ratio, 3),
      z         = round(cont_df$z.ratio, 3),
      p_value   = round(cont_df$p.value, 4),
      eff_pct   = round(abs(eff_pct), 1),
      direction = if_else(eff_pct >= 0, "reduction", "increase")
    )
  }) |>
  bind_rows()

print(guild_glmm_results)

# Combined model to test whether the road effect differs across guilds
m_combined <- glmmTMB(
  n_sp ~ diet * treatment + (1 | pair_id),
  data = guild_richness_cam, family = poisson
)
print(summary(m_combined))

# Tukey-corrected contrasts within and between guilds
emm_within   <- emmeans(m_combined, ~treatment | diet, type = "response")
emm_between  <- emmeans(m_combined, ~diet | treatment, type = "response")

cat("\nTreatment contrasts within each guild (Tukey):\n")
print(pairs(emm_within, reverse = TRUE, adjust = "tukey"))

cat("\nGuild contrasts within each treatment (Tukey):\n")
print(pairs(emm_between, adjust = "tukey"))

# Interaction contrasts: does the road effect size differ between guilds?
cat("\nInteraction contrasts (Tukey):\n")
print(contrast(
  emmeans(m_combined, ~treatment * diet, type = "response"),
  interaction = "pairwise", adjust = "tukey"
))

# Final summary joining per-guild p-values with Tukey-corrected estimates
tukey_df <- as_tibble(pairs(emm_within, reverse = TRUE, adjust = "tukey")) |>
  select(diet, ratio_tukey = ratio, z_tukey = z.ratio, p_tukey = p.value) |>
  mutate(across(where(is.numeric), round, 3))

guild_final <- guild_glmm_results |>
  left_join(tukey_df, by = "diet") |>
  mutate(tukey_sig = case_when(
    p_tukey < 0.001 ~ "***",
    p_tukey < 0.01  ~ "**",
    p_tukey < 0.05  ~ "*",
    p_tukey < 0.10  ~ ".",
    TRUE            ~ "ns"
  ))

print(guild_final)


# Body mass × treatment GLMM ----------------------------------------------
# A significant negative interaction (log10_mass_kg:treatmentRoadside) means
# larger species show disproportionate avoidance at roadside cameras.

sp_detection <- animals_site |>
  group_by(cam_label, pair_id, treatment, scientificName) |>
  summarise(detected = 1L, .groups = "drop") |>
  complete(nesting(cam_label, pair_id, treatment), scientificName,
           fill = list(detected = 0L)) |>
  left_join(guild_lookup |> select(scientificName, diet, log10_mass_kg),
            by = "scientificName") |>
  filter(!is.na(log10_mass_kg)) |>
  mutate(treatment = relevel(factor(treatment), ref = "Interior"))

m_mass <- glmmTMB(
  detected ~ log10_mass_kg * treatment + (1 | pair_id) + (1 | scientificName),
  data = sp_detection, family = poisson
)
print(summary(m_mass))

pred_mass <- expand_grid(
  log10_mass_kg = seq(min(sp_detection$log10_mass_kg),
                      max(sp_detection$log10_mass_kg), length.out = 200),
  treatment = factor(c("Interior", "Roadside"), levels = c("Interior", "Roadside"))
) |>
  mutate(pred = predict(m_mass, newdata = pick(everything()),
                        re.form = NA, type = "response"))

p_mass_effect <- pred_mass |>
  mutate(body_mass_kg = 10^log10_mass_kg) |>
  ggplot(aes(body_mass_kg, pred, colour = treatment)) +
  geom_line(linewidth = 1.2) +
  scale_colour_treatment() +
  scale_x_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000),
                labels = scales::comma) +
  labs(
    title    = "H1: Detection Probability vs Body Mass",
    subtitle = "Steeper roadside decline = larger species more road-sensitive",
    x = "Body mass (kg)", y = "Detection probability", colour = NULL
  ) +
  theme_bwindi()

print(p_mass_effect)
bwindi_save(p_mass_effect, "H1_bodymass_detection", width = "double")


# Species-level contribution check ----------------------------------------

sp_contribution <- animals_site |>
  group_by(scientificName, treatment, diet, log10_mass_kg) |>
  summarise(total_det = n(), .groups = "drop") |>
  pivot_wider(names_from = treatment, values_from = total_det, values_fill = 0L) |>
  mutate(
    total          = Roadside + Interior,
    roadside_share = Roadside / total,
    complete_avoid = Roadside == 0
  ) |>
  filter(total >= 1) |>
  arrange(diet, roadside_share)

print(sp_contribution)

# Carnivore camera-level check: is the roadside signal spatially concentrated?
carnivore_cam <- animals_site |>
  filter(diet == "Carnivore") |>
  count(cam_label, treatment, name = "n_det") |>
  complete(nesting(cam_label, treatment), fill = list(n_det = 0))

carnivore_cam |>
  filter(treatment == "Roadside") |>
  arrange(desc(n_det)) |>
  print(n = Inf)


# Plots -------------------------------------------------------------------

# A: Paired dot-and-line — per-camera richness
p_rich_paired <- richness_cam |>
  mutate(treatment = factor(treatment, levels = c("Roadside", "Interior"))) |>
  ggplot(aes(treatment, n_species)) +
  geom_line(aes(group = pair_id), colour = "grey70", alpha = 0.6) +
  geom_point(aes(colour = treatment), size = 3) +
  stat_summary(fun = mean, geom = "point", size = 5, shape = 18, colour = "#1a1a1a") +
  stat_summary(fun = mean, geom = "line", aes(group = 1),
               linewidth = 1, colour = "#1a1a1a", linetype = "dashed") +
  scale_colour_treatment(guide = "none") +
  labs(
    title    = "H1: Per-Camera Species Richness by Treatment",
    subtitle = "Lines connect paired cameras — diamond = treatment mean",
    x = NULL, y = "Species detected"
  ) +
  theme_bwindi()

print(p_rich_paired)
bwindi_save(p_rich_paired, "H1_richness_paired", width = "single")


# B: Species roadside share
pal_diet_binary <- c(Herbivore = "#6B8E23", Carnivore = "#D7263D", Omnivore = "#d98b19")

p_species_road_share <- sp_contribution |>
  mutate(
    label = glue::glue("{scientificName}  (n={total})"),
    label = fct_reorder(label, roadside_share)
  ) |>
  ggplot(aes(roadside_share, label, fill = diet)) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30", linewidth = 0.6) +
  geom_text(
    data  = \(x) filter(x, complete_avoid),
    aes(x = 0.03, label = "zero roadside"),
    hjust = 0, size = 3, colour = "grey20", fontface = "italic"
  ) +
  annotate("text", x = 0.48, y = 0.5, label = "Equal split",
           hjust = 1, size = 3, colour = "grey40", angle = 90) +
  scale_fill_manual(values = pal_diet_binary) +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1),
                     expand = c(0, 0.01)) +
  labs(
    title    = "Proportion of Detections at Roadside per Species",
    subtitle = "Roadside TN = 682, Interior TN = 712",
    x = "Proportion of detections at roadside cameras", y = NULL, fill = "Diet guild"
  ) +
  theme_bwindi()

print(p_species_road_share)
bwindi_save(p_species_road_share, "H1_species_roadside_share", width = "double")


# C & D: Guild richness — boxplot and violin
guild_richness_cam <- animals_site |>
  group_by(cam_label, treatment, diet) |>
  summarise(n_sp = n_distinct(scientificName), .groups = "drop") |>
  complete(nesting(cam_label, treatment), diet, fill = list(n_sp = 0)) |>
  mutate(treatment = factor(treatment, levels = c("Roadside", "Interior")))

p_guild_rich <- rai_boxplot(
  guild_richness_cam, x = treatment, y = n_sp, facet = diet,
  ylab     = "Species detected",
  title    = "H1: Diet Guild Richness by Treatment",
  subtitle = "Each point = one camera"
)

p_guild_rich_violin <- guild_richness_cam |>
  ggplot(aes(treatment, n_sp, fill = treatment, colour = treatment)) +
  geom_violin(width = 0.5, alpha = 0.18, linewidth = 0.7, trim = TRUE) +
  geom_quasirandom(width = 0.15, size = 2.5, alpha = 0.85,
                   shape = 21, stroke = 0.6, fill = "white") +
  stat_summary(fun = median, geom = "crossbar", width = 0.3,
               colour = "#000000", linewidth = 0.2, alpha = 0.8) +
  scale_fill_treatment() +
  scale_colour_treatment() +
  facet_wrap(~diet, scales = "free_y", nrow = 1) +
  labs(
    title    = "H1: Diet Guild Richness by Treatment",
    subtitle = "Violin = distribution | Bar = median | Points = cameras (n = 18 per treatment)",
    y = "Species detected per camera", fill = NULL, colour = NULL
  ) +
  theme_bwindi() +
  theme(legend.position = "none")

print(p_guild_rich)
print(p_guild_rich_violin)
bwindi_save(p_guild_rich,        "H1_guild_richness",        width = "double")
bwindi_save(p_guild_rich_violin, "H1_guild_richness_violin", width = "double")

plot_grid(p_guild_rich_violin, p_mass_effect, labels = "AUTO", ncol = 2)