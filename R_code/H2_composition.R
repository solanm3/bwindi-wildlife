# 04_H2_composition.R
#
# H2: Community composition will differ significantly between roadside and
# interior, with large-bodied herbivores underrepresented and
# carnivores/omnivores overrepresented at roadside cameras.
#
# Approach: PERMANOVA on Hellinger-transformed RAI, betadisper for
# dispersion, NMDS ordination, SIMPER to identify driving species.
# Body size treated as continuous log10(kg) throughout.
# --------------------------------------------------------------------------


# Camera × species RAI matrix ---------------------------------------------
# cam_effort joined before complete() so trap_nights survives the expansion.

rai_matrix_long <- animals_site |>
  group_by(cam_label, treatment, pair_id, scientificName) |>
  summarise(n_det = n(), .groups = "drop") |>
  left_join(cam_effort |> select(cam_label, trap_nights), by = "cam_label") |>
  mutate(RAI = (n_det / trap_nights) * 100) |>
  complete(
    nesting(cam_label, treatment, pair_id, trap_nights),
    scientificName,
    fill = list(n_det = 0L, RAI = 0)
  )

rai_wide <- rai_matrix_long |>
  select(cam_label, treatment, pair_id, scientificName, RAI) |>
  pivot_wider(names_from = scientificName, values_from = RAI, values_fill = 0)

meta_h2 <- rai_wide |> select(cam_label, treatment, pair_id)
sp_mat  <- rai_wide |> select(-cam_label, -treatment, -pair_id) |> as.matrix()
rownames(sp_mat) <- meta_h2$cam_label

sp_hell <- decostand(sp_mat, method = "hellinger")


# PERMANOVA ---------------------------------------------------------------

set.seed(42)
perm_result <- adonis2(
  sp_hell ~ treatment,
  data = meta_h2, method = "euclidean",
  strata = meta_h2$pair_id, permutations = 999
)
print(perm_result)


# Betadisper --------------------------------------------------------------
# Tests whether treatments differ in within-group spread, which would
# confound the PERMANOVA centroid test.

dist_hell <- dist(sp_hell, method = "euclidean")
betadisp  <- betadisper(dist_hell, meta_h2$treatment)
beta_perm <- permutest(betadisp, permutations = 999)
print(beta_perm)


# NMDS --------------------------------------------------------------------

set.seed(42)
nmds <- metaMDS(sp_hell, distance = "euclidean", k = 2, trymax = 100, trace = FALSE)

cat(sprintf("NMDS stress = %.3f", nmds$stress))
if      (nmds$stress < 0.10) {
  cat(" (excellent)\n")
} else if (nmds$stress < 0.20) {
  cat(" (acceptable)\n")
}else                          {
      cat(" (poor — consider k = 3)\n")
  }

nmds_scores <- as_tibble(scores(nmds, display = "sites")) |> bind_cols(meta_h2)

sp_scores <- as_tibble(scores(nmds, display = "species"), rownames = "scientificName") |>
  left_join(guild_lookup |> select(scientificName, diet, log10_mass_kg),
            by = "scientificName")


# Plots -------------------------------------------------------------------

# Betadisper ordination
betadisp_pts <- as_tibble(betadisp$vectors) |>
  bind_cols(meta_h2) |>
  rename(Dim1 = PCoA1, Dim2 = PCoA2)

betadisp_ctr <- as_tibble(betadisp$centroids) |>
  mutate(treatment = rownames(betadisp$centroids)) |>
  rename(Dim1 = PCoA1, Dim2 = PCoA2)

p_betadisp <- ggplot() +
  geom_segment(
    data = left_join(betadisp_pts,
                     betadisp_ctr |> rename(cDim1 = Dim1, cDim2 = Dim2),
                     by = "treatment"),
    aes(x = Dim1, y = Dim2, xend = cDim1, yend = cDim2, colour = treatment),
    alpha = 0.3, linewidth = 0.4
  ) +
  geom_point(data = betadisp_pts, aes(Dim1, Dim2, colour = treatment),
             size = 2.5, alpha = 0.8) +
  geom_point(data = betadisp_ctr, aes(Dim1, Dim2, colour = treatment),
             size = 5, shape = 18) +
  scale_colour_treatment() +
  labs(
    title    = "Within-Group Dispersion (betadisper)",
    subtitle = sprintf("betadisper p = %.3f, df = 34", beta_perm$tab$`Pr(>F)`[1]),
    x = "PCoA Axis 1", y = "PCoA Axis 2", colour = NULL
  ) +
  theme_bwindi()

print(p_betadisp)
bwindi_save(p_betadisp, "H2_betadisp", width = "single")


# Distance to centroid boxplot
p_betadisp_box <- tibble(distance = betadisp$distances, treatment = meta_h2$treatment) |>
  ggplot(aes(treatment, distance, fill = treatment, colour = treatment)) +
  geom_boxplot(width = 0.4, alpha = 0.75, outlier.shape = 21, linewidth = 0.45) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.65, fill = "white",
              shape = 21, stroke = 0.8) +
  scale_fill_treatment(guide = "none") +
  scale_colour_treatment(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, NA)) +
  labs(
    title    = "H2: Distance to Treatment Centroid",
    subtitle = "Larger spread = more variable community composition",
    x = NULL, y = "Distance to centroid"
  ) +
  theme_bwindi()

print(p_betadisp_box)
bwindi_save(p_betadisp_box, "H2_betadisp_box", width = "single")


# NMDS ordination — camera points + species overlay sized by body mass
p_nmds <- ggplot() +
  stat_ellipse(
    data = nmds_scores,
    aes(NMDS1, NMDS2, colour = treatment, fill = treatment),
    geom = "polygon", alpha = 0.12, level = 0.95, linetype = "dashed"
  ) +
  geom_point(data = nmds_scores, aes(NMDS1, NMDS2, colour = treatment),
             size = 3, alpha = 0.9) +
  geom_text(data = nmds_scores, aes(NMDS1, NMDS2, label = cam_label, colour = treatment),
            size = 2.5, vjust = -0.8) +
  geom_point(data = sp_scores,
             aes(NMDS1, NMDS2, colour = diet, size = log10_mass_kg),
             alpha = 0.75, shape = 17) +
  geom_text(data = sp_scores, aes(NMDS1, NMDS2, label = scientificName),
            size = 2.5, colour = "grey35", fontface = "italic", vjust = -0.9) +
  annotate("text", x = Inf, y = -Inf,
           label = sprintf("Stress = %.3f", nmds$stress),
           hjust = 1.1, vjust = -0.5, size = 3.5, colour = "grey50") +
  scale_colour_manual(values = c(pal_treatment, pal_diet), na.value = "grey70") +
  scale_fill_treatment() +
  scale_size_continuous(
    range  = c(1.5, 5), name = "Body mass (kg)",
    breaks = log10(c(8, 40, 100, 4000)), labels = c("8", "40", "100", "4,000")
  ) +
  labs(
    title    = "NMDS Ordination of Camera Communities",
    subtitle = "Hellinger-transformed RAI; 95% confidence ellipses\nTriangles = species, sized by body mass (kg)",
    colour = NULL, fill = NULL
  ) +
  theme_bwindi(legend_pos = "right")

print(p_nmds)
bwindi_save(p_nmds, "H2_nmds", width = "double")


# Mean RAI per species, faceted by body mass tertile
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

sp_rai_summary <- rai_matrix_long |>
  left_join(guild_lookup |> select(scientificName, log10_mass_kg, diet),
            by = "scientificName") |>
  filter(!is.na(log10_mass_kg)) |>
  group_by(scientificName, treatment, log10_mass_kg, diet) |>
  summarise(mean_RAI = mean(RAI), se_RAI = sd(RAI) / sqrt(n()), .groups = "drop") |>
  mutate(mass_tertile = make_mass_tertile(log10_mass_kg))

p_sp_rai <- sp_rai_summary |>
  mutate(scientificName = fct_reorder(scientificName, mean_RAI, sum)) |>
  ggplot(aes(mean_RAI, scientificName, colour = treatment)) +
  geom_linerange(aes(xmin = mean_RAI - se_RAI, xmax = mean_RAI + se_RAI),
                 position = position_dodge(0.5), linewidth = 0.8) +
  geom_point(aes(size = log10_mass_kg), position = position_dodge(0.5), alpha = 0.85) +
  scale_colour_treatment() +
  scale_size_body(range = c(2, 6)) +
  facet_wrap(~mass_tertile, scales = "free_y", ncol = 1) +
  labs(
    title    = "H2: Mean RAI per Species — Roadside vs Interior",
    subtitle = "Mean ± SE across cameras; faceted by body mass tertile\nPoint size = body mass (kg)",
    x = "Mean RAI (detections per 100 trap-nights)", y = NULL, colour = NULL
  ) +
  theme_bwindi(legend_pos = "right") +
  theme(axis.text.y = element_text(face = "italic"))

print(p_sp_rai)
bwindi_save(p_sp_rai, "H2_species_rai", width = "single", height = 160)


# Community composition by diet guild
p_comp_guild <- animals_site |>
  filter(!is.na(diet)) |>
  count(treatment, diet) |>
  group_by(treatment) |>
  mutate(prop = n / sum(n), diet = factor(diet, levels = c("Herbivore", "Omnivore", "Carnivore"))) |>
  ungroup() |>
  ggplot(aes(treatment, prop, fill = diet)) +
  geom_col(width = 0.5) +
  scale_fill_diet() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "H2: Community Composition by Diet Guild",
    subtitle = "Proportion of all independent detections",
    x = NULL, y = "Proportion of detections", fill = NULL
  ) +
  theme_bwindi()

print(p_comp_guild)
bwindi_save(p_comp_guild, "H2_comp_diet", width = "single")

glimpse(animals_site)
# Community composition by body mass tertile
comp_size_data <- animals_site |>
  filter(!is.na(log10_mass_kg)) |>
  mutate(mass_tertile = make_mass_tertile(log10_mass_kg)) |>
  count(treatment, mass_tertile) |>
  group_by(treatment) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

p_comp_size <- comp_size_data |>
  ggplot(aes(treatment, prop, fill = mass_tertile)) +
  geom_col(width = 0.5, alpha = 0.85) +
  scale_fill_manual(
    values = c("Small-bodied" = "#008080", "Mid-bodied" = "#D7263D",
               "Large-bodied" = "#F4B400"),
    name = "Body mass tertile"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = "H2: Community Composition by Body Mass",
    subtitle = "Proportion of detections; tertiles from log\u2081\u2080(kg) distribution",
    x = NULL, y = "Proportion of detections", fill = NULL
  ) +
  theme_bwindi()

print(p_comp_size)
bwindi_save(p_comp_size, "H2_comp_size", width = "single")


# SIMPER — species driving the roadside vs interior separation
simp    <- simper(sp_mat, meta_h2$treatment, permutations = 999)
simp_df <- summary(simp, ordered = TRUE)[["Interior_Roadside"]] |>
  as_tibble(rownames = "scientificName") |>
  filter(cumsum(average) <= 0.70 | row_number() <= 5) |>
  left_join(guild_lookup |> select(scientificName, diet, log10_mass_kg),
            by = "scientificName") |>
  mutate(scientificName = fct_reorder(scientificName, average))

p_simper <- simp_df |>
  ggplot(aes(average, scientificName, fill = diet, size = log10_mass_kg)) +
  geom_col(width = 0.65, alpha = 0.85, show.legend = c(fill = TRUE, size = FALSE)) +
  geom_point(aes(x = average + 0.001, size = log10_mass_kg),
             colour = "grey30", alpha = 0.7, shape = 21) +
  scale_fill_diet() +
  scale_size_body(range = c(2, 7)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "Species Contributions to Community Dissimilarity (SIMPER)",
    subtitle = "Top species explaining \u226470% of Interior \u2013 Roadside Bray-Curtis dissimilarity\nPoint size = body mass (kg)",
    x = "Average contribution to dissimilarity", y = NULL, fill = "Diet guild"
  ) +
  theme_bwindi(legend_pos = "right") +
  theme(axis.text.y = element_text(face = "italic"))

print(p_simper)
bwindi_save(p_simper, "H2_simper", width = "double")
