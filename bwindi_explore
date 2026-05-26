library(tidyverse)
library(lubridate)
library(readxl)

# ── Paths ─────────────────────────────────────────────────────────────────────

AGOUTI   <- file.path(ONEDRIVE, "03_Processed_Data/Agouti_Exports/annotations-unvalidated")
SITE_XLS <- file.path(ONEDRIVE, "03_Processed_Data/Site_Info/CT_Site_Info.xlsx")
OUT_DIR  <- file.path(ONEDRIVE, "05_R_Analysis/Outputs")


# ══════════════════════════════════════════════════════════════════════════════
# 1. LOAD & CLEAN AGOUTI EXPORTS
# ══════════════════════════════════════════════════════════════════════════════

deployments  <- read_csv(file.path(AGOUTI, "deployments_updated.csv"), show_col_types = FALSE)
observations <- read_csv(file.path(AGOUTI, "observations.csv"),        show_col_types = FALSE)

deps_clean <- deployments %>%
  mutate(
    deploymentStart = ymd_hms(deploymentStart, quiet = TRUE),
    deploymentEnd   = ymd_hms(deploymentEnd,   quiet = TRUE),
    trap_nights = {
      tn <- as.numeric(difftime(deploymentEnd, deploymentStart, units = "days"))
      if_else(tn > 0, tn, NA_real_)
    },
    raw_label = str_remove(locationName, "^Bwindi\\s*-\\s*"),
    cam_label = str_sub(str_trim(raw_label), 1, 3),
    treatment = case_when(
      str_starts(cam_label, "F") ~ "Roadside",
      str_starts(cam_label, "C") ~ "Interior",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(treatment))


# ── Camera-level deployment summary ───────────────────────────────────────────

cam_meta <- deps_clean %>%
  group_by(cam_label, treatment) %>%
  summarise(
    n_batches      = n(),
    trap_nights    = sum(trap_nights, na.rm = TRUE),
    latitude       = first(na.omit(latitude)),
    longitude      = first(na.omit(longitude)),
    study_start    = min(deploymentStart, na.rm = TRUE),
    study_end      = max(deploymentEnd,   na.rm = TRUE),
    deployment_ids = paste(str_sub(deploymentID, 1, 8), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(treatment, cam_label)


# ── Verification ───────────────────────────────────────────────────────────────

message("── Deployment checks ────────────────────────────────────────")

# Camera counts
cam_meta %>% count(treatment) %>% print()

# Multi-batch cameras
deps_clean %>%
  filter(cam_label %in% filter(cam_meta, n_batches > 1)$cam_label) %>%
  select(cam_label, raw_label, deploymentStart, deploymentEnd, trap_nights) %>%
  arrange(cam_label, deploymentStart) %>%
  print(n = Inf)

# Invalid trap nights
bad_tn <- deps_clean %>% filter(is.na(trap_nights) | trap_nights == 0)
if (nrow(bad_tn) == 0) {
  message("  All trap nights valid ✓")
} else {
  print(bad_tn %>% select(raw_label, cam_label, deploymentStart, deploymentEnd, trap_nights))
}


# ══════════════════════════════════════════════════════════════════════════════
# 2. OBSERVATIONS — JOIN & FILTER
# ══════════════════════════════════════════════════════════════════════════════

obs_drop <- c("mediaID", "cameraSetupType", "lifeStage", "sex", "behavior",
              "individualID", "individualPositionRadius", "individualPositionAngle",
              "individualSpeed", "bboxX", "bboxY", "bboxHeight", "bboxWidth",
              "classificationMethod", "classifiedBy", "classificationTimestamp",
              "classificationProbability", "observationTags", "observationComments")

obs_joined <- observations %>%
  mutate(
    eventStart = ymd_hms(eventStart, quiet = TRUE),
    eventEnd   = ymd_hms(eventEnd,   quiet = TRUE),
    hour       = hour(eventStart),
    date       = as_date(eventStart)
  ) %>%
  left_join(
    deps_clean %>% select(deploymentID, raw_label, cam_label, treatment),
    by = "deploymentID"
  ) %>%
  select(-any_of(obs_drop))

n_unmatched <- sum(is.na(obs_joined$cam_label))
message(sprintf("  Unmatched observations (0 is goood!): %d / %d (%.1f%%)",
                n_unmatched, nrow(obs_joined), 100 * n_unmatched / nrow(obs_joined)))


# ── Species exclusions ─────────────────────────────────────────────────────────
# Excluded: birds (Aves), rodents (Rodentia + specific spp.), arboreal primates
# that rarely trigger ground-level cameras and inflate detection counts.

exclude_spp <- c(
  "Aves", "Rodentia",
  "Cricetomys gambianus", "Rattus rattus", "Protoxerus stangeri",
  "Colobus guereza", "Cercopithecus mitis"
)

animals_raw <- obs_joined %>%
  filter(
    observationType == "animal",
    !is.na(scientificName),
    scientificName != "Homo sapiens",
    !is.na(cam_label),
    !scientificName %in% exclude_spp
  )


# ══════════════════════════════════════════════════════════════════════════════
# 3. GUILD LOOKUP & BODY MASS
# ══════════════════════════════════════════════════════════════════════════════
# Body mass (kg) from PanTHERIA / MammalBase via GBIF Backbone.
# log10 transform used in models — mass spans ~4 orders of magnitude.

guild_lookup <- tribble(
  ~scientificName,           ~diet,        ~activity,    ~specialist, ~body_mass_kg,
  
  "Genetta victoriae",       "Carnivore",  "Nocturnal",  "No",         3.0,
  "Leptailurus serval",      "Carnivore",  "Nocturnal",  "No",        10.5,
  "Cephalophus nigrifrons",  "Herbivore",  "Cathemeral", "Yes",       14.5,
  "Canis adustus",           "Omnivore",   "Nocturnal",  "No",         7.5,
  "Civettictis civetta",     "Omnivore",   "Nocturnal",  "No",        13.5,
  "Caracal aurata",          "Carnivore",  "Nocturnal",  "Yes",       10.0,
  "Mellivora capensis",      "Omnivore",   "Cathemeral", "No",         9.75,
  "Allochrocebus lhoesti",   "Omnivore",   "Diurnal",    "No",         6.75,
  "Cephalophus silvicultor", "Herbivore",  "Cathemeral", "No",        62.5,
  "Pan troglodytes",         "Omnivore",   "Diurnal",    "No",        45.0,
  "Potamochoerus larvatus",  "Omnivore",   "Nocturnal",  "Yes",       82.5,
  "Tragelaphus sylvaticus",  "Herbivore",  "Cathemeral", "Yes",       52.0,
  "Gorilla beringei",        "Herbivore",  "Diurnal",    "Yes",      135.0,
  "Loxodonta cyclotis",      "Herbivore",  "Cathemeral", "Yes",     4350.0
) %>%
  mutate(
    log10_mass_kg = log10(body_mass_kg),
    diet       = factor(diet,       levels = c("Herbivore", "Carnivore", "Omnivore")),
    activity   = factor(activity,   levels = c("Diurnal", "Nocturnal", "Cathemeral")),
    specialist = factor(specialist, levels = c("Yes", "No"))
  )

animals_raw <- animals_raw %>%
  left_join(guild_lookup, by = "scientificName")


# ══════════════════════════════════════════════════════════════════════════════
# 4. EVENT COLLAPSING — 30-MIN INDEPENDENCE THRESHOLD
# ══════════════════════════════════════════════════════════════════════════════

collapse_events <- function(df, gap_minutes = 30) {
  df %>%
    arrange(cam_label, scientificName, eventStart) %>%
    group_by(cam_label, scientificName) %>%
    mutate(
      last_kept = accumulate(eventStart, function(anchor, current) {
        if (as.numeric(difftime(current, anchor, units = "mins")) > gap_minutes)
          current else anchor
      }),
      keep = eventStart == last_kept
    ) %>%
    filter(keep) %>%
    select(-last_kept, -keep) %>%
    ungroup()
}

animals <- collapse_events(animals_raw)

# Collapsing summary
before <- nrow(animals_raw)
after  <- nrow(animals)
message(sprintf(
  "── Event collapsing: %d → %d detections (%d removed, %.1f%%)",
  before, after, before - after, 100 * (before - after) / before
))

left_join(
  animals_raw %>% count(scientificName, name = "before"),
  animals     %>% count(scientificName, name = "after"),
  by = "scientificName"
) %>%
  mutate(removed = before - after, pct = round(100 * removed / before, 1)) %>%
  arrange(desc(removed)) %>%
  print(n = Inf)


# ── Independence validation ────────────────────────────────────────────────────

violations <- animals %>%
  arrange(cam_label, scientificName, eventStart) %>%
  group_by(cam_label, scientificName) %>%
  mutate(gap = as.numeric(difftime(eventStart, lag(eventStart), units = "mins"))) %>%
  filter(!is.na(gap), gap < 30) %>%
  ungroup()

if (nrow(violations) == 0) {
  message("  All retained events ≥30 min apart ✓")
} else {
  message(sprintf("  ✗ %d independence violation(s)", nrow(violations)))
  print(violations %>% select(cam_label, scientificName, eventStart, gap))
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. CAMERA-LEVEL DETECTION SUMMARY
# ══════════════════════════════════════════════════════════════════════════════

det_summary <- obs_joined %>%
  filter(!is.na(cam_label)) %>%
  group_by(cam_label, treatment) %>%
  summarise(
    n_events        = n(),
    n_animal_events = sum(observationType == "animal" &
                            !is.na(scientificName) &
                            scientificName != "Homo sapiens"),
    n_blank         = sum(observationType == "blank"),
    n_species       = n_distinct(scientificName[
      observationType == "animal" &
        !is.na(scientificName) &
        scientificName != "Homo sapiens"
    ]),
    .groups = "drop"
  ) %>%
  left_join(cam_meta %>% select(cam_label, trap_nights), by = "cam_label") %>%
  mutate(RAI = (n_animal_events / trap_nights) * 100)

cam_effort <- det_summary %>% select(cam_label, treatment, trap_nights)

species_summary <- animals %>%
  group_by(scientificName) %>%
  summarise(
    total_detections = n(),
    n_cameras        = n_distinct(cam_label),
    roadside_det     = sum(treatment == "Roadside"),
    interior_det     = sum(treatment == "Interior"),
    .groups = "drop"
  ) %>%
  arrange(desc(total_detections))


# ══════════════════════════════════════════════════════════════════════════════
# 6. HUMAN ACTIVITY
# ══════════════════════════════════════════════════════════════════════════════
# Vehicle detections log both a "vehicle" and "human" row under the same
# eventID — deduplicate to one row per event before summarising.

human_activity <- obs_joined %>%
  filter(observationType %in% c("vehicle", "human"), !is.na(cam_label)) %>%
  distinct(eventID, cam_label, treatment, hour, date, eventStart)

human_summary <- human_activity %>%
  count(cam_label, treatment, name = "n_events") %>%
  complete(nesting(cam_label, treatment), fill = list(n_events = 0)) %>%
  left_join(cam_effort, by = c("cam_label", "treatment")) %>%
  mutate(RAI = (n_events / trap_nights) * 100)


# ══════════════════════════════════════════════════════════════════════════════
# 7. SITE INFORMATION
# ══════════════════════════════════════════════════════════════════════════════

site_raw <- read_excel(SITE_XLS)

site <- site_raw %>%
  rename(
    cam_label              = `Camera ID`,
    sd_card_id             = `SD Card ID`,
    form_start             = start,
    form_end               = end,
    deployment_date        = `Deployment date`,
    deployment_time        = `Deployment time (24hr)`,
    lat_gps                = `_GPS location_latitude`,
    lon_gps                = `_GPS location_longitude`,
    alt_gps                = `_GPS location_altitude`,
    elevation_m            = `Elevation (m.a.s.l.)`,
    dist_road_m            = `Distance to road edge (m)`,
    dist_park_boundary_m   = `Distance to park boundary (m)`,
    dist_water_m           = `Distance to nearest water (m0`,
    cam_height_cm          = `Camera height above ground (cm)`,
    dist_detection_zone_m  = `Distance to detection zone (m)`,
    cam_bearing_deg        = `Camera bearing (degrees)`,
    tree_dbh_cm            = `Tree DBH (cm)`,
    canopy_pct_1           = `Canopy cover reading 1 (%)`,
    canopy_pct_2           = `Canopy cover reading 2 (%)`,
    canopy_pct_3           = `Canopy cover reading 3 (%)`,
    canopy_mean_pct        = canopy_mean,
    understory_density_pct = `Understory density (%)`,
    detection_range_m      = `Detection range (m)`,
    veg_1                  = `Dominant vegetation 1`,
    veg_2                  = `Dominant vegetation 2`,
    veg_3                  = `Dominant vegetation 3`,
    game_trail_present     = `Game trail present`,
    trail_dist_m           = `Trail distance from camera (m)`,
    trail_width_m          = `Trail width (m)`,
    downed_logs_present    = `Downed logs present`,
    n_logs                 = `Number of logs`,
    largest_log_diam_cm    = `Largest log diameter (cm)`,
    log_dist_m             = `Log distance from camera`,
    road_width_m           = `Road width (m)`,
    road_surface_condition = `Road surface condition`,
    road_surface_type      = `Surface type`,
    vehicles_observed      = `Vehicles observed`,
    vehicle_types          = `Vehicle types observed`,
    veh_tourist            = `Vehicle types observed/Tourist`,
    veh_ranger             = `Vehicle types observed/Ranger`,
    veh_local              = `Vehicle types observed/Local`,
    veh_other              = `Vehicle types observed/Other`,
    veg_clearance_width_m  = `Vegetation clearance width (m)`,
    road_edge_condition    = `Road edge condition`,
    road_surface_photo     = `Road surface photo`,
    notes                  = `Additional notes`,
    form_index             = `_index`
  ) %>%
  mutate(
    treatment = case_when(
      str_starts(cam_label, "F") ~ "Roadside",
      str_starts(cam_label, "C") ~ "Interior",
      TRUE ~ NA_character_
    ),
    game_trail_present  = game_trail_present == "Yes",
    downed_logs_present = downed_logs_present == "Yes",
    canopy_mean_pct     = as.numeric(canopy_mean_pct),
    canopy_mean_recalc  = rowMeans(cbind(canopy_pct_1, canopy_pct_2, canopy_pct_3), na.rm = TRUE)
  )


# ── Pair IDs ──────────────────────────────────────────────────────────────────

pair_lookup <- tribble(
  ~cam_label, ~partner,
  "F01", "C04",  "F02", "C12",  "F03", "C05",  "F04", "C03",
  "F05", "C01",  "F06", "C07",  "F07", "C06",  "F08", "C08",
  "F09", "C10",  "F10", "C09",  "F11", "C11",  "F12", "C14",
  "F13", "C02",  "F14", "C18",  "F15", "C15",  "F16", "C16",
  "F17", "C17",  "F18", "C13"
) %>%
  mutate(pair_id = as.integer(str_extract(cam_label, "\\d+"))) %>%
  pivot_longer(c(cam_label, partner), values_to = "cam_label") %>%
  select(cam_label, pair_id) %>%
  arrange(pair_id, cam_label)

site <- left_join(site, pair_lookup, by = "cam_label")


# ── Site checks ───────────────────────────────────────────────────────────────

message("── Site checks ──────────────────────────────────────────────")
message(sprintf("  Rows: %d (expect 36)", nrow(site)))
site %>% count(treatment) %>% print()

# Pair completeness
pair_counts <- site %>% count(pair_id)
if (all(pair_counts$n == 2)) {
  message("  All 18 pairs have 2 cameras ✓")
} else {
  print(pair_counts %>% filter(n != 2))
}

# GPS coverage
missing_gps <- site %>% filter(is.na(lat_gps) | is.na(lon_gps))
if (nrow(missing_gps) == 0) message("  GPS complete ✓") else print(missing_gps %>% select(cam_label, lat_gps, lon_gps))

# Canopy mean cross-check
canopy_mismatch <- site %>%
  mutate(diff = abs(canopy_mean_pct - canopy_mean_recalc)) %>%
  filter(diff > 0.1)
if (nrow(canopy_mismatch) == 0) message("  Canopy means match ✓") else print(canopy_mismatch %>% select(cam_label, canopy_mean_pct, canopy_mean_recalc))

# Road variables should only appear on roadside cameras
interior_road <- site %>% filter(treatment == "Interior", !is.na(road_width_m))
if (nrow(interior_road) == 0) message("  Road variables clean ✓") else print(interior_road %>% select(cam_label, road_width_m))


# ── Sub-tables ────────────────────────────────────────────────────────────────

site_core <- site %>%
  select(
    cam_label, treatment, pair_id,
    lat_gps, lon_gps, elevation_m,
    dist_road_m, dist_park_boundary_m, dist_water_m,
    cam_height_cm, cam_bearing_deg, detection_range_m,
    canopy_mean_pct, canopy_pct_1, canopy_pct_2, canopy_pct_3,
    understory_density_pct, tree_dbh_cm,
    veg_1, veg_2, veg_3,
    game_trail_present, trail_dist_m, trail_width_m,
    downed_logs_present, n_logs, largest_log_diam_cm, log_dist_m,
    deployment_date, notes
  ) %>%
  arrange(treatment, cam_label)

road_chars <- site %>%
  filter(treatment == "Roadside") %>%
  select(
    cam_label,
    road_width_m, road_surface_condition, road_surface_type,
    road_edge_condition, veg_clearance_width_m,
    vehicles_observed, vehicle_types,
    veh_tourist, veh_ranger, veh_local, veh_other
  ) %>%
  arrange(cam_label)


# ── Joined tables ─────────────────────────────────────────────────────────────

cam_site <- left_join(cam_meta, site_core, by = c("cam_label", "treatment"))

det_site <- det_summary %>%
  left_join(site_core, by = c("cam_label", "treatment")) %>%
  left_join(road_chars, by = "cam_label")

animals_site <- animals %>%
  left_join(
    site_core %>% select(cam_label, pair_id, elevation_m, dist_road_m,
                         canopy_mean_pct, understory_density_pct,
                         dist_water_m, game_trail_present),
    by = "cam_label"
  )

message(sprintf("  animals_site: %d rows, %d missing elevation",
                nrow(animals_site), sum(is.na(animals_site$elevation_m))))


# ══════════════════════════════════════════════════════════════════════════════
# 8. PLOTS
# ══════════════════════════════════════════════════════════════════════════════

pal <- c("Roadside" = "#E05C3A", "Interior" = "#3A7DBE")

# P1: Trap nights — effort QC
p1 <- cam_meta %>%
  mutate(cam_label = fct_reorder(cam_label, trap_nights)) %>%
  ggplot(aes(trap_nights, cam_label, fill = treatment)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(title = "Trap nights per camera",
       subtitle = "Summed across all Agouti batches",
       x = "Total trap nights", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# P2: RAI per camera
p2 <- det_summary %>%
  mutate(cam_label = fct_reorder(cam_label, RAI)) %>%
  ggplot(aes(RAI, cam_label, fill = treatment)) +
  geom_col() +
  scale_fill_manual(values = pal) +
  labs(title = "RAI per camera",
       subtitle = "Animal detections per 100 trap-nights",
       x = "RAI", y = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# P3: Total detections by species
p3 <- species_summary %>%
  mutate(scientificName = fct_reorder(scientificName, total_detections)) %>%
  ggplot(aes(total_detections, scientificName)) +
  geom_col(fill = "#5A9E6F") +
  labs(title = "Total detections by species",
       x = "Detections", y = NULL) +
  theme_minimal(base_size = 12)

# P4: Diel activity
p4 <- animals %>%
  count(hour, treatment) %>%
  ggplot(aes(hour, n, fill = treatment)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = 0:23) +
  labs(title = "Diel activity — all species",
       x = "Hour of day", y = "Detections", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# P5: Camera locations
p5 <- cam_meta %>%
  ggplot(aes(longitude, latitude, colour = treatment, label = cam_label)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 2.8) +
  scale_colour_manual(values = pal) +
  coord_fixed() +
  labs(title = "Camera trap locations",
       x = "Longitude", y = "Latitude", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

walk(list(p1, p2, p3, p4, p5), print)


# ══════════════════════════════════════════════════════════════════════════════
# 9. EXPORT
# ══════════════════════════════════════════════════════════════════════════════

save(
  cam_meta, cam_effort, cam_site, pair_lookup,
  site_core, road_chars,
  det_summary, det_site,
  animals, animals_site,
  guild_lookup, species_summary,
  human_activity, human_summary,
  file = file.path(OUT_DIR, "bwindi_analysis_data.RData")
)

message("Saved: bwindi_analysis_data.RData")
message(sprintf("  %-20s %d × %d", "cam_meta",      nrow(cam_meta),      ncol(cam_meta)))
message(sprintf("  %-20s %d × %d", "det_summary",   nrow(det_summary),   ncol(det_summary)))
message(sprintf("  %-20s %d × %d", "animals",       nrow(animals),       ncol(animals)))
message(sprintf("  %-20s %d × %d", "animals_site",  nrow(animals_site),  ncol(animals_site)))
message(sprintf("  %-20s %d × %d", "guild_lookup",  nrow(guild_lookup),  ncol(guild_lookup)))
