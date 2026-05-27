# 01_load_data.R
# Load pre-processed RData objects and join weather covariates.
# DATA_DIR must be set before sourcing, e.g. in .Rprofile or a project setup script.

load(file.path(DATA_DIR, "bwindi_analysis_data.RData"))
load(file.path(DATA_DIR, "bwindi_weather_data.RData"))


# Join deployment-period weather onto the camera-level summary
det_site <- det_site |>
  left_join(
    weather_deployment |>
      select(cam_label, mean_temp_c, mean_rh_pct, mean_wind_ms,
             mean_soil_moisture, sd_temp_c, max_temp_c, min_temp_c, n_wet_days),
    by = "cam_label"
  ) |>
  select(-any_of(c(
    "road_width_m", "road_surface_condition", "road_surface_type",
    "road_edge_condition", "veg_clearance_width_m", "vehicles_observed",
    "vehicle_types", "veh_tourist", "veh_ranger", "veh_local", "veh_other"
  )))

# Join daily weather onto individual detections
animals_site <- animals_site |>
  mutate(date = as_date(eventStart)) |>
  left_join(
    weather_daily |> select(date, temp_c, rh_pct, wind_speed_ms, soil_moisture),
    by = "date"
  )

message(sprintf("det_site:     %d × %d", nrow(det_site), ncol(det_site)))
message(sprintf("animals_site: %d × %d", nrow(animals_site), ncol(animals_site)))
message(sprintf("Missing temp_c in animals_site: %d", sum(is.na(animals_site$temp_c))))


# Sampling effort summary -------------------------------------------------
# Reported before hypothesis tests; RAI (detections / 100 trap-nights)
# standardises comparisons across unequal effort.

effort_summary <- cam_effort |>
  group_by(treatment) |>
  summarise(
    n_cameras = n(),
    total_tn  = sum(trap_nights),
    mean_tn   = round(mean(trap_nights), 1),
    sd_tn     = round(sd(trap_nights), 1),
    .groups   = "drop"
  )

cat("Sampling effort by treatment\n")
print(effort_summary)
cat(sprintf("\nTotal independent detections (30-min filter): %d\n", nrow(animals_site)))
cat(sprintf("Species detected: %d\n", n_distinct(animals_site$scientificName)))
