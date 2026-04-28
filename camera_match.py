from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsGeometry
from qgis.core import QgsField, QgsFields
from PyQt5.QtCore import QVariant

# -----------------------------------------------------------
# Config
# -----------------------------------------------------------

ELEV_TOLERANCE       = 50
WATER_TOLERANCE      = 2000
MIN_ROAD_DIST        = 400
MIN_CONTROL_SEP      = 400
MIN_CANDIDATE_ELEV   = 2225
MAX_CANDIDATE_ELEV   = 2525
TOP_N                = 3

WEIGHT_ELEV  = 8
WEIGHT_WATER = 2

LAYER_FEATURES   = 'possible_feature_total'
LAYER_CANDIDATES = 'control_candidates_500m'
LAYER_ROAD       = 'ruhija_road_in_park'

FIELD_ELEV  = 'elevation1'
FIELD_WATER = 'water_dist'


# -----------------------------------------------------------
# Helpers
# -----------------------------------------------------------

def get_layer(name):
    layers = QgsProject.instance().mapLayersByName(name)
    if not layers:
        raise RuntimeError(f"Could not find layer: '{name}'")
    return layers[0]


def distance_to_road(geom, road_layer):
    min_dist = float('inf')
    for road_feat in road_layer.getFeatures():
        dist = geom.distance(road_feat.geometry())
        if dist < min_dist:
            min_dist = dist
    return min_dist


def is_too_close_to_placed(geom, placed_controls):
    for placed in placed_controls:
        if geom.distance(placed) < MIN_CONTROL_SEP:
            return True
    return False


def calc_score(elev_diff, water_diff):
    return (elev_diff * WEIGHT_ELEV) + (water_diff * WEIGHT_WATER)


def build_output_layer(crs_id):
    fields = QgsFields()
    fields.append(QgsField('pair_id',       QVariant.Int))
    fields.append(QgsField('option_rank',   QVariant.Int))
    fields.append(QgsField('feature_elev',  QVariant.Double))
    fields.append(QgsField('control_elev',  QVariant.Double))
    fields.append(QgsField('elev_diff',     QVariant.Double))
    fields.append(QgsField('feature_water', QVariant.Double))
    fields.append(QgsField('control_water', QVariant.Double))
    fields.append(QgsField('water_diff',    QVariant.Double))
    fields.append(QgsField('road_dist',     QVariant.Double))
    fields.append(QgsField('pair_dist',     QVariant.Double))
    fields.append(QgsField('match_score',   QVariant.Double))
    fields.append(QgsField('recommended',   QVariant.String))

    out = QgsVectorLayer(f'Point?crs={crs_id}', 'Control_Cameras_Interior_Forest', 'memory')
    out.dataProvider().addAttributes(fields)
    out.updateFields()
    return out


# -----------------------------------------------------------
# Load layers
# -----------------------------------------------------------

feat_layer = get_layer(LAYER_FEATURES)
cand_layer = get_layer(LAYER_CANDIDATES)
road_layer = get_layer(LAYER_ROAD)

print(f"Features loaded: {feat_layer.featureCount()}")
print(f"Candidates loaded: {cand_layer.featureCount()}")

out_layer = build_output_layer(feat_layer.crs().authid())
print("Output layer created")

# -----------------------------------------------------------
# Matching
# -----------------------------------------------------------

print("\nFinding matches...")
print(f"  Elevation tolerance: ±{ELEV_TOLERANCE}m")
print(f"  Min road distance: {MIN_ROAD_DIST}m")
print(f"  Min control separation: {MIN_CONTROL_SEP}m")
print()

all_features    = []
placed_controls = []
no_match_ids    = []
few_match_ids   = []

for feat in feat_layer.getFeatures():
    feat_id    = feat.id()
    feat_geom  = feat.geometry()
    feat_elev  = feat[FIELD_ELEV]
    feat_water = feat[FIELD_WATER]

    valid_candidates = []

    for cand in cand_layer.getFeatures():
        cand_elev  = cand[FIELD_ELEV]
        cand_water = cand[FIELD_WATER]
        cand_geom  = cand.geometry()

        # Pre-filter by elevation range
        if not (MIN_CANDIDATE_ELEV <= cand_elev <= MAX_CANDIDATE_ELEV):
            continue

        # Must be far enough from road
        road_dist = distance_to_road(cand_geom, road_layer)
        if road_dist < MIN_ROAD_DIST:
            continue

        # Elevation and water must be similar to the feature camera
        elev_diff  = abs(feat_elev - cand_elev)
        water_diff = abs(feat_water - cand_water)

        if elev_diff > ELEV_TOLERANCE or water_diff > WATER_TOLERANCE:
            continue

        # Don't place controls too close to each other
        if is_too_close_to_placed(cand_geom, placed_controls):
            continue

        pair_dist   = feat_geom.distance(cand_geom)
        match_score = calc_score(elev_diff, water_diff)

        valid_candidates.append({
            'geometry':   cand_geom,
            'elev':       cand_elev,
            'water':      cand_water,
            'elev_diff':  elev_diff,
            'water_diff': water_diff,
            'road_dist':  road_dist,
            'pair_dist':  pair_dist,
            'score':      match_score,
        })

    # Sort by score (lower is better) and take top N
    valid_candidates.sort(key=lambda x: x['score'])
    top_matches = valid_candidates[:TOP_N]

    # Reserve the primary match location so other pairs don't overlap it
    if top_matches:
        placed_controls.append(top_matches[0]['geometry'])

    # Build output features
    for rank, match in enumerate(top_matches, start=1):
        label = 'PRIMARY' if rank == 1 else f'BACKUP_{rank - 1}'

        out_feat = QgsFeature()
        out_feat.setGeometry(match['geometry'])
        out_feat.setAttributes([
            feat_id,
            rank,
            feat_elev,
            match['elev'],
            match['elev_diff'],
            feat_water,
            match['water'],
            match['water_diff'],
            match['road_dist'],
            match['pair_dist'],
            match['score'],
            label,
        ])
        all_features.append(out_feat)


# -----------------------------------------------------------
# Save and add to map
# -----------------------------------------------------------

out_layer.dataProvider().addFeatures(all_features)
out_layer.updateExtents()
QgsProject.instance().addMapLayer(out_layer)

# -----------------------------------------------------------
# Summary
# -----------------------------------------------------------

total      = feat_layer.featureCount()
full_count = total - len(no_match_ids) - len(few_match_ids)

print(f"\nDone.")
print(f"  Total cameras:         {total}")
print(f"  Total control options: {len(all_features)}")
print(f"  Full ({TOP_N} options): {full_count}")
print(f"  Partial:               {len(few_match_ids)}")
print(f"  No matches:            {len(no_match_ids)}")

if no_match_ids:
    print(f"\n  Cameras with no matches: {no_match_ids}")
