import geopandas as gpd
import pandas as pd
from pathlib import Path


from .general_functions import (
    remove_z_dims,
    snap_unconnected_endpoints_to_endpoint_or_line,
    split_waterways_by_endpoints,
    check_duplicate_codes,
)


def preprocess_hydroobjecten(hydroobjecten, snapping_distance=0.05):
    # explode
    hydroobjecten = hydroobjecten.explode()

    # Setup hydroobjecten correctly
    hydroobjecten = remove_z_dims(hydroobjecten)

    # check duplicates
    hydroobjecten = check_duplicate_codes(hydroobjecten, "code")

    # check for invalid or duplicate geometries (linestrings forming a ring)
    hydroobjecten_old = hydroobjecten.copy()
    hydroobjecten = hydroobjecten[~hydroobjecten['geometry'].apply(lambda geom: geom.is_closed)]
    hydroobjecten = hydroobjecten.loc[~hydroobjecten["geometry"].duplicated(keep="first")]
    hydroobjecten_removed = hydroobjecten_old.loc[~hydroobjecten_old.index.isin(hydroobjecten.index)]

    # Snap hydroobjecten
    hydroobjecten = snap_unconnected_endpoints_to_endpoint_or_line(
        hydroobjecten, snapping_distance=snapping_distance
    )

    hydroobjecten_snapped = hydroobjecten.copy()

    # Split_hydroobjecten
    hydroobjecten = split_waterways_by_endpoints(hydroobjecten, hydroobjecten)

    return hydroobjecten, hydroobjecten_snapped, hydroobjecten_removed
