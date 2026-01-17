import geopandas as gpd
import pandas as pd
from pathlib import Path


from .general_functions import (
    remove_z_dims,
    snap_unconnected_endpoints_to_endpoint_or_line,
    split_waterways_by_endpoints,
    check_duplicate_codes,
)


def preprocess_hydroobject(hydroobject, snapping_distance=0.05):
    # explode
    hydroobject = hydroobject.explode()

    # Setup hydroobject correctly
    hydroobject = remove_z_dims(hydroobject)

    # check duplicates
    hydroobject = check_duplicate_codes(hydroobject, "code")

    # Snap hydroobject
    hydroobject = snap_unconnected_endpoints_to_endpoint_or_line(
        hydroobject, snapping_distance=snapping_distance
    )

    hydroobject_snapped = hydroobject.copy()

    # Split_hydroobject
    hydroobject = split_waterways_by_endpoints(hydroobject, hydroobject)

    return hydroobject, hydroobject_snapped
