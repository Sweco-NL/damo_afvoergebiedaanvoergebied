import geopandas as gpd
import pandas as pd
from pathlib import Path


from .general_functions import (
    remove_z_dims,
    snap_unconnected_endpoints_to_endpoint_or_line,
    split_waterways_by_endpoints,
)


def preprocess_hydroobjecten(hydroobjecten):
    # Setup hydroobjecten correctly
    hydroobjecten = remove_z_dims(hydroobjecten)

    # Snap hydroobjecten
    hydroobjecten = snap_unconnected_endpoints_to_endpoint_or_line(
        hydroobjecten, snapping_distance=0.05
    )

    # Split_hydroobjecten
    hydroobjecten = split_waterways_by_endpoints(hydroobjecten, hydroobjecten)

    return hydroobjecten
