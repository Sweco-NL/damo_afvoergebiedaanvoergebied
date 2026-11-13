import logging
import warnings
import time
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pydantic import ConfigDict
from shapely.geometry import LineString, Point
from itertools import combinations

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.folium_map import generate_folium_map
from ..utils.general_functions import (
    calculate_angle_difference,
    calculate_angle_end,
    calculate_angle_reverse,
    calculate_angle_start,
    check_and_flip,
    line_to_vertices,
    split_waterways_by_endpoints,
)

def _format_duration(duration_seconds: float) -> str:
    minutes, seconds = divmod(int(duration_seconds), 60)
    return f"{minutes:02d}:{seconds:02d}"


# Suppress specific warnings
warnings.filterwarnings("ignore", message="Geometry column does not contain geometry")


class GeneratorDuikers(GeneratorBasis):
    """ "Module to guess (best-guess) the locations of culverts
    based on existing water network, other water bodies (c-watergangen),
    roads and level areas (peilgebied)."""

    generator: str = "generator_duikers"
    path: Path

    case_id: str = None
    case_name: str = None
    base_dir: Path = None
    waterschap: str = None

    hydroobject: gpd.GeoDataFrame = None
    overige_watergang: gpd.GeoDataFrame = None
    rws_water: gpd.GeoDataFrame = None
    outflow_nodes: gpd.GeoDataFrame = None

    bebouwing: gpd.GeoDataFrame = None
    kering: gpd.GeoDataFrame = None
    nwb: gpd.GeoDataFrame = None
    peilgebied: gpd.GeoDataFrame = None
    snelweg: gpd.GeoDataFrame = None
    spoorweg: gpd.GeoDataFrame = None

    read_results: bool = False
    write_results: bool = False

    snapping_distance: float = 0.05
    distance_vertices: float = 10.0
    max_culvert_length: float = None

    water_line_pnts: gpd.GeoDataFrame = None
    duplicates: gpd.GeoDataFrame = None

    # potential culverts
    potential_culverts_0: gpd.GeoDataFrame = None           # alle binnen 40m
    potential_culverts_1: gpd.GeoDataFrame = None           # filter kruizingen
    potential_culverts_2: gpd.GeoDataFrame = None           # scores
    potential_culverts_2a: gpd.GeoDataFrame = None          # fictieve lengtes
    potential_culverts_2b: gpd.GeoDataFrame = None          # filter op score en lengte
    potential_culverts_3: gpd.GeoDataFrame = None           # eerste resultaat
    potential_culverts_pre_filter: gpd.GeoDataFrame = None  # eerste resultaat
    potential_culverts_4: gpd.GeoDataFrame = None           # resultaat met nabewerking
    potential_culverts_5: gpd.GeoDataFrame = None           # resultaat met flips

    # hydroobject including splits by culverts
    hydroobject_processed_0: gpd.GeoDataFrame = None
    # overige_watergang including splits by culverts
    overige_watergang_processed_0: gpd.GeoDataFrame = None
    # overige_watergang including splits by culverts and post process
    overige_watergang_processed_1: gpd.GeoDataFrame = None
    # overige_watergang including splits by culverts
    overige_watergang_processed_2: gpd.GeoDataFrame = None

    # overige_watergang including outflow_nodes_hydro to hydroobjects and redirected
    overige_watergang_processed_3: gpd.GeoDataFrame = None
    overige_watergang_processed_3_nodes: gpd.GeoDataFrame = None

    # outflow points from overige watergangen to hydroobject
    outflow_nodes_overig: gpd.GeoDataFrame = None

    folium_map: folium.Map = None
    folium_html_path: Path = None


    def generate_vertices_along_waterlines(
        self,
        distance_vertices: float = 10.0,
        waterlines: list[str] = ["hydroobject", "overige_watergang"],
        read_results: bool = None,
        write_results: bool = None,
    ) -> gpd.GeoDataFrame:
        """Generate vertices (water_line_pnts) along waterlines every x meters

        Parameters
        ----------
        distance_vertices : float, optional
            distance between vertices, by default 10.0
        waterlines : list[str], optional
            list of attributes to be included, by default ["hydroobject", "overige_watergang"]
        read_results : bool, optional
            option (True/False) to read previous results from gpkg, by default None
        write_results : bool, optional
            option (True/False) to write results to case folder in gpkg, by default None

        Returns
        -------
        vertices (water_line_pnts): gpd.GeoDataFrame
            All points within a geodataframe
        """
        if isinstance(read_results, bool):
            self.read_results = read_results

        if self.read_results and self.water_line_pnts is not None:
            logging.info(
                f"   x {len(self.water_line_pnts)} vertices for {len(waterlines)} waterlines already generated"
            )

            # Identify duplicates among the "dangling" points
            start_end_points = self.water_line_pnts[
                self.water_line_pnts["line_type"].isin(
                    ["dangling_start", "dangling_end"]
                )
            ]

            # Find duplicate indices based on geometry
            duplicate_indices = start_end_points[
                start_end_points.duplicated(subset="geometry", keep=False)
            ].index
            self.duplicates = self.water_line_pnts.loc[duplicate_indices].copy()
            return self.water_line_pnts

        gdf_waterlines = None
        for waterline_name in waterlines:
            waterline = getattr(self, waterline_name)
            waterline["waterline_type"] = waterline_name

            if gdf_waterlines is None:
                gdf_waterlines = waterline.copy()
            else:
                gdf_waterlines = pd.concat([gdf_waterlines, waterline])

        logging.info(f"   x generate vertices for {len(gdf_waterlines)} waterlines")
        self.water_line_pnts = line_to_vertices(
            gdf_waterlines, distance=distance_vertices
        )
        self.water_line_pnts["unique_id"] = self.water_line_pnts.reset_index(
            drop=True
        ).index

        # Identify duplicates among the "dangling" points
        start_end_points = self.water_line_pnts[
            self.water_line_pnts["line_type"].isin(["dangling_start", "dangling_end"])
        ]

        # Find duplicate indices based on geometry
        duplicate_indices = start_end_points[
            start_end_points.duplicated(subset="geometry", keep=False)
        ].index
        self.duplicates = self.water_line_pnts.loc[duplicate_indices].copy()

        # Update all duplicated "dangling" points to "other"
        self.water_line_pnts.loc[duplicate_indices, "line_type"] = "other"

        if isinstance(write_results, bool):
            self.write_results = write_results

        return self.water_line_pnts, self.duplicates

    def find_potential_culvert_locations(
        self,
        water_line_pnts: gpd.GeoDataFrame = None,
        max_culvert_length: float = 40.0,
        read_results: bool = None,
        write_results: bool = None,
    ) -> gpd.GeoDataFrame:
        """Find potential culvert locations based on water_line_pnts.
        The connections between two points from different waterlines
        with a maximum distance of x m (max_culvert_length)

        Parameters
        ----------
        water_line_pnts : _type_, optional
            points every x m along the waterlines; includes a water_line_id, by default None
        max_culvert_length : int, optional
            maximum culvert length: in case of larger distance between points, connections are not made, by default 40
        read_results : bool, optional
            option (True/False) to read previous results from gpkg, by default None
        write_results : bool, optional
            option (True/False) to write results to case folder in gpkg, by default None

        Returns
        -------
        potential_culverts_0: gpd.GeoDataFrame
            Locations potential culverts (between different water_lines)
        """
        # check read_results and write_results
        if read_results is not None and isinstance(read_results, bool):
            self.read_results = read_results
        if write_results is not None and isinstance(write_results, bool):
            self.write_results = write_results

        if (
            not (
                isinstance(max_culvert_length, int)
                or isinstance(max_culvert_length, float)
            )
            or max_culvert_length <= 0.0
        ):
            raise ValueError(" x max_culvert_length is not a correct value")
        else:
            self.max_culvert_length = max_culvert_length

        # check if water_line_pnts are given. if not known, use function to generate.
        if water_line_pnts is None:
            if self.water_line_pnts is None:
                self.generate_water_line_pnts_along_waterlines()
        else:
            self.water_line_pnts = water_line_pnts

        # check if potential culverts already exists and should be read
        if self.read_results and self.potential_culverts_0 is not None:
            logging.info(
                f"   x {len(self.potential_culverts_0)} potential culverts already generated"
            )
            return self.potential_culverts_0

        logging.info("   x find potential culvert locations")

        # Filter for end points (only overige watergangen and not when connected)
        end_pnts = self.water_line_pnts[
            (self.water_line_pnts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (self.water_line_pnts["waterline_type"] == "overige_watergang")
        ]
        end_pnts = end_pnts.drop_duplicates(subset="geometry", keep=False)
        end_pnts = end_pnts.rename(
            columns={
                "code": "dangling_code",
                "unique_id": "dangling_id",
                "line_type": "line_type_start_end",
            }
        )
        end_pnts_orig_geometry = end_pnts.geometry

        # end points have buffer as geometry, perform spatial join with points
        logging.info("     - spatial join vertices")
        end_pnts["geometry"] = end_pnts.buffer(self.max_culvert_length)
        end_pnts = gpd.sjoin(
            end_pnts,
            self.water_line_pnts[["geometry", "unique_id", "code", "line_type"]],
            how="left",
            predicate="intersects",
        ).drop(columns="index_right")

        # remove connections with same hydroobject
        end_pnts = end_pnts[end_pnts["dangling_code"] != end_pnts["code"]]

        # make water_line_pnts id into unique_id again,
        # restore dangling geometry to points,
        # merge with original water_line_pnts,
        logging.info("     - add data to potential culverts")
        end_pnts["geometry"] = end_pnts_orig_geometry
        end_pnts = end_pnts.rename(columns={"geometry": "geometry2"})

        end_pnts = pd.merge(
            self.water_line_pnts,
            end_pnts[
                [
                    "unique_id",
                    "dangling_id",
                    "dangling_code",
                    "geometry2",
                    "line_type_start_end",
                ]
            ],
            on="unique_id",
            how="inner",
        )

        # create potential culverts (linestrings)
        potential_culverts_0 = end_pnts.copy()
        potential_culverts_0.dropna(subset=["dangling_id"], inplace=True)
        potential_culverts_0["geometry"] = potential_culverts_0.apply(
            lambda x: LineString([x["geometry2"], x["geometry"]]), axis=1
        )
        self.potential_culverts_0 = potential_culverts_0.drop(columns="geometry2")

        if self.write_results:
            logging.info(
                f"     - potential_culverts_0 not saved to save space"
            )
            # self.export_results_to_gpkg_or_nc(list_layers=[
            #     "potential_culverts_0",
            # ])

        logging.info(
            f"     - {len(self.potential_culverts_0)} potential culverts generated"
        )
        return self.potential_culverts_0


    def check_intersections_potential_culverts(self, read_results=None) -> gpd.GeoDataFrame:
        """Check intersections of culverts with other objects like roads, etc"""

        if isinstance(read_results, bool):
            self.read_results = read_results
        if self.read_results and self.potential_culverts_1 is not None:
            logging.info(f"   x {len(self.potential_culverts_1)} potential culverts already generated")
            return self.potential_culverts_1

        crossing_objects = [
            "hydroobject",
            "overige_watergang",
            "kering",
            "nwb",
            "peilgebied",
            "snelweg",
            "spoorweg",
        ]

        # Work directly on the original culverts
        culverts = self.potential_culverts_0.copy()

        logging.info("   x check intersections culverts with objects")
        for crossing_object in crossing_objects:
            logging.info(f"     - check for crossings with {crossing_object}")
            crossing_gdf = getattr(self, crossing_object)

            if crossing_gdf is None:
                culverts[f"crossings_{crossing_object}"] = np.nan
                continue
            
            # Perform spatial join
            joined = gpd.sjoin(
                culverts,
                crossing_gdf[["code", "geometry"]].rename(columns={"code": f"{crossing_object}_code"}),
                predicate="crosses",
            )

            # Aggregate to avoid duplicates
            joined = joined.groupby(joined.index)[f"{crossing_object}_code"].first()

            # Merge without creating duplicates
            culverts = culverts.merge(joined, how="left", left_index=True, right_index=True)
            
            if crossing_object in ["hydroobject", "overige_watergang"]:
                # Step 1: Create a copy with only the relevant columns
                culverts_copy = culverts[
                    [f"{crossing_object}_code", "dangling_code", "code", "unique_id"]
                ].copy()

                # Step 2: Ensure columns are of the same type (string) for comparison
                culverts_copy[f"{crossing_object}_code"] = culverts_copy[
                    f"{crossing_object}_code"
                ].astype("string")
                culverts_copy["dangling_code"] = culverts_copy["dangling_code"].astype(
                    "string"
                )
                culverts_copy["code"] = culverts_copy["code"].astype("string")

                # Step 3: Check where f"{crossing_object}_code" matches either dangling_code or code
                # Create boolean mask for each comparison
                mask_dangling = (
                    culverts_copy[f"{crossing_object}_code"]
                    == culverts_copy["dangling_code"]
                )
                mask_code = (
                    culverts_copy[f"{crossing_object}_code"] == culverts_copy["code"]
                )

                # Combine the masks: identify rows where f"{crossing_object}_code" matches either dangling_code or code
                mask_combined = mask_dangling | mask_code

                culverts.loc[
                    culverts_copy[mask_combined].index, f"{crossing_object}_code"
                ] = np.nan

            culverts[f"crossings_{crossing_object}"] = culverts[f"{crossing_object}_code"].notna()

            logging.info(f"     - {crossing_object} ({culverts[f'crossings_{crossing_object}'].sum()} crossings)")

        # Remove unwanted culverts directly
        culverts = culverts.loc[
            culverts[["crossings_hydroobject", "crossings_overige_watergang",
                    "crossings_kering", "crossings_snelweg", "crossings_spoorweg"]].any(axis=1) == False
        ]

        logging.info(f"     - {len(culverts)} potential culverts remaining")

        # Avoid unnecessary file writes if not needed
        self.potential_culverts_1 = culverts.copy()
        if self.write_results:
            logging.info(
                f"     - potential_culverts_1 not saved to save space"
            )
            # self.export_results_to_gpkg_or_nc(list_layers=[
            #     "potential_culverts_1",
            # ])
            
        return self.potential_culverts_1


    def assign_scores_to_potential_culverts(
        self, read_results=None
    ) -> gpd.GeoDataFrame:
        """Assign scores to all potential culverts based on the connected vertice
        and crossings with roads and peilgebied borders.

        Returns:
            gpd.GeoDataFrame: potential culverts with scores
        """
        if isinstance(read_results, bool):
            self.read_results = read_results
        if self.read_results and self.potential_culverts_2 is not None:
            logging.info(
                f"   x {len(self.potential_culverts_2)} potential culverts already generated"
            )
            return self.potential_culverts_2

        culverts = self.potential_culverts_1.copy()
        
        logging.info("   x assigning scores to potential culverts")
        # 1e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 1

        # 2e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 2

        # 3e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 3

        # 4e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 4

        # 5e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 5

        # 6e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 6

        # 7e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 7

        # 8e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 8

        # 9e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 9

        # 10e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 10

        # 11e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 11

        # 12e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 12

        # 13e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 13

        # 14e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "hydroobject")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 14

        # 15e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 15

        # 16e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["waterline_type"] == "overige_watergang")
            & (culverts["crossings_peilgebied"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 16

        logging.info("   x removing culverts without score")

        # Drop culverts that don't score 1-16
        culverts = culverts[culverts["score"].notna()]

        # Set copy and save data
        self.potential_culverts_2 = culverts.copy()
        if self.write_results:
            logging.info(
                f"     - potential_culverts_2 not saved to save space"
            )
            # self.export_results_to_gpkg_or_nc(list_layers=[
            #     "potential_culverts_2",
            # ])
        logging.info(
            f"     - {len(self.potential_culverts_2)} potential culverts remaining"
        )
        return self.potential_culverts_2

    def calculate_length_and_angles_potential_culverts_2a(
        self,
        factor_angle_on_length=2,
    ):
        """Calculate lengths, angles and fictive lengths.

        Input
        -----
        self.potential_culverts_2

        Output
        ------
        self.potential_culverts_2a
        """
        step_start_time = time.perf_counter()

        culverts = self.potential_culverts_2.copy()

        if culverts.empty:
            self.potential_culverts_2a = culverts.copy()
            return self.potential_culverts_2a

        culverts["length"] = culverts.geometry.length

        is_dangling_start = (
            culverts["line_type_start_end"] == "dangling_start"
        )

        culverts["angle_culvert"] = None

        culverts.loc[is_dangling_start, "angle_culvert"] = (
            culverts.loc[is_dangling_start, "geometry"]
            .map(calculate_angle_reverse)
        )

        culverts.loc[~is_dangling_start, "angle_culvert"] = (
            culverts.loc[~is_dangling_start, "geometry"]
            .map(calculate_angle_start)
        )

        overige_watergang = (
            self.overige_watergang
            .explode(index_parts=False)
            .copy()
        )

        overige_watergang["angle_start"] = (
            overige_watergang["geometry"]
            .map(calculate_angle_start)
        )

        overige_watergang["angle_end"] = (
            overige_watergang["geometry"]
            .map(calculate_angle_end)
        )

        if overige_watergang["code"].duplicated().any():
            duplicate_codes = (
                overige_watergang.loc[
                    overige_watergang["code"].duplicated(keep=False),
                    "code",
                ]
                .unique()
            )

            logging.warning(
                f"{len(duplicate_codes)} waterline codes have duplicate geometries"
            )

        overige_watergang_lookup = (
            overige_watergang
            .drop_duplicates(subset="code", keep="last")
        )

        angle_lookup = (
            overige_watergang_lookup
            .set_index("code")[["angle_start", "angle_end"]]
            .to_dict("index")
        )
        def get_waterline_angle(code, line_type):
            values = angle_lookup.get(code)

            if values is None:
                return None

            if line_type == "dangling_start":
                return values["angle_start"]

            if line_type == "dangling_end":
                return values["angle_end"]

            return None

        culverts["angle_waterline"] = [
            get_waterline_angle(code, line_type)
            for code, line_type in zip(
                culverts["dangling_code"],
                culverts["line_type_start_end"],
            )
        ]

        culverts["angle_difference"] = [
            calculate_angle_difference(
                angle_culvert,
                angle_waterline,
            )
            for angle_culvert, angle_waterline in zip(
                culverts["angle_culvert"],
                culverts["angle_waterline"],
            )
        ]

        angle_factor = (
            1
            + (factor_angle_on_length - 1)
            * culverts["angle_difference"].abs()
            / 90
        )

        culverts["fictive_length"] = (
            culverts["length"] * angle_factor
        )

        self.potential_culverts_2a = culverts.copy()

        logging.info(
            "   x Calculated angles and fictive lengths in "
            f"{_format_duration(time.perf_counter() - step_start_time)}"
        )

        return self.potential_culverts_2a

    def filter_potential_culverts_2b(self):
        """Filter potential culverts by length and score groups.

        Input
        -----
        self.potential_culverts_2a

        Output
        ------
        self.potential_culverts_2b
        """
        step_start_time = time.perf_counter()

        culverts = self.potential_culverts_2a.copy()

        def _log_score_counts(self, label, culverts):
            """Log the number of culverts per score."""

            if culverts.empty or "score" not in culverts.columns:
                logging.info(f"{label}: no scores available")
                return

            score_counts = (
                culverts["score"]
                .value_counts()
                .sort_index()
                .to_dict()
            )

            logging.info(
                f"{label}: score counts = {score_counts}"
            )

        _log_score_counts(self, "2a input", culverts)

        self.potential_culverts_pre_filter = culverts.copy()

        if culverts.empty:
            self.potential_culverts_2b = culverts.copy()
            return self.potential_culverts_2b

        # Remove candidates exceeding the maximum fictive length.
        culverts = culverts.loc[
            culverts["fictive_length"] <= self.max_culvert_length
        ].copy()

        _log_score_counts(self, "After max length filter", culverts)

        if culverts.empty:
            self.potential_culverts_2b = culverts.copy()
            return self.potential_culverts_2b

        # Keep the shortest complete row per dangling_id and score.
        culverts = (
            culverts
            .sort_values(
                by=[
                    "dangling_id",
                    "score",
                    "fictive_length",
                ],
                kind="stable",
            )
            .drop_duplicates(
                subset=["dangling_id", "score"],
                keep="first",
            )
            .copy()
        )

        _log_score_counts(self, "After duplicate removal", culverts)

        # Scores 1-4: group 1
        # Scores 5-8: group 2
        # Scores 9-12: group 3
        # Scores 13-16: group 4
        culverts["score_group"] = (
            (culverts["score"] - 1) // 4
        ) + 1

        # Keep the best available score group per dangling_id.
        best_score_group = (
            culverts
            .groupby("dangling_id")["score_group"]
            .transform("min")
        )
        _log_score_counts(self, "After score-group filter", culverts)

        culverts = culverts.loc[
            culverts["score_group"] == best_score_group
        ].drop(columns=["score_group"])

        # Remove scores 9 or higher if the same dangling_code
        # has a score 8 or lower.
        has_low_score = (
            culverts
            .groupby("dangling_code")["score"]
            .transform("min")
            <= 8
        )

        remove_high_score = (
            has_low_score
            & (culverts["score"] >= 9)
        )

        culverts = culverts.loc[
            ~remove_high_score
        ].copy()

        _log_score_counts(self, "After dangling-code filter", culverts)

        self.potential_culverts_2b = culverts

        logging.info(
            "   x Filtered potential culverts in "
            f"{_format_duration(time.perf_counter() - step_start_time)}"
        )
        logging.info(
            f"     - {len(culverts)} potential culverts remaining"
        )

        return self.potential_culverts_2b

    def select_best_potential_culverts_2c(
        self,
        offset_1_2=0.75,
        offset_12_34=0.6,
        offset_3_4=0.75,
    ):
        """Select the best score within each dangling_id.

        Input
        -----
        self.potential_culverts_2b

        Output
        ------
        self.potential_culverts_3
        """
        step_start_time = time.perf_counter()

        culverts = self.potential_culverts_2b.copy()

        if culverts.empty:
            self.potential_culverts_3 = culverts.copy()
            return self.potential_culverts_3

        def select_best_score(lengths_by_score):
            for base_score in (1, 5, 9, 13):
                score_a = base_score
                score_b = base_score + 1
                score_c = base_score + 2
                score_d = base_score + 3

                length_a = lengths_by_score.get(score_a)
                length_b = lengths_by_score.get(score_b)
                length_c = lengths_by_score.get(score_c)
                length_d = lengths_by_score.get(score_d)

                selected_score = None
                selected_length = None

                # Score A
                if length_a is not None:
                    selected_score = score_a
                    selected_length = length_a

                # Score B
                if (
                    length_b is not None
                    and (
                        selected_length is None
                        or length_b <= offset_1_2 * selected_length
                    )
                ):
                    selected_score = score_b
                    selected_length = length_b

                # Score C
                if selected_length is None and length_c is not None:
                    selected_score = score_c
                    selected_length = length_c

                if (
                    length_c is not None
                    and selected_score in (score_a, score_b)
                    and length_c <= offset_12_34 * selected_length
                ):
                    selected_score = score_c
                    selected_length = length_c

                # Score D
                if length_d is not None:
                    if (
                        selected_score == score_c
                        and length_d <= offset_3_4 * selected_length
                    ):
                        selected_score = score_d
                        selected_length = length_d

                    elif (
                        selected_score in (score_a, score_b)
                        and length_d <= offset_12_34 * selected_length
                        and (
                            length_c is None
                            or length_d <= offset_3_4 * length_c
                        )
                    ):
                        selected_score = score_d
                        selected_length = length_d

                # Preserve the original fallback behavior.
                if selected_length is None and length_d is not None:
                    selected_score = score_d

                if selected_score is not None:
                    return selected_score

            return None

        selected_scores = {}

        for dangling_id, group in culverts.groupby(
            "dangling_id",
            sort=False,
        ):
            lengths_by_score = dict(
                zip(
                    group["score"].astype(int),
                    group["fictive_length"],
                )
            )

            selected_scores[dangling_id] = select_best_score(
                lengths_by_score
            )

        culverts["selected_score"] = culverts["dangling_id"].map(
            selected_scores
        )

        culverts = culverts.loc[
            culverts["selected_score"] == culverts["score"]
        ].copy()

        self.potential_culverts_3 = culverts

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=["potential_culverts_3"]
            )

        logging.info(
            "   x Selected best potential culverts in "
            f"{_format_duration(time.perf_counter() - step_start_time)}"
        )
        logging.info(
            f"     - {len(culverts)} potential culverts remaining"
        )

        return self.potential_culverts_3

    def post_process_potential_culverts(self):
        """Postprocess culverts. removing double culverts, or culverts going to one point from the same waterline. 

        Returns
        -------
        self.potential_culverts_4: gpd.GeoDataFrame
            Geodataframe with culverts after postprocessing 
        """
        culverts = self.potential_culverts_3.copy()

        # Check for already connected watergangen and remove potential connections
        # Check and create a dictionary for duplicate groups with non-empty code lists
        duplicate_code_groups = (
            self.duplicates.groupby("geometry")["code"]
            .apply(
                lambda codes: list(codes) if len(codes) > 1 else []
            )  # Only create list if more than one code
            .to_dict()
        )

        # Chose shortest culvert when two culverts go to the same point from both dangling_points of a watergang
        culverts["code_pair"] = culverts.apply(
            lambda row: frozenset([row["dangling_code"], row["code"]]), axis=1
        )

        culverts = culverts.loc[culverts.groupby("code_pair")["fictive_length"].idxmin()]

        # # Define the filter function
        # def should_remove(row):
        #     dangling_code = row["dangling_code"]
        #     code = row["code"]

        #     # Check if both dangling_code and code appear in any of the lists in duplicate_code_groups
        #     for codes_list in duplicate_code_groups.values():
        #         if dangling_code in codes_list and code in codes_list:
        #             return True  # Mark for removal if both are found in the same list

        #     return False  # Otherwise, don't remove

        # # Filter culverts using the modified function
        # culverts = culverts[~culverts.apply(should_remove, axis=1)].copy()
        
        forbidden_pairs = set()

        for codes in duplicate_code_groups.values():
            for code_a, code_b in combinations(codes, 2):
                forbidden_pairs.add(
                    frozenset((code_a, code_b))
                )

        culverts["code_pair"] = [
            frozenset((dangling_code, code))
            for dangling_code, code in zip(
                culverts["dangling_code"],
                culverts["code"],
            )
        ]

        culverts = culverts.loc[
            ~culverts["code_pair"].isin(forbidden_pairs)
        ].copy()

        culverts = culverts.drop(columns="code_pair")

        # check if overige watergang already has another connection and choose shortest:
        duplicate_culverts = culverts.loc[culverts["waterline_type"]=="hydroobject"]
        duplicate_culverts = duplicate_culverts.sort_values(
            [
                "dangling_code", 
                "selected_score", 
                "code", 
                "fictive_length"
            ], ascending=[True, True, True, True]
        )
        duplicate_culverts = duplicate_culverts[
            duplicate_culverts.duplicated(subset=["dangling_code"], keep=False)
        ]
        duplicate_culverts_to_drop = duplicate_culverts.loc[
            duplicate_culverts.duplicated(subset=['dangling_code'], keep="first")
        ]
        logging.info(f"   - found {len(duplicate_culverts)} duplicate_culverts to drop")
        culverts = culverts[~culverts.index.isin(duplicate_culverts_to_drop.index)]

        self.potential_culverts_4 = culverts.copy()
        logging.info(
            f"     - {len(self.potential_culverts_4)} potential culverts remaining"
        )
        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "potential_culverts_4",
            ])
        return self.potential_culverts_4



    def splits_hydroobject_by_endpoints_of_culverts_and_combine(self):
        """Splits hydroobjects and overige_watergangt the location where culverts are connected. 
        This is done to create a complete and connected network.
        The function also generates the outflow point of the overige_watergang in hydroobjects.

        Returns
        -------
        self.overige_watergang_processed_0: gpd.GeoDataFrame
            overige_watergang after the splits from the culverts
        self.hydroobject_processed_0: gpd.GeoDataFrame
            hydroobject after the splits from the culverts
        self.outflow_nodes_overig: gpd.GeoDataFrame
            Outflow points of overige_watergang in hydroobjects
           """
        step_start = time.perf_counter()
        
        potential = self.potential_culverts_4

        # Create the common mask once
        other_mask = potential["line_type"].eq("other")

        other_culverts_hydro = potential.loc[
            other_mask
            & potential["waterline_type"].eq("hydroobject")
        ].copy()

        other_culverts_other = potential.loc[
            other_mask
            & potential["waterline_type"].eq("overige_watergang")
        ].copy()

        logging.info(
            "     - selected %d hydroobject culverts and %d overige-watergang culverts "
            "in %.2f seconds",
            len(other_culverts_hydro),
            len(other_culverts_other),
            time.perf_counter() - step_start,
        )

        # These calls are probably the main bottleneck
        step_start = time.perf_counter()

        self.hydroobject_processed_0 = split_waterways_by_endpoints(
            self.hydroobject,
            other_culverts_hydro,
        )

        logging.info(
            "     - hydroobject gesplit in %.2f seconds",
            time.perf_counter() - step_start,
        )

        step_start = time.perf_counter()

        self.overige_watergang_processed_0 = split_waterways_by_endpoints(
            self.overige_watergang,
            other_culverts_other,
        )

        logging.info(
            "     - overige watergangen gesplit in %.2f seconds",
            time.perf_counter() - step_start,
        )

        if self.write_results:
            step_start = time.perf_counter()

            self.export_results_to_gpkg_or_nc(
                list_layers=[
                    "hydroobject_processed_0",
                    "overige_watergang_processed_0",
                ]
            )

            logging.info(
                "     - split results written in %.2f seconds",
                time.perf_counter() - step_start,
            )

        # Select hydroobject-connected culverts.
        # .copy() is important because a new column is created below.
        culverts_hydro = potential.loc[
            potential["waterline_type"].eq("hydroobject")
        ].copy()

        # Keep only valid LineStrings
        valid_lines = (
            culverts_hydro.geometry.notna()
            & culverts_hydro.geometry.geom_type.eq("LineString")
        )

        culverts_hydro = culverts_hydro.loc[valid_lines].copy()

        # A list comprehension is generally faster than apply(axis=1)
        endpoint_geometries = [
            Point(geometry.coords[-1])
            for geometry in culverts_hydro.geometry
        ]

        # Build the GeoDataFrame directly; no temporary 'end_point' column
        self.outflow_nodes_overig = gpd.GeoDataFrame(
            culverts_hydro[
                ["unique_id", "dangling_code"]
            ].copy(),
            geometry=endpoint_geometries,
            crs=culverts_hydro.crs,
        )

        logging.info(
            "     - created %d outflow nodes",
            len(self.outflow_nodes_overig),
        )

        if self.write_results:
            step_start = time.perf_counter()

            self.export_results_to_gpkg_or_nc(
                list_layers=["outflow_nodes_overig"]
            )

            logging.info(
                "     - outflow nodes written in %.2f seconds",
                time.perf_counter() - step_start,
            )

        return (
            self.overige_watergang_processed_0,
            self.hydroobject_processed_0,
            self.outflow_nodes_overig,
        )


    # def check_culverts_direction(self):
    #     """Culverts can be generated in the wrong direction compared to the waterlines, if this is the case, they are flipped.

    #     Returns
    #     -------
    #     self.potential_culverts_5: gpd.GeoDataFrame
    #         Geodataframe containing corrected culverts
    #     """
    #     culvert = self.potential_culverts_4.copy()
    #     lines = pd.concat(
    #         [self.hydroobject_processed_0, self.overige_watergang_processed_0],
    #         ignore_index=True,
    #     )

    #     # Extract starting and ending points from gdf2
    #     gdf2_start_points = lines.geometry.apply(lambda geom: geom.coords[0])
    #     gdf2_end_points = lines.geometry.apply(lambda geom: geom.coords[-1])

    #     # Initialize the 'flipped' column if it doesn't exist
    #     if "flipped" not in culvert.columns:
    #         culvert["flipped"] = 0

    #     # Apply the check_and_flip function to each line in gdf1
    #     results = culvert.geometry.apply(
    #         lambda line: check_and_flip(
    #             line, gdf2_start_points.tolist(), gdf2_end_points.tolist()
    #         )
    #     )

    #     # Update the geometry and flipped columns
    #     culvert["geometry"] = results.apply(lambda x: x[0])
    #     culvert["flipped"] += results.apply(lambda x: 1 if x[1] else 0)
    #     logging.info("     - culvert direction checked")

    #     self.potential_culverts_5 = culvert.copy()

    #     if self.write_results:
    #         self.export_results_to_gpkg_or_nc(list_layers=[
    #             "potential_culverts_5",
    #         ])
    #     return self.potential_culverts_5


    def add_culvert_to_overige_watergang_as_seperate_object(self):
        """Add culverts as separate objects to overige_watergang_processed_0."""

        lines = self.overige_watergang_processed_0.copy()
        culvert = self.potential_culverts_4.copy()

        # Keep only valid line geometries
        culvert = culvert.loc[
            culvert.geometry.notna()
            & culvert.geometry.geom_type.eq("LineString")
        ].copy()

        # Set required attributes
        culvert["code"] = culvert["dangling_code"].astype(str) + "_culvert"
        culvert = culvert.rename(columns={"length": "geometry_len"})
        culvert = culvert[["code", "geometry_len", "geometry"]]

        # Concatenate and fillna
        combined = pd.concat([lines, culvert], ignore_index=True)
        combined["status"] = combined["status"].fillna("unchanged")
        combined["waterline_type"] = combined["waterline_type"].fillna("culvert")

        # Recreate GeoDataFrame
        self.overige_watergang_processed_1 = gpd.GeoDataFrame(
            combined,
            geometry="geometry",
            crs=lines.crs,
        )

        logging.info("     - culverts added to watergangen as separate objects")

        if self.write_results:
            self.export_results_to_gpkg_or_nc(
                list_layers=["overige_watergang_processed_1"]
            )

        return self.overige_watergang_processed_1

    def combine_culvert_with_line(self):
        """Culverts are combined with the one of the waterlines they are connected to, to form a single linestring. After this funtion the 'check_culverts_direction' function be run again.

        Returns
        -------
        self.overige_watergang_processed_1: gpd.GeoDataFrame
            Geodataframe containing processed overige_watergang where the waterline is combined with culverts
        """
        culvert = self.potential_culverts_5.copy()
        culvert_dict = (
            culvert.groupby("dangling_code")["geometry"].apply(list).to_dict()
        )
        lines = self.overige_watergang_processed_0.copy()

        def get_base_code(code):
            return code.split("-")[0]

        # Create base code column for merging
        lines["base_code"] = lines["code"].apply(get_base_code)

        # Define the combine_lines function
        def combine_lines(line, culvert_dict):
            base_code = line["base_code"]
            line_geom = line["geometry"]

            if base_code in culvert_dict:
                for culv in culvert_dict[base_code]:
                    if (
                        culv.coords[-1] == line_geom.coords[0]
                    ):  # Culvert end to line start
                        line_geom = LineString(
                            list(culv.coords) + list(line_geom.coords)
                        )
                    elif (
                        culv.coords[0] == line_geom.coords[-1]
                    ):  # Line end to culvert start
                        line_geom = LineString(
                            list(line_geom.coords) + list(culv.coords)
                        )

            return line_geom

        lines["geometry"] = lines.apply(
            lambda line: combine_lines(line, culvert_dict), axis=1
        )

        logging.info("     - culverts combined with watergangen")

        self.overige_watergang_processed_1 = lines.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "overige_watergang_processed_1",
            ])
        return self.overige_watergang_processed_1


    def splits_hydroobject_by_endpoind_of_culverts_and_combine_2(
        self, write_results=True
    ):
        """After a second check of the culvert direction, the hydroobjects and overige_watergang must be split again for some instances and the remaining flipped culverts are combined with the overige_watergang.

        Returns
        -------
        self.overige_watergang_processed_1: gpd.GeoDataFrame
            Geodataframe containing processed overige_waterganghere the waterline is combined with culverts
        """
        # split overige watergangen opnieuw
        overige_watergang = self.overige_watergang_processed_1.copy()
        overige_watergang = split_waterways_by_endpoints(
            overige_watergang, overige_watergang
        )
        logging.info("     - overige watergangen weer gesplit")

        self.overige_watergang_processed_2 = overige_watergang.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "overige_watergang_processed_2",
            ])
        return self.overige_watergang_processed_2


    def get_shortest_path_from_overige_watergang_to_hydroobjects(
        self, write_results=False
    ):
        """Calculate the shortest path to a hydroobject for the overige_watergang. This is done to identify overige_watergang with an incorrect direction, these are then flipped.

        Parameters
        ----------
        write_results : bool, optional
            Option (True/False) to write results to case folder in gpkg, by default False

        Returns
        -------
        outflow_nodes_hydro: gpd.GeoDataFrame
            Outflow nodes of the overige_watergang in hydroobject
        overige_watergang_nodes: gpd.GeoDataFrame
            All outflow nodes of overige_watergang
        edges: gpd.GeoDataFrame
            All overige_watergang as edges with information about previous and subsequent overige_watergang
        edges_cleaned: gpd.GeoDataFrame
            overige_watergang cleaned through removal of unconnected waterlines and flipped waterlines when required
        """
        logging.info(f"   x redirect 'overige watergangen' based on shortest path")
        # create networkx graph
        logging.info(f"     - create sub networks")
        nodes, edges, _ = create_graph_from_edges(
            self.overige_watergang_processed_1, directed=True
        )
        _, _, G = create_graph_from_edges(
            self.overige_watergang_processed_1, directed=False
        )

        # get outflow points and add nodes information
        outflow_nodes_hydro = (
            self.outflow_nodes_overig[
                ["unique_id", "dangling_code", "geometry"]
            ]
            .sjoin(nodes, how="inner")
            .reset_index(drop=True)
            .drop(columns=["index_right"], errors="ignore")
        )
        self.outflow_nodes_overig = outflow_nodes_hydro.copy()
        
        if outflow_nodes_hydro.empty:
            overige_watergang_nodes = None
            edges_cleaned = None
        else:
            # get shortest path including length from all nodes to outflow points
            logging.info(f"     - find shortest path")
            len_outflow_nodes_hydro, matrix = nx.multi_source_dijkstra(
                G, [int(n) for n in outflow_nodes_hydro["nodeID"].values], weight="geometry_len"
            )
            node_to_outflow_nodes_hydro = pd.DataFrame(
                {
                    "nodeID": matrix.keys(),
                    "outflow_nodes_hydro": [v[0] for v in matrix.values()],
                    "outflow_nodes_hydro_dist": [len_outflow_nodes_hydro[node] for node in matrix.keys()],
                }
            )

            # Merge nodes with node_to_outflow_nodes_hydro
            overige_watergang_nodes = nodes.merge(
                node_to_outflow_nodes_hydro, how="outer", left_on="nodeID", right_on="nodeID"
            )
            overige_watergang_nodes["outflow_nodes_hydro"] = (
                overige_watergang_nodes["outflow_nodes_hydro"].fillna(-999).astype(int)
            )

            # To get length to outflow point at start and end: merge edges with nodes, first on node_start, then on node_end
            edges = (
                edges.merge(
                    overige_watergang_nodes[
                        ["nodeID", "outflow_nodes_hydro", "outflow_nodes_hydro_dist"]
                    ].rename(
                        columns={
                            "outflow_nodes_hydro": "outflow_start",
                            "outflow_nodes_hydro_dist": "outflow_start_dist",
                        }
                    ),
                    how="left",
                    left_on="node_start",
                    right_on="nodeID",
                )
                .merge(
                    overige_watergang_nodes[
                        ["nodeID", "outflow_nodes_hydro", "outflow_nodes_hydro_dist"]
                    ].rename(
                        columns={
                            "outflow_nodes_hydro": "outflow_end",
                            "outflow_nodes_hydro_dist": "outflow_end_dist",
                        }
                    ),
                    how="left",
                    left_on="node_end",
                    right_on="nodeID",
                )
                .drop(columns=["nodeID_x", "nodeID_y"])
            )

            # select shortest paths for each each
            logging.info(f"     - select direction with shortest path")

            def select_shortest_direction(edge):
                if edge.outflow_end_dist <= edge.outflow_start_dist:
                    edge.outflow_nodes_hydro = edge.outflow_end
                    edge.outflow_nodes_hydro_dist = edge.outflow_end_dist
                else:
                    edge.outflow_nodes_hydro = edge.outflow_start
                    edge.outflow_nodes_hydro_dist = edge.outflow_start_dist
                    edge.reversed_direction = True
                return edge

            edges["outflow_nodes_hydro"] = -999
            edges["outflow_nodes_hydro_dist"] = -999.0
            edges["reversed_direction"] = False
            edges["outflow_end_dist"] = edges["outflow_end_dist"].fillna(99999.9)
            edges["outflow_start_dist"] = edges["outflow_start_dist"].fillna(99999.9)
            edges = edges.apply(lambda x: select_shortest_direction(x), axis=1)

            # clean edges by removing unconnected edges and reversing direction if required
            edges_cleaned = edges[edges["outflow_nodes_hydro"] != -999]
            edges_cleaned.loc[edges_cleaned["reversed_direction"], "geometry"] = (
                edges_cleaned.loc[edges_cleaned["reversed_direction"], "geometry"].reverse()
            )

            self.overige_watergang_processed_3_nodes = overige_watergang_nodes.copy()
            self.overige_watergang_processed_3 = edges_cleaned.copy()

        if self.write_results:
            self.export_results_to_gpkg_or_nc(list_layers=[
                "outflow_nodes_overig",
                "overige_watergang_processed_3_nodes",
                "overige_watergang_processed_3",
            ])
        return outflow_nodes_hydro, overige_watergang_nodes, edges, edges_cleaned


    def generate_folium_map(
        self, 
        html_file_name:str="", 
        show_other_waterways_culverts: bool = True,
        **kwargs
    ):
        if html_file_name == '':
            html_file_name = self.case_id + "_culvert_locations"
        self.folium_map = generate_folium_map(
            self, 
            html_file_name=html_file_name, 
            show_other_waterways_culverts=show_other_waterways_culverts,
            **kwargs
        )
        return self.folium_map
