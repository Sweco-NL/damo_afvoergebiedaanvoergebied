import logging
import warnings
import webbrowser
from pathlib import Path

import folium
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from pydantic import ConfigDict
from shapely.geometry import LineString, Point

from ..generator_basis import GeneratorBasis
from ..utils.create_graph import create_graph_from_edges
from ..utils.folium_utils import (
    add_basemaps_to_folium_map,
    add_categorized_lines_to_map,
    add_labels_to_points_lines_polygons,
)
from ..utils.general_functions import (
    calculate_angle_difference,
    calculate_angle_end,
    calculate_angle_reverse,
    calculate_angle_start,
    check_and_flip,
    line_to_vertices,
    split_waterways_by_endpoints,
)
from ..utils.preprocess import preprocess_hydroobjecten

# Suppress specific warnings
warnings.filterwarnings("ignore", message="Geometry column does not contain geometry")


class GeneratorCulvertLocations(GeneratorBasis):
    """ "Module to guess (best-guess) the locations of culverts
    based on existing water network, other water bodies (c-watergangen),
    roads and level areas (peilgebieden)."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None

    hydroobjecten: gpd.GeoDataFrame = None
    overige_watergangen: gpd.GeoDataFrame = None
    bebouwing: gpd.GeoDataFrame = None
    keringen: gpd.GeoDataFrame = None
    nwb: gpd.GeoDataFrame = None
    peilgebieden: gpd.GeoDataFrame = None
    snelwegen: gpd.GeoDataFrame = None
    spoorwegen: gpd.GeoDataFrame = None

    preprocess_hydroobjecten: bool = True
    read_results: bool = False
    write_results: bool = False

    max_culvert_length: float = None

    water_line_pnts: gpd.GeoDataFrame = None
    duplicates: gpd.GeoDataFrame = None
    potential_culverts_0: gpd.GeoDataFrame = None  # alle binnen 40m
    potential_culverts_1: gpd.GeoDataFrame = None  # filter kruizingen
    potential_culverts_2: gpd.GeoDataFrame = None  # scores
    potential_culverts_3: gpd.GeoDataFrame = None  # eerste resultaat
    potential_culverts_pre_filter: gpd.GeoDataFrame = None  # eerste resultaat
    potential_culverts_4: gpd.GeoDataFrame = None  # resultaat met nabewerking
    potential_culverts_5: gpd.GeoDataFrame = None  # resultaat met flips

    # hydroobjecten including splits by culverts
    hydroobjecten_processed_0: gpd.GeoDataFrame = None
    # overige_watergangen including splits by culverts
    overige_watergangen_processed_0: gpd.GeoDataFrame = None
    # overige_watergangen including splits by culverts and post process
    overige_watergangen_processed_1: gpd.GeoDataFrame = None
    # overige_watergangen including splits by culverts
    overige_watergangen_processed_2: gpd.GeoDataFrame = None

    # overige_watergangen including outflow_node to hydroobjects and redirected
    overige_watergangen_processed_3: gpd.GeoDataFrame = None
    overige_watergangen_processed_3_nodes: gpd.GeoDataFrame = None

    # outflow points from overige watergangen to hydroobjecten
    outflow_nodes_overige_watergangen: gpd.GeoDataFrame = None

    folium_map: folium.Map = None
    folium_html_path: Path = None

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.generate_or_use_preprocessed_hydroobjecten()

    def use_processed_hydroobjecten(self):
        logging.info("     - culvert generator will generate processed hydroobjecten")

    def generate_or_use_preprocessed_hydroobjecten(
        self, preprocessed_file="preprocessed"
    ):
        hydroobjecten_preprocessed_file = None
        files_in_dir = self.dir_results.glob("**/*")
        for f in files_in_dir:
            if f"hydroobjecten_{preprocessed_file}" == f.stem:
                hydroobjecten_preprocessed_file = f

        if hydroobjecten_preprocessed_file is None:
            logging.info(
                "     - hydroobjecten_snap_split.gpkg not in directory, preprocessing hydroobjecten"
            )
            self.hydroobjecten, hydroobjecten_snapped = preprocess_hydroobjecten(
                self.hydroobjecten
            )
            hydroobjecten_snapped.to_file(
                Path(self.dir_results, "hydroobjecten_snapped.gpkg")
            )
            self.hydroobjecten.to_file(
                Path(self.dir_results, "hydroobjecten_preprocessed.gpkg")
            )
        else:
            logging.info("     - get dataset preprocessed hydroobjecten")
            self.hydroobjecten = gpd.read_file(hydroobjecten_preprocessed_file)

    def generate_vertices_along_waterlines(
        self,
        distance_vertices: float = 10.0,
        waterlines: list[str] = ["hydroobjecten", "overige_watergangen"],
        read_results: bool = None,
        write_results: bool = None,
    ) -> gpd.GeoDataFrame:
        """Generate vertices (water_line_pnts) along waterlines every x meters

        Parameters
        ----------
        distance_vertices : float, optional
            distacne between vertices, by default 10.0
        waterlines : list[str], optional
            list of attributes to be included, by default ["hydroobjecten", "overige_watergangen"]
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
            waterline["WaterLineType"] = waterline_name
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

        if self.write_results:
            self.water_line_pnts.to_file(Path(self.dir_results, "water_line_pnts.gpkg"))
        return self.water_line_pnts, self.duplicates

    def find_potential_culvert_locations(
        self,
        water_line_pnts: gpd.GeoDataFrame = None,
        max_culvert_length: float = 40.0,
        read_results: bool = False,
        write_results: bool = False,
    ) -> gpd.GeoDataFrame:
        """Find potential culvert locations based on water_line_pnts.
        THe connections between two points from different waterlines
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
        if isinstance(read_results, bool):
            self.read_results = read_results
        if isinstance(write_results, bool):
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
            & (self.water_line_pnts["WaterLineType"] == "overige_watergangen")
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

        if write_results:
            self.potential_culverts_0.to_file(
                Path(self.dir_results, "potential_culverts_0.gpkg")
            )

        logging.info(
            f"     - {len(self.potential_culverts_0)} potential culverts generated"
        )
        return self.potential_culverts_0

    def check_intersections_potential_culverts(
        self,
        read_results=None,
    ) -> gpd.GeoDataFrame:
        """Check intersections of culverts with other objects like roads, etc

        Args:
            shorten_line_offset (float, optional): shorten lines in analysis. Defaults to 0.01.

        Returns:
            gpd.GeoDataFrame: potential culverts without impossible intersections
        """
        # Crossing objects
        if isinstance(read_results, bool):
            self.read_results = read_results
        if self.read_results and self.potential_culverts_1 is not None:
            logging.info(
                f"   x {len(self.potential_culverts_1)} potential culverts already generated"
            )
            return self.potential_culverts_1

        crossing_objects = [
            "hydroobjecten",
            "overige_watergangen",
            "keringen",
            "nwb",
            "peilgebieden",
            "snelwegen",
            "spoorwegen",
        ]

        # Set culverts gdf
        culverts = self.potential_culverts_0.copy()

        logging.info("   x check intersections culverts with objects")
        for crossing_object in crossing_objects:
            crossing_gdf = getattr(self, crossing_object)

            # Merge operation
            culverts = culverts.merge(
                gpd.sjoin(
                    culverts,
                    crossing_gdf[["code", "geometry"]].rename(
                        columns={"code": f"{crossing_object}_code"}
                    ),
                    predicate="crosses",
                )[f"{crossing_object}_code"],
                how="left",
                left_index=True,
                right_index=True,
            )
            if crossing_object in ["hydroobjecten", "overige_watergangen"]:
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

            # Set true in new column for features with crossing objects.
            culverts[f"crossings_{crossing_object}"] = culverts[
                f"{crossing_object}_code"
            ].notna()

            no_intersections = len(
                culverts[culverts[f"crossings_{crossing_object}"] == True]
            )
            logging.info(f"     - {crossing_object} ({no_intersections} crossings)")

        # Remove potential culverts crossing: main road, railways, barriers,
        # hydroobjects and other waterways. (Drop columns for these objects).
        for crossing_object in crossing_objects:
            if crossing_object in [
                "hydroobjecten",
                "overige_watergangen",
                "keringen",
                "snelwegen",
                "spoorwegen",
            ]:
                culverts = culverts[culverts[f"crossings_{crossing_object}"] == False]

        # Write data
        self.potential_culverts_1 = culverts.copy()
        self.potential_culverts_1.to_file(
            Path(self.dir_results, "potential_culverts_1.gpkg")
        )
        logging.info(
            f"     - {len(self.potential_culverts_1)} potential culverts remaining"
        )
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
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 1

        # 2e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 2

        # 3e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 3

        # 4e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 4

        # 5e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 5

        # 6e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 6

        # 7e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 7

        # 8e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == False)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 8

        # 9e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 9

        # 10e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 10

        # 11e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 11

        # 12e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 12

        # 13e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == False),
            "score",
        ] = 13

        # 14e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "hydroobjecten")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 14

        # 15e voorkeur
        culverts.loc[
            (culverts["line_type"].isin(["dangling_start", "dangling_end"]))
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 15

        # 16e voorkeur
        culverts.loc[
            (culverts["line_type"] == "other")
            & (culverts["WaterLineType"] == "overige_watergangen")
            & (culverts["crossings_peilgebieden"] == True)
            & (culverts["crossings_nwb"] == True),
            "score",
        ] = 16

        logging.info("   x removing culverts without score")

        # Drop culverts that don't score 1-16
        culverts = culverts[culverts["score"].notna()]

        # Set copy and save data
        self.potential_culverts_2 = culverts.copy()
        self.potential_culverts_2.to_file(
            Path(self.dir_results, "potential_culverts_2.gpkg")
        )
        logging.info(
            f"     - {len(self.potential_culverts_2)} potential culverts remaining"
        )
        return self.potential_culverts_2

    def select_correct_score_based_on_score_and_length(
        self,
        read_results=None,
        factor_angle_on_length=3,
    ) -> gpd.GeoDataFrame:
        if isinstance(read_results, bool):
            self.read_results = read_results
        if self.read_results and self.potential_culverts_3 is not None:
            logging.info(
                f"   x {len(self.potential_culverts_3)} potential culverts already generated"
            )
            return self.potential_culverts_3

        # Create copy of potential culverts and calculate length
        culverts = self.potential_culverts_2.copy()
        culverts["length"] = culverts.geometry.length
        culverts["angle_culvert"] = culverts.apply(
            lambda row: calculate_angle_reverse(row["geometry"])
            if row["line_type_start_end"] == "dangling_start"
            else calculate_angle_start(row["geometry"]),
            axis=1,
        )
        culverts["angle_waterline"] = None

        self.overige_watergangen = self.overige_watergangen.explode()

        self.overige_watergangen["angle_start"] = self.overige_watergangen[
            "geometry"
        ].apply(calculate_angle_start)
        self.overige_watergangen["angle_end"] = self.overige_watergangen[
            "geometry"
        ].apply(calculate_angle_end)

        code_to_angle_start = self.overige_watergangen.set_index("code")[
            "angle_start"
        ].to_dict()
        code_to_angle_end = self.overige_watergangen.set_index("code")[
            "angle_end"
        ].to_dict()

        culverts["angle_waterline"] = culverts.apply(
            lambda row: code_to_angle_start.get(row["dangling_code"])
            if row["line_type_start_end"] == "dangling_start"
            else code_to_angle_end.get(row["dangling_code"])
            if row["line_type_start_end"] == "dangling_end"
            else None,
            axis=1,
        )

        culverts["angle_difference"] = culverts.apply(
            lambda row: calculate_angle_difference(
                row["angle_culvert"], row["angle_waterline"]
            ),
            axis=1,
        )

        def calculate_fictive_length(
            length, angle, factor_angle_on_length=factor_angle_on_length
        ):
            # Calculate the factor based on the angle
            factor = 1 + (factor_angle_on_length - 1) * abs(angle) / 90
            # Calculate the fictive length
            fictive_length = length * factor
            return fictive_length

        # Calculate fictive length
        culverts["fictive_length"] = culverts.apply(
            lambda row: calculate_fictive_length(
                row["length"], row["angle_difference"]
            ),
            axis=1,
        )

        self.potential_culverts_pre_filter = culverts.copy()
        # Drop features with a fictive length larger than 40
        culverts = culverts[culverts["fictive_length"] <= self.max_culvert_length]

        # Keep only shortest potential culverts when the score is equal
        culverts = culverts.sort_values(by=["dangling_id", "score", "fictive_length"])
        culverts = culverts.groupby(["dangling_id", "score"]).first().reset_index()

        logging.info(f"     - {len(culverts)} potential culverts remaining")

        # Remove potential culverts when they are not in the highest score group of 4 that is available
        # Define score groups
        def get_score_group(score):
            if score <= 4:
                return 1
            elif score <= 8:
                return 2
            elif score <= 12:
                return 3
            else:
                return 4

        # Create a new column for score group
        culverts["score_group"] = culverts["score"].apply(get_score_group)
        # Find the highest score group for each dangling_id
        highest_group = (
            culverts.groupby("dangling_id")["score_group"].min().reset_index()
        )
        highest_group.columns = ["dangling_id", "highest_score_group"]
        # Merge back to the original DataFrame to filter based on highest score group
        culverts = culverts.merge(highest_group, on="dangling_id")
        # Filter to keep only scores in the highest score group
        culverts = culverts[culverts["score_group"] == culverts["highest_score_group"]]
        # Drop the helper column
        culverts = culverts.drop(columns=["score_group", "highest_score_group"])

        logging.info(f"     - {len(culverts)} potential culverts remaining")

        # Remove potential culverts of score 9 and lower when the other side of the watergang that a culvert with score 8 or less.
        # Find dangling IDs that have a corresponding higher score in the same dangling_code
        def filter_low_scores(group):
            # Find the lowest score in the group
            low_score_ids = group[group["score"] >= 9]["dangling_id"].unique()

            # Check for any scores >= 8 in the group
            if any(group["score"] >= 8):
                # Drop all lines with low scores from dangling IDs
                group = group[~group["dangling_id"].isin(low_score_ids)]

            return group

        # Apply the filtering function to each group
        culverts = (
            culverts.groupby("dangling_code")
            .apply(filter_low_scores)
            .reset_index(drop=True)
        )

        logging.info(f"     - {len(culverts)} potential culverts remaining")

        # Select correct score within group of 4 based on defined logic between score and lenght.
        def select_best_lines(df, offset_1_2=0.75, offset_12_34=0.6, offset_3_4=0.75):
            dangling_ids = df.dangling_id.unique()

            best_lines = []

            for dangling_id in dangling_ids:
                df_dangling_id = df[df.dangling_id == dangling_id]

                selected_score = None
                selected_length = None
                # Loop through the defined groups
                for base_score in [1, 5, 9, 13]:
                    # Find the corresponding scores in the group
                    score_a = base_score
                    score_b = base_score + 1
                    score_c = base_score + 2
                    score_d = base_score + 3

                    # Get lines for the score group
                    line_a = df_dangling_id[df_dangling_id["score"] == score_a]
                    line_b = df_dangling_id[df_dangling_id["score"] == score_b]
                    line_c = df_dangling_id[df_dangling_id["score"] == score_c]
                    line_d = df_dangling_id[df_dangling_id["score"] == score_d]

                    # Select from score_a and score_b
                    if not line_a.empty:
                        selected_length = line_a["fictive_length"].values[0]
                        selected_score = score_a

                    if not line_b.empty and (
                        selected_length is None
                        or line_b["fictive_length"].values[0]
                        <= offset_1_2 * selected_length
                    ):
                        selected_length = line_b["fictive_length"].values[0]
                        selected_score = score_b

                    # Check conditions for selecting score_c
                    if selected_length is None and not line_c.empty:
                        selected_length = line_c["fictive_length"].values[0]
                        selected_score = score_c

                    if (
                        not line_c.empty
                        and (selected_score in [score_a, score_b])
                        and (
                            line_c["fictive_length"].values[0]
                            <= offset_12_34 * selected_length
                        )
                    ):
                        selected_length = line_c["fictive_length"].values[0]
                        selected_score = score_c

                    # Check conditions for selecting score_d
                    if not line_d.empty:
                        if selected_score == score_c and (
                            line_d["fictive_length"].values[0]
                            <= offset_3_4 * selected_length
                        ):
                            selected_length = line_d["fictive_length"].values[0]
                            selected_score = score_d
                        elif selected_score in [score_a, score_b] and (
                            line_d["fictive_length"].values[0]
                            <= offset_12_34 * selected_length
                        ):
                            if line_c.empty or (
                                line_d["fictive_length"].values[0]
                                <= offset_3_4 * line_c["fictive_length"].values[0]
                            ):
                                selected_length = line_d["fictive_length"].values[0]
                                selected_score = score_d

                    # If selected_length is still None, check for line_d
                    if selected_length is None and not line_d.empty:
                        selected_length = line_d["fictive_length"].values[0]
                        selected_score = score_d

                if selected_score is not None:
                    best_lines.append(
                        {"dangling_id": dangling_id, "selected_score": selected_score}
                    )

            return pd.DataFrame(best_lines)

        culverts = culverts.merge(
            select_best_lines(culverts), how="left", on="dangling_id"
        )
        culverts = culverts[culverts["selected_score"] == culverts["score"]]
        culverts = culverts.set_crs(28992)

        # Set copy and save data
        self.potential_culverts_3 = culverts.copy()
        logging.info(
            f"     - {len(self.potential_culverts_3)} potential culverts remaining"
        )
        self.potential_culverts_3.to_file(
            Path(self.dir_results, "potential_culverts_3.gpkg")
        )
        self.potential_culverts_pre_filter.to_file(
            Path(self.dir_results, "potential_culverts_pre_filter.gpkg"),
        )

        return self.potential_culverts_3, self.potential_culverts_pre_filter

    def post_process_potential_culverts(self):
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

        culverts = culverts.loc[culverts.groupby("code_pair")["length"].idxmin()]

        # Define the filter function
        def should_remove(row):
            dangling_code = row["dangling_code"]
            code = row["code"]

            # Check if both dangling_code and code appear in any of the lists in duplicate_code_groups
            for codes_list in duplicate_code_groups.values():
                if dangling_code in codes_list and code in codes_list:
                    return True  # Mark for removal if both are found in the same list

            return False  # Otherwise, don't remove

        # Filter culverts using the modified function
        culverts = culverts[~culverts.apply(should_remove, axis=1)].copy()

        self.potential_culverts_4 = culverts.copy()
        logging.info(
            f"     - {len(self.potential_culverts_4)} potential culverts remaining"
        )
        self.potential_culverts_4.to_file(
            Path(self.dir_results, "potential_culverts_4.gpkg"),
        )
        return self.potential_culverts_4

    def splits_hydroobjecten_by_endpoints_of_culverts_and_combine(self):
        other_culverts_hydro = self.potential_culverts_4[
            (self.potential_culverts_4["line_type"] == "other")
            & (self.potential_culverts_4["WaterLineType"] == "hydroobjecten")
        ].copy()
        other_culverts_other = self.potential_culverts_4[
            (self.potential_culverts_4["line_type"] == "other")
            & (self.potential_culverts_4["WaterLineType"] == "overige_watergangen")
        ].copy()

        self.hydroobjecten_processed_0 = split_waterways_by_endpoints(
            self.hydroobjecten, other_culverts_hydro
        )
        logging.info("     - hydroobjecten gesplit")
        self.overige_watergangen_processed_0 = split_waterways_by_endpoints(
            self.overige_watergangen, other_culverts_other
        )

        logging.info("     - overige watergangen gesplit")

        self.hydroobjecten_processed_0.to_file(
            Path(self.dir_results, "hydroobjecten_processed_0.gpkg")
        )

        self.overige_watergangen_processed_0.to_file(
            Path(self.dir_results, "overige_watergangen_processed_0.gpkg")
        )

        culverts_hydro = self.potential_culverts_4[
            (self.potential_culverts_4["WaterLineType"] == "hydroobjecten")
        ]

        # Extract end points and retain 'unique_id' and 'dangling_code'
        culverts_hydro["end_point"] = culverts_hydro.geometry.apply(
            lambda line: Point(line.coords[-1])
            if isinstance(line, LineString)
            else None
        )

        # Create GeoDataFrame with end points
        end_points_gdf = gpd.GeoDataFrame(
            culverts_hydro[["end_point", "unique_id", "dangling_code"]].dropna(
                subset=["end_point"]
            ),
            geometry="end_point",
            crs=culverts_hydro.crs,
        )

        # Rename the geometry column to 'geometry'
        end_points_gdf = end_points_gdf.rename(
            columns={"end_point": "geometry"}
        ).set_geometry("geometry")

        self.outflow_nodes_overige_watergangen = end_points_gdf.copy()

        # Save to file
        self.outflow_nodes_overige_watergangen.to_file(
            Path(self.dir_results, "outflow_nodes_overige_watergangen.gpkg")
        )

        return (
            self.overige_watergangen_processed_0,
            self.hydroobjecten_processed_0,
            self.outflow_nodes_overige_watergangen,
        )

    def check_culverts_direction(self):
        culvert = self.potential_culverts_4.copy()
        lines = pd.concat(
            [self.hydroobjecten_processed_0, self.overige_watergangen_processed_0],
            ignore_index=True,
        )

        # Extract starting and ending points from gdf2
        gdf2_start_points = lines.geometry.apply(lambda geom: geom.coords[0])
        gdf2_end_points = lines.geometry.apply(lambda geom: geom.coords[-1])

        # Initialize the 'flipped' column if it doesn't exist
        if "flipped" not in culvert.columns:
            culvert["flipped"] = 0

        # Apply the check_and_flip function to each line in gdf1
        results = culvert.geometry.apply(
            lambda line: check_and_flip(
                line, gdf2_start_points.tolist(), gdf2_end_points.tolist()
            )
        )

        # Update the geometry and flipped columns
        culvert["geometry"] = results.apply(lambda x: x[0])
        culvert["flipped"] += results.apply(lambda x: 1 if x[1] else 0)
        logging.info("     - culvert direction checked")

        self.potential_culverts_5 = culvert.copy()

        self.potential_culverts_5.to_file(
            Path(self.dir_results, "potential_culverts_5.gpkg")
        )
        return self.potential_culverts_5

    def combine_culvert_with_line(self):
        culvert = self.potential_culverts_5.copy()
        culvert_dict = (
            culvert.groupby("dangling_code")["geometry"].apply(list).to_dict()
        )
        lines = self.overige_watergangen_processed_0.copy()

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

        self.overige_watergangen_processed_1 = lines.copy()

        self.overige_watergangen_processed_1.to_file(
            Path(self.dir_results, "overige_watergangen_processed_1.gpkg")
        )
        return self.overige_watergangen_processed_1

    def splits_hydroobjecten_by_endpoind_of_culverts_and_combine_2(
        self, write_results=True
    ):
        # split overige watergangen opnieuw
        overige_watergangen = self.overige_watergangen_processed_1.copy()
        overige_watergangen = split_waterways_by_endpoints(
            overige_watergangen, overige_watergangen
        )
        logging.info("     - overige watergangen weer gesplit")

        self.overige_watergangen_processed_2 = overige_watergangen.copy()

        if write_results:
            self.overige_watergangen_processed_2.to_file(
                Path(self.dir_results, "overige_watergangen_processed_2.gpkg"),
            )

        return self.overige_watergangen_processed_2

    def get_shortest_path_from_overige_watergangen_to_hydroobjects(
        self, write_results=False
    ):
        logging.info(f"   x redirect 'overige watergangen' based on shortest path")
        # create networkx graph
        logging.info(f"     - create sub networks")
        nodes, edges, _ = create_graph_from_edges(
            self.overige_watergangen_processed_2, directed=True
        )
        _, _, G = create_graph_from_edges(
            self.overige_watergangen_processed_2, directed=False
        )

        # get outflow points and add nodes information
        outflow_nodes = (
            self.outflow_nodes_overige_watergangen[
                ["unique_id", "dangling_code", "geometry"]
            ]
            .sjoin(nodes, how="inner")
            .reset_index(drop=True)
            .drop(columns=["index_right"], errors="ignore")
        )
        self.outflow_nodes_overige_watergangen = outflow_nodes.copy()

        # get shortest path including length from all nodes to outflow points
        logging.info(f"     - find shortest path")
        len_outflow_node, matrix = nx.multi_source_dijkstra(
            G, [int(n) for n in outflow_nodes["nodeID"].values], weight="geometry_len"
        )
        node_to_outflow_node = pd.DataFrame(
            {
                "nodeID": matrix.keys(),
                "outflow_node": [v[0] for v in matrix.values()],
                "outflow_node_dist": [len_outflow_node[node] for node in matrix.keys()],
            }
        )

        # Merge nodes with node_to_outflow_node
        overige_watergangen_nodes = nodes.merge(
            node_to_outflow_node, how="outer", left_on="nodeID", right_on="nodeID"
        )
        overige_watergangen_nodes["outflow_node"] = (
            overige_watergangen_nodes["outflow_node"].fillna(-999).astype(int)
        )

        # To get length to outflow point at start and end: merge edges with nodes, first on node_start, then on node_end
        edges = (
            edges.merge(
                overige_watergangen_nodes[
                    ["nodeID", "outflow_node", "outflow_node_dist"]
                ].rename(
                    columns={
                        "outflow_node": "outflow_start",
                        "outflow_node_dist": "outflow_start_dist",
                    }
                ),
                how="left",
                left_on="node_start",
                right_on="nodeID",
            )
            .merge(
                overige_watergangen_nodes[
                    ["nodeID", "outflow_node", "outflow_node_dist"]
                ].rename(
                    columns={
                        "outflow_node": "outflow_end",
                        "outflow_node_dist": "outflow_end_dist",
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
                edge.outflow_node = edge.outflow_end
                edge.outflow_node_dist = edge.outflow_end_dist
            else:
                edge.outflow_node = edge.outflow_start
                edge.outflow_node_dist = edge.outflow_start_dist
                edge.reversed_direction = True
            return edge

        edges["outflow_node"] = -999
        edges["outflow_node_dist"] = -999.0
        edges["reversed_direction"] = False
        edges["outflow_end_dist"] = edges["outflow_end_dist"].fillna(99999.9)
        edges["outflow_start_dist"] = edges["outflow_start_dist"].fillna(99999.9)
        edges = edges.apply(lambda x: select_shortest_direction(x), axis=1)

        # clean edges by removing unconnected edges and reversing direction if required
        edges_cleaned = edges[edges["outflow_node"] != -999]
        edges_cleaned.loc[edges_cleaned["reversed_direction"], "geometry"] = (
            edges_cleaned.loc[edges_cleaned["reversed_direction"], "geometry"].reverse()
        )

        self.overige_watergangen_processed_3_nodes = overige_watergangen_nodes.copy()
        self.overige_watergangen_processed_3 = edges_cleaned.copy()

        if write_results:
            outflow_nodes.to_file(
                Path(self.dir_results, "outflow_nodes_overige_watergangen.gpkg")
            )
            self.overige_watergangen_processed_3_nodes.to_file(
                Path(
                    self.dir_results, "overige_watergangen_processed_3_nodes.gpkg"
                )
            )
            self.overige_watergangen_processed_3.to_file(
                Path(self.dir_results, "overige_watergangen_processed_3.gpkg")
            )
        return outflow_nodes, overige_watergangen_nodes, edges, edges_cleaned

    def generate_folium_map(
        self, html_file_name=None, base_map="Light Mode", open_html=False, zoom_start=12
    ):
        # Make figure

        hydro_4326 = self.hydroobjecten_processed_0.to_crs(4326)
        # Calculate the extent (bounding box) of your GeoDataFrame
        bounds = hydro_4326.total_bounds  # returns (minx, miny, maxx, maxy)

        # Center the map around the mean coordinates of the bounds
        center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]
        m = folium.Map(
            location=center,
            zoom_start=zoom_start,
            tiles=None,
        )

        folium.GeoJson(
            self.hydroobjecten_processed_0.geometry,
            name="A/B-Watergangen",
            color="blue",
            fill_color="blue",
            zoom_on_click=True,
            z_index=2,
        ).add_to(m)

        folium.GeoJson(
            self.overige_watergangen.geometry,
            name="C-Watergangen - Zonder duikers",
            color="lightblue",
            fill_color="blue",
            zoom_on_click=True,
            z_index=0,
        ).add_to(m)

        folium.GeoJson(
            self.potential_culverts_5.geometry,
            name="C-Watergangen - Gevonden Duikers",
            color="red",
            fill_color="blue",
            zoom_on_click=True,
            z_index=1,
        ).add_to(m)

        if self.outflow_nodes_overige_watergangen is not None:
            folium.GeoJson(
                self.outflow_nodes_overige_watergangen.geometry,
                name="C-Watergangen - Uitstroompunten",
                marker=folium.Circle(
                    radius=3,
                    fill_color="orange",
                    fill_opacity=1.0,
                    color="orange",
                    weight=1,
                    z_index=3,
                ),
                show=False,
            ).add_to(m)

        if self.overige_watergangen_processed_3 is not None:
            add_categorized_lines_to_map(
                m=m,
                lines_gdf=self.overige_watergangen_processed_3,
                layer_name=f"C-Watergangen - Gegroepeerd per uitstroompunt",
                control=True,
                lines=True,
                line_color_column="outflow_node",
                line_color_cmap=None,
                show=False,
                z_index=2,
            )

            if "outflow_node" in self.overige_watergangen_processed_3.columns:
                fg = folium.FeatureGroup(
                    name="C-Watergangen - Labels uitstroompunten",
                    control=True,
                    show=False,
                ).add_to(m)

                add_labels_to_points_lines_polygons(
                    gdf=self.overige_watergangen_processed_3,
                    column="outflow_node",
                    label_fontsize=7,
                    label_decimals=0,
                    fg=fg,
                )

                add_labels_to_points_lines_polygons(
                    gdf=self.outflow_nodes_overige_watergangen,
                    column="nodeID",
                    label_fontsize=8,
                    label_decimals=0,
                    fg=fg,
                )

        m = add_basemaps_to_folium_map(m=m, base_map=base_map)

        folium.LayerControl(collapsed=False).add_to(m)
        m.add_child(folium.plugins.MeasureControl())

        self.folium_map = m
        if html_file_name is None:
            html_file_name = self.name + "_culvert_locations"

        self.folium_html_path = Path(self.path, f"{html_file_name}.html")
        m.save(self.folium_html_path)

        logging.info(f"   x html file saved: {html_file_name}.html")

        if open_html:
            webbrowser.open(Path(self.path, f"{html_file_name}.html"))
        return m
