import logging
from pathlib import Path

from pydantic import BaseModel, ConfigDict
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point
import numpy as np
import folium
from folium.features import DivIcon

from ..generator_basis import GeneratorBasis
from ..utils.general_functions import (
    shorten_line_two_vertices,
    line_to_vertices,
    split_waterways_by_endpoints,
    check_and_flip,
)
from ..utils.folium_utils import add_basemaps_to_folium_map
from ..utils.preprocess import preprocess_hydroobjecten

import warnings

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

    read_results: bool = False
    write_results: bool = False

    water_line_pnts: gpd.GeoDataFrame = None
    duplicates: gpd.GeoDataFrame = None
    potential_culverts_0: gpd.GeoDataFrame = None  # alle binnen 40m
    potential_culverts_1: gpd.GeoDataFrame = None  # filter kruizingen
    potential_culverts_2: gpd.GeoDataFrame = None  # scores
    potential_culverts_3: gpd.GeoDataFrame = None  # eerste resultaat
    potential_culverts_4: gpd.GeoDataFrame = None  # resultaat met nabewerking

    hydroobjecten_processed: gpd.GeoDataFrame = (
        None  # hydroobjecten including splits by culverts
    )
    overige_watergangen_processed: gpd.GeoDataFrame = (
        None  # overige_watergangen including splits by culverts
    )
    combined_hydroobjecten: gpd.GeoDataFrame = (
        None  # combined A, B, en C watergangen including splits
    )
    combined_end_product: gpd.GeoDataFrame = (
        None  # combined A, B, en C watergangen, and culverts. including splits
    )
    outflow_points_overig_to_hydro: gpd.GeoDataFrame = (
        None  # outflow points from overige watergangen to hydroobjecten
    )

    folium_map: folium.Map = None

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
                self.water_line_pnts["line_type"] == "dangling"
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
            self.water_line_pnts["line_type"] == "dangling"
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
            dir_results = Path(self.path, "1_tussenresultaat")
            self.water_line_pnts.to_file(
                Path(dir_results, "water_line_pnts.gpkg"), layer="water_line_pnts"
            )
        return self.water_line_pnts, self.duplicates

    def find_potential_culvert_locations(
        self,
        water_line_pnts=None,
        max_culvert_length=40,
        read_results=None,
        write_results=None,
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

        logging.info("  x find potential culvert locations")

        # Filter for end points (only overige watergangen and not when connected)
        end_pnts = self.water_line_pnts[
            (self.water_line_pnts["line_type"] == "dangling")
            & (self.water_line_pnts["WaterLineType"] == "overige_watergangen")
        ]
        end_pnts = end_pnts.drop_duplicates(subset="geometry", keep=False)
        end_pnts = end_pnts.rename(
            columns={"code": "dangling_code", "unique_id": "dangling_id"}
        )
        end_pnts_orig_geometry = end_pnts.geometry

        # end points have buffer as geometry, perform spatial join with points
        logging.debug("    - spatial join vertices")
        end_pnts["geometry"] = end_pnts.buffer(max_culvert_length)
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
        logging.debug("    - add data to potential culverts")
        end_pnts["geometry"] = end_pnts_orig_geometry
        end_pnts = end_pnts.rename(columns={"geometry": "geometry2"})

        end_pnts = pd.merge(
            self.water_line_pnts,
            end_pnts[["unique_id", "dangling_id", "dangling_code", "geometry2"]],
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
            dir_results = Path(self.path, "1_tussenresultaat")
            self.potential_culverts_0.to_file(
                Path(dir_results, "potential_culverts_0.gpkg"),
                layer="potential_culverts_0",
            )

        logging.debug(
            f"    - {len(self.potential_culverts_0)} potential culverts generated"
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
            logging.debug(f"    - {crossing_object} ({no_intersections} crossings)")

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
            Path(self.path, "1_tussenresultaat", "potential_culverts_1.gpkg"),
            layer="potential_culverts_1",
        )
        logging.debug(
            f"    - {len(self.potential_culverts_1)} potential culverts remaining"
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            (culverts["line_type"] == "dangling")
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
            Path(self.path, "1_tussenresultaat", "potential_culverts_2.gpkg"),
            layer="potential_culverts_2",
        )
        logging.debug(
            f"    - {len(self.potential_culverts_2)} potential culverts remaining"
        )
        return self.potential_culverts_2

    def select_correct_score_based_on_score_and_length(
        self, read_results=None
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

        # Keep only shortest potential culverts when the score is equal
        culverts = culverts.sort_values(by=["dangling_id", "score", "length"])
        culverts = culverts.groupby(["dangling_id", "score"]).first().reset_index()

        logging.debug(f"    - {len(culverts)} potential culverts remaining")

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

        logging.debug(f"    - {len(culverts)} potential culverts remaining")

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

        logging.debug(f"    - {len(culverts)} potential culverts remaining")

        # Select correct score within group of 4 based on defined logic between score and lenght.
        def select_best_lines(df, offset_1_2=0.75, offset_12_34=0.5, offset_3_4=0.5):
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
                        selected_length = line_a["length"].values[0]
                        selected_score = score_a

                    if not line_b.empty and (
                        selected_length is None
                        or line_b["length"].values[0] <= offset_1_2 * selected_length
                    ):
                        selected_length = line_b["length"].values[0]
                        selected_score = score_b

                    # Check conditions for selecting score_c
                    if selected_length is None and not line_c.empty:
                        selected_length = line_c["length"].values[0]
                        selected_score = score_c

                    if (
                        not line_c.empty
                        and (selected_score in [score_a, score_b])
                        and (
                            line_c["length"].values[0] <= offset_12_34 * selected_length
                        )
                    ):
                        selected_length = line_c["length"].values[0]
                        selected_score = score_c

                    # Check conditions for selecting score_d
                    if not line_d.empty:
                        if selected_score == score_c and (
                            line_d["length"].values[0] <= offset_3_4 * selected_length
                        ):
                            selected_length = line_d["length"].values[0]
                            selected_score = score_d
                        elif selected_score in [score_a, score_b] and (
                            line_d["length"].values[0] <= offset_12_34 * selected_length
                        ):
                            if line_c.empty or (
                                line_d["length"].values[0]
                                <= offset_3_4 * line_c["length"].values[0]
                            ):
                                selected_length = line_d["length"].values[0]
                                selected_score = score_d

                    # If selected_length is still None, check for line_d
                    if selected_length is None and not line_d.empty:
                        selected_length = line_d["length"].values[0]
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
        logging.debug(
            f"    - {len(self.potential_culverts_3)} potential culverts remaining"
        )
        self.potential_culverts_3.to_file(
            Path(self.path, "1_tussenresultaat", "potential_culverts_3.gpkg"),
            layer="potential_culverts_3",
        )

        return self.potential_culverts_3

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
        logging.debug(
            f"    - {len(self.potential_culverts_4)} potential culverts remaining"
        )
        self.potential_culverts_4.to_file(
            Path(self.path, "1_tussenresultaat", "potential_culverts_4.gpkg"),
            layer="potential_culverts_4",
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

        self.hydroobjecten_processed = split_waterways_by_endpoints(
            self.hydroobjecten, other_culverts_hydro
        )
        logging.debug("hydroobjecten gesplit")
        self.overige_watergangen_processed = split_waterways_by_endpoints(
            self.overige_watergangen, other_culverts_other
        )

        logging.debug("overige watergangen gesplit")

        self.combined_hydroobjecten = pd.concat(
            [self.hydroobjecten_processed, self.overige_watergangen_processed],
            ignore_index=True,
        )

        logging.debug("hydroobjecten en overige watergangen gecombineerd")

        self.hydroobjecten_processed.to_file(
            Path(self.path, "1_tussenresultaat", "hydroobjecten_processed.gpkg"),
            layer="hydroobjecten_processed",
        )

        self.overige_watergangen_processed.to_file(
            Path(self.path, "1_tussenresultaat", "overige_watergangen_processed.gpkg"),
            layer="overige_watergangen_processed",
        )

        self.combined_hydroobjecten.to_file(
            Path(self.path, "1_tussenresultaat", "combined_hydroobjecten.gpkg"),
            layer="combined_hydroobjecten",
        )

        return (
            self.combined_hydroobjecten,
            self.overige_watergangen_processed,
            self.hydroobjecten_processed,
        )

    def check_culverts_direction(self):
        culvert = self.potential_culverts_4.copy()
        lines = self.combined_hydroobjecten.copy()

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
        logging.debug("culvert direction checked")

        self.potential_culverts_4 = culvert.copy()

        self.potential_culverts_4.to_file(
            Path(self.path, "1_tussenresultaat", "potential_culverts_4.gpkg"),
            layer="potential_culverts_4",
        )
        return self.potential_culverts_4

    def combine_culvert_with_line(self):
        culvert = self.potential_culverts_4.copy()
        culvert_dict = (
            culvert.groupby("dangling_code")["geometry"].apply(list).to_dict()
        )
        lines = self.overige_watergangen_processed.copy()

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

        logging.debug("culverts combined with watergangen")

        self.overige_watergangen_processed = lines.copy()

        self.overige_watergangen_processed.to_file(
            Path(self.path, "1_tussenresultaat", "overige_watergangen_processed.gpkg"),
            layer="overige_watergangen_processed",
        )
        return self.overige_watergangen_processed

    def splits_hydroobjecten_by_endpoind_of_culverts_and_combine_2(self):
        # split overige watergangen opnieuw
        overige_watergangen = self.overige_watergangen_processed.copy()
        overige_watergangen = split_waterways_by_endpoints(
            overige_watergangen, overige_watergangen
        )
        logging.debug("overige watergangen weer gesplit")

        self.overige_watergangen_processed = overige_watergangen.copy()
        self.combined_hydroobjecten = pd.concat(
            [self.hydroobjecten_processed, self.overige_watergangen_processed],
            ignore_index=True,
        )

        logging.debug("hydroobjecten en overige watergangen gecombineerd")

        self.overige_watergangen_processed.to_file(
            Path(self.path, "1_tussenresultaat", "overige_watergangen_processed.gpkg"),
            layer="overige_watergangen_processed",
        )

        self.combined_hydroobjecten.to_file(
            Path(self.path, "1_tussenresultaat", "combined_hydroobjecten.gpkg"),
            layer="combined_hydroobjecten",
        )

        return self.overige_watergangen_processed, self.combined_hydroobjecten

    def generate_outflow_points_overige_watergangen_to_hydroobjecten(self):
        culverts_hydro = self.potential_culverts_4[
            (self.potential_culverts_4["WaterLineType"] == "hydroobjecten")
            & (self.potential_culverts_4["flipped"] % 2 == 0)
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

        self.outflow_points_overig_to_hydro = end_points_gdf.copy()

        # Save to file
        self.outflow_points_overig_to_hydro.to_file(
            Path(self.path, "1_tussenresultaat", "outflow_points_overig_to_hydro.gpkg"),
            layer="outflow_points_overig_to_hydro",
        )

        return self.outflow_points_overig_to_hydro

    def generate_folium_map(self, base_map="OpenStreetMap"):
        # Make figure

        hydro_4326 = self.hydroobjecten_processed.to_crs(4326)
        # Calculate the extent (bounding box) of your GeoDataFrame
        bounds = hydro_4326.total_bounds  # returns (minx, miny, maxx, maxy)

        # Center the map around the mean coordinates of the bounds
        center = [(bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2]
        m = folium.Map(
            location=center,
            zoom_start=14,
            tiles=None,
        )

        folium.GeoJson(
            self.hydroobjecten_processed.geometry,
            name="A/B-Watergangen",
            color="blue",
            fill_color="blue",
            zoom_on_click=True,
            z_index=1,
        ).add_to(m)

        folium.GeoJson(
            self.overige_watergangen_processed.geometry,
            name="C-Watergangen",
            color="cadetblue",
            fill_color="blue",
            zoom_on_click=True,
            z_index=1,
        ).add_to(m)

        folium.GeoJson(
            self.potential_culverts_4.geometry,
            name="Duikers",
            color="red",
            fill_color="blue",
            zoom_on_click=True,
            z_index=1,
        ).add_to(m)

        m = add_basemaps_to_folium_map(m=m, base_map=base_map)

        folium.LayerControl(collapsed=False).add_to(m)

        self.folium_map = m
        return m
