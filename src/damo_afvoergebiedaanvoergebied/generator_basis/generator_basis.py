import logging
from pathlib import Path

import folium
import geopandas as gpd
import pandas as pd
import rioxarray
import xarray as xr
from pydantic import BaseModel, ConfigDict
from ..utils.preprocess import preprocess_hydroobject
from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import (
    calculate_angles_of_edges_at_nodes,
    define_list_upstream_downstream_edges_ids,
    find_node_edge_ids_in_directed_graph,
    calculate_discharges_of_edges_at_nodes,
    select_downstream_upstream_edges_angle,
    select_downstream_upstream_edges_discharge,
)


def drop_band_dim(da):
    if "band" not in da.dims:
        return da
    if da.sizes["band"] == 1:
        return da.squeeze("band", drop=True)
    return da.isel(band=0)


def select_downstream_upstream_edges(
    nodes, 
    min_difference_angle=10.0, 
    min_difference_discharge_factor=2.0
):
    logging.info("   x find downstream upstream edges")
    if "specifieke_afvoer" not in nodes:
        logging.info(f"     - use angle using min_difference_angle [{min_difference_angle}deg]")
        nodes = select_downstream_upstream_edges_angle(
            nodes, min_difference_angle=min_difference_angle
        )
    else:
        logging.info(f"     - use discharge distribution [factor{min_difference_discharge_factor:.3f}]")
        nodes = select_downstream_upstream_edges_discharge(
            nodes, min_difference_discharge_factor=min_difference_discharge_factor
        )
    return nodes



def analyse_netwerk_add_information_to_nodes_edges(
    edges,
    nodes, 
    min_difference_angle=10.0, 
    min_difference_discharge_factor=2.0
):
    logging.info("   x get upstream downstream edges ids")
    nodes = define_list_upstream_downstream_edges_ids(
        node_ids=nodes.nodeID.values, nodes=nodes, edges=edges
    )
    logging.info("   x calculate angles of edges to nodes")
    nodes, edges = calculate_angles_of_edges_at_nodes(
        nodes=nodes, edges=edges
    )
    logging.info("   x calculate discharges of edges to nodes")
    nodes, edges = calculate_discharges_of_edges_at_nodes(
        nodes=nodes, edges=edges
    )
    logging.info("   x select downstream upstream edges")
    nodes = select_downstream_upstream_edges(
        nodes=nodes,
        min_difference_angle=min_difference_angle,
        min_difference_discharge_factor=min_difference_discharge_factor
    )
    return nodes, edges


class GeneratorBasis(BaseModel):
    """Basis class for all Generators
    
    Basis class for reading all basis datasets (based on attributes) 
    from the subdirectory basisdata (dir_basisdata) and optionally read results

    Parameters
    ----------
    path : pathlib.Path
        windowspath to analysis folder
    name : str
        String representing name of case (equal to folder name)
    dir_basisdata : str | pathlib.Path
        String representing subfolder with basisdata
    dir_results : str | pathlib.Path
        String representing subfolder with results
    read_results : bool
        setting to know whether results in dir_results should be read
    write_results : bool
        setting to know whether results should be written in dir_results
    required_results : list[str]
        attributes required as input (in the results directory)
    folium_map : folium.Map
        folium map
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    path: Path = None
    name: str = None
    base_dir: Path = None
    waterschap: str = None

    dir_basisdata: str | Path = "0_basisdata"
    dir_results: str | Path | None = "1_resultaat"

    read_results: bool = False
    write_results: bool = False

    required_results: list[str] = []

    folium_map: folium.Map = None


    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        if self.path is not None:
            self.check_case_path_directory(path=self.path)
            self.read_data_from_case()
            self.read_required_data_from_case()


    def check_case_path_directory(self, path: Path):
        """Check on existence case directory and structure

        Parameters
        ----------
        path : Path
            path to case directory. name of directory is used as case name.
            self.path and self.name are set

        Raises ValueErrors in case directory and 0_basisdata directory not exist
        """
        if not path.exists() and path.is_dir():
            raise ValueError(
                f"provided path [{path}] does not exist or is not a directory"
            )
        self.path = path
        self.name = self.path.name
        logging.info(f' ### Waterschap "{self.waterschap.capitalize()}" Case "{self.name.capitalize()}" ###')

        # check if directories 0_basisdata and 1_resultaat exist
        if isinstance(self.dir_basisdata, str):
            self.dir_basisdata = Path(self.path, self.dir_basisdata)
        if not isinstance(self.dir_basisdata, Path) or not self.dir_basisdata.exists():
            raise ValueError(
                f"provided [{self.dir_basisdata}] is not a path or does not exist"
            )

        if self.dir_results is not None:
            if isinstance(self.dir_results, str):
                self.dir_results = Path(self.path, self.dir_results)
            if isinstance(self.dir_results, Path):
                if not self.dir_results.exists():
                    self.dir_results.mkdir(parents=True, exist_ok=True)

        logging.info(f"     - dir basisdata    = {self.dir_basisdata}")
        logging.info(f"     - dir results      = {self.dir_results}")


    def read_data_from_case(self, path: Path = None, read_results: bool = None):
        """Read data from case: including basis data and intermediate results

        Parameters
        ----------
        path : Path, optional
            Path to the case directory including directories 0_basisdata and
            1_resultaat. Directory name is used as name for the case,
            by default None
        read_results : bool, optional
            if True, it reads already all resulst from, by default None
        """
        if path is not None and path.exists():
            self.check_case_path_directory(path=path)

        def read_attributes_from_folder(path_dir: Path):
            for f in path_dir.glob("**/*"):
                if hasattr(self, f.stem):
                    logging.info(f"     - load dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        with rioxarray.open_rasterio(f) as raster:
                            ds = drop_band_dim(raster.load())
                            ds.name = f.stem
                            setattr(self, f.stem, ds)

        logging.info(f"   x read basisdata")
        if self.dir_basisdata is not None and self.dir_basisdata.exists():
            read_attributes_from_folder(self.dir_basisdata)
        
        if self.read_results:
            logging.info(f"   x read results")
            if self.dir_results is not None and self.dir_results.exists():
                read_attributes_from_folder(self.dir_results)


    def read_required_data_from_case(self):
        """Check if required results (from previous analyses) is available

        This function check if all required datasets are imported.

        Raises
        ------
        ValueError
            if required dataset is not available
        """
        for required_dataset in self.required_results:
            for f in self.dir_results.glob("**/*"):
                if f.stem != required_dataset:
                    continue
                if hasattr(self, f.stem) and getattr(self, f.stem) is None:
                    logging.info(f"     - load dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        with rioxarray.open_rasterio(f) as raster:
                            ds = drop_band_dim(raster.load())
                            ds.name = f.stem
                            setattr(self, f.stem, ds)
            if getattr(self, required_dataset) is None:
                logging.info(f" * dataset {required_dataset} is missing - check if absolutely required")


    def use_processed_hydroobject(self, processed_file="processed", force_preprocess=False, snapping_distance=0.05):
        """actualize hydroobject and overige_watergang

        replaces hydroobject and overige_watergang with the newest processed attributes

        Parameters
        ----------
        processed_file : str, optional
            suffix of processed files, by default "processed"
        """
        if self.snapping_distance is not None:
            snapping_distance = self.snapping_distance

        for watergang in ["hydroobject", "overige_watergang"]:
            if getattr(self, watergang, None) is None:
                logging.info(f"     - attribute {watergang} does not exist")
                continue
            
            watergang_processed_file_name = None
            attributes = dir(self)
            files_in_dirs = list(self.dir_basisdata.glob("**/*")) + list(self.dir_results.glob("**/*"))

            for file in files_in_dirs:
                if (
                    f"{watergang}_{processed_file}" in file.stem
                    and "nodes" not in file.stem
                    and file.stem in self.required_results
                ):
                    watergang_processed_file_name = file
           
            if force_preprocess or watergang_processed_file_name is None:
                logging.info(f"     - preprocessing dataset {watergang}")
                waterline = self.generate_or_use_preprocessed_hydroobject(
                    waterline=watergang,
                    snapping_distance=snapping_distance if watergang == "hydroobject" else None
                )
                setattr(self, watergang, waterline)
            else:
                logging.info(
                    f"     - use processed dataset {watergang}: {watergang_processed_file_name.name}"
                )
                setattr(self, watergang, gpd.read_file(watergang_processed_file_name))


    def generate_or_use_preprocessed_hydroobject(
        self, waterline, preprocessed_file="preprocessed", snapping_distance=0.05
    ):
        files_in_dir = self.dir_results.glob("**/*")
        waterline_preprocessed_file = Path(self.dir_results, f"{waterline}_{preprocessed_file}.gpkg")
        
        if waterline_preprocessed_file in files_in_dir:
            logging.info(f"     - load dataset preprocessed {waterline}")
            gdf_waterline = gpd.read_file(waterline_preprocessed_file)

            return gdf_waterline

        else:
            logging.info(
                f"     - no {waterline}_preprocessed.gpkg, preprocessing {waterline}"
            )
            gdf_waterline = getattr(self, waterline)
            len_gdf_waterline = len(gdf_waterline)
            if snapping_distance is not None:
                gdf_waterline, gdf_waterline_snapped = preprocess_hydroobject(
                    gdf_waterline, snapping_distance=snapping_distance
                )
                if self.write_results:
                    gdf_waterline_snapped.to_file(Path(self.dir_results, f"{waterline}_snapped.gpkg"))
                    gdf_waterline.to_file(Path(self.dir_results, f"{waterline}_preprocessed.gpkg"))
                logging.info(f"     - preprocessing done: {waterline}")
            else:
                logging.info(f"     - no preprocessing: {waterline}")

            # check for invalid or duplicate geometries (linestrings forming a ring)
            gdf_waterline_old = gdf_waterline.copy()
            gdf_waterline = gdf_waterline[~gdf_waterline['geometry'].apply(lambda geom: geom.is_closed)]
            gdf_waterline = gdf_waterline.loc[~gdf_waterline["geometry"].duplicated(keep="first")]
            logging.info(f"     - removed {len_gdf_waterline-len(gdf_waterline)} waterlines [{waterline}]")

            return gdf_waterline
        

    def create_graph_from_network(self, processed="processed"):
        """Turns a linestring layer containing waterlines into a graph of edges and nodes. 

        Returns
        -------
        self.nodes: gpd.GeoDataFrame
            Geodataframe containing nodes between waterlines
        self.edges: gpd.GeoDataFrame
            Geodataframe containing edges (waterlines)
        self.graph: nx.DiGraph
            Networkx graph containing the edges and nodes
        """
        edges = None

        if self.edges_hydro is not None:
            self.edges_hydro["source"] = "hydroobject"
            if self.edges_overig is not None:
                edges = pd.concat([
                    self.edges_hydro,
                    self.edges_overig
                ])
            else:
                edges = self.edges_hydro.copy()
        else:
            for line_type in ["hydroobject", "overige_watergang"]:
                gdf = getattr(self, line_type)
                for i in range(10):
                    if not hasattr(self, f"{line_type}_{processed}_{i}"):
                        break
                    gdf_processed = getattr(self, f"{line_type}_{processed}_{i}")
                    if gdf_processed is None:
                        break
                    else:
                        gdf = gdf_processed.copy()
                if gdf is None:
                    continue
                gdf["source"] = line_type
                if edges is None:
                    edges = gdf.copy()
                else:
                    edges = pd.concat([edges, gdf.explode()])

        for abbr, line_type in zip(["hydro", "overig"], ["hydroobject", "overige_watergang"]):
            edges_abbr = edges[edges["source"] == line_type]
            if edges_abbr.empty:
                setattr(self, f"edges_{abbr}", None)
            else:
                setattr(self, f"edges_{abbr}", edges_abbr)

        # only hydroobject
        self.nodes_hydro, self.edges_hydro, self.graph_hydro = create_graph_from_edges(self.edges_hydro)
        self.nodes_hydro, self.edges_hydro = analyse_netwerk_add_information_to_nodes_edges(
            self.edges_hydro,
            self.nodes_hydro, 
            min_difference_angle=10.0, 
            min_difference_discharge_factor=2.0
        )

        # all edges
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        self.nodes, self.edges = analyse_netwerk_add_information_to_nodes_edges(
            self.edges,
            self.nodes, 
            min_difference_angle=10.0, 
            min_difference_discharge_factor=2.0
        )

        logging.info(
            f"   x create network graph ({len(self.edges)} edges, {len(self.nodes)} nodes)"
        )
        return self.nodes, self.edges, self.graph


    def export_results_to_gpkg_or_nc(self, list_layers: list[str] = None, dir_output: str | Path = None):
        """Export results to geopackages in folder 1_resultaat"""
        if dir_output is None:
            dir_output = self.dir_results

        logging.info(f"   x export results")
        if list_layers is None:
            return
        for layer in list_layers:
            result = getattr(self, layer)
            if result is None:
                logging.info(f"     - {layer} not available")
            elif isinstance(result, gpd.GeoDataFrame):
                logging.info(f"     - {layer} ({len(result)})")
                result.to_file(Path(dir_output, f"{layer}.gpkg"))
            elif isinstance(result, xr.DataArray) or isinstance(result, xr.Dataset):
                logging.info(f"     - {layer} (netcdf)")
                netcdf_file_path = Path(dir_output, f"{layer}.nc")
                if netcdf_file_path.exists():
                   netcdf_file_path.unlink()
                encoding = {
                    layer: {
                        'dtype': str(result.dtype),
                        'zlib': True,
                        'complevel': 9,
                    },
                }
                result.to_netcdf(netcdf_file_path, mode='w', encoding=encoding)
            else:
                raise ValueError("type not exportable")
