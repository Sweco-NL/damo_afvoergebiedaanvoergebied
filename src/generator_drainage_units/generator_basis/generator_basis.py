import logging
from pathlib import Path

import folium
import geopandas as gpd
import rioxarray
import xarray
from pydantic import BaseModel, ConfigDict
from ..utils.preprocess import preprocess_hydroobjecten
from ..utils.create_graph import create_graph_from_edges
from ..utils.network_functions import (
    calculate_angles_of_edges_at_nodes,
    define_list_upstream_downstream_edges_ids,
    find_node_edge_ids_in_directed_graph,
    select_downstream_upstream_edges,
)


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
        logging.info(f' ### Case "{self.name.capitalize()}" ###')

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
                    logging.info(f"     - get dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        with rioxarray.open_rasterio(f) as raster:
                            setattr(self, f.stem, raster.load())

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
                    logging.info(f"     - get dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        setattr(self, f.stem, rioxarray.open_rasterio(f))
            if getattr(self, required_dataset) is None:
                logging.info(f" * dataset {required_dataset} is missing - check if absolutely required")


    def use_processed_hydroobjecten(self, processed_file="processed", force_preprocess=False):
        """actualize hydroobjecten and overige_watergangen

        replaces hydroobjecten and overige_watergangen with the newest processed attributes

        Parameters
        ----------
        processed_file : str, optional
            suffix of processed files, by default "processed"
        """
        for watergang in ["hydroobjecten", "overige_watergangen"]:
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
                if watergang == "hydroobjecten":
                    logging.info(f"     - preprocessing dataset {watergang}")
                    waterline = self.generate_or_use_preprocessed_hydroobjecten(waterline=watergang)
                    setattr(self, watergang, waterline[0])
            else:
                logging.info(
                    f"     - use processed dataset {watergang}: {watergang_processed_file_name.name}"
                )
                setattr(self, watergang, gpd.read_file(watergang_processed_file_name))


    def generate_or_use_preprocessed_hydroobjecten(
        self, waterline, preprocessed_file="preprocessed"
    ):
        files_in_dir = self.dir_results.glob("**/*")
        waterline_preprocessed_file = Path(self.dir_results, f"{waterline}_{preprocessed_file}.gpkg")
        
        if waterline_preprocessed_file in files_in_dir:
            logging.info(f"     - get dataset preprocessed {waterline}")
            gdf_waterline = gpd.read_file(waterline_preprocessed_file)
            return gdf_waterline, gdf_waterline.copy()

        else:
            logging.info(
                f"     - {waterline}_preprocessed.gpkg not in directory, preprocessing hydroobjecten"
            )
            gdf_waterline, gdf_waterline_snapped = preprocess_hydroobjecten(
                getattr(self, waterline)
            )
            if self.write_results:
                gdf_waterline_snapped.to_file(Path(self.dir_results, f"{waterline}_snapped.gpkg"))
                gdf_waterline.to_file(Path(self.dir_results, f"{waterline}_preprocessed.gpkg"))
            logging.info(f"     - preprocessing {waterline}: done")
            return gdf_waterline, gdf_waterline_snapped
        

    def create_graph_from_network(self, water_lines=["hydroobjecten"], processed="processed"):
        """Turns a linestring layer containing waterlines into a graph of edges and nodes. 

        Parameters
        ----------
        water_lines : list, optional
            List of waterline files names used to create graph, must refer to geopackages containing linestrings, by default ["hydroobjecten"]

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
        for water_line in water_lines:
            gdf_water_line = getattr(self, water_line)
            for i in range(10):
                if not hasattr(self, f"{water_line}_{processed}_{i}"):
                    break
                gdf_water_line_processed = getattr(self, f"{water_line}_{processed}_{i}")
                if gdf_water_line_processed is None:
                    break
                else:
                    gdf_water_line = gdf_water_line_processed.copy()
            if gdf_water_line is None:
                continue
            if edges is None:
                edges = gdf_water_line.explode()
            else:
                edges = pd.concat([edges, gdf_water_line.explode()])
        self.nodes, self.edges, self.graph = create_graph_from_edges(edges)
        logging.info(
            f"   x create network graph ({len(self.edges)} edges, {len(self.nodes)} nodes)"
        )
        return self.nodes, self.edges, self.graph


    def analyse_netwerk_add_information_to_nodes_edges(self, min_difference_angle=20.0):
        self.nodes = define_list_upstream_downstream_edges_ids(
            node_ids=self.nodes.nodeID.values, nodes=self.nodes, edges=self.edges
        )
        logging.info("   x calculate angles of edges to nodes")
        self.nodes, self.edges = calculate_angles_of_edges_at_nodes(
            nodes=self.nodes, edges=self.edges
        )
        logging.info("   x find downstream upstream edges")
        self.nodes = select_downstream_upstream_edges(
            self.nodes, min_difference_angle=min_difference_angle
        )
        return self.nodes, self.edges


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
            elif isinstance(result, xarray.DataArray) or isinstance(result, xarray.Dataset):
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
