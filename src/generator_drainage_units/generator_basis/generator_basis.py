import logging
from pathlib import Path

import folium
import geopandas as gpd
import rioxarray
from pydantic import BaseModel, ConfigDict


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
            # self.use_processed_hydroobjecten(processed_file="processed")


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
                if hasattr(self, f.stem) and getattr(self, f.stem) is None:
                    logging.info(f"     - get dataset {f.stem.upper()}")
                    if f.suffix == ".gpkg":
                        setattr(self, f.stem, gpd.read_file(f, layer=f.stem))
                    if f.suffix in [".nc", ".NC"]:
                        setattr(self, f.stem, rioxarray.open_rasterio(f))
            if getattr(self, required_dataset) is None:
                raise ValueError(f" * dataset {required_dataset} is missing")


    def use_processed_hydroobjecten(self, processed_file="processed"):
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

            watergang_processed_file = None
            attributes = dir(self)
            files_in_dirs = list(self.dir_basisdata.glob("**/*")) + list(self.dir_results.glob("**/*"))

            for file in files_in_dirs:
                if (
                    f"{watergang}_{processed_file}" in file.stem
                    and "nodes" not in file.stem
                    and file.stem in self.required_results
                ):
                    watergang_processed_file = file

            if watergang_processed_file is not None:
                logging.info(
                    f"     - use dataset processed {watergang}: {watergang_processed_file.name}"
                )
                setattr(self, watergang, gpd.read_file(watergang_processed_file))

