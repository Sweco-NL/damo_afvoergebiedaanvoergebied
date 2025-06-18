import pyflwdir
import logging
import xarray as xr
import numpy as np
import time


def run_pyflwdir(dem: xr.Dataset, waterways: xr.Dataset, iterations: int = 2000, iteration_group: int = 100) -> xr.Dataset:
    # create pyflwdir object
    logging.info("     - create pyflwdir object to calculate downstream direction")

    flw_pyflwdir = pyflwdir.from_dem(
        data=dem.data[0],
        nodata=dem._FillValue,
        transform=dem.rio.transform(),
        latlon=False,
    )

    # get upstream values
    def get_upstream_values(
        flw_mask: np.ndarray, 
        flw_idxs_ds: np.ndarray, 
        waterways_flat: np.ndarray, 
        iterations: int,
        iteration_start: int
    ):
        for i in range(iterations):
            time_start = time.time()
            upstream_values = waterways_flat.copy()
            upstream_values[flw_mask] = upstream_values[flw_idxs_ds[flw_mask]]
            
            new_filled_cells = (waterways_flat == -1.0) & (upstream_values != -1.0)
            waterways_flat = np.where(new_filled_cells, upstream_values, waterways_flat)

            number_new_filled_cells = new_filled_cells.sum()

            if number_new_filled_cells == 0:
                print("     * break at iteration: ", i + iteration_start)
                break
            print1 = f"  * iteration: {i + iteration_start}"
            print2 = f" | number new cells: {number_new_filled_cells}"
            print3 = f"({round(time.time()-time_start, 2)} seconds)"
            print(print1 + print2 + print3, end="\r")
        return waterways_flat, number_new_filled_cells
    
    logging.info(f"     - get upstream area of each waterway: {iterations} iterations")
    waterways_flat_new = flw_pyflwdir._check_data(waterways.data, "data")

    time_start_groups = time.time()
    for i in range(0, iterations, iteration_group):
        time_start_group = time.time()
        waterways_flat_new, number_new_filled_cells = get_upstream_values(
            flw_mask=flw_pyflwdir.mask, 
            flw_idxs_ds=flw_pyflwdir.idxs_ds, 
            waterways_flat=waterways_flat_new, 
            iterations=min(iterations-i, iteration_group),
            iteration_start=i
        )
        print("")
        print(f"iteration ({i+iteration_group}/{iterations}): {round(time.time()-time_start_group, 2)}s/{round(time.time()-time_start_groups, 2)}s")
        if number_new_filled_cells == 0:
            break

    waterways.data = waterways_flat_new.reshape(
        waterways.data.shape
    )
    return waterways, flw_pyflwdir
