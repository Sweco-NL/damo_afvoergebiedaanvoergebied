import pyflwdir
import logging
import xarray as xr
import numpy as np
import time
import logging
import numba as nb


def define_flowdirection_raster_d16(raster):
    list_dxdy_factorxfactory = [
        [1, 1, 1/2**0.5, -1/2**0.5],
        [0, 1, 0.0, -1.0],
        [-1, 1, -1/2**0.5, -1/2**0.5],
        [-1, 0, -1.0, 0.0],
        [-1, -1, -1/2**0.5, 1/2**0.5],
        [0, -1, 0.0, 1.0],
        [1, -1, 1/2**0.5, 1/2**0.5],
        [1, 0, 1.0, 0.0],
    ]

    flow = raster.copy()
    flow.data[flow.data<-10.0] = np.nan
    flow_x = None
    flow_y = None

    for [dx, dy, factorx, factory] in list_dxdy_factorxfactory:
        dflow = flow.shift(x=dx, y=dy)
        if flow_x is None:
            flow_x = dflow * factorx
        else:
            flow_x += dflow * factorx
        if flow_y is None:
            flow_y = dflow * factory
        else:
            flow_y += dflow * factory

    flow_direction = flow.copy()
    flow_direction.data = np.arctan2(flow_x.data, flow_y.data) * 180 / np.pi
    flow_direction.data = flow_direction.data % 360
    return flow_direction, flow_x, flow_y            


def find_flow_direction_indices(flow_direction, level=1):
    if level == 1:
        min_angles = [360.0 / 8.0 * (-1.5 + float(i)) for i in range(8)]
        max_angles = [360.0 / 8.0 * (-0.5 + float(i)) for i in range(8)]
        dxs = [-1, 0, 1, 1, 1, 0, -1, -1]
        dys = [-1, -1, -1, 0, 1, 1, 1, 0]
        d_infinities = [int(i)+1 for i in range(8)]
    elif level == 2:
        min_angles = [360.0 / 16.0 * (-2.5 + float(i)) for i in range(16)]
        max_angles = [360.0 / 16.0 * (-1.5 + float(i)) for i in range(16)]
        dxs = [-2, -1, 0, 1, 2, 2, 2, 2, 2, 1, 0, -1, -2, -2, -2, -2]
        dys = [-2, -2, -2, -2, -2, -1, 0, 1, 2, 2, 2, 2, 2, 1, 0, -1]
        d_infinities = [int(i)+1 for i in range(16)]
    else:
        raise ValueError(f"   x level {level} not implemented")

    dindices = [dy*len(flow_direction.x) + dx for dx, dy in zip(dxs, dys)]

    flow_direction_dinf = flow_direction.copy()
    flow_direction_ind = flow_direction.copy()

    for min_angle, max_angle, dindex, dinf in zip(min_angles, max_angles, dindices, d_infinities):
        mask = (flow_direction.data >= min_angle) & (flow_direction.data < max_angle)
        flow_direction_ind.data[mask] = dindex
        flow_direction_dinf.data[mask] = dinf
    
    return flow_direction_dinf.fillna(-9999).astype(np.int32), flow_direction_ind.fillna(-9999).astype(np.int32)


# get upstream values
@nb.njit()
def get_upstream_values_d16(
    flw_mask, 
    flw_idxs_ds, 
    drainage_units_flat, 
    iterations,
    iteration_start
):
    for i in range(iterations):
        # time_start = time.time()
        upstream_values = drainage_units_flat.copy()
        upstream_values[flw_mask] = upstream_values[flw_idxs_ds[flw_mask]]
        
        new_filled_cells = (drainage_units_flat == -1) & (upstream_values != -1)
        drainage_units_flat = np.where(new_filled_cells, upstream_values, drainage_units_flat)

        number_new_filled_cells = new_filled_cells.sum()

        if number_new_filled_cells == 0:
            # print("     * break at iteration: ", i + iteration_start)
            break
        # print1 = f"  * iteration: {i + iteration_start}"
        # print2 = f" | number new cells: {number_new_filled_cells}"
        # print3 = f"({round(time.time()-time_start, 2)} seconds)"
        # logging.debug(print1 + print2 + print3, end="\r")
    return drainage_units_flat, number_new_filled_cells


# get upstream values
def get_upstream_values_d8(
    flw_mask: np.ndarray, 
    flw_idxs_ds: np.ndarray, 
    waterways_flat: np.ndarray, 
    iterations: int,
    iteration_start: int
):
    for i in range(iterations):
        # time_start = time.time()
        upstream_values = waterways_flat.copy()
        upstream_values[flw_mask] = upstream_values[flw_idxs_ds[flw_mask]]
        
        new_filled_cells = (waterways_flat == -1) & (upstream_values != -1)
        waterways_flat = np.where(new_filled_cells, upstream_values, waterways_flat)

        number_new_filled_cells = new_filled_cells.sum()

        if number_new_filled_cells == 0:
            # print("     * break at iteration: ", i + iteration_start)
            break
        # print1 = f"  * iteration: {i + iteration_start}"
        # print2 = f" | number new cells: {number_new_filled_cells}"
        # print3 = f"({round(time.time()-time_start, 2)} seconds)"
        # logging.debug(print1 + print2 + print3, end="\r")
    return waterways_flat, number_new_filled_cells


def run_pyflwdir(
    dem: xr.Dataset, 
    waterways: xr.Dataset, 
    iterations: int = 2000, 
    iteration_group: int = 100, 
    flow_method="d8", 
    flow_direction_d16_fills=None
) -> xr.Dataset:

    if flow_method == "d8":
        iterations_d8 = iterations
    else:
        iterations_d8 = 2
    
    # create pyflwdir object
    logging.info("     - create pyflwdir object to calculate downstream direction (d8)")
    flw_pyflwdir = pyflwdir.from_dem(
        data=dem.data[0],
        nodata=dem._FillValue,
        transform=dem.rio.transform(),
        latlon=False,
    )
    
    logging.info(f"     - get upstream area of each waterway: {iterations_d8} iterations (d8)")
    waterways_flat_new = flw_pyflwdir._check_data(waterways.data, "data").astype(np.int32)

    time_start_groups = time.time()
    for i in range(0, iterations_d8, iteration_group):
        time_start_group = time.time()
        waterways_flat_new, number_new_filled_cells = get_upstream_values_d8(
            flw_mask=flw_pyflwdir.mask, 
            flw_idxs_ds=flw_pyflwdir.idxs_ds, 
            waterways_flat=waterways_flat_new, 
            iterations=min(iterations_d8-i, iteration_group),
            iteration_start=i
        )
        logging.info(f"       o iteration ({i+iteration_group}/{iterations_d8}): {round(time.time()-time_start_group, 2)}s/{round(time.time()-time_start_groups, 2)}s [{number_new_filled_cells} new cells]")
        if number_new_filled_cells == 0:
            break
    
    flow_direction = flw_pyflwdir.idxs_ds.copy()

    if flow_method == "d16":
        # create dinf object
        logging.info("     - create d-infinity object to calculate downstream direction")

        # use fill depressions from pyflwdir
        dem.data[0] = pyflwdir.fill_depressions(
            dem.data[0],
            nodata=dem._FillValue,
        )[0]
        dem.data[dem.data<-9999.0] = np.nan

        flow_direction, _, _ = define_flowdirection_raster_d16(dem)
        
        # flow_direction.data[flow_direction > 360 - 360.0 / 8.0 * 1.5] = \
        #     flow_direction.data[flow_direction > 360 - 360.0 / 8.0 * 1.5] - 360.0
        # flow_direction_d8, flow_direction_d8_ind = find_flow_direction_indices(flow_direction, level=1)
        flow_direction.data[flow_direction > 360 - 360.0 / 16.0 * 2.5] = \
            flow_direction.data[flow_direction > 360 - 360.0 / 16.0 * 2.5] - 360.0
        
        flow_direction_d16_fills_deg = flow_direction_d16_fills.copy()
        for i, j in zip(range(16), np.arange(-45.0, 315.0, 22.5)):
            flow_direction_d16_fills_deg.data[flow_direction_d16_fills.data==i+1] = j%360

        flow_direction = flow_direction_d16_fills_deg.where(flow_direction_d16_fills_deg > -1.0, flow_direction)
        flow_direction_d16, flow_direction_d16_ind = find_flow_direction_indices(flow_direction, level=2)
        
        # waterways_flat_d8_index = np.arange(waterways_flat_new.size) + flow_direction_d8_ind.data.ravel()
        # waterways_flat_d8_index[np.isnan(waterways_flat_d8_index)] = -1.0
        # waterways_flat_d8_index = waterways_flat_d8_index.astype(np.int32)

        waterways_flat_d16_index = np.arange(waterways_flat_new.size) + flow_direction_d16_ind.data.ravel()
        waterways_flat_d16_index[np.isnan(waterways_flat_d16_index)] = -1.0
        waterways_flat_d16_index = waterways_flat_d16_index.astype(np.int32)
        
        logging.info(f"     - get upstream area of each waterway: {iterations} iterations (d16)")

        # time_start_groups = time.time()
        for i in range(0, iterations, iteration_group):
            time_start_group = time.time()
            for j in range(iteration_group):
                waterways_flat_new, number_new_filled_cells = get_upstream_values_d16(
                    flw_mask=(dem.data[0] > -100.0).ravel(), 
                    flw_idxs_ds=waterways_flat_d16_index, 
                    drainage_units_flat=waterways_flat_new, 
                    iterations=1,
                    iteration_start=0
                )
                if number_new_filled_cells == 0:
                    # print("     * break at iteration: ", i + j)
                    break
                # print1 = f"  * iteration: {i + j}"
                # print2 = f" | number new cells: {number_new_filled_cells}"
                # print3 = f"({round(time.time()-time_start_group, 2)} seconds)"
                # print(print1 + print2 + print3, end="\r")
                # logging.info(f"       o iteration {i+j}/{iterations} [{number_new_filled_cells} new cells]")
            logging.info(f"       o iteration {i+iteration_group}/{iterations}: {round(time.time()-time_start_group, 2)}s/{round(time.time()-time_start_groups, 2)}s [{number_new_filled_cells} new cells]")
            if number_new_filled_cells == 0:
                break
        
        flow_direction = flow_direction_d16.copy()

    waterways.data = waterways_flat_new.reshape(
        waterways.data.shape
    )
    return waterways, dem, flow_direction
