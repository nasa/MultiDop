import numpy as np
import datetime as dt
import xarray
from copy import deepcopy


def make_new_grid(grid_list, filen):
    """
    This function takes the DDA input grids and output analysis file, and
    creates a new CF- and Py-ART-compliant Grid object. Wind information
    outside the mutual radar coverage region is masked.

    Parameters
    ----------
    grid_list : 2-3 element list of Py-ART Grid objects
        The Py-ART Grids are the original input grids to the DDA engine.
    filen : str
        Name of DDA output file.

    Returns
    -------
    new_grid : Py-ART Grid object
        CF- and Py-ART-compliant Grid object
    """
    new_grid = deepcopy(grid_list[0])

    if len(grid_list) == 2:
        new_grid.radar_latitude['data'] = np.array(
            [grid_list[0].radar_latitude['data'][0],
             grid_list[1].radar_latitude['data'][0]])
        new_grid.radar_longitude['data'] = np.array(
            [grid_list[0].radar_longitude['data'][0],
             grid_list[1].radar_longitude['data'][0]])
        new_grid.radar_altitude['data'] = np.array(
            [grid_list[0].radar_altitude['data'][0],
             grid_list[1].radar_altitude['data'][0]])
        new_grid.radar_name['data'] = np.array(
            [''.join(grid_list[0].radar_name['data'][0]),
             ''.join(grid_list[1].radar_name['data'][0])])
        dt1 = dt.datetime.strptime(
            grid_list[0].radar_time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
        dt2 = dt.datetime.strptime(
            grid_list[1].radar_time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
        dtlist = [dt1, dt2]
        imin = np.argmin(dtlist)
        new_grid.radar_time['data'] = np.array(
            [(dtlist[0]-dtlist[imin]).total_seconds(),
             (dtlist[1]-dtlist[imin]).total_seconds()])
        new_grid.radar_time['units'] = grid_list[imin].radar_time['units']
        new_grid.nradar = 2

    elif len(grid_list) == 3:
        new_grid.radar_latitude['data'] = np.array(
            [grid_list[0].radar_latitude['data'][0],
             grid_list[1].radar_latitude['data'][0],
             grid_list[2].radar_latitude['data'][0]])
        new_grid.radar_longitude['data'] = np.array(
            [grid_list[0].radar_longitude['data'][0],
             grid_list[1].radar_longitude['data'][0],
             grid_list[2].radar_longitude['data'][0]])
        new_grid.radar_altitude['data'] = np.array(
            [grid_list[0].radar_altitude['data'][0],
             grid_list[1].radar_altitude['data'][0],
             grid_list[2].radar_altitude['data'][0]])
        new_grid.radar_name['data'] = np.array(
            [''.join(grid_list[0].radar_name['data'][0]),
             ''.join(grid_list[1].radar_name['data'][0]),
             ''.join(grid_list[2].radar_name['data'][0])])
        dt1 = dt.datetime.strptime(
            grid_list[0].radar_time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
        dt2 = dt.datetime.strptime(
            grid_list[1].radar_time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
        dt3 = dt.datetime.strptime(
            grid_list[2].radar_time['units'][14:], '%Y-%m-%dT%H:%M:%SZ')
        dtlist = [dt1, dt2, dt3]
        imin = np.argmin(dtlist)
        new_grid.radar_time['data'] = np.array(
            [(dtlist[0]-dtlist[imin]).total_seconds(),
             (dtlist[1]-dtlist[imin]).total_seconds(),
             (dtlist[2]-dtlist[imin]).total_seconds()])
        new_grid.radar_time['units'] = grid_list[imin].radar_time['units']
        new_grid.nradar = 3

    else:
        raise ValueError('Expecting 2 or 3 Doppler radars!')

    # Add new fields to grid, and remove old ones
    orig_keylist = [key for key in grid_list[0].fields.keys()]
    ds = xarray.open_dataset(filen)
    dz = np.array(ds.MAXDBZ.T)
    dz = np.ma.masked_where(np.logical_or(dz <= 0, dz > 100), dz)
    cvg = np.array(ds.CVG.T)
    u = np.ma.asanyarray(ds.U.T)
    u.mask = np.logical_or(dz.mask, cvg != 1)
    v = np.ma.asanyarray(ds.V.T)
    v.mask = np.logical_or(dz.mask, cvg != 1)
    w = np.ma.asanyarray(ds.W.T)
    w.mask = np.logical_or(dz.mask, cvg != 1)
    new_grid = _add_field_to_grid(
        dz, new_grid, 'reflectivity', 'dBZ', 'Merged Reflectivity from Radars',
        'Merged Reflectivity', -32768)
    new_grid = _add_field_to_grid(
        u, new_grid, 'eastward_wind', 'meters_per_second',
        'Eastward Component of the Wind', 'Eastward Wind', -32768)
    new_grid = _add_field_to_grid(
        v, new_grid, 'northward_wind', 'meters_per_second',
        'Northward Component of the Wind', 'Northward Wind', -32768)
    new_grid = _add_field_to_grid(
        w, new_grid, 'upward_air_velocity', 'meters_per_second',
        'Vertical Component of the Wind', 'Vertical Wind', -32768)
    for key in orig_keylist:
        del(new_grid.fields[key])

    return new_grid


def _add_field_to_grid(
        field, grid, field_name, units, long_name, standard_name, fill_value):
    """
    Adds a field to the Py-ART Grid object.
    """
    field_dict = {'data': field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    grid.add_field(field_name, field_dict, replace_existing=True)
    return grid
