import xarray as xr


def calculate_med_water(ds):
    # ds must be theta >= 5
    # Process Mediterranean Water influence
    mw_mask = (
            (ds.salinity - 35 > (0.8 / 9) * (ds.theta - 4)) &
            (ds.theta < 12.2) &
            (ds.latitude > 15) &
            (ds.pressure > 400)
    )

    # Calculate Mediterranean Water parameters where where mask is True
    where_mw = mw_mask.where(mw_mask)
    if where_mw.any():
        ds['salinityCorrectedByMW'] = xr.where(
            mw_mask,
            (4 - 11.74 + 36.5 * (ds.theta - 11.74) / (ds.salinity - 36.5) - 34.9 * 9 / 0.8) /
            ((ds.theta - 11.74) / (ds.salinity - 36.5) - 9 / 0.8),
            ds.salinity
        )

        ds['percentMW'] = xr.where(
            mw_mask,
            1 - (ds.salinity - 36.5) / (ds['salinityCorrectedByMW'] - 36.5),
            0
        )

        # Update corrected parameters
        ds['thetaCorrectedByMW'] = xr.where(
            mw_mask,
            (ds.theta - 11.74) / (1 - ds.percentMW) + 11.74,
            ds.theta
        )
        ds['noConservativeCorrectedByMW'] = xr.where(
            mw_mask,
            (ds.no_conservative - 307) / (1 - ds.percentMW) + 307,
            ds.no_conservative
        )
        ds['poConservativeCorrectedByMW'] = xr.where(
            mw_mask,
            (ds.po_conservative - 304) / (1 - ds.percentMW) + 304,
            ds.po_conservative
        )
    