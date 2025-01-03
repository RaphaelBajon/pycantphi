import xarray as xr


def calculate_delta_cdiseq(theta, salinity, percent_mw, no_conservative, po_conservative, si_teo, lat):

    if percent_mw == 0:
        if theta < 5:
            delta_cdiseq = (
                    -15 + 3.44 * (theta-10) + 20.37*(salinity-35) +
                    0.109*(no_conservative-300) +
                    0.139*(po_conservative-300) -
                    0.406*si_teo
            )
        elif theta < 8:
            dummy_correction = 0
            if 0 < lat < 50:
                dummy_correction = 0.72 * (theta-10)
            delta_cdiseq = (
                    dummy_correction - 5.1 - 0.75*1.72*(theta-10) - 2 + 0.75 * 14.8 * (salinity-35)
            )
            # Anomalous DCdis higher than Zero for T<8C
            if delta_cdiseq > 0:
                delta_cdiseq = 0
        elif theta < 18:
            if lat > 20:
                delta_cdiseq = (
                        -13.4 + 1.56*(theta-10) + 0.166*(no_conservative-300))
            elif lat < -20:
                delta_cdiseq = (
                        4-14.7
                        - 1.82*(theta-10)
                        + 9.68*(salinity-35)
                        - 0.117*(no_conservative-300)
                        + 0.134*(po_conservative-300)
                )
            else:  # lats between -20 and 20
                delta_cdiseq = -4.4 + 0.9*(theta-10)
        else:  # temperatures >= 18C
            if -20 < lat < 20:  # lats between -20 and 20
                delta_cdiseq = (
                        -3.93 + 2.55*(theta-10) - 11.3*(salinity-35)
                )
            else:
                delta_cdiseq = (
                        -38.6 + 1.67*(theta-10)
                        + 16.3*(salinity-35)
                        - 0.32*(no_conservative-300)
                        + 0.524*(po_conservative-300)
                )
    else:
        delta_cdiseq = - 6 * (1 - percent_mw) + (-11) * percent_mw

    return delta_cdiseq


def calculate_delta_cdiseq_dataset(ds, condition_theta, northern_sills):
    """
    Apply delta CDiseq calculation to a xarray.Dataset.

    Parameters:
    - ds (xarray.Dataset): Dataset containing variables:
        - `thetaCorrectedByMW`
        - `salinityCorrectedByMW`
        - `percentMW`
        - `noConservativeCorrectedByMW`
        - `poConservativeCorrectedByMW`
        - `silicate`
        - `latitude`

    Returns:
    - xarray.DataArray: Delta CDiseq values.
    """
    delta_cdiseq_if = xr.where(
        condition_theta & northern_sills,
        2.5 * (ds["thetaCorrectedByMW"] - 5) - 6,
        ds.deltaCDiseq
    ).rename('deltaCDiseq')

    ds_prime = ds.where(condition_theta & ~northern_sills)

    delta_cdiseq_else = xr.apply_ufunc(
        calculate_delta_cdiseq,
        ds_prime["thetaCorrectedByMW"], ds_prime["salinityCorrectedByMW"], ds_prime["percentMW"],
        ds_prime["noConservativeCorrectedByMW"], ds_prime["poConservativeCorrectedByMW"],
        ds_prime["silicate_teorical"], ds_prime["latitude"],
        input_core_dims=[[], [], [], [], [], [], []],  # No dimension reduction
        vectorize=True,  # Apply element-wise
        dask="parallelized",  # Support for dask arrays
        output_dtypes=[float],  # Specify output type
    ).rename('deltaCDiseq')

    delta_cdiseq = xr.merge([delta_cdiseq_if, delta_cdiseq_else])

    return delta_cdiseq['deltaCDiseq']
