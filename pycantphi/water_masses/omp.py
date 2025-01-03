import numpy as np
import xarray as xr
from scipy.optimize import nnls
import gsw


def calculate_omp_single_sample(
        salinity, theta, pressure, latitude, nitrate_preformed, phosphate_preformed, silicate,
        sal_types, tpot_types, no_preformed_types, po_preformed_types, sio2_types,
        talk_preformed_types, cdis_types, params_means, params_stds, tipos_values_nor
):
    """
    Perform OMP analysis for a single sample.
    """
    if (
            np.isnan(salinity) or np.isnan(theta) or np.isnan(pressure) or
            np.isnan(latitude) or np.isnan(nitrate_preformed) or
            np.isnan(phosphate_preformed) or np.isnan(silicate)
    ):
        return None, None, None, np.array([None] * 6)

    # Calculate density (sigma)
    sigma = gsw.density.rho(salinity, theta, pressure) - 1000

    # Calculate density at sample depth for water masses
    sigma_wm = gsw.density.rho(sal_types, tpot_types, pressure) - 1000

    # Pressure factor
    Mprs = 0.01 * (-0.64781 + 0.64781e-4 * pressure + 0.0162e-84 * pressure ** 2)

    # Initialize xiso array
    xiso = np.zeros(10)

    # First pair (N5-NF)
    if sigma < sigma_wm[0]:
        xiso[1] = 0
    elif sigma > sigma_wm[1]:
        xiso[1] = 1
    else:
        a = Mprs * (tpot_types[1] - tpot_types[0]) ** 2
        b = sigma_wm[1] - sigma_wm[0] - a
        c = sigma - sigma_wm[0]
        xiso[1] = (np.sqrt(b ** 2 + 4 * a * c) - b) / (2 * a)
    xiso[0] = 1 - xiso[1]

    # Second pair (S5-SF)
    if sigma < sigma_wm[2]:
        xiso[3] = 0
    elif sigma > sigma_wm[3]:
        xiso[3] = 1
    else:
        a = Mprs * (tpot_types[3] - tpot_types[2]) ** 2
        b = sigma_wm[3] - sigma_wm[2] - a
        c = sigma - sigma_wm[2]
        xiso[3] = (np.sqrt(b ** 2 + 4 * a * c) - b) / (2 * a)
    xiso[2] = 1 - xiso[3]

    # Third pair (CDW-CDW2)
    if sigma < sigma_wm[4]:
        xiso[5] = 0
    elif sigma > sigma_wm[5]:
        xiso[5] = 1
    else:
        a = Mprs * (tpot_types[5] - tpot_types[4]) ** 2
        b = sigma_wm[5] - sigma_wm[4] - a
        c = sigma - sigma_wm[4]
        xiso[5] = (np.sqrt(b ** 2 + 4 * a * c) - b) / (2 * a)
    xiso[4] = 1 - xiso[5]

    # Remaining xiso values
    xiso[7] = 0
    xiso[6] = 1 - xiso[7]
    xiso[9] = 0
    xiso[8] = 1 - xiso[9]

    # Valtip matrix (5x5)
    Valtip = np.zeros((5, 5))
    for i, types in enumerate([tpot_types, sal_types, no_preformed_types, po_preformed_types, sio2_types]):
        Valtip[0, i] = types[0] + (types[1] - types[0]) * xiso[1]
        Valtip[1, i] = types[2] + (types[3] - types[2]) * xiso[3]
        Valtip[2, i] = types[4] + (types[5] - types[4]) * xiso[5]
        Valtip[3, i] = types[2] + (types[0] - types[2]) * xiso[7]
        Valtip[4, i] = types[2] + (types[4] - types[2]) * xiso[9]

    # Normalize
    Mediso = np.mean(Valtip, axis=0)
    STDiso = np.std(Valtip, axis=0, ddof=1)
    ValtipNor = (Valtip - Mediso) / STDiso

    # Prepare optimization variables
    b = np.array([theta, salinity, nitrate_preformed, phosphate_preformed, silicate])
    bNor = (b - Mediso) / STDiso
    bNorO = np.concatenate(([1], bNor))
    TiposSelNor = np.vstack((np.ones(len(ValtipNor)), ValtipNor.T))

    # Solve optimization
    weights = np.array([100, 8, 3, 1, 1, 1])
    W = np.diag(weights)
    TiposSelNorW = W @ TiposSelNor
    bNorW = W @ bNorO
    Xabc = nnls(TiposSelNorW, bNorW)[0]

    # Calculate mixing ratios
    X = np.zeros(6)
    X[0] = Xabc[0] * xiso[0] + Xabc[3] * xiso[7]
    X[1] = Xabc[0] * xiso[1]
    X[2] = Xabc[1] * xiso[2] + Xabc[3] * xiso[6] + Xabc[4] * xiso[8]
    X[3] = Xabc[1] * xiso[3]
    X[4] = Xabc[2] * xiso[4] + Xabc[4] * xiso[9]
    X[5] = Xabc[2] * xiso[5]

    mass_weight = 100
    weights = np.array([10, 4, 2, 2, .1, .1, .1, .1, .1, .1, 3])

    if (pressure < 300) and (latitude < -45):
        weights[1] = weights[1] / 10

    IndTipos = np.array([0, 1, 2, 3, 4, 5])

    b = np.array([
        theta,
        salinity,
        nitrate_preformed,
        phosphate_preformed,
        X[0],
        X[1],
        X[2],
        X[3],
        X[4],
        X[5],
        silicate,
    ])

    bNor = (b - params_means) / params_stds
    bNorO = np.concatenate(([1], bNor))

    if latitude > 48:
        IndTipos = np.array([0, 1, 2, 3, 4, 0])

    if latitude > 59:
        IndTipos = np.array([0, 1, 2, 0, 0, 0])

    y = IndTipos
    y = np.concatenate((y[0], y[y != 0]), axis=None)

    NIndTipos = y.size
    TiposSelNor = tipos_values_nor[y, :]

    # Add mass conservation equation
    TiposSelNor = np.vstack((np.ones(NIndTipos), TiposSelNor.T))

    # Conservation of mass
    PesosSel = np.concatenate(([mass_weight], weights))

    W = np.diag(PesosSel)  # Diagonal matrix

    # Solve optimization
    TiposSelNorW = W @ TiposSelNor
    bNorW = W @ bNorO

    # Use NNLS as in MATLAB
    Xx = nnls(TiposSelNorW, bNorW)[0]

    X1 = np.zeros(6)
    X1[y] = Xx[:NIndTipos]

    # Final parameters
    alk_preformed = np.sum(X1 * talk_preformed_types)
    delta_cdiseq = np.sum(X1 * cdis_types)
    si_teo = min(np.sum(X1 * sio2_types), silicate)
    si_teo = max(si_teo, 0)

    return alk_preformed, delta_cdiseq, si_teo, X1


def calculate_omp_dataset(ds, sal_types, tpot_types, no_preformed_types, po_preformed_types,
                          sio2_types, talk_preformed_types, cdis_types, params_means, params_stds, tipos_values_nor):
    """
    Apply OMP calculation to an xarray.Dataset, processing only non-NaN values.

    Parameters:
    - ds (xarray.Dataset): Dataset with input variables.
    - sal_types, tpot_types, no_preformed_types, po_preformed_types, sio2_types,
      talk_preformed_types, cdis_types: Parameters for OMP calculation.
    - params_means, params_stds, TiposValuesNor: Normalization parameters.

    Returns:
    - xarray.Dataset: Dataset with the calculated OMP results.
    """

    # Apply the single-sample function with apply_ufunc
    results = xr.apply_ufunc(
        calculate_omp_single_sample,
        ds["salinity"], ds["theta"], ds["pressure"], ds["latitude"],
        ds["nitrate_preformed"], ds["phosphate_preformed"], ds["silicate"],
        input_core_dims=[[], [], [], [], [], [], []],
        output_core_dims=[[], [], [], ["mixing_ratios"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float, float, float, float],
        kwargs={
            "sal_types": sal_types,
            "tpot_types": tpot_types,
            "no_preformed_types": no_preformed_types,
            "po_preformed_types": po_preformed_types,
            "sio2_types": sio2_types,
            "talk_preformed_types": talk_preformed_types,
            "cdis_types": cdis_types,
            "params_means": params_means,
            "params_stds": params_stds,
            "tipos_values_nor": tipos_values_nor,
        },
    )

    # Unpack results
    alk_preformed, delta_cdiseq, si_teo, mixing_ratios = results

    # Get the actual dimensions from input data
    dims = ds["salinity"].dims
    mixing_dims = [*ds.theta.dims, 'mix_dim']

    # Create output dataset with NaN values where input was invalid
    ds_out = xr.Dataset(coords=ds.coords)
    ds_out = ds_out.assign(
        alkalinity_preformed_parametric=(dims, alk_preformed.data),
        deltaCDiseq=(dims, delta_cdiseq.data),
        silicate_teorical=(dims, si_teo.data),
        mixing_ratios=(mixing_dims, mixing_ratios.data),
    )

    return ds_out
