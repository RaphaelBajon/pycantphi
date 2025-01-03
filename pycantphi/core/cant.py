import xarray as xr
import numpy as np


class CantCalculator:

    @staticmethod
    def calculate_cstar(ct, aou, ct_preformed, delta_ca_corrected):
        return ct - aou / 1.45 - ct_preformed - delta_ca_corrected

    @staticmethod
    def calculate_cphi(ds):
        ah8 = 0.55
        # Initialize cAntPhiCt0ML with NaN
        ds['cAntPhiCt0ML'] = xr.full_like(ds.theta, np.nan)

        # Calculate where conditions are met (AOU > -9 & pressure > 90)
        phi_mask = (ds.aou > -9) & (ds.pressure > 90)
        ds['cAntPhiCt0ML'] = xr.where(
            phi_mask,
            (ds.cStar - ds.deltaCDiseq) /
            (1 + ah8 / (ds.salinity / 35 * (46 + 0.85 * ds.theta)) *
             np.abs(ds.deltaCDiseq)),
            ds['cAntPhiCt0ML']
        )

        # Fill surface values (pressure < 25)
        surface_mask = (ds.pressure < 25) & ds.cAntPhiCt0ML.isnull()
        ds['cAntPhiCt0ML'] = xr.where(
            surface_mask,
            ds.cAntSatML,
            ds['cAntPhiCt0ML']
        )

        # Apply saturation limit
        ds['cAntSatLimitML'] = 1.05 * ds.cAntSatML
        ds['cAntPhiCt0ML'] = xr.where(
            ds.cAntPhiCt0ML > ds.cAntSatLimitML,
            ds.cAntSatLimitML,
            ds.cAntPhiCt0ML
        )
        return ds

    @staticmethod
    def calculate_ctroca(ct, at, oxygen, theta):
        """Calculate anthropogenic CO2 using TrOCA method"""
        # TrOCA equation parameters
        a = 1.279
        b = 7.511
        c = 0.01087
        d = 781000
        
        # Calculate TrOCA
        troca = (oxygen + a * (ct - at / 2) - np.exp(b - c * theta - d / at**2)) / a
        
        return troca