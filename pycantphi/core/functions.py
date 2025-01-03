import numpy as np
from dataclasses import dataclass
import pandas as pd
import xarray as xr
from ..data.sills import read_sills_data


@dataclass
class Functions:

    @staticmethod
    def filter_non_nan(ds, variables):
        """
        Filter an xarray.Dataset to include only samples where all specified variables are non-NaN.

        Parameters:
        - ds (xarray.Dataset): The dataset to filter.
        - variables (list of str): List of variable names to check for non-NaN values.

        Returns:
        - xarray.Dataset: A new dataset containing only samples with non-NaN values for all specified variables.
        """
        # Start with a mask of all True
        valid_mask = xr.full_like(ds[variables[0]], True, dtype=bool)

        # Combine masks for all variables
        for var in variables:
            valid_mask &= ds[var].notnull()

        # Apply the mask to all variables in the dataset
        return ds.where(valid_mask, drop=True)

    @staticmethod
    def wrap_to_180lon(longitude):
        """
        Wraps longitude values to the range [-180, 180].

        Parameters:
        longitude (float or array-like): Longitude value(s) to be wrapped.

        Returns:
        float or np.ndarray: Longitude value(s) wrapped to the range [-180, 180].
        """
        return ((longitude + 180) % 360) - 180

    @staticmethod
    def celsiustokelvin(temp_celsius):  # Okay
        return 273.15 + temp_celsius
    
    @staticmethod
    def _calculate_vapor_pressure(temp_kelvin, salinity):  # Okay
        return 24.4543 - 67.4509 * (100 / temp_kelvin) - 4.8489 * np.log(temp_kelvin / 100) - 0.000544 * salinity
    
    @staticmethod
    def _calculate_virial_coef(temp_kelvin):  # Okay
        return -1636.75 + 12.0408 * temp_kelvin - 0.0327957 * temp_kelvin**2 + 3.16528e-5 * temp_kelvin**3
    
    @staticmethod
    def _calculate_cross_virial_coef(temp_kelvin): 
        return 57.7 - 0.118 * temp_kelvin

    def calculate_fco2_hum(self, theta, salinity, co2_molar_fraction, atm_pressure):
        """Calculate fCO2 for human-influenced conditions"""
        temp_kelvin = self.celsiustokelvin(theta)
        watervapor_pressure = np.exp(self._calculate_vapor_pressure(temp_kelvin, salinity))
        virial = self._calculate_virial_coef(temp_kelvin)
        sig = self._calculate_cross_virial_coef(temp_kelvin)

        fco2 = co2_molar_fraction * (atm_pressure - watervapor_pressure) * np.exp((atm_pressure * (virial + 2 * sig)) /
                                                                                  (0.082 * 1000 * temp_kelvin))
        
        return fco2

    @staticmethod
    def calculate_o2_saturation(salinity, theta):  # Okay
        """Calculate oxygen saturation using Gordon and Garcia equations"""
        ts = np.log((298.15 - theta) / (273.15 + theta))

        # Constants for Gordon and Garcia
        # in mumol/kg
        a0 = 5.80871
        a1 = 3.20291
        a2 = 4.17887
        a3 = 5.10006
        a4 = - 9.86643e-2
        a5 = 3.80369
        b0 = - 7.01577e-3
        b1 = - 7.70028e-3
        b2 = - 1.13864e-2
        b3 = - 9.51519e-3
        c0 = - 2.75915e-7

        # Calculate oxygen saturation
        o2sat = np.exp(a0 + a1 * ts + a2 * ts ** 2 + a3 * ts ** 3 + a4 * ts ** 4 + a5 * ts ** 5 +
                       salinity * (b0 + b1 * ts + b2 * ts ** 2 + b3 * ts ** 3) +
                       c0 * salinity ** 2)

        return o2sat

    @staticmethod
    def calculate_aou(o2sat, oxygen):  # Okay
        return o2sat - oxygen

    @staticmethod
    def calculate_delta_ca(alkalinity, alk_preformed, aou, constants):  # Okay
        """Calculate delta Ca"""
        delta_ca = (alkalinity - alk_preformed + aou * (1 / constants['RN'] + 1 / constants['RP'])) / 2
        return delta_ca

    @staticmethod
    def sample_is_at_north_of_sills(ds):
        """
        Determine whether samples in an xarray Dataset are north of the sills boundary.

        Parameters:
        - ds (xarray.Dataset): The dataset containing `latitude` and `longitude` variables.
        - sills (xarray.Dataset): Dataset containing `longitude` and `latitude` of the sill boundary.

        Returns:
        - xarray.DataArray: A boolean array indicating True where samples are north of the sill boundary.
        """
        sills = read_sills_data()
        # Interpolate the sill latitude at each sample's longitude
        limit_lat = np.interp(ds["longitude"], sills["longitude"], sills["latitude"])

        # Compare sample latitudes with interpolated sill latitudes
        mask = ds["latitude"] > limit_lat

        return mask

    @staticmethod
    def calculate_preformed_alkalinity(
            mw_results, nitrate_preformed, phosphate_preformed, si_teo, northern_sills
    ):
        """Calculate preformed alkalinity"""
        mw_results['alkalinity_preformed_parametric'] = mw_results.alkalinity_preformed_parametric.where(
            ~northern_sills, mw_results["alkalinity"]
        )

        # Constants
        as8 = 0.1
        ar8 = 1

        mw_results['alkalinity_preformed_parametric'] = mw_results.alkalinity_preformed_parametric.where(
            northern_sills,
            (
                585.73 +
                46.2234 * mw_results['salinityCorrectedByMW'] +
                3.26643 * mw_results['thetaCorrectedByMW'] +
                0.2397 * mw_results['noConservativeCorrectedByMW'] +
                0.72779 * si_teo -
                nitrate_preformed -
                phosphate_preformed
            ) * (1 - mw_results['percentMW']) +
            mw_results['percentMW'] * 2414
            - (45.4 + 0.93 * mw_results['thetaCorrectedByMW']) * as8  # Acidification
            - 4 * ar8  # Global warming
        )

        return mw_results['alkalinity_preformed_parametric']
