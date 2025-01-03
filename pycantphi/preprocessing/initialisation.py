import xarray as xr
import numpy as np
import gsw
from typing import Union

from ..core.constants import MethodConstants, ProcessingConstants, ChemicalConstants
from ..core.functions import Functions
from ..water_masses.definitions import Waters

from .validator import DataValidator

class Initialization(Waters):
    def __init__(self, ds: Union[str, xr.Dataset]):
        super().__init__()
        # Define base constants
        self.CONSTANTS = {
            'ATMOSPHERIC_PRESSURE': 1,
            'PREINDUSTRIAL_CO2_MOLAR_FRACTION': 278,
            'K1K2': 4,  # previously 10: Lueker, 2000. according best practices guide in 2007 - NOTE: Previously it was option 4, MEHRBACH refit BY DICKSON AND MILLERO.
            'KSO4': 1,  # Dickson's KSO4 (constant of choice in 2007)
            'RN': 9,
            'RP': 135,
            'MAUNA_LOA_OFFSET_YEARS': 5
        }

        # parameters needed to compute the method 
        self.needed_parameters = MethodConstants.REQUIRED_PARAMETERS

        # Output years for calculations
        self.output_years = ProcessingConstants.OUTPUT_YEARS

        # Initialize based and additional parameters
        self.base_params_to_return = ProcessingConstants.base_params_to_return
        self.additional_params = ProcessingConstants.additional_params

        # load properly ds
        self.ds = ds
        self.ds = self._process_input_data()

        # check if missing params & replace -9, -999 values
        DataValidator(ds)
        self._clean_data()

        # move longitude & calculate general parameters
        self.ds['longitude'] = Functions.wrap_to_180lon(self.ds['longitude'])
        self.ds = self._calculate_common_parameters()

        # create empty datasets for all possible output variables
        self._emptys()

    def _process_input_data(self):  # Okay
        """Process and validate input data"""
        if isinstance(self.ds, xr.Dataset):
            return self.ds
        elif isinstance(self.ds, str):
            try:
                return xr.open_dataset(self.ds)
            except Exception as e:
                raise ValueError(f"Error reading netcdf file: {e}")
        else:
            raise ValueError("Input must be xarray Dataset or path to netcdf file")

    def _clean_data(self):  # Okay
        """Clean data by removing invalid values"""
        # Replace -9 and -999 with NaN
        self.ds = self.ds.where((self.ds != -9) & (self.ds != -999), np.nan)

        # Remove rows with NaN in needed parameters
        self.ds = self.ds.dropna(
            dim='time', subset=self.needed_parameters
        ).dropna(dim='location').dropna(dim='pressure')

    def _calculate_common_parameters(self):
        """Calculate common parameters using xarray operations"""
        # Calculate O2 saturation
        ds = self.ds.copy()
        ds['o2_sat'] = Functions.calculate_o2_saturation(ds.salinity, ds.theta)

        # Calculate AOU
        ds['aou'] = Functions.calculate_aou(ds.o2_sat, ds.oxygen)

        # Calculate preformed parameters
        ds['nitrate_preformed'] = ds['nitrate'] - ds['aou'] / self.CONSTANTS['RN']
        ds['phosphate_preformed'] = ds['phosphate'] - ds['aou'] / self.CONSTANTS['RP']

        # Calculate conservative parameters
        ds['no_conservative'] = ds['nitrate'] * self.CONSTANTS['RN'] + ds['oxygen']
        ds['po_conservative'] = ds['phosphate'] * self.CONSTANTS['RP'] + ds['oxygen']

        # Calculate sigma
        ds['sigma'] = gsw.density.rho(
            ds['salinity'],
            ds['theta'],
            ds['pressure'],
        ) - 1000

        return ds

    def _calculate_sigma_wm(self):
        # Calculate density at sample depth for water masses
        self.sigma_wm = gsw.density.rho(
            xr.DataArray(self.sal_types),
            xr.DataArray(self.tpot_types),
            self.ds['pressure'],
        ) - 1000

    def _emptys(self):
        for param in self.additional_params:
            self.ds[param] = xr.full_like(self.ds.theta, np.nan)
