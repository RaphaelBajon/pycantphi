import xarray as xr
from scipy.interpolate import pchip_interpolate
import PyCO2SYS as pyco2

from ..preprocessing.initialisation import Initialization

from .cant import CantCalculator
from .functions import Functions

from ..data.mauna_loa import load_maunaloa_pco2
from ..water_masses.omp import calculate_omp_dataset
from ..utils.cdiseq import calculate_delta_cdiseq_dataset
from ..water_masses.definitions import Waters
from ..water_masses.mediterranean import calculate_med_water


class CantPhiCt0(Initialization):
    '''
    CantPhiCt0 is a class to calculate the PhiCt method. 
    '''
    def __init__(self, ds):
        # Initialisation of the different parameters and also the water masses characteristics
        super(CantPhiCt0, self).__init__(ds=ds)
        self.ds_exit = None

    def process(self):
        ds = self.ds.copy()
        northern_sills = (ds["latitude"] > 62) & (Functions.sample_is_at_north_of_sills(ds))
        #
        ds_1 = self._calculate_minus5(ds, northern_sills)
        ds_2 = self._calculate_plus5(ds, northern_sills)
        #
        ds = xr.merge([ds_1, ds_2])
        #
        ds['carbon_preformed'] = self.calculate_ct_preformed(ds.theta, ds.salinity,
                                                             ds.alkalinity_preformed_parametric_revised,
                                                             ds.silicate_teorical, ds.phosphate_preformed)

        ds['deltaCa_corrected'] = xr.where(ds.deltaCa > 0, ds.deltaCa, 0)

        ds['cAntSatML'] = self._calculate_cant_saturation_ml(
            ds.theta, ds.salinity, ds.alkalinity_preformed_parametric_revised,
            ds.silicate, ds.phosphate_preformed, ds.year
        )

        # Put all cant values
        self._calculate_all_cantvalues(ds)

        for year in self.output_years:
            ds[f'cAntSatML{year}'] = self._calculate_cant_saturation_ml(
                theta=ds['theta'],
                salinity=ds['salinity'],
                at_preformed=ds['alkalinity_preformed_parametric_revised'],
                silicate=ds.silicate,
                po4_preformed=ds.phosphate_preformed,
                year=year
            )

            ds[f'cAntPhiCt0ML{year}'] = ds['cAntPhiCt0ML'] * ds[f'cAntSatML{year}'] / ds['cAntSatML']

            ds[f'cAntTrocaML{year}'] = ds['cAntTroca'] * ds[f'cAntSatML{year}'] / ds['cAntSatML']

        self.ds_exit = ds

        return ds

    def _apply_alkalinity_corrections(self, ds: xr.Dataset, condition_theta: xr.DataArray) -> xr.Dataset:
        ds = ds.where(condition_theta)
        """Apply corrections to alkalinity parameters."""
        ds['deltaCa'] = Functions.calculate_delta_ca(ds.alkalinity, ds.alkalinity_preformed_parametric, ds.aou,
                                                     self.CONSTANTS)
        # Correction over preformed alkalinity
        ds['alkalinity_preformed_parametric_revised'] = xr.where(
            (ds.deltaCa < 0) & condition_theta,
            ds.alkalinity_preformed_parametric + 2 * ds.deltaCa,
            ds.alkalinity_preformed_parametric
        )
        # Correction over deltaCa
        ds['deltaCa'] = Functions.calculate_delta_ca(ds.alkalinity, ds.alkalinity_preformed_parametric_revised, ds.aou,
                                                     self.CONSTANTS)
        return ds

    def _calculate_minus5(self, ds, northern_sills):
        # Take only theta < 5
        condition_theta = ds["theta"] < 5
        ds = ds.where(condition_theta)

        # Case 1: Sample is northern than Sills
        case1 = northern_sills & condition_theta
        ds['alkalinity_preformed_parametric'] = ds.alkalinity_preformed_parametric.where(~case1, ds["alkalinity"])
        ds['deltaCDiseq'] = ds.deltaCDiseq.where(~case1, 2.5 * (ds["theta"] - 5) - 6)
        ds['silicate_teorical'] = ds.silicate_teorical.where(~case1, ds["silicate"])

        # Case 2: Perform OMP (Placeholder for additional logic)
        case2 = condition_theta & ~northern_sills
        # Here calculate OMP !!!
        new_ds = ds.where(case2)
        new_ds = calculate_omp_dataset(
            new_ds,
            self.sal_types,
            self.tpot_types,
            self.no_preformed_types,
            self.po_preformed_types,
            self.sio2_types,
            self.talk_preformed_types,
            self.cdis_types,
            self.params_means,
            self.params_stds,
            self.tipos_values_nor
        )
        # need to merge to not overwrite the first values for alkalinity_preformed_parametric,
        # deltaCDiseq and silicate_teorical
        ds = xr.merge([new_ds, ds])
        # now I can revise the different variables
        ds = self._apply_alkalinity_corrections(ds, condition_theta)
        return ds

    def _calculate_plus5(self, ds, northern_sills):
        # Take only theta >= 5
        condition_theta = (ds["theta"] >= 5)
        ds = ds.where(condition_theta)

        # put 0 for all values where theta >= 5
        ds['silicate_teorical'] = ds.silicate_teorical.where(~condition_theta, 0)
        # Case 1: Sample is northern than Sills
        case1 = northern_sills & condition_theta
        ds['silicate_teorical'] = ds.silicate_teorical.where(~case1, ds["silicate"])
        # Case 2:
        case2 = condition_theta & ~northern_sills & (ds["theta"] < 14)
        ds['silicate_teorical'] = ds.silicate_teorical.where(~case2, -0.7 * ds.theta + 9.9)

        # Med Waters here
        calculate_med_water(ds)
        # Calculate preformed alkalinity
        ds['alkalinity_preformed_parametric'] = Functions.calculate_preformed_alkalinity(
            ds, ds.nitrate_preformed, ds.phosphate_preformed, ds.silicate_teorical, northern_sills
        )
        # Calculate deltaCDiseq using
        ds['deltaCDiseq'] = calculate_delta_cdiseq_dataset(ds, condition_theta, northern_sills)
        # Calculate deltaCa
        ds = self._apply_alkalinity_corrections(ds, condition_theta)
        return ds

    @staticmethod
    def _calculate_all_cantvalues(ds):
        """Calculate final parameters including cAnt values"""
        # Calculate TrOCA
        ds['cAntTroca'] = CantCalculator.calculate_ctroca(ds['carbon'], ds['alkalinity'], ds['oxygen'], ds['theta'])
        # Calculate cStar
        ds['cStar'] = CantCalculator.calculate_cstar(ds['carbon'], ds['aou'], ds['carbon_preformed'], ds['deltaCa_corrected'])
        # change existing values (Todo: Changes in Cdis fiz 20161027)
        # Create a mask for the conditions
        mask = ((ds.latitude > 40) &
                (ds.pressure < 400) &
                (ds.salinity < 34.8))
        # Use xr.where() to calculate deltaCDiseq conditionally
        ds['deltaCDiseq'] = xr.where(
            mask,
            45.9 * (ds.salinity - 34.8) - 12,
            ds.deltaCDiseq  # Keep existing values where mask is False
        )
        # Final step: calculate cphi (which is using deltaCDiseq)
        CantCalculator.calculate_cphi(ds)
        return ds

    def _calculate_co2_system(self, alk, fco2, salinity, temp, pressure, si, po4):
        """Calculate CO2 system parameters using PyCO2SYS"""
        # Define CO2SYS parameters
        par1_type = 1  # alkalinity
        par2_type = 5  # fCO2
        ph_scale_in = 1  # total scale
        # use CO2SYS
        results = pyco2.sys(
            par1=alk,
            par2=fco2,
            par1_type=par1_type,
            par2_type=par2_type,
            salinity=salinity,
            temperature=temp,
            temperature_out=temp,
            pressure=pressure,
            pressure_out=pressure,
            total_silicate=si,
            total_phosphate=po4,
            opt_pH_scale=ph_scale_in,
            opt_k_carbonic=self.CONSTANTS['K1K2'],
            opt_k_bisulfate=self.CONSTANTS['KSO4'],
        )
        return results

    def calculate_ct_preformed(self, theta, salinity, at_preformed_parametric_revised, silicate, po4_preformed):  

        fco2_preind = Functions().calculate_fco2_hum(
            theta=theta,
            salinity=salinity,
            co2_molar_fraction=self.CONSTANTS['PREINDUSTRIAL_CO2_MOLAR_FRACTION'],  # in ppm
            atm_pressure=self.CONSTANTS['ATMOSPHERIC_PRESSURE'],
        )

        result_preind_pref = self._calculate_co2_system(
            alk=at_preformed_parametric_revised,
            fco2=fco2_preind,
            salinity=salinity,
            temp=theta,
            pressure=0,  # pressure at surface
            si=silicate,
            po4=po4_preformed
        )

        return xr.DataArray(
            name='carbon_preformed',
            data=result_preind_pref['dic'],
            dims=theta.dims, coords=theta.coords
        )

    def _calculate_cant_saturation_ml(self, theta, salinity, at_preformed, silicate, po4_preformed, year):
        """
        Calculate delta in CT by comparing with preindustrial pCO2 state using Mauna Loa data
        Original MATLAB function: cant_saturation_ml.m
        IIM-CSIC 2016 (avelo, ffperez)

        Args:
            theta: Temperature array
            salinity: Salinity array
            at_preformed: Preformed alkalinity array
            silicate: Silicate concentration array
            po4_preformed: Preformed phosphate array
            year: Year for calculation

        Returns:
            Array: Difference between current and preindustrial CT
        """

        try:
            # Load Mauna Loa dataset
            mauna_loa = load_maunaloa_pco2()  
        except Exception as ex:
            raise RuntimeError(f"Error loading Mauna Loa pCO2 data: {str(ex)}")

        # Calculate preindustrial CT
        fco2_preind = Functions().calculate_fco2_hum(
            theta=theta,
            salinity=salinity,
            co2_molar_fraction=self.CONSTANTS['PREINDUSTRIAL_CO2_MOLAR_FRACTION'],  # in ppm
            atm_pressure=self.CONSTANTS['ATMOSPHERIC_PRESSURE'],
        )

        result_preind = self._calculate_co2_system(
            alk=at_preformed,
            fco2=fco2_preind,
            salinity=salinity,
            temp=theta,
            pressure=0,  # pressure at surface
            si=silicate,
            po4=po4_preformed
        )

        ct_preind = result_preind['dic']  # Get DIC from results

        ml_pco2 = pchip_interpolate(
            mauna_loa['year'],
            mauna_loa['pco2'],
            year - self.CONSTANTS['MAUNA_LOA_OFFSET_YEARS'],
        )

        # Calculate current year CT
        fco2_current = Functions().calculate_fco2_hum(
            theta=theta,
            salinity=salinity,
            co2_molar_fraction=ml_pco2,  # ppm
            atm_pressure=self.CONSTANTS['ATMOSPHERIC_PRESSURE']
        )

        result_current = self._calculate_co2_system(
            alk=at_preformed,
            fco2=fco2_current,
            salinity=salinity,
            temp=theta,
            pressure=0,  # pressure at surface
            si=silicate,
            po4=po4_preformed,
        )

        ct_year = result_current['dic']  # Get DIC from results

        return xr.DataArray(
            name='cAntSatML',
            data=(ct_year - ct_preind),
            dims=theta.dims,
            coords=theta.coords,
        )


# To implement it after like this
@xr.register_dataset_accessor('cantaccessor')
class CantAccessor:
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
