import numpy as np 


class Waters:
    def __init__(self, *args, **kwargs):
        self._initialize_water_masses()

    def _initialize_water_masses(self):
        """Initialize water masses types and constants"""
        # Predefined Water Masses Types
        self.water_masses_names = ['N5', 'NF', 'S5', 'SF', 'CDW', 'CDW2']

        # Initialize water masses parameters
        # maybe xr.DataArray(
        # np.array([5, -1.1, 5, -1.72, 1.5, -0.7]), coords={'name': self.water_masses_names},
        # name='sal_types'
        # )
        self.tpot_types = np.array([5, -1.1, 5, -1.72, 1.5, -0.7])
        self.sal_types = np.array([35.20, 34.88, 34.1, 34.35, 34.7, 34.65])
        self.no_preformed_types = np.array([10.5, 9.8, 20.6, 37, 15.1, 20.5])
        self.po_preformed_types = np.array([0.63, 0.67, 1.43, 2.6, 1.07, 1.48])
        self.x1_iso_types = np.array([1, 0, 0, 0, 0, 0])
        self.x2_iso_types = np.array([0, 1, 0, 0, 0, 0])
        self.x3_iso_types = np.array([0, 0, 1, 0, 0, 0])
        self.x4_iso_types = np.array([0, 0, 0, 1, 0, 0])
        self.x5_iso_types = np.array([0, 0, 0, 0, 1, 0])
        self.x6_iso_types = np.array([0, 0, 0, 0, 0, 1])
        self.sio2_types = np.array([11, 6, 16, 72, 106, 162])

        # From excel file
        self.cdis_types = np.array([-6.9, -21.5, 0.0, 4.5, -20.0, -29.8])
        self.talk_preformed_types = np.array([2310.077578, 2285.655901, 2281.588882,
                                              2343.657338, 2319.485664, 2347.211967])
        self.sio2_preformed_types = np.array([6, 6, 7, 70, 39.25, 80])

        # Calculate means and stds for normalization
        independent_vars = np.vstack([
            self.tpot_types, self.sal_types, self.no_preformed_types,
            self.po_preformed_types, self.x1_iso_types, self.x2_iso_types,
            self.x3_iso_types, self.x4_iso_types, self.x5_iso_types,
            self.x6_iso_types, self.sio2_types,

        ]).T

        self.params_means = np.mean(independent_vars, axis=0)
        self.params_stds = np.std(independent_vars, axis=0, ddof=1)

        self.tipos_values_nor = (independent_vars - self.params_means) / self.params_stds
