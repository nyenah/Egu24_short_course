# -*- coding: utf-8 -*-
"""
Priestly-Taylor potential evapotranspiration

Example values are not realistic
"""


class get_evapotranspiration:
    """ Compute potential evapotranspiration."""

    def __init__(self):
        self.pt_coeff_arid = 1.76
        self.pt_coeff_humid = 1.26

    def calculate_pet(self, slope_svp, psy_cons, net_rad, region):
        """
        Compute Priestly-Taylor potential evapotranspiration.

        Parameters
        ----------
        slope_svp : TYPE
            DESCRIPTION.
        psy_cons : TYPE
            DESCRIPTION.
        net_rad : TYPE
            DESCRIPTION.
        region : TYPE
            DESCRIPTION.

        Returns
        -------
        pet : TYPE
            DESCRIPTION.

        """

        if region == "humid":
            pet = self.pt_coeff_humid * (slope_svp * net_rad / (slope_svp + psy_cons))
        else:
            pet = self.pt_coeff_arid * (slope_svp * net_rad / (slope_svp + psy_cons))
        return pet

evap_init = get_evapotranspiration()
PET = evap_init.calculate_pet(2, 3, 4, "arid")
print(f"Result: {PET}")