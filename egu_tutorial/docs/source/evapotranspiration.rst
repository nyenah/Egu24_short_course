Evapotranspiration
==================

.. autoclass:: evapotranspiration.get_evapotranspiration
   :members: calculate_pet
   :special-members: __init__

The potential evapotranspiration :math:`{E}_{pot}` :math:`[mm/d]` is calculated with the **Priestley–Taylor** equation according to Shuttleworth (1993) [1]_, as:


.. math::
   {E}_{pot} = \alpha\Big(\frac{S_a R}{S_a + g}\Big)

:math:`\alpha` is set to 1.26 in humid and to 1.74 in (semi)arid cells (see Appendix B in Müller et al. [2]_). :math:`R` is the net radiation :math:`[mm/d]` that depends on land cover (Table C2, Müller et al. [2]_). :math:`{S_a}` is the slope of the saturation vapor pressure–temperature relationship, and :math:`g` is the psychrometric constant :math:`[{\frac{kPa}{°C}}]`. 


References 
----------
.. [1] Shuttleworth, W.: Evaporation, in: Handbook of Hydrology, edited by: Maidment, D., McGraw-Hill, New York, 1–4, 1993
.. [2]  Müller Schmied, H., Müller, R., Sanchez-Lorenzo, A., Ahrens, B., and Wild, M.: Evaluation of radiation components in a global freshwater model with station-based observations, Water, 8, 450, https://doi.org/10.3390/w8100450, 2016b

