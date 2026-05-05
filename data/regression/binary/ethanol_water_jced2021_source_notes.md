# Ethanol + Water VLE Source Notes

This folder contains a compact transcription of the 100 kPa slice from Table 6 of:

Susial, P.; Garcia Vera, D.; Montesdeoca, I.; Santiago, D. E.; Lopez Beltran, J. "High-Pressure Phase Equilibria in an Ethanol/Water Binary System: Experimental Data and Modeling." Journal of Chemical & Engineering Data 2021, 66, 928-946. DOI: 10.1021/acs.jced.0c00686.

The source reports isobaric VLE data for ethanol (1) + water (2) at 100, 1500, and 2000 kPa. The CSV here stores the 100 kPa records as package regression inputs with SI pressure in Pa and full `x_`/`y_` mole-fraction columns for `fit_binary_pair(...)`.

The same paper reports a PC-SAFT literature/reference value `k_12 = k_21 = -0.0269` for the 100 kPa ethanol/water case using the Gross-Sadowski associating-system reference parameterization. The focused regression test uses that value as a broad sanity band, not as an exact identity target, because this repository's `2012_Held` water and ethanol pure parameters are not byte-identical to every Table 12 parameter set in the paper.
