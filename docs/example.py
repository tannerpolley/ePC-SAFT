import numpy as np
from epcsaft import create_parameter_template, ePCSAFTMixture

# Toluene
mixture = ePCSAFTMixture.from_params(
    {
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    },
    species=["Toluene"],
)

t = 320  # K
p = 101325  # Pa
state = mixture.state(T=t, x=np.asarray([1.0]), P=p)
den = state.density()
den_mass = state.density(units="mass")
den_molar = state.molar_density()
print("Density of toluene at {} K:".format(t), den, "mol m^-3")
print("Mass density of toluene at {} K:".format(t), den_mass, "kg m^-3")
print("Explicit molar density of toluene at {} K:".format(t), den_molar, "mol m^-3")

# Water using a user-owned external parameter folder
template_root = create_parameter_template(
    location=r"C:\Users\Tanner\Documents\my_epcsaft_data",
    folder_name="water_salt_case",
    species=["H2O", "Na+", "Cl-"],
)
mixture = ePCSAFTMixture.from_dataset(
    template_root,
    ["H2O", "Na+", "Cl-"],
    np.asarray([0.9998, 1e-4, 1e-4]),
    274.0,
)
t = 274
p = 101325
state = mixture.state(T=t, x=np.asarray([0.9998, 1e-4, 1e-4]), P=p)
den = state.density()
den_mass = state.density(units="mass")
den_molar = state.molar_density()
print("Density of water at {} K:".format(t), den, "mol m^-3")
print("Mass density of water at {} K:".format(t), den_mass, "kg m^-3")
print("Explicit molar density of water at {} K:".format(t), den_molar, "mol m^-3")
