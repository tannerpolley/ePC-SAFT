import numpy as np
from pcsaft import PCSAFTMixture

# Toluene
mixture = PCSAFTMixture.from_params(
    ["Toluene"],
    {
        "m": np.asarray([2.8149]),
        "s": np.asarray([3.7169]),
        "e": np.asarray([285.69]),
    },
)

t = 320 # K
p = 101325 # Pa
state = mixture.state(T=t, x=np.asarray([1.0]), P=p)
den = state.density()
print('Density of toluene at {} K:'.format(t), den, 'mol m^-3')

# Water using the packaged dataset
mixture = PCSAFTMixture.from_dataset("2012_Held", ["Water"])
t = 274
p = 101325
state = mixture.state(T=t, x=np.asarray([1.0]), P=p)
den = state.density()
print('Density of water at {} K:'.format(t), den, 'mol m^-3')
