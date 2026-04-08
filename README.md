# PC-SAFT

`pcsaft` is a Cython/C++ implementation of the PC-SAFT equation of state with dipole, association, and electrolyte terms.

## Install

```bash
pip install .
```

For editable development:

```bash
pip install -e . --no-build-isolation
```

## Example

```python
import numpy as np
from pcsaft import PCSAFTMixture

mixture = PCSAFTMixture.from_params(
    {"m": np.asarray([2.8149]), "s": np.asarray([3.7169]), "e": np.asarray([285.69])},
    species=["Toluene"],
)
state = mixture.state(T=320.0, x=np.asarray([1.0]), P=101325.0)
print(state.density())
```

## Package Layout

- `src/pcsaft/`: installable runtime package, Cython/C++ sources, and packaged parameter datasets
- `data/`: datasets and figures that are not required at runtime
- `docs/`: user documentation and reference material

## Public API

The main entry points are `PCSAFTMixture`, `PCSAFTState`, and the structured result objects returned by solver-style methods.

## Author

Tanner Polley

## License

GNU General Public License v3.0
