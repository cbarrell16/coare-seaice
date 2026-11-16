# Draft bug report for pycoare

*For Chris to post as a GitHub issue on https://github.com/pyCOARE/coare/issues. Not
sent. Checked against pycoare 0.4.3 (PyPI) and the development branch as of 16 Nov 2026.*

**Title:** `qair` modifies the caller's `rh` array in place, so reusing it gives wrong fluxes

**Body:**

When `rh` is passed as a float numpy array, `pycoare.util.qair` divides the caller's
array by 100 in place. A second call with the same array then uses humidities 100 times
too small, with no warning.

```python
import numpy as np
from pycoare import coare_35

kw = dict(t=np.array([-10.0, 5.0]), ts=np.array([-1.0, 8.0]), zu=10, zt=10, zq=10,
          p=np.array([1000.0, 1010.0]), lat=np.array([75.0, 60.0]))
rh = np.array([80.0, 90.0])

a = coare_35(np.array([8.0, 6.0]), rh=rh, **kw)
print(a.fluxes.hlb)   # [71.90 40.23]  correct
print(rh)             # [0.8 0.9]      caller's array has changed
b = coare_35(np.array([8.0, 6.0]), rh=rh, **kw)
print(b.fluxes.hlb)   # [122.32 161.28]  wrong
```

The same happens when the array is a view into an xarray Dataset or pandas DataFrame
(for example `ds["rh"].values`), so the user's dataset is changed too.

**Cause:** in `pycoare/util.py`, `qair` does

```python
rh = np.asarray(rh, dtype=float)
rh /= 100.0
```

`np.asarray` returns the same object when the input is already a float64 array, so the
in-place division changes it.

**Possible fix:** `rh = np.asarray(rh, dtype=float) / 100.0`, which makes a new array.
A test could check that the input array is unchanged after a call.

Environment: pycoare 0.4.3, numpy 2.x, Python 3.10–3.13.
