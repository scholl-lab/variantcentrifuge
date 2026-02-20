# File: variantcentrifuge/association/weights.py
# Location: variantcentrifuge/variantcentrifuge/association/weights.py
"""
Variant weight schemes for rare variant association burden tests.

Provides three public functions:

- ``beta_maf_weights``: Beta(MAF; a, b) weights (default a=1, b=25) matching
  the SKAT R package convention — rare variants receive higher weights.
- ``uniform_weights``: All-ones weight vector for unweighted burden tests.
- ``get_weights``: String-spec parser that dispatches to either function.

Weight spec string format
-------------------------
``"beta:a,b"``
    Beta distribution weights. ``a`` and ``b`` are float parameters.
    Example: ``"beta:1,25"`` (default, SKAT convention).

``"uniform"``
    Uniform weights; all variants receive weight 1.0.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.stats import beta as _beta_dist

logger = logging.getLogger("variantcentrifuge")


def beta_maf_weights(
    mafs: np.ndarray,
    a: float = 1.0,
    b: float = 25.0,
) -> np.ndarray:
    """
    Compute Beta(MAF; a, b) density weights for each variant.

    Follows the SKAT R package convention (default a=1, b=25): rare variants
    (low MAF) receive substantially higher weights than common variants.

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies in [0, 1].
    a : float
        Beta distribution shape parameter alpha. Default: 1.0.
    b : float
        Beta distribution shape parameter beta. Default: 25.0.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Weight for each variant. Values are strictly positive.

    Notes
    -----
    MAF values are clipped to ``[1e-8, 1 - 1e-8]`` before evaluation to
    avoid numerical edge issues at 0 and 1.
    """
    maf_clipped = np.clip(np.asarray(mafs, dtype=np.float64), 1e-8, 1 - 1e-8)
    return _beta_dist.pdf(maf_clipped, a=a, b=b)


def uniform_weights(n_variants: int) -> np.ndarray:
    """
    Return a uniform (all-ones) weight vector.

    Parameters
    ----------
    n_variants : int
        Number of variants.

    Returns
    -------
    np.ndarray, shape (n_variants,), float64
        Vector of ones.
    """
    return np.ones(n_variants, dtype=np.float64)


def get_weights(mafs: np.ndarray, weight_spec: str) -> np.ndarray:
    """
    Parse a weight specification string and return the corresponding weights.

    Parameters
    ----------
    mafs : np.ndarray, shape (n_variants,)
        Minor allele frequencies, required for Beta weight computation.
        Passed to ``beta_maf_weights``; ignored for ``"uniform"``.
    weight_spec : str
        Specification string. Supported formats:
        - ``"beta:a,b"`` — calls ``beta_maf_weights(mafs, a, b)``
        - ``"uniform"`` — calls ``uniform_weights(len(mafs))``

    Returns
    -------
    np.ndarray, shape (n_variants,), float64

    Raises
    ------
    ValueError
        If ``weight_spec`` does not match a known format.

    Examples
    --------
    >>> import numpy as np
    >>> mafs = np.array([0.01, 0.05, 0.10])
    >>> get_weights(mafs, "beta:1,25")
    array([...])
    >>> get_weights(mafs, "uniform")
    array([1., 1., 1.])
    """
    mafs_arr = np.asarray(mafs, dtype=np.float64)

    if weight_spec == "uniform":
        return uniform_weights(len(mafs_arr))

    if weight_spec.startswith("beta:"):
        params_str = weight_spec[len("beta:") :]
        try:
            parts = params_str.split(",")
            if len(parts) != 2:
                raise ValueError()
            a = float(parts[0].strip())
            b = float(parts[1].strip())
        except (ValueError, IndexError) as err:
            raise ValueError(
                f"Invalid beta weight spec '{weight_spec}'. "
                "Expected format: 'beta:a,b' where a and b are floats "
                "(e.g. 'beta:1,25')."
            ) from err
        return beta_maf_weights(mafs_arr, a=a, b=b)

    raise ValueError(
        f"Unknown weight spec '{weight_spec}'. Supported specs: 'beta:a,b' or 'uniform'."
    )
