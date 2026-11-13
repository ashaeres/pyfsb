import healpy as hp
import numpy as np

from pyfsb import FSB, get_filters


def test_get_filters():
    nside, nbands = 16, 4
    fls = get_filters(nbands, nside)
    assert fls.shape == (nbands, 3 * nside)
    # filters are binary and do not overlap
    assert set(np.unique(fls)) == {0.0, 1.0}
    assert fls.sum(axis=0).max() == 1.0
    # monopole and dipole are excluded
    assert not fls[:, :2].any()


def test_fsb_smoke():
    nside, nbands = 16, 3
    npix = hp.nside2npix(nside)
    rng = np.random.default_rng(0)
    delta = rng.normal(size=npix)
    mask = np.ones(npix)
    fls = get_filters(nbands, nside)

    fsb = FSB(delta, mask1=mask, filters=fls, ells_per_bin=8)
    binned = fsb.fsb_binned
    assert binned.shape[0] == nbands
    assert np.all(np.isfinite(binned))
