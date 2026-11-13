# `pyfsb`: a `pymaster` Extension for Filtered Squared Bispectra

PyFSB is a Python module to compute Filtered-Squared Bispectra (FSBs) from the existing pseudo-Cl framework implemented in [NaMaster](https://github.com/LSSTDESC/NaMaster.git) ([Alonso et al. 2019](https://arxiv.org/abs/1809.09603)).


### Installation

Install directly from the repository with pip:

```
pip install git+https://github.com/ashaeres/pyfsb.git
```

or, to work on the code, clone it and install in editable mode:

```
git clone https://github.com/ashaeres/pyfsb.git
cd pyfsb
pip install -e ".[test]"
pytest
```

This installs the dependencies `numpy`, `healpy` and `pymaster` (NaMaster). If the `pymaster` installation gives you trouble, see the [NaMaster installation instructions](https://namaster.readthedocs.io/en/latest/source/installation.html); installing it with conda first (`conda install -c conda-forge namaster`) is often easiest.

### Usage

```python
from pyfsb import FSB, get_filters

filters = get_filters(nbands, nside)
fsb = FSB(delta, mask1=mask, filters=filters, ells_per_bin=17)
fsb.fsb_binned       # the binned FSBs, one per filter
fsb.get_full_cov()   # the analytical covariance
```

### Examples

You can familiarize yourself with the code by checking out the [`howtouse`](examples/howtouse.ipynb) jupyter notebook. The extra packages it needs can be installed with `pip install "pyfsb[examples]"`.

### Licensing

This code open source and freely available under terms consistent with GPLv3 licensing.

If you use PyFSB for any scientific publication, we ask you to cite this github repository and the companion paper https://arxiv.org/inpreparation. 

Please contact the author of this repository via github issues or email (lea.harscouet@physics.ox.ac.uk) -- any feedback is welcome!
