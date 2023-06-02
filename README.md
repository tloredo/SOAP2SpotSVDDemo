# SOAP2SpotSVDDemo
SOAP2 sunspot spectrum SVD demo

This script and notebook run in an `eprv10` Conda environment created as follows:

```bash
conda create --name eprv10 -c conda-forge -c defaults python=3.10 ipython jupyter scipy matplotlib ipympl \
  pooch tqdm h5py beautifulsoup4 html5lib bleach pandas sortedcontainers \
  pytz setuptools mpmath bottleneck jplephem asdf pyarrow colorcet hypothesis astropy copier gsl
```

Note that this adds some packages to our previous `eprv10` env. Those may be added manually from within the activated env via:

```
conda install -c conda-forge -c defaults ipympl pooch tqdm
```

