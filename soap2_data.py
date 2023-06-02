"""
A module providing local, cached access to archived SOAP2 HDF5 data files
stored on a remote server.

This implementation accesses files stored on Dropbox, used as an
interim archive.  For long-term storage, another server will be used,
such as Dataverse or Zenodo.

Created 2023-05-22 by Tom Loredo
"""
from pathlib import Path

import pooch
import h5py


def prep_data(dcdir=None):
    """
    Prepare for data access, returning a Pooch data fetching object.

    `dcdir` may specify a path to a directory to use to cache the data;
    if unspecified, an OS-specific location appropriate for data caches
    will be used.  This location is stored in the `path` attribute of
    the returned data fetcher.
    """

    if dcdir is None:
        data_cache_dir = pooch.os_cache("SOAP2-1Spot")
    else:
        data_cache_dir = Path(dcdir)

    fetcher = pooch.create(
        path=data_cache_dir,
        base_url="",
        # We are not currently versioning this data.
        # version=version,
        # version_dev="main",
        registry={
        "lambda-3923-4010-phases-100.h5" : None,
        "lambda-3923-6664-phases-4.h5" : None
        },
        # Now specify custom URLs for some of the files in the registry.
        urls={
            "lambda-3923-4010-phases-100.h5" : "https://www.dropbox.com/s/5b9m1pq5qif5obf/lambda-3923-4010-phases-100.h5?dl=1",
            "lambda-3923-6664-phases-4.h5" : "https://www.dropbox.com/s/pyeapovhk4q6az0/lambda-3923-6664-phases-4.h5?dl=1"
        },
    )

    return fetcher


class DynamicSpectrum:
    """
    Load and access a SOAP dynamic spectrum stored in an HDF5 file.
    """
 
    def __init__(self, h5_path):
        """
        Load a SOAP dynamic spectrum stored in an HDF5 file
        """
        self.h5_path = h5_path
        store = h5py.File(h5_path, 'r')
        print('Loaded HDF5 store {}:'.format(h5_path))

        # Load all items in the HDF5 store into the namespace.
        for name, value in store.attrs.items():
            print('  {}: {}'.format(name, value))
            setattr(self, name, value)
        # Store the four key arrays as NumPy arrays.
        names = ['lambdas', 'quiet', 'phases', 'active']
        for name in names:
            print(store[name])
            setattr(self, name, store[name][()])
        self.nphases = len(self.phases)
        store.close()


def fetch_full_spec(fetcher):
    """
    Load the full spectrum, computed at just 4 phases.
    """
    h5_path = fetcher.fetch("lambda-3923-6664-phases-4.h5", progressbar=True)
    dynspec = DynamicSpectrum(h5_path)
    return dynspec

def fetch_ca_spec(fetcher):
    """
    Load the spectrum for the Ca H & K line region, computed at 100 phases.
    """
    h5_path = fetcher.fetch("lambda-3923-4010-phases-100.h5", progressbar=True)
    dynspec = DynamicSpectrum(h5_path)
    return dynspec


if __name__ == '__main__':
    
    fetcher = prep_data('SOAP2-1Spot')

    # Full spectrum, just 4 phases:
    full_spec = fetch_full_spec(fetcher)

    # Calcium H & K line region, 100 phases:
    ca_spec = fetch_ca_spec(fetcher)
