"""Class for neutrino flux models"""

import numpy as np
import os.path

from . import units

MODEL_DIR = os.path.join(os.path.dirname(__file__), 'default_models')


class Model:
    """
    Class for storing details of a flux model.

    Parameters
    ----------
    filename : str
        Name of the model's data file. Can be a simple name for default models,
        but must be the file path for any custom models.

    Attributes
    ----------
    name : str
        The name of the model, from the "NAME" field of the file.
    description : str
        A description of the model, from the "SOURCE" field of the file.
    energies : ndarray
        An array of energies at which the fluxes are given.
    fluxes, band_min, band_max : ndarray
        The flux data for the model. Typically either `fluxes` is non-zero or
        the pair `band_min` and `band_max` are non-zero, though sometimes there
        is data for all three.
    metadata : dict
        A dictionary containing other metadata about the model from the file.
        This includes entries for the energy units, flux units, energy power,
        the model's data type, and the energy binning.

    """
    def __init__(self, filename):
        smart_filename = filename
        if not smart_filename.endswith('.txt'):
            smart_filename += '.txt'
        smart_filename = os.path.join(MODEL_DIR, smart_filename)
        if os.path.isfile(smart_filename):
            self.load_from_file(smart_filename)
        else:
            self.load_from_file(filename)

    def load_from_file(self, filename):
        """
        Read in data to the model from the given `filename`.

        Paramters
        ---------
        filename : str
            The file path of the model's data file.

        """
        name = None
        source = None
        energy_col = None
        flux_col = None
        min_col = None
        max_col = None
        e_power = None
        bins_per_decade = None
        data_type = None
        # Interpret metadata using comments at start of model file
        with open(filename, 'r') as f:
            column_number = -1
            for line in f:
                line = line.rstrip()
                if line=='':
                    continue
                if not line.startswith('#'):
                    break
                line = line.strip('#')
                words = line.split()
                column_number += 1
                if 'name' in line.lower():
                    name = " ".join(words[1:])
                    continue
                if 'source' in line.lower():
                    source = " ".join(words[1:])
                    continue
                if 'data type:' in line.lower():
                    data_type = " ".join(words[2:]).lower()
                    continue
                if 'bins per decade' in line.lower():
                    bins_per_decade = float(words[-1])
                    continue
                if words[0].lower().startswith("column"):
                    column_number = -1
                    continue
                if words[0].lower()=='energy':
                    unit = words[1].strip('[').rstrip(']')
                    if not hasattr(units, unit):
                        raise ValueError("Unable to interpret unit "+words[2])
                    energy_unit = getattr(units, unit)
                    energy_col = column_number
                elif words[0].lower()=='flux' and words[1].lower()!='band':
                    flux_unit = 1
                    for word in words[1:]:
                        word = word.strip('[').rstrip(']')
                        bits = word.rsplit("^", 1)
                        unit = bits[0]
                        power = float(bits[1]) if len(bits)>1 else 1
                        if not hasattr(units, unit):
                            raise ValueError("Unable to interpret unit ["
                                             +word+"]")
                        if unit.endswith('eV'):
                            e_power = power
                        elif unit.endswith('m'):
                            if power!=-2:
                                raise ValueError("Expected unit ["+word+
                                                 "] to be to the -2 power")
                        else:
                            if power!=-1:
                                raise ValueError("Expected unit ["+word+
                                                 "] to be to the -1 power")
                        flux_unit *= getattr(units, unit)**power
                    flux_col = column_number
                elif 'minimum' in words:
                    min_col = column_number
                elif 'maximum' in words:
                    max_col = column_number

        # Store data from file in standardized units
        data = np.loadtxt(filename)
        self.name = name
        self.description = source
        self.energies = data[:, energy_col] * energy_unit
        self.fluxes = data[:, flux_col] * flux_unit
        self.band_min = data[:, min_col] * flux_unit
        self.band_max = data[:, max_col] * flux_unit
        # Ensure raw flux units are used (no E^2 scaling)
        self.fluxes *= self.energies**(-1-e_power)
        self.band_min *= self.energies**(-1-e_power)
        self.band_max *= self.energies**(-1-e_power)
        # Ensure 1 bin per decade is used for upper limits
        if 'limit' in data_type:
            self.fluxes /= bins_per_decade
            self.band_min /= bins_per_decade
            self.band_max /= bins_per_decade
        # Store metadata from file so operations can be traced back
        self.metadata = {
            "energy_unit": energy_unit,
            "flux_unit": flux_unit,
            "energy_power": e_power,
            "data_type": data_type,
            "bins_per_decade": bins_per_decade
        }
