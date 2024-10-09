import numpy as np
from scipy import constants
from scipy.optimize import fsolve
from collections import OrderedDict
import json
import os

# Speed of light in m/s
c = 299792458  

class BoltzmannDistribution:
    def __init__(self):
        self.representative_point = 4.7849648084645400e-20  # 50% point at 5000K
        self.reference_temperature = 5000
        self.cache_size = 10000
        self.cache = OrderedDict()
        self.cache_filename = "boltzmann_cache.json"
        self.load_cache()

    def scale_to_temperature(self, T):
        return T / self.reference_temperature

    def point_energy(self, T):
        return self.representative_point * self.scale_to_temperature(T)

    def peak_energy(self, T):
        return self.point_energy(T) / 0.2457097414250071665

    def area_under_curve(self, T):
        return self.point_energy(T) * T * 10e-22/50.20444590190353665093425661

    def peak_frequency(self, T):
        return self.point_energy(T) * 10e32 / .1627836661598892

    def fwhm(self, T):
        return self.point_energy(T)/.162935865000977884

    def wavelength_from_frequency(self, frequency):
        return c / frequency

    def frequency_from_wavelength(self, wavelength):
        return c / wavelength

    def energy_at_percentage(self, T, percentage):
        k = constants.Boltzmann
        def equation(E):
            return np.exp(-E / (k * T)) - percentage
        return fsolve(equation, k*T)[0]

    def calculate_energy_ratio(self, percentage):
        # Calculate energy at the reference temperature (5000K)
        return self.energy_at_percentage(self.reference_temperature, percentage)

    def get_energy_ratio(self, percentage, temperature):
        # Check if percentage is in cache
        if percentage in self.cache:
            reference_ratio = self.cache[percentage]
            self.cache.move_to_end(percentage)  # Move to end to mark as recently used
        else:
            reference_ratio = self.calculate_energy_ratio(percentage)
            self.cache[percentage] = reference_ratio
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)  # Remove least recently used item

        # Scale to requested temperature and return
        return reference_ratio * self.scale_to_temperature(temperature)

    def load_cache(self):
        if os.path.exists(self.cache_filename):
            with open(self.cache_filename, 'r') as f:
                loaded_cache = json.load(f)
                self.cache = OrderedDict((float(k), v) for k, v in loaded_cache.items())

    def save_cache(self):
        with open(self.cache_filename, 'w') as f:
            json.dump({str(k): v for k, v in self.cache.items()}, f)

    def __del__(self):
        self.save_cache()


