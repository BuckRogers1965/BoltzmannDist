import numpy as np
from scipy import constants
from scipy.optimize import fsolve
from collections import OrderedDict
import json
import os
from typing import Dict, Union
from pathlib import Path

class BoltzmannDistribution:
    """
    A class representing the Boltzmann distribution through a single representative point.
    This class provides methods to calculate various properties of the distribution at
    different temperatures.
    
    Attributes:
        SPEED_OF_LIGHT (float): Speed of light in m/s
        representative_point (float): The 50% point at 5000K reference temperature
        reference_temperature (float): Reference temperature in Kelvin
        cache_size (int): Maximum number of entries in the energy ratio cache
        cache (OrderedDict): LRU cache for energy ratios
        cache_path (Path): Path to the cache file
    """
    
    SPEED_OF_LIGHT = 299792458  # m/s

    def __init__(self, 
                 reference_temperature: float = 5000,
                 cache_size: int = 10000,
                 cache_filename: str = "boltzmann_cache.json"):
        """
        Initialize the BoltzmannDistribution instance.

        Args:
            reference_temperature (float): Reference temperature in Kelvin
            cache_size (int): Maximum number of entries in the cache
            cache_filename (str): Name of the file to store the cache
        """
        self.representative_point = 4.7849648084645400e-20  # 50% point at 5000K
        self.reference_temperature = reference_temperature
        self.cache_size = cache_size
        self.cache: OrderedDict = OrderedDict()
        self.cache_path = Path(cache_filename)
        self._load_cache()

    def scale_to_temperature(self, temperature: float) -> float:
        """
        Scale a value from reference temperature to target temperature.

        Args:
            temperature (float): Target temperature in Kelvin

        Returns:
            float: Scaling factor
        """
        return temperature / self.reference_temperature

    def point_energy(self, temperature: float) -> float:
        """
        Calculate the point energy at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: Point energy
        """
        return self.representative_point * self.scale_to_temperature(temperature)

    def peak_energy(self, temperature: float) -> float:
        """
        Calculate the peak energy of the distribution at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: Peak energy
        """
        return self.point_energy(temperature) / 0.245671510374651163
    
    def area_under_curve_blackbody_energy(self, T):
        """
        Calculate the black body energy at a temperature using standard formula converted to using these same scaling factors.
        I did it this way to put it in line with the other formulas I am decomposing in this way to scaling factors
        Same as sigma * T**4

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: Radient black body energy
        """
        return  T**4 * 10e-8/1.763551974009425578

    def area_under_curve(self, temperature: float) -> float:
        """
        Calculate the area under the distribution curve at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: Area under the curve
        """
        return (self.point_energy(temperature) * temperature * 
                10e-22 / 50.20444590190353665093425661)

    def peak_frequency(self, temperature: float) -> float:
        """
        Calculate the peak frequency of the distribution at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: Peak frequency
        """
        return self.point_energy(temperature) * 10e32 / 0.1627836661598892

    def fwhm(self, temperature: float) -> float:
        """
        Calculate the Full Width at Half Maximum (FWHM) at a given temperature.

        Args:
            temperature (float): Temperature in Kelvin

        Returns:
            float: FWHM value
        """
        return self.point_energy(temperature) / 0.162935865000977884

    @staticmethod
    def wavelength_from_frequency(frequency: float) -> float:
        """
        Convert frequency to wavelength using the speed of light.

        Args:
            frequency (float): Frequency in Hz

        Returns:
            float: Wavelength in meters
        """
        return BoltzmannDistribution.SPEED_OF_LIGHT / frequency

    @staticmethod
    def frequency_from_wavelength(wavelength: float) -> float:
        """
        Convert wavelength to frequency using the speed of light.

        Args:
            wavelength (float): Wavelength in meters

        Returns:
            float: Frequency in Hz
        """
        return BoltzmannDistribution.SPEED_OF_LIGHT / wavelength

    def energy_at_percentage(self, temperature: float, percentage: float) -> float:
        """
        Calculate the energy at a given percentage of the distribution.

        Args:
            temperature (float): Temperature in Kelvin
            percentage (float): Percentage point to calculate (between 0 and 1)

        Returns:
            float: Energy at the specified percentage
        """
        k = constants.Boltzmann
        
        def equation(E):
            return np.exp(-E / (k * temperature)) - percentage
        
        initial_guess = k * temperature
        return float(fsolve(equation, initial_guess)[0])

    def get_energy_ratio(self, percentage: float, temperature: float) -> float:
        """
        Get the energy ratio for a given percentage and temperature, using caching.

        Args:
            percentage (float): Percentage point to calculate
            temperature (float): Temperature in Kelvin

        Returns:
            float: Energy ratio
        """
        if percentage in self.cache:
            reference_ratio = self.cache[percentage]
            self.cache.move_to_end(percentage)
        else:
            reference_ratio = self.energy_at_percentage(
                self.reference_temperature, percentage)
            self.cache[percentage] = reference_ratio
            if len(self.cache) > self.cache_size:
                self.cache.popitem(last=False)

        return reference_ratio * self.scale_to_temperature(temperature)

    def _load_cache(self) -> None:
        """Load the cache from disk if it exists."""
        if self.cache_path.exists():
            with self.cache_path.open('r') as f:
                loaded_cache = json.load(f)
                self.cache = OrderedDict(
                    (float(k), v) for k, v in loaded_cache.items())

    def save_cache(self) -> None:
        """Save the current cache to disk."""
        with self.cache_path.open('w') as f:
            json.dump({str(k): v for k, v in self.cache.items()}, f)

    def __del__(self):
        """Save the cache when the object is destroyed."""
        self.save_cache()
