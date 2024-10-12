import numpy as np
from scipy import constants
from scipy.optimize import fsolve
import importlib.util
import sys
from pathlib import Path
import numpy as np
from scipy import constants
from scipy.integrate import quad
from scipy.optimize import fsolve
from BoltzmannDistribution import BoltzmannDistribution

# Traditional methods for comparison
def planck_law(freq, T):
    h = constants.h
    c = constants.c
    k = constants.k
    return 2 * h * freq**3 / (c**2 * (np.exp(h * freq / (k * T)) - 1))

def calculate_fwhm_blackbody(T):
    h = constants.h
    c = constants.c
    k = constants.k
    def find_peak_freq(T):
        return 2.821439 * k * T / h  # Wien's displacement law

    peak_freq = find_peak_freq(T)
    peak_intensity = planck_law(peak_freq, T)
    half_max = peak_intensity / 2

    def half_max_equation(freq):
        return planck_law(freq, T) - half_max

    # Find the two frequencies at half maximum
    f1 = fsolve(half_max_equation, peak_freq / 2)[0]
    f2 = fsolve(half_max_equation, peak_freq * 2)[0]

    # Calculate FWHM in frequency
    fwhm_freq = f2 - f1

    # Convert FWHM to energy (E = hf)
    fwhm_energy = constants.h * fwhm_freq

    return fwhm_energy

def calculate_peak_energy_blackbody(T):
    k = constants.k
    return 2.821439 * k * T


def planck_peak_frequency(T):
    return 2.821439 * constants.k * T / constants.h


def boltzmann(E, T):
    return np.exp(-E / (constants.k * T))

def auc_traditional(T):
    # Integrate the Boltzmann distribution from 0 to a sufficiently large value
    # We use 100kT as the upper limit, which should cover most of the significant area
    k = constants.Boltzmann
    return quad(lambda E: boltzmann(E, T), 0, 100*k*T)[0] * k * T


def calculate_area_under_curve_blackbody(T):
    return constants.sigma * T**4

def calculate_peak_frequency_blackbody(T):
    h = constants.h
    k = constants.k
    return 2.821439 * k * T / h

# Test harness
def run_tests():
    temperatures = [4, 500, 1000, 3000, 5000, 7000, 10000, 15000, 20000]
    bd = BoltzmannDistribution()

    print("1. FWHM Test")
    print("Temperature(K) | FWHM (Black Body) | FWHM (Your Method) | Ratio %")
    print("-" * 75)
    for T in temperatures:
        fwhm_blackbody = calculate_fwhm_blackbody(T)
        fwhm_your_method = bd.fwhm(T)
        ratio = fwhm_your_method / fwhm_blackbody * 100
        print(f"{T:14d} | {fwhm_blackbody:.6e} | {fwhm_your_method:.6e} | {ratio:.16f}")

    print("\n2. Peak Energy Test")
    print("Temperature(K) | Peak Energy (Black Body) | Peak Energy (Your Method) | Ratio %")
    print("-" * 85)
    for T in temperatures:
        peak_energy_blackbody = calculate_peak_energy_blackbody(T)
        peak_energy_your_method = bd.peak_energy(T)
        ratio = peak_energy_your_method / peak_energy_blackbody * 100
        print(f"{T:14d} | {peak_energy_blackbody:.6e} | {peak_energy_your_method:.6e} | {ratio:.16f}")

    print("\n3. Area Under Boltzman Curve Test")
    print("Temperature(K) | Area (Black Body) | Area (Your Method) | Ratio %")
    print("-" * 75)
    for T in temperatures:
        area_blackbody = auc_traditional(T)
        area_your_method = bd.area_under_curve(T)
        ratio = area_your_method / area_blackbody * 100
        print(f"{T:14d} | {area_blackbody:.6e} | {area_your_method:.6e} | {ratio:.16f}")

    print("\n3. Area Under Curve Blackbody energy Test")
    print("Temperature(K) | Area (Black Body) | Area (Your Method) | Ratio %")
    print("-" * 75)
    for T in temperatures:
        area_blackbody = calculate_area_under_curve_blackbody(T)
        area_your_method = bd.area_under_curve_blackbody_energy(T)
        ratio = area_your_method / area_blackbody * 100
        print(f"{T:14d} | {area_blackbody:.6e} | {area_your_method:.6e} | {ratio:.16f}")

    print("\n4. Peak Frequency Test")
    print("Temperature(K) | Peak Freq (Black Body) | Peak Freq (Your Method) | Ratio %")
    print("-" * 85)
    for T in temperatures:
        peak_freq_blackbody = calculate_peak_frequency_blackbody(T)
        peak_freq_your_method = bd.peak_frequency(T)
        ratio = peak_freq_your_method / peak_freq_blackbody * 100
        print(f"{T:14d} | {peak_freq_blackbody:.6e} | {peak_freq_your_method:.6e} | {ratio:.16f}")
    
    percentages = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]

    print("Energy Ratio Test\n")
    print("Temperature(K) | Percentage | Energy Ratio (Old) | Energy Ratio (New) | Ratio %")
    print("-" * 80)
    for T in temperatures:
        for P in percentages:
            # Energy ratio with old method (direct calculation)
            old_ratio = bd.energy_at_percentage(T, P)
            # Energy ratio with new method (caching and scaling)
            new_ratio = bd.get_energy_ratio(P, T)
            # Calculate ratio of results
            ratio = (new_ratio / old_ratio) * 100
            print(f"{T:14d} | {P:.2f}        | {old_ratio:.6e} | {new_ratio:.6e} | {ratio:.6f}")
        print("-" * 80)  # Separator between temperature sets

if __name__ == "__main__":
    run_tests()
