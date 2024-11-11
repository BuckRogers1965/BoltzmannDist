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


# Traditional methods for comparison

def traditional_entropy(T):
    k = constants.Boltzmann
    def integrand(x):
        # Use dimensionless energy x = E / (kT)
        return x * np.exp(-x) * (1 - np.log(x) - np.log(k*T))
    # Integrate from 0 to a large number (e.g., 100) instead of infinity
    result, _ = quad(integrand, 0, 100)
    return k * result

def traditional_partition_function(T):
    k = constants.Boltzmann
    return k * T

def traditional_average_energy(T):
    return 1.5 * constants.Boltzmann * T

def traditional_heat_capacity(T):
    return 1.5 * constants.Boltzmann

def traditional_free_energy(T):
    k = constants.Boltzmann
    return -k * T * np.log(k * T)

#def traditional_energy_distribution(T, E):
    #k = constants.Boltzmann
    #return np.exp(-E / (k * T)) / (k * T)

def traditional_most_probable_energy(T):
    return 0.5 * constants.Boltzmann * T

# Test harness
def run_tests():
    temperatures = [4, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]
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

    print("\n5. Entropy Test")
    print("Temperature(K) | Entropy (Traditional) | Entropy (Your Method) | Ratio %")
    print("-" * 80)
    prev = 0
    for T in temperatures:
        entropy_traditional = traditional_entropy(T)
        entropy_your_method = bd.entropy(T)
        ratio = entropy_your_method / entropy_traditional * 100
        print(f"{T:14d} | {entropy_traditional:.6e} | {entropy_your_method:.6e} | {ratio:.16f} {ratio-prev:.16f}")
        prev = ratio

    print("\n6. Partition Function Test")
    print("Temperature(K) | Partition Func (Traditional) | Partition Func (Your Method) | Ratio %")
    print("-" * 90)
    for T in temperatures:
        pf_traditional = traditional_partition_function(T)
        pf_your_method = bd.partition_function(T)
        ratio = pf_your_method / pf_traditional * 100
        print(f"{T:14d} | {pf_traditional:.6e} | {pf_your_method:.6e} | {ratio:.16f}")

    print("\n7. Average Energy Test")
    print("Temperature(K) | Avg Energy (Traditional) | Avg Energy (Your Method) | Ratio %")
    print("-" * 85)
    for T in temperatures:
        avg_energy_traditional = traditional_average_energy(T)
        avg_energy_your_method = bd.average_energy(T)
        ratio = avg_energy_your_method / avg_energy_traditional * 100
        print(f"{T:14d} | {avg_energy_traditional:.6e} | {avg_energy_your_method:.6e} | {ratio:.16f}")

    print("\n8. Heat Capacity Test")
    print("Temperature(K) | Heat Capacity (Traditional) | Heat Capacity (Your Method) | Ratio %")
    print("-" * 90)
    for T in temperatures:
        hc_traditional = traditional_heat_capacity(T)
        hc_your_method = bd.heat_capacity(T)
        ratio = hc_your_method / hc_traditional * 100
        print(f"{T:14d} | {hc_traditional:.6e} | {hc_your_method:.6e} | {ratio:.16f}")

    print("\n9. Free Energy Test")
    print("Temperature(K) | Free Energy (Traditional) | Free Energy (Your Method) | Ratio %")
    print("-" * 85)
    prev = 0
    for T in temperatures:
        fe_traditional = traditional_free_energy(T)
        fe_your_method = bd.free_energy(T)
        ratio = fe_your_method / fe_traditional * 100
        print(f"{T:14d} | {fe_traditional:.16e} | {fe_your_method:.16e} | {ratio:.16f} | {ratio - prev}")
        prev = ratio

    print("\n10. Most Probable Energy Test")
    print("Temperature(K) | MPE (Traditional) | MPE (Your Method) | Ratio %")
    print("-" * 75)
    for T in temperatures:
        mpe_traditional = traditional_most_probable_energy(T)
        mpe_your_method = bd.most_probable_energy(T)
        ratio = mpe_your_method / mpe_traditional * 100
        print(f"{T:14d} | {mpe_traditional:.6e} | {mpe_your_method:.6e} | {ratio:.16f}")

    
'''
    percentages = [0.01, 0.05, 0.1, 0.2, 0.3, 0.333, 0.4, 0.5, 0.53, 0.6, 0.7, 0.8, 0.9, 0.99]

    print("\nEnergy Ratio Test\n")
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
            print(f"{T:14d} | {P:.2f}        | {old_ratio:.6e} | {new_ratio:.6e} | {ratio:.16f}")
        print("-" * 80)  # Separator between temperature sets
'''

if __name__ == "__main__":
    run_tests()
