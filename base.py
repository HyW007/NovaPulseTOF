import itertools as it
import os
from abc import ABC, abstractmethod
from warnings import warn

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.units import UnitTypeError, get_physical_type
from astropy.units.quantity import Quantity
from scipy.special import loggamma
from snewpy._model_downloader import LocalFileLoader

from snewpy.neutrino import Flavor
from snewpy.flavor_transformation import NoTransformation
from functools import wraps

from snewpy.flux import Flux
from pathlib import Path


class SupernovaModel(ABC, LocalFileLoader):

    def get_delayed_flux(self, distance, arrival_time, energies, masses, flavor_xform):
        """
        Calculate the delayed flux of neutrinos accounting for time-of-flight delays.
        Neutrino masses can be a single Quantity or an array of Quantities.

        Parameters
        ----------
        arrival_time : astropy.Quantity
            Arrival time of neutrinos in the detector (same for all neutrinos, independent of energy).
        energies : array-like
            List of neutrino energies in MeV.
        masses : Quantity or array-like Quantities
            Neutrino mass in eV.
        flavor_xform : FlavorTransformation

        Returns
        -------
        dict
           Dictionary of delayed spectra, keyed by neutrino flavor. The dictionary structure is:
        delayed_spectrum[index in the array 'masses'][flavor][index in the array 'energies'] == flux
        """

        D = distance.to(u.kpc)
        c = 299792458 * u.m / u.s
        tof_light = D.to(u.m) / c # Time of flight for light 
        # print(t)
        # print(D)
        AT = arrival_time.to(u.ms)
        # print(T)
        E = energies.to(u.MeV)
        M = masses.to (u.eV)
         
        if M.isscalar:   # If masses is a single value, convert it to a list
            M = [M]

        reslut = {}
        emission_time={}
        delayed_spectrum={}
        for i, mass in enumerate(M):
            emission_time[i] = {}
            for j, energy in enumerate(E):
                # TOF delay formula from arXiv:1006.1889
                tof = 0.57 * (D.value / 10) ** 2 * (mass.value ** 2) * (30 / energy.value) ** 2 * u.ms
                # Calculate the emission time
                emission_time[i][j] = AT - tof

        # Calculate the transformed spectra at the emission time
        for i, mass in enumerate(M):
            delayed_spectrum[i] = {}
            for j, energy in enumerate(E):
                reslut = self.get_transformed_spectra(emission_time[i][j], E, flavor_xform)
                for flavor in Flavor:
                    delayed_spectrum[i][flavor] = reslut[flavor]

        return delayed_spectrum


    def get_delayed_flux(self, distance, arrival_time, energies, masses, flavor_xform):
            """
            Calculate the delayed flux of neutrinos accounting for time-of-flight delays.
            Neutrino masses can only be a single Quantity.

            Parameters
            ----------
            arrival_time : astropy.Quantity
              Arrival time of neutrinos in the detector (same for all neutrinos, independent of energy).
            energies : array-like
              List of neutrino energies in MeV.
            masses : float
              Neutrino mass in eV.
            flavor_xform : FlavorTransformation
              An instance from the flavor_transformation module.

            Returns
            delayed_spectrum[flavor][index in the array 'energies']== flux
            -------
            dict
                Dictionary of delayed spectra, keyed by neutrino flavor.
            """
            D = distance.to(u.kpc)
            c = 299792458 * u.m / u.s
            tof_light = D.to(u.m) / c
            # print(t)
            # print(D)
            AT = arrival_time.to(u.ms)
            # print(T)
            M = masses.to (u.eV)

            reslut = {}
            emission_time={}
            delayed_spectrum={}
    
            for j, energy in enumerate(energies):
                # TOF delay formula from arXiv:1006.1889
                tof = 0.57 * (D.value / 10) ** 2 * (M.value ** 2) * (30 / energy.value) ** 2 * u.ms
                  # Calculate the emission time
                emission_time[j] = AT - tof
                  
                   # Calculate the transformed spectra at the emission time

            for j, energy in enumerate(energies):
                reslut = self.get_transformed_spectra(emission_time[j], energies, flavor_xform)
                for flavor in Flavor:
                    delayed_spectrum[flavor] = reslut[flavor]

            return delayed_spectrum



