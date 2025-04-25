---
title: 'pyROX: Rapid Opacity X-sections'
tags:
  - Python
  - astronomy
  - exoplanets
  - planets
  - brown dwarfs
  - stars
  - atmospheres
  - opacities
  - cross-sections
  - spectroscopy
authors:
  - name:
      given-names: Sam
      non-dropping-particle: de
      surname: Regt
    orcid: 0000-0003-4760-6168
    corresponding: true
    affiliation: 1
  - name: Siddharth Gandhi
    orcid: 0000-0001-9552-3709
    affiliation: "2, 3"
  - name: Louis Siebenaler
    orcid: 0009-0005-3389-8819
    affiliation: 1
  - name: Darío González Picos
    orcid: 0000-0001-9282-9462
    affiliation: 1
affiliations:
  - name: Leiden Observatory, Leiden University, P.O. Box 9513, 2300 RA, Leiden, The Netherlands
    index: 1
  - name: Department of Physics, University of Warwick, Coventry CV4 7AL, UK
    index: 2
  - name: Centre for Exoplanets and Habitability, University of Warwick, Gibbet Hill Road, Coventry CV4 7AL, UK
    index: 3
date: 24 April 2025
bibliography: paper.bib

---

# Summary
In recent years, significant advances have been made in exoplanet and brown dwarf observations. By using state-of-the-art modelling software, astronomers can infer atmospheric properties, such as chemical and cloud compositions. Accurate and up-to-date opacities are critical to avoid inconclusive or biased results, but it can be challenging to compute opacity cross-sections from the so-called line lists provided by various online databases. 

We introduce `pyROX`, a user-friendly Python package to calculate molecular and atomic cross-sections. ... <!-- Summarise features (CPU, custom wave/P/T-grids) -->

# Statement of need
The advent of a new generation of telescopes and instruments has led to a dramatic increase in the quality of observations of exoplanets and brown dwarfs. For example, spectroscopic observations with JWST cover sub-stellar objects over a wide wavelength-range (1-20 um) that was previously difficult to access with ground- or space-based facilities [e.g. @Matthews_ea_2025; @Miles_ea_2023; @Rustamkulov_ea_2023]. Similarly, developments in direct-imaging instrumentation allow astronomers to measure the emission of exoplanet companions at closer separations to their host stars [e.g. @Landman_ea_2024; @Xuan_ea_2024]. Significant progress has also been made in atmospheric modelling with various software used for radiative transfer, chemistry, cloud-condensation, etc. Recently, these tools are coupled with sampling algorithms to infer properties of the observed objects [e.g. @Brogi_ea_2019; @Gibson_ea_2020]. 

Opacity cross-sections are a critical component to accurately modelling sub-stellar atmospheres. ... <!-- Energy transport --> Furthermore, studies at high spectral resolution require the frequencies of transition lines to be well-determined. Inaccuracies in line-list data can result in biased abundance constraints [e.g. @Brogi_ea_2019; @de_Regt_ea_2024] or suspected non-detections of certain molecules [e.g. @de_Regt_ea_2022]. For these reasons, it is important that the most up-to-date and complete opacity data are used when studying sub-stellar atmospheres. However, it can be challenging to calculate opacity cross-sections in a computationally efficient manner for line lists from different databases which can sometimes consist of billions of transitions. 

To help resolve this challenge, we present `pyROX` an easy-to-use Python package to calculate molecular and atomic cross-sections from the ExoMol, HITRAN/HITEMP, and Kurucz databases. As another important opacity source in sub-stellar atmospheres, `pyROX` supports calculations of Collision-Induced Absorption (CIA) coefficients from the HITRAN and Borysow databases. 

# Functionality of `pyROX`
Documentation for `pyROX` is available at [https://py-rox.readthedocs.io/en/latest/](https://py-rox.readthedocs.io/en/latest/) and includes tutorial examples to running the code. Here, we outline the main functionality of `pyROX`:

- **Download and read files**: The necessary input files (line lists, partition functions, broadening parameters or CIA files) can be downloaded with a simple command. When reading the relevant parameters, `pyROX` handles the different data structures of the supported databases (ExoMol, HITRAN/HITEMP, Kurucz, Borysow). 
- **Compute line-strengths and -widths**: For line-opacity calculations, `pyROX` calculates the strength and broadening-widths for each line transition at the user-provided pressure and temperature. Support is offered for various [pressure-broadening descriptions](https://py-rox.readthedocs.io/en/latest/notebooks/pressure_broadening.html).
- **Compute line profiles**: Next, `pyROX` computes the Voigt profiles as the real part of the Faddeeva function [Eq. 12 of @Gandhi_ea_2020], using the [`scipy.special.wofz`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wofz.html) implementation. Notably, the line profiles can be calculated on any user-provided wavelength or wavenumber grid. This enables cross-sections to be calculated at fixed resolutions, for instance, which can simplify comparisons with observations at high spectral resolutions.
- **Combine and save**: The line profiles are summed into wavelength-dependent cross-sections for each temperature-pressure point. These cross-sections are saved into an efficient HDF5 output file. For CIA calculations, `pyROX` restructures the coefficients read from the input files into a wavelength- and temperature-dependent grid and also saves these data to an HDF5 file.

# Similar tools
Existing open source codes, such as [`Cthulhu`](https://github.com/MartianColonist/Cthulhu) [@Agrawal_ea_2024], [`ExoCross`](https://github.com/Trovemaster/exocross) [@Yurchenko_ea_2018; @Zhang_ea_2024] and [`HELIOS-K`](https://github.com/exoclime/HELIOS-K) [@Grimm_ea_2015; @Grimm_ea_2021], can calculate cross-sections at a similar performance to `pyROX`. 
<!-- Python + CPUs-only + custom wavelength-grid -->

<!-- 
# Future developments 
VALD?
Output conversions to other radiative transfer codes (besides pRT)?
We welcome suggestions...
-->

<!-- 
# Citations
...
 -->

# Documentation
Documentation for pyROX is available at [https://py-rox.readthedocs.io/en/latest/](https://py-rox.readthedocs.io/en/latest/).

<!-- # Acknowledgements -->

# References

<!-- Summary: Has a clear description of the high-level functionality and purpose of the software for a diverse, non-specialist audience been provided? -->
<!-- A statement of need: Does the paper have a section titled ‘Statement of need’ that clearly states what problems the software is designed to solve, who the target audience is, and its relation to other work? -->
<!-- State of the field: Do the authors describe how this software compares to other commonly-used packages? -->
<!-- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. -->
