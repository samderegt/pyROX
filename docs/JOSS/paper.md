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
date: 28 April 2025
bibliography: paper.bib

---

# Summary
In recent years, significant advances have been made in exoplanet and brown dwarf observations. By using state-of-the-art models, astronomers can determine properties of their atmospheres, such as temperatures, the presence of clouds, or the chemical abundances of molecules and atoms. Accurate and up-to-date opacities are crucial to avoid inconclusive or biased results, but it can be challenging to compute opacity cross-sections from the line lists provided by various online databases. 

We introduce `pyROX`, an easy-to-use Python package to calculate molecular and atomic cross-sections. Since `pyROX` works on CPUs, it can compute a small line list on a regular workstation, but it is also easily parallelised on a cluster for larger line lists. In addition to line opacities, `pyROX` also supports calculations of collision-induced absorption. Tutorials are provided in the online documentation, explaining the configuration parameters and different functionalities of `pyROX` (e.g. custom wavelength-, pressure-, or temperature-grid, pressure-broadening descriptions, etc.). 


# Statement of need
The advent of a new generation of telescopes and instruments has led to a dramatically increased quality in observations of exoplanets and brown dwarfs. Such sub-stellar objects are now observed over a wide wavelength range (1-20 µm) with JWST spectra [e.g. @August_ea_2023; @Carter_ea_2024; @Matthews_ea_2025; @Miles_ea_2023], for instance, which was previously difficult to access. Developments in ground-based instrumentation allow astronomers to measure young exoplanet companions at closer separations to their host stars [e.g. @Landman_ea_2024; @Xuan_ea_2024b] and at high spectral resolutions [e.g. @Nortmann_ea_2025; @Xuan_ea_2024a]. At the same time, progress has also been made in atmospheric modelling using software for radiative transfer, chemistry, climate models, etc. [e.g. @Molliere_ea_2019; @Stock_ea_2018; @Wardenier_ea_2021]. Recently, these observations and software are coupled with sampling algorithms to characterise the atmospheres of the sub-stellar objects [e.g. @Barrado_ea_2023; @Brogi_ea_2019; @Gibson_ea_2020; @Line_ea_2015]. 

Opacity cross-sections play a key role in accurately modelling sub-stellar atmospheres. Opacity governs the dominant energy transport mechanism (i.e. radiative or convective) which affects the thermal structure of the atmosphere [@Marley_ea_2021]. <!-- Explain how/when this can be an issue --> Furthermore, high-resolution studies require well-determined frequencies for the transition lines. Inaccuracies in line-list data can result in biased abundance constraints [e.g. @Brogi_ea_2019; @de_Regt_ea_2024] or ambiguous (non)-detections of certain molecules [e.g. @de_Regt_ea_2022; @Merritt_ea_2020; @Serindag_ea_2021]. As such, it is important that the most up-to-date and complete opacity data are used. However, it can be difficult to efficiently calculate opacity cross-sections from line lists that sometimes consist of billions of transitions. 

To help resolve this challenge, we present `pyROX` a user-friendly Python package to calculate molecular and atomic cross-sections for applications in models of sub-stellar atmospheres. `pyROX` supports line opacity calculations from the ExoMol [@Tennyson_ea_2024], HITRAN [@Gordon_ea_2020], HITEMP [@Rothman_ea_2010], and Kurucz[^1] databases. Collision-Induced Absorption (CIA) coefficients can also be calculated from the HITRAN and Borysow[^2] databases.

[^1]: http://kurucz.harvard.edu/
[^2]: https://www.astro.ku.dk/~aborysow/programs/index.html

<!-- Work pyROX or a predecessor has already been used for? [@Gonzales_Picos_ea_2024?; @Siebenaler_ea_2025; @de_Regt_ea_2025] -->


# Functionality of `pyROX`
Documentation for `pyROX` is available at [https://py-rox.readthedocs.io/en/latest/](https://py-rox.readthedocs.io/en/latest/) and includes tutorial examples to running the code. Here, we outline the main functionality of `pyROX`:

- **Download and read files**: The necessary input files (line lists, partition functions, broadening parameters or CIA files) can be downloaded with a simple command. When reading the relevant parameters, `pyROX` handles the different data structures of the supported databases (ExoMol, HITRAN/HITEMP, Kurucz, Borysow). 
- **Compute line-strengths and -widths**: For line-opacity calculations, `pyROX` calculates the strength and broadening-widths for each line transition at the user-provided pressure and temperature. Support is offered for various [pressure-broadening descriptions](https://py-rox.readthedocs.io/en/latest/notebooks/pressure_broadening.html).
- **Compute line profiles**: Next, `pyROX` computes the Voigt profiles as the real part of the Faddeeva function [Eq. 12 of @Gandhi_ea_2020], using the [`scipy.special.wofz`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.wofz.html) implementation. 
- **Combine and save**: The line profiles are summed into wavelength-dependent cross-sections for each temperature-pressure point. These cross-sections are saved into an efficient HDF5 output file. For CIA calculations, `pyROX` restructures the coefficients read from the input files into a wavelength- and temperature-dependent grid and also saves these data to an HDF5 file.

Currently, `pyROX` offers built-in support for converting its output into the high-resolution opacities used by petitRADTRANS [@Molliere_ea_2019]. In future releases, we plan to add conversions for other radiative transfer codes that are popular in the exoplanet and brown dwarf community. We welcome suggestions for new features, which can be done by [opening an issue](https://github.com/samderegt/pyROX/issues) on GitHub. If you want to contribute to `pyROX`, please read the [documented guidelines](https://py-rox.readthedocs.io/en/latest/contributing.html).


# Similar tools
Existing open source codes, such as [`Cthulhu`](https://github.com/MartianColonist/Cthulhu) [@Agrawal_ea_2024], [`ExoCross`](https://github.com/Trovemaster/exocross) [@Yurchenko_ea_2018; @Zhang_ea_2024] and [`HELIOS-K`](https://github.com/exoclime/HELIOS-K) [@Grimm_ea_2015; @Grimm_ea_2021], can calculate cross-sections at comparable performances to `pyROX`. However, `ExoCross` is written in Fortran and `HELIOS-K` utilises GPU-acceleration which can limit their use to experts with the appropriate hardware. `pyROX` is a Python code that runs only on CPUs which should make it accessible for the opacity needs of most astronomers. Notably, `pyROX` supports cross-section calculations on any user-provided wavelength or wavenumber grid. This enables the user to fix the spectral resolution ($\mathcal{R}=\lambda/\Delta\lambda$) which cannot be achieved with equal wavelength- or wavenumber-spacing.

<!-- 
# Citations
...
 -->

# Acknowledgements
S.d.R. and D.G.P. acknowledge funding from NWO grant OCENW.M.21.010.

# References

<!-- Summary: Has a clear description of the high-level functionality and purpose of the software for a diverse, non-specialist audience been provided? -->
<!-- A statement of need: Does the paper have a section titled ‘Statement of need’ that clearly states what problems the software is designed to solve, who the target audience is, and its relation to other work? -->
<!-- State of the field: Do the authors describe how this software compares to other commonly-used packages? -->
<!-- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it. -->
