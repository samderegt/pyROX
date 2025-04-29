==========================
Configuration per Database
==========================

pyROX supports calculations from several online databases. Molecular and atomic line opacities can be computed from the `ExoMol <https://www.exomol.com/>`_, `HITRAN <https://hitran.org/>`_, `HITEMP <https://hitran.org/hitemp/>`_, and `Kurucz <http://kurucz.harvard.edu/>`_ line lists. Collision-Induced Absorption (CIA) can be computed from the `HITRAN <https://hitran.org/cia/>`_ and `Borysow <https://www.astro.ku.dk/~aborysow/programs/index.html>`_ databases. 

The different formats of these databases sometimes require or allow different configuration parameters to be given to pyROX. Here, we describe some notable differences between the databases. Examples of the configuration files are available on the `GitHub repository <https://github.com/samderegt/pyROX/blob/main/examples/>`_.

Downloading
-----------
pyROX can automatically download (most) required data files from the online databases. The ``urls`` list in the configuration file should contain the following urls:

- **ExoMol**: The line list's definition file (``*.def.json``) should be given in the ``urls`` list. The ``.def.json``-file contains information to download the line list's states, transitions, and partition function files (``.states``, ``.trans``, ``.pf``). Adding the broadening files (``*.broad``) to the ``urls`` list is optional, but recommended for ExoMol line lists (if available). 
- **HITRAN/HITEMP**: The partition function file (``q*.txt``) should be given in the ``urls`` list. You can find the url of the respective isotopologue on the `HITRAN website <https://hitran.org/docs/iso-meta/>`_ under "Documentation" > "Isotopologues". We also recommend adding ExoMol's ``.broad``-files to the ``urls`` list. The ``.par``-file contains all transitions and should be downloaded as well. However, it requires a login, so we recommend downloading this manually from the `HITRAN/HITEMP website <https://hitran.org/>`_ using your account login. 
- **Kurucz**: If the ``species`` name is provided, pyROX can automatically download the transitions (``gf*.pos``) and states (``*.tsv``, from `NIST <https://physics.nist.gov/PhysRefData/ASD/levels_form.html>`_) files.
- **CIA_HITRAN/CIA_Borysow**: For the CIA calculations, you should specify the exact urls from the `HITRAN website <https://hitran.org/cia/>`_ or Aleksandra `Borysow's website <https://www.astro.ku.dk/~aborysow/programs/index.html>`_.

Input-data files
----------------
pyROX should be run with a ``files`` dictionary pointing to the downloaded input data in the ``input_data_dir`` directory. The following files are required for each database:

 - **ExoMol**: 

   - ``transitions`` (*list* or *str*): Line transitions in ``*.trans``-file (or compressed as ``bz2``).
   - ``states`` (*str*): State energies in ``*.states``-file (or compressed as ``bz2``).
   - ``partition_function`` (*str*): Partition function in ``*.pf``-file.
 - **HITRAN/HITEMP**: 

   - ``transitions`` (*list* or *str*): Line transitions in ``*.par``-file (or compressed as ``bz2``).
   - ``partition_function`` (*str*): Partition function in ``q*.txt``-file.
 - **Kurucz**: 

   - ``transitions`` (*str*): Line transitions in ``gf*.pos``-file.
   - ``states`` (*str*): ``*.tsv``-file (from NIST), which is used for calculating the partition function.
 - **CIA_HITRAN/CIA_Borysow** (see also `"Collision-Induced Absorption" <./notebooks/collision_induced_absorption.html>`_): 

   - ``cia`` (*list of tuples* or *list of strings*): If list of tuples, the expected structure is (``cia_file``, ``temperature_mask``, ``wavenumber_mask``). If list of strings, only filenames are expected. 

Pressure-broadening description
-------------------------------
The line list databases can use different pressure-broadening descriptions (see also `"Pressure Broadening" <./notebooks/pressure_broadening.html>`_). 

 * **ExoMol**: The broadening coefficients can be made to depend on the rotational quantum numbers J. This is handled for the "a0" and "m0" broadening diets if ``.broad``-files are given (see `"Coefficients from a file" <./notebooks/pressure_broadening.html#Coefficients-from-a-file>`_). Alternatively, you can provide a callable function that takes the lower state J-quantum number as input and returns the broadening coefficients ``gamma`` and ``n`` (see `"Parameterised coefficients" <./notebooks/pressure_broadening.html#Parameterised-coefficients>`_). 
 * **HITRAN/HITEMP**: Currently only handles constant broadening coefficients, that is the same for each transition. If a ``.broad``-file is given, pyROX will use the average values.
 - **Kurucz**: can use Schweitzer et al. 1996 for alkalis, otherwise uses Sharp & Burrows 2007 (can give polarisability alpha/C)?