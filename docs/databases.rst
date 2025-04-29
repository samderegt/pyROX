======================
Database Configuration
======================

pyROX supports calculations from several online databases. Molecular and atomic line opacities can be computed from the `ExoMol <https://www.exomol.com/>`_, `HITRAN <https://hitran.org/>`_, `HITEMP <https://hitran.org/hitemp/>`_, and `Kurucz <http://kurucz.harvard.edu/>`_ line lists. Collision-Induced Absorption (CIA) can be computed from the `HITRAN <https://hitran.org/cia/>`_ and `Borysow <https://hitran.org/cia/>`_ databases. 

The different formats of these databases sometimes require or allow different configuration parameters to be given to pyROX. Here, we describe some notable differences between the databases. Examples of the configuration files are available on the `GitHub page <https://github.com/samderegt/pyROX/blob/main/examples/>`_.

Downloading
-----------
pyROX can automatically download (most) required data files from the online databases. The ``urls`` list in the configuration file should contain the following urls:

 - **ExoMol**: The line list's definition file (``*.def.json``) should be given in the ``urls`` list. The ``.def.json``-file contains information to download the line list's states, transitions, and partition function files (``.states``, ``.trans``, ``.pf``). Adding the broadening files (``*.broad``) to the ``urls`` list is optional, but recommended for ExoMol line lists (if available). 
 - **HITRAN/HITEMP**: The partition function file (``q*.txt``) should be given in the ``urls`` list. You can find the url of the respective isotopologue on the `HITRAN website <https://hitran.org/docs/iso-meta/>`_ under "Documentation" > "Isotopologues". We also recommend adding ExoMol's ``.broad``-files to the ``urls`` list. The ``.par``-file contains all transitions and should be downloaded as well. However, it requires a login, so we recommend downloading this manually from the `HITRAN/HITEMP website <https://hitran.org/>`_ using your account login. 
 - **Kurucz**: For Kurucz line lists, pyROX can automatically download the transitions (``*.pos``) and states (``*.tsv``, from `NIST <https://physics.nist.gov/PhysRefData/ASD/levels_form.html>`_) files. Providing the ``species`` name is sufficient to download the required files.
 - **CIA_HITRAN**: ...
 - **CIA_Borysow**: ...

Input data files
----------------
...

 - **ExoMol**: .trans(.bz2), .states(.bz2), .pf
 - **HITRAN/HITEMP**: .par(.bz2), Q<aa>.txt
 - **Kurucz**: .pos, .tsv (NIST, used for pf)
 - **CIA_HITRAN**: ...
 - **CIA_Borysow**: ...

Pressure-broadening information
-------------------------------
...

 - **ExoMol**: J-dependent broadening (.broad diets or callable function)
 - **HITRAN/HITEMP**: can use ExoMol's .broad files, but assumes a constant coefficient
 - **Kurucz**: can use Schweitzer et al. 1996 for alkalis, otherwise uses Sharp & Burrows 2007 (can give polarisability alpha/C)?