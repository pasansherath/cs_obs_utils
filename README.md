# cs_obs_utils (Controlled-source Ocean Bottom Seismograph Utilities)

### A set of Python utilities to process controlled source seismic data recorded by ocean bottom seismographs (OBSs).

## orient.py
orient.py is a Python3 script to estimate the angle between the airgun shot line and the horizontal components of OBSs. 
It uses a hodogram (particle motion) analysis and a polarisation analysis of the first arriving direct water wave from
the airgun shot recorded by the ocean bottom seismograph.

This code is a part from my PhD research work [Herath, 2021], inspired by the work of Flinn [1965], Perelberg [1994] and Eccles [2008].

#### Required libraries:
1. numpy
2. obspy
3. matplotlib
4. mpl_toolkits

It is recommended that you used a conda environmet with these libraries installed.

#### Required data:
1. SEGY format gathers for the two horizontal components of an OBS (H1, H2)

#### Other required parameters:
1. A trace header that can identify traces within the direct water arrival (first arrival) and its start and end values
2. Sliding window parameters for particle motion plotting and polarisation angle calculation

#### Outputs:
Plots of particle motion (hodogram) and polarisation angles on the horizontal plane of the direct water arrival to estimate 
the angle between the shot line and the horizontal components. 

#### Example: 
Two segy files (H1.segy and H2.segy) are provided to test and visualise the functionality.
To try this, clone the repository to your local computer and run the following command.
`python3 orient.py`

#### References:
1. Flinn, E. A. (1965). Signal Analysis Using Rectilinearity and Direction of Particle Motion. Proceedings of the IEEE, 53(12), 1874–1876. https://doi.org/10.1109/PROC.1965.4462
2. Perelberg, A. I., & Hornbostel, S. C. (1994). Applications of seismic polarization analysis. GEOPHYSICS, 59(1), 119–130. https://doi.org/10.1190/1.1443522
3. Eccles, J. D. (2008). Shear Wave Analysis of Volcanic Rifted Continental Margins in the North Atlantic. University of Cambridge.
4. Herath, P. (2021). Lithospheric Structure of the Hikurangi Subduction Margin [Victoria University of Wellington]. https://doi.org/10.26686/wgtn.14653491.v2
