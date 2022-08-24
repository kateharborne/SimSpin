---
layout: default
title: FITS observation format
parent: Examples
nav_order: 2
last_modified_date: "Wed, 3 Aug 2022 14:08:00 AWST"
---

# Working with SimSpin FITS files 

SimSpin can generate astronomy standard FITS files through the [`build_datacube`](/SimSpin/docs/build_datacube) and [`write_simspin_FITS`](/SimSpin/docs/write_simspin_FITS) functions. In this example, we take you through how to work with these files. 
{: .fs-5 .fw-300 .pb-2 }
---




| EXT 	| Name     	| Description                                                                                                                                                                                                                            	|
|-----	|----------	|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| 1   	|          	| Header containing general information relating to the individual <br>galaxy observed, the code version run, and placeholder values to <br>maintain typical FITS layout.                                                                	|
| 2   	| DATA     	| The kinematic cube with axes x, y, v_los. Each z-axis bin <br>corresponds to a given LOS velocity, given by the axis labels. <br>Values within each bin correspond to the amount of r-band flux at <br>a given x-y projected location. 	|
| 3   	| OB_TABLE 	| A table that contains all of the SimSpin run information such that <br>a specific data cube can be recreated. Contains three columns (Name, <br>Value, Units).                                                                         	|
| 4   	| OB_FLUX  	| An image of the *observed* flux within the **r-band** in CGS units.                                                                                                                                                                    	|
| 5   	| OB_VEL   	| An image of the *observed* line-of-sight (LOS) velocity in units of km/s.                                                                                                                                                              	|
| 6   	| OB_DISP  	| An image of the *observed* LOS velocity dispersion in units of km/s.                                                                                                                                                                   	|
| 7   	| OB_H3    	| An image of the *observed* LOSVD higher order kinematic parameter, h3.                                                                                                                                                                 	|
| 8   	| OB_H4    	| An image of the *observed* LOSVD higher order kinematic parameter, h4.                                                                                                                                                                 	|
| 9   	| RAW_FLUX 	| An image of the raw particle **r-band** flux in CGS.                                                                                                                                                                                   	|
| 10  	| RAW_VEL  	| An image of the raw particle LOS velocities in units of km/s.                                                                                                                                                                          	|
| 11  	| RAW_DISP 	| An image of the raw particle LOS velocity dispersions in units of km/s.                                                                                                                                                                	|
| 12  	| RAW_AGE  	| An image of the raw particle stellar ages in units of Gyr.                                                                                                                                                                             	|
| 13  	| RAW_Z    	| An image of the raw particle metallicities in units of Z_sol.                                                                                                                                                                          	|
| 14  	| NPART    	| An image of the raw number of particles per pixel.                                                                                                                                                                                     	|