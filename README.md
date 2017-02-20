This repository contains the 'grid version' of SAGE (https://github.com/darrencroton/sage). 

post_processed_SAGE: This is functionally identical to the original SAGE with the output being the main difference.  In this version, each galaxy has a number of grids that tracks properties (such as Stellar Mass) through redshift.  Due to the increase in the output files, it is suggested (and coded with this in mind) that only a single snapshot (the smallest redshift generally) should be saved. Furthermore, due to the fact that base SAGE disregards merged galaxies for future evolution and saving, a separate output with the tag 'MergedGalaxies' is produced for all the merged galaxies. This is important because these galaxies exist at high redshift but are not present at low redshift (where the save routine is run).

post_processed_grid: This codebase takes the modified output of post_processed_SAGE and creates ionizing photon grids. This step was separated to the running of SAGE as the user can now alter the ionizing photon prescriptions without having to run the entirety of SAGE again.  The output of this code is named in such a way to be used with Anne Hutter's semi-numerical reionization code (https://github.com/annehutter/grid-model).

Any issues/questions/ranting can be sent to Jacob Seiler: jseiler@swin.edu.au 
