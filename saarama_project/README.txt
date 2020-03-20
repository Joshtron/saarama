Please make sure you use Python 3 and that MDAnalysis is installed. This can be done with pip or conda.

>> pip install MDAnalysis 

>> conda install -c conda-forge mdanalysis  

Now run setup.py to set up everything.

>> python setup.py develop

You are now able to run saarama on the command line.

>> saarama --topology /path/to/saarama_project/example_files/md_0_1.gro --xtc /path/to/saarama_project/example_files/md_0_1_noPBC.xtc plot
