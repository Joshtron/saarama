# SAARAMA (Single Amino Acid Ramachandran Tool)

SAARAMA helps you with calculating and plotting of φ (Phi) and ψ (Psi) angles in single amino acid simulations. So far SAARAMA 
only accepts GROMACS .gro topology files and .xtc trajectories but I'm happy to expend it to different trajectory file types
if needed. Instead of calculating torsion angles around peptide bonds the tool uses an approximation by using the H-Atom for
φ and the O-Atom for ψ. 

* φ: H1-N-CA-C
* ψ: N-CA-C-OT1

## Installation

First you need to clone the repository:

```
git clone https://github.com/Joshtron/saarama
```

Please make sure to install MDAnalysis. pip and conda commands are listed below:

```
pip install mdanalysis
```

```
conda install -c conda-forge mdanalysis
```

After successful MDAnalysis installation, call the setup.py script which will enable you to use the SAARAMA command line 
tool. It will take care of dependencies such us matplotlib, click and MDAnalysis.

```
python /path/to/saarama_project/setup.py develop
```

## Command line tool

Like mentioned above, SAARAMA uses a topology and a trajectory file to calculate the torsion angles over the different time 
frames. It can either print the angles to stdout [terminal] or plot the angles in a Ramachandran plot fashion [plot] with 
the help of Matplotlib. Commands should look like this:

```
saarama --topology /path/to/topology --xtc /path/to/trajectory [COMMAND]
```

## Example

Example files are provided in the example_files folder and can be used to check if the tool is working. The resulting plot
should look like this:

![](https://github.com/Joshtron/saarama/blob/master/saarama_project/example_plot.png)
