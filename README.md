# SAARAMA (Single Amino Acid Ramachandran Tool)

SAARAMA helps you with calculating and plotting of φ (Phi) and ψ (Psi) angles in single amino acid simulations. So far SAARAMA 
only accepts GROMACS .gro topology files and .xtc trajectories but I'm happy to expend it to different trajectory file types
if needed. Instead of calculating torsion angles around peptide bonds the tool uses an approximation by using the H-atom for
φ and the O-atom for ψ. 

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
the help of Matplotlib. It takes 'naked' as well as N-terminal acetyl capped amino acids (ACE-X-NME). The script will automatically detect if the first residue is a acetyl residue and will treat it accordingly.
Commands should look like this:

```
saarama --topology /path/to/topology --xtc /path/to/trajectory [COMMAND]
```

## Example

Example files are provided in the example_files folder and can be used to check if the tool is working. The resulting plot coming from the Alanine data should look like the image below. It is visible that there's less range of motion around ψ because the O-atom has a double bound character and is therefore more stiff. The H-atom has a lot of space to wiggle around which explains φ range of motion.

![](https://github.com/Joshtron/saarama/blob/master/saarama_project/example_files/alanine.png)

Proline on the other hand looks different than other amino acids. The nitrogen is bound in a ring which results in a constrained φ angle. This is the reason Proline is usually excluded from standard Ramachandran plots. 

![](https://github.com/Joshtron/saarama/blob/master/saarama_project/example_files/proline.png)
