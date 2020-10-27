# MESP based optimizer written by Anmol Kumar

This program performs orientational optimization of molecule/s
in effect of electric field of a parent molecule. This helps in building molecular clusters and 
provides a good starting point for QM optimization of molecular clusters. 
QM optimization would rather be computationally time consuming, if initial geometry of the cluster is far away from minimum.

The primiary requirement for this program is .damqt and .dmqtv files of the
parent molecule. These files are binaries produced by DAMQT program upon initialization
of atomic density calculation and used for MED and MESP evaluation at any given point.
Please check usage of DAMQT to create these files. 


# Input File Structure for creating molecular cluster
"The molecules to be optimized can be provided in two different ways."

# === Input file ====

$OPTIONS

preprocfile="filename.xyz" # Provide geometry of all the molecules to be optimized in one file. The coordinates of these molecules will be used as is.

nocharge=.false. # if you do have charges in xyz file. Default is .true. If you do not have charges then the program will use obabel to obtain qtpie charges,    which are very similar to MESP charges.

tssize=0.5   #  Step size of translation

rssize=20    # Step size of rotation  This is an integer

$END

projectname  # Basename of parent molecule for which .damqt and .dmqtv files are present in the directory.

# === End of Input file ====

"==========OR============="
preprocfile argument can be replaced by following two arguments.

templatefile="filename.xyz"

insertlocfile="filen_containing_locations_where_template_coordinates_is_to_be_placed.xyz"

'''
"Molecule in the template file will be used and moved to location specified in insertlocfile.
Only "x" symbol in insertlocfile will be used to place the template molecule.
This allows easy use of cps-v.xyz file as insertlocfile.
You may or may not have charges as fifth column in xyz file of molecules to be optimized."
'''

# To install obabel "sudo apt-get install openbabel"

Run the program
DAMESPIMIZER.exe < mespimizer.inp
