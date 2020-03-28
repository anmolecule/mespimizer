####################################################
MESP based optimizer written by Anmol Kumar
####################################################

This program performs orientation optimization of molecule/s
in effect of MESP of a parent molecule. This provides a good starting point
for QM optimization of molecular clusters, which would rather be computationally 
time consuming, if initial geometry of the clsuter is far away from minimum.

The primiary requirement for this program is .damqt and .dmqtv files related to
parent molecule. These files are binaries produced by DAMQT upon initialization
of atomic density calculation and used for MED and MESP evaluation at any given point.

Input File Structure
$OPTIONS
#The molecules to be moptimized can be provided in two different ways.

######a. 
preprocfile="filename.xyz"
#   Provide geometry of all the molecules to be optimized in one file. 
#   The coordinates of these molecules will be used as is.
#
#    =========OR=============
#

######b. 
templatefile="filename.xyz"
insertlocfile="filen_containing_locations_where_template_coordinates_is_to_be_placed.xyz"

#   Molecule in the template file will be used and moved to location specified in insertlocfile.
#   Only "x" symbol in insertlocfile will be used to place the template molecule.
#   This allows easy use of cps-v.xyz file as insertlocfile.

#You may or may not have charges as fifth column in xyz file of molecules to be moptimized.

nocharge=.false. # if you do have charges in xyz file. Default is .true.

#If you do not have charges then the program will use obabel to obtain qtpie charges,
#which are very similar to MESP charges.
#To install obabel "sudo apt-get install openbabel"

tssize=0.5   #  Step size of translation
rssize=20    # Step size of rotation  This is an integer
$END
projectname  # Basename of parent molecule for which .damqt and .dmqtv files are present in the directory.


Run the program
../../DAMESPIMIZER.exe < mespimizer.inp
