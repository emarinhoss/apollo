##
# Configure variables for WarpM build environment.  All of these could alteranitively
# be passed at the SCons command line.  Command line arguments of the same name will
# overide what is set here.
##

# Settings are host specific
import os
import sys

########################################
###########  TEMPLATES  ################
########################################

# Generic Ubuntu Install following the Wiki page
mpi_base = "/usr"
petsc_base = "/home/sousae/software/petsc-moab"
boost_base = "/usr"
#---------------------------------------------------
