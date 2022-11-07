############################################################################
# VESDES (Volcanic Event Spatial Density Estimation Suite)
#
# VESDES includes a number of kernel functions and bandwidth 
# estimation techniques used to estimate spatial density of 
# volcanic events. VESDES is written in the C programming
# language and includes calls to the R statistical package.
#
#    Copyright (C) 2021-2022  
#    Laura Connor (ljconnor@gmail.com)
#    Charles Connor (chuck.connor@gmail.com)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################    

#Makefile for compiling spd_2022 logL-inversion

# This makefile file compiles spatial_density.c 
# into the executable: spd_2022
#
# Linking and compiling variables
# Alter as needed for your system.
#
export CC = gcc
export INSTALLPATH = .

all clean install uninstall spd_2022:
	$(MAKE) -C src $@
