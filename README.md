# PROTISTcode
Fortran code for the primary production module PROTIST of the aquatic ecosystem software Delft3D-WAQ. PROTIST describes the growth and motality of phytoplankton, protozooplankton and mixoplankton. The following is a short description of each file.

## PROTISTprocesses_pdf.txt
Contains the input parameters and state variables, the output parameters and the flux and their direction for each protist functional type.

## protistFunctions.f90
Code that contains all the functions for the physiological processes. These functions are combined to form the code for each protist functional type. 

## protistDiat.f90
Code required to simulate the growth and mortality of diatoms.

## protistDiatSedi.f90
Code required to simulate the sedimentation of diatoms.

## protistGreen. f90
Code required to simulate the growth and mortality of green algae (non-diatom phytoplankton).

## protistCM.f90
Code required to simulate the growth and mortality of consitutive mixoplankton. 

## protistNCM.f90
Code required to simulate the growth and mortality of non-consitutive mixoplankton.

## protistZoo.f90
Code required to simulate the growth and mortality of protozooplankton.

## protistAtten.f90
Code descrbing the proces of light attenuation by the phototrophic protist functional types.

