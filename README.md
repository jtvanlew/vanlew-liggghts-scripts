# vanlew-liggghts-scripts
'modularized' scripts for running liggghts sims

This script file (potentially this will be plural in teh future) contains some functions to make liggghts input scripts more compact. They are rather static right now and mostly deal with the 'box' configurations I deal with in ceramic breeder blankets for fusion. There's not much documentation for the script itself but I try to make the function names self-explanatory so they can be read through and adapted to future simulations.

in your .bashrc file (or whatever you have), add this directory so you can import vl_dem from python. For instance, on my machine I have:

export PYTHONPATH=$PYTHONPATH:/home/jon/LIGGGHTS/vanlew-liggghts-scripts

then you can do 

  . ~/.bashrc

at the terminal to source the bash file and then from python you can say 

  from vl_dem import *