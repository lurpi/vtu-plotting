# vtu-plotting
Plotting OGS6 results from vtu with python scripts

# 2d-slice-sigma.py
Script to plot two sets of graphs , one with principal stresses 1&3, one with principal stresses 1&2

To change the section plotted: go to line 187 - 192: first the filter is along the y-axis (distance from target "y_section" smaller than 0.2 meters), then exclude poins inside the tunnel radius (assumed to be 1.55 m).

For visualization:
- the scale of x and y is set at lines 294-295 and 397-398
- the size of the arrows is defined by the value "scale_par"
- the frequency of the arrows can be changed by multipling the value "ji" in the lines 267 and 374 (now is plotting all the values in the plot for 1&3, only 1 point every 3 for the plot 1&2).
