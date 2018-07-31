import numpy as np
import matplotlib.pyplot as plt
import twutils.plot as twplot

# Simple example of how to use the plotter

file_to_plot = '/Users/gordon/Run/air.dvdat'

# Make it buffered if you want to have time available as a plotting axis
# For huge files make it unbuffered
plotter = twplot.plotter(file_to_plot,buffered=False)
plotter.display_info()

# 2D Figure Example
plt.figure(1,figsize=(7,5),dpi=75)
# slicing_spec : first two characters are axes to plot, second two are the slice axes
# slice_to_plot : select the slices to plot (order is whatever is in slicing_spec)
data_slice,plot_dict = plotter.falsecolor2d(slicing_spec='xyzt',slice_to_plot=(64,-1),dyn_range=2.0)
plt.imshow(data_slice,origin='lower',extent=plot_dict['extent'],vmin=plot_dict['vmin'],vmax=plot_dict['vmax'],cmap='jet')
plt.colorbar(label=file_to_plot.split('/')[-1])
plt.xlabel(plot_dict['xlabel'],fontsize=18)
plt.ylabel(plot_dict['ylabel'],fontsize=18)

# Lineout Example
plt.figure(2,figsize=(7,5),dpi=75)
abcissa,ordinate,plot_dict = plotter.lineout(slicing_spec='zxyt',slice_to_plot=(64,64,-1))
plt.plot(abcissa,ordinate)
plt.xlabel(plot_dict['xlabel'],fontsize=18)
plt.ylabel(plot_dict['ylabel'],fontsize=18)

plt.show()
