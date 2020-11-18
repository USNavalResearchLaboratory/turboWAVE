import pathlib
import numpy as np
import matplotlib.pyplot as plt
import twutils.plot

# Simple example of how to use the plotter class
# to develop your own visualization script

file_to_plot = pathlib.Path.home() / 'Run' / 'Ez.npy'

# First make a plotter object tied to a specific file.
# Make it buffered if you want to load all data at once.
# For huge files make it unbuffered.
plotter = twutils.plot.plotter(file_to_plot,units='cgs',buffered=False)
plotter.display_info()
plt.figure(1,figsize=(10,4),dpi=75)

# 2D Figure Example
plt.subplot(121)
# slicing_spec : first two characters are axes to plot, second two are the slice axes
# slice_to_plot : select the slices to plot (order is whatever is in slicing_spec)
data_slice,plot_dict = plotter.falsecolor2d(slicing_spec='zxyt',slice_to_plot=(0,-1))
plt.imshow(data_slice,aspect='auto',origin='lower',extent=plot_dict['extent'],vmin=plot_dict['vmin'],vmax=plot_dict['vmax'],cmap='jet')
bar = plt.colorbar()
plt.xlabel(plot_dict['xlabel'],fontsize=18)
plt.ylabel(plot_dict['ylabel'],fontsize=18)
bar.set_label(plot_dict['blabel'],fontsize=18)

# Lineout Example
plt.subplot(122)
abcissa,ordinate,plot_dict = plotter.lineout(slicing_spec='zxyt',slice_to_plot=(int(plot_dict['dims'][1]/2),0,-1))
plt.plot(abcissa,ordinate)
plt.xlabel(plot_dict['xlabel'],fontsize=18)
plt.ylabel(plot_dict['ylabel'],fontsize=18)

plt.tight_layout()
plt.show()
