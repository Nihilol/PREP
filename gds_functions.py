import numpy as np
import gdspy
import gdstk
import matplotlib as plt


def gds_cut(read_file, x_lims, y_lims, out_file, view = False):
    loaded_lib = gdspy.GdsLibrary(infile = read_file)
    new_lib = gdspy.GdsLibrary()
    main_cell = loaded_lib.top_level()[0] 
    cell1 = new_lib.new_cell('Sliced')
    for i in main_cell.polygons:
        i1 = gdspy.slice(i, x_lims, axis = 0, layer = i.layers[0])
        i2 = gdspy.slice(i1[1], y_lims, axis = 1, layer = i.layers[0])
        cell1.add(i2[1])
    
    new_lib.write_gds(out_file)
    if view:
        new_gds = gdspy.GdsLibrary(infile = out_file)
        gdspy.LayoutViewer(new_gds)

class gate:
    def __init__(self, layer, x_placement, y_placement, voltage):
        self.state = 0
        self.layer = layer
        self.outline = np.array([x_placement, y_placement])
        self.voltage = voltage
        
    def plot_outline(self, axis):
        axis.plot(self.outline[0], self.outline[1])