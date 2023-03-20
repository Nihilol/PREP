import numpy
import gdspy
from matplotlib import pyplot as plt


gdsii = gdspy.GdsLibrary(infile="1x4_and_2x3.gds")
main_cell = gdsii.top_level()[0]  # Assume a single top level cell
print(main_cell)
points = main_cell.polygons[0].polygons[0]
for p in points:
    print("Points: {}".format(p)) 

gdspy.LayoutViewer(gdsii)