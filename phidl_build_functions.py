from phidl import Device, Layer, LayerSet, make_device
from phidl import quickplot as qp # Rename "quickplot()" to the easier "qp()"
import phidl.geometry as pg
import phidl.routing as pr
import phidl.utilities as pu
import phidl
import numpy as np
import phidl.path as pp
from phidl import Path, CrossSection, Device
import phidl.path as pp
from phidl import Group

lys = LayerSet()
lys.add_layer('die', gds_layer = 0, gds_datatype = 0, alpha = 0.15)
lys.add_layer('marker', gds_layer = 1, gds_datatype = 0)
lys.add_layer('SiO2', gds_layer = 2, gds_datatype = 0)
lys.add_layer('ohmic Al', gds_layer = 3, gds_datatype = 0)
lys.add_layer('GL1 Pd/Ti', gds_layer = 4, gds_datatype = 0)
lys.add_layer('GL2 Pd/Ti', gds_layer = 5, gds_datatype = 0)
lys.add_layer('wiring ohmic Al', gds_layer = 6, gds_datatype = 0)
lys.add_layer('wiring GL1 Pd/Ti', gds_layer = 7, gds_datatype = 0)
lys.add_layer('wiring GL2 Pd/Ti', gds_layer = 8, gds_datatype = 0)
lys.add_layer('text', gds_layer = 9, gds_datatype = 0)
lys.add_layer('screening', gds_layer = 10, gds_datatype = 0)

scale_factor = 1

radius_dot = 120/2*scale_factor #nm

width_external_box = 50000*scale_factor
height_external_box = 50000*scale_factor
N_ports_per_edge = 12
width_external_port = 1000*scale_factor
spacing_h_external_port = (width_external_box - N_ports_per_edge*width_external_port)/(N_ports_per_edge+1)
spacing_v_external_port = (height_external_box - N_ports_per_edge*width_external_port)/(N_ports_per_edge+1)

xmax_ext_box = width_external_box/2
xmin_ext_box = -width_external_box/2

ymax_ext_box = 2*radius_dot + height_external_box/2 #- height_external_box/20
ymin_ext_box = 2*radius_dot - height_external_box/2 #- height_external_box/20

width_inner_marker = 500 #nm
length_inner_marker = 10e3 #nm

width_outer_marker = 5e3 #nm
length_outer_marker = 180e3 #nm

edge_box_cut_outer_marker = 12e3 #nm

width_local_marker = 1000 #nm
length_local_marker = 20e3 #nm
length_to_local = 270e3
layer_marker =  lys['marker']

width_bonding_pad = 130e3 #nm
length_bonding_pad = 130e3

width_bonding_box = 2.1e6*scale_factor #it's actually the inner box, not the outer one
height_bonding_box = 2.1e6*scale_factor

layer_dot_bond_pad = lys['wiring GL2 Pd/Ti']

layer_barrier_gates_bond_pad = lys['wiring GL1 Pd/Ti']
layer_screening_gates_bond_pad = layer_barrier_gates_bond_pad
layer_ohmic_bond_pad = lys['wiring ohmic Al']

overlap_ext_ports_and_bond_fanout = 1e3 #nm edge of square to be partly subtracted

width_port_bonding_pad = width_bonding_pad/3*1
width_port_ohmic_bonding_pad = width_bonding_pad*2/3

width_external_port = 1000*scale_factor

width_middle_port = width_external_port - 600

closing_distance = 600

protection_extra = 20e3
length_extension_protection_bonding_pad = 15e3 #nm
width_extension_protection_bonding_pad = 15e3 #nm
layer_protection_bonding_pad = lys['SiO2']

def x_midpoint_ext_port_h(N):
    x0 = spacing_h_external_port*N+width_external_port*(N-1)+ width_external_port/2 + xmin_ext_box
    return x0
def y_midpoint_ext_port_v(N):
    y0 = spacing_v_external_port*N+width_external_port*(N-1)+ width_external_port/2 + ymin_ext_box
    return y0

def inner_marker_generator(width, length, layer):

    cross_inner_marker = pg.cross(length = length, width = width, layer = layer)

    inner_marker_joined = pg.union(cross_inner_marker, by_layer = True)

    return inner_marker_joined

def outer_marker_generator():
    cross_outer_marker =  pg.cross(length = length_outer_marker, width = width_outer_marker, layer = layer_marker)

    outer_marker_device_joined = pg.union(cross_outer_marker, by_layer = True)

    rect_box_cut_outer_marker = pg.rectangle(size = (edge_box_cut_outer_marker, edge_box_cut_outer_marker), layer = layer_marker)
    rect_box_cut_outer_marker.center = [0,0]

    outer_marker_joined_cut = pg.boolean(outer_marker_device_joined,rect_box_cut_outer_marker,'A-B', layer = layer_marker)

    return outer_marker_joined_cut

def place_sem_structure(device, marker_structure, sem_structure, ref_structure):

    ref_dummy_devices_1 = device.add_ref(marker_structure)
    ref_dummy_devices_1.xmax = ref_structure.xmax
    ref_dummy_devices_1.ymin = ref_structure.ymin

    ref_dummy_devices_2 = device.add_ref(sem_structure)

    ref_dummy_devices_2.xmin = ref_structure.xmin
    ref_dummy_devices_2.ymin = ref_structure.ymin

    ref_dummy_devices_3 = device.add_ref(sem_structure)

    ref_dummy_devices_3.xmax = ref_structure.xmax
    ref_dummy_devices_3.ymax = ref_structure.ymax

    ref_dummy_devices_4 = device.add_ref(sem_structure)

    ref_dummy_devices_4.xmin = ref_structure.xmin
    ref_dummy_devices_4.ymax = ref_structure.ymax



def marker_device_generator(x, y):

    layer = lys['marker']

    amount_of_markers = 10
    marker_spacing = 200000

    marker_device = Device()

    inner_markers, outer_markers = [], []

    for i in range(amount_of_markers):
        inner_markers.append("inner_left_" + str(i))
        outer_markers.append("outer_left_" + str(i))
        inner_markers[i] = marker_device << inner_marker_generator(width = width_inner_marker, length = length_inner_marker, layer = layer)
        outer_markers[i] = marker_device << outer_marker_generator()

    def marker_mover(list1, list2, i, x, y):
        list1[i].move(destination = (x*marker_spacing, y *marker_spacing))
        list2[i].move(destination = (x*marker_spacing, y *marker_spacing))

    for j in range(amount_of_markers):
        marker_mover(inner_markers, outer_markers, j, x[j], y[j])

    return marker_device


def protection_pads(device, width, length, destination):
    protection_pad = device << pg.rectangle(size = (width, length), layer = layer_protection_bonding_pad)
    protection_pad.move(destination = destination)
    return protection_pad

def create_and_place_dummy_bond_pads(device, amount_dummy_bonding_pads, width, length, spacing, direction, destination):
    layer_dummy_bonds = [lys['wiring ohmic Al'], lys['wiring GL2 Pd/Ti'], lys['wiring GL1 Pd/Ti']]
    dummy_bonding_pads = []
    j = 0
    for i in range(amount_dummy_bonding_pads):
        if j == 3:
            j = 0
        dummy_bonding_pads.append(device << pg.rectangle(size = (length, width), layer = layer_dummy_bonds[j]))
        j += 1
    dummy_bonding_pads_group = Group(dummy_bonding_pads)
    dummy_bonding_pads_group.distribute(direction = direction, spacing = spacing)
    dummy_bonding_pads_group.move(destination = destination)
    return dummy_bonding_pads

def build_2x3structure(layer, layer_dot, layer_barrier, layer_screening, layer_ohmic):
    D = pg.preview_layerset(lys, size = 100, spacing = 130)

    overlap = 10

    scale_factor = 1
    # inner dots parameters
    radius_dot = 120/2*scale_factor #nm
    #layer_dot = lys['GL2 Pd/Ti']
    angle_resolution_dot = 2.5 # unclear
    spacing_dots = 50*scale_factor
    # sensor dot parameters
    radius_sensor_dot = 180/2*scale_factor #nm

    angle_resolution_sensor_dot = 2.5 # unclear

    tilt_ports = 25*scale_factor

    small_tilt = 15*scale_factor


    # plunger gate inner dots
    width_rect_plunger_inner_dot = 40*scale_factor #nm
    length_rect_plunger_inner_dot = 200*scale_factor #nm minimum length

    # plunger gate sensor dots
    width_rect_plunger_sensor_dot = 40*scale_factor #nm
    length_rect_plunger_sensor_dot = 200*scale_factor #nm


    layer_spacing = spacing_dots + 2*radius_dot


    extra_spacing_dot_6_5 = 70*scale_factor #nm

    D = Device()
    dot1 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])
    dot2 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])
    dot3 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])
    dot4 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])
    dot5 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])
    dot6 = D << pg.circle(radius = radius_dot, angle_resolution = angle_resolution_dot, layer = layer_dot[layer])

    sensordot1 = D << pg.circle(radius = radius_sensor_dot, angle_resolution = angle_resolution_sensor_dot, layer = layer_dot[layer])
    sensordot2 = D << pg.circle(radius = radius_sensor_dot, angle_resolution = angle_resolution_sensor_dot, layer = layer_dot[layer])

    #creating a group of the vertical array of dots
    vertical_array_dots_1 = Group([dot1,dot2,dot3])

    vertical_array_dots_2 = Group([dot4,dot5,dot6])

    #distribute vertical array of dots
    vertical_array_dots_1.distribute(direction = 'x', spacing = spacing_dots)
    vertical_array_dots_2.distribute(direction = 'x', spacing = spacing_dots)
    vertical_array_dots_1.move(destination = (0, 0))
    vertical_array_dots_2.move(destination = (0, layer_spacing))

    sensordot1.move(destination = (dot1.center[0] - radius_sensor_dot - radius_dot - spacing_dots + (width_rect_plunger_inner_dot - 10), dot1.center[1]))
    sensordot2.move(destination = (dot6.center[0] + radius_sensor_dot + radius_dot + spacing_dots - (width_rect_plunger_inner_dot - 10), dot6.center[1]))

    # #Adding rectangular plunger gates of inner dots
    rect_plunger_dot_1 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot), layer = layer_dot[layer])
    rect_plunger_dot_2 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot + 75), layer = layer_dot[layer])
    rect_plunger_dot_3 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot), layer = layer_dot[layer])
    rect_plunger_dot_4 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot), layer = layer_dot[layer])
    rect_plunger_dot_5 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot + 75), layer = layer_dot[layer])
    rect_plunger_dot_6 = D << pg.rectangle(size = (width_rect_plunger_inner_dot, length_rect_plunger_inner_dot), layer = layer_dot[layer])

    rect_plunger_sensordot_1 = D << pg.rectangle(size = (length_rect_plunger_sensor_dot, width_rect_plunger_sensor_dot), layer = layer_dot[layer])
    rect_plunger_sensordot_2 = D << pg.rectangle(size = (length_rect_plunger_sensor_dot, width_rect_plunger_sensor_dot), layer = layer_dot[layer])

    horizontal_plunger_of_array_dots_1 = Group([rect_plunger_dot_1, rect_plunger_dot_3])
    horizontal_plunger_of_array_dots_1.distribute(direction = 'x', spacing = 4*radius_dot + 2*(dot2.xmin - dot1.xmax) - width_rect_plunger_inner_dot)
    horizontal_plunger_of_array_dots_1.move(destination = (dot1.center[0] - 0.5*width_rect_plunger_inner_dot, dot1.center[1] - length_rect_plunger_inner_dot - radius_dot/1.1))

    rect_plunger_sensordot_1.move(destination = (sensordot1.center[0] - length_rect_plunger_sensor_dot + overlap - radius_sensor_dot, sensordot1.center[1] - 0.5*width_rect_plunger_sensor_dot))
    rect_plunger_sensordot_2.move(destination = (sensordot2.center[0] + radius_sensor_dot - overlap, sensordot2.center[1] - 0.5*width_rect_plunger_sensor_dot))

    vertical_array_plunger_dots_1 = Group([rect_plunger_dot_2, rect_plunger_dot_5])
    vertical_array_plunger_dots_1.distribute(direction = 'y', spacing = 4*radius_dot + (dot5.ymin - dot2.ymax)/1.22)
    vertical_array_plunger_dots_1.move(destination = (dot2.center[0] - 0.5*width_rect_plunger_inner_dot, dot2.center[1] - radius_dot/1.1 - (length_rect_plunger_inner_dot+75)))

    horizontal_plunger_of_array_dots_2 = Group([rect_plunger_dot_4, rect_plunger_dot_6])
    horizontal_plunger_of_array_dots_2.distribute(direction = 'x', spacing = 4*radius_dot + 2*(dot5.xmin - dot4.xmax) - width_rect_plunger_inner_dot)
    horizontal_plunger_of_array_dots_2.move(destination = (dot4.center[0] - 0.5*width_rect_plunger_inner_dot, dot4.center[1] + radius_dot/1.1))

    sensordot_1_port = D << pg.polygon_ports(xpts = [rect_plunger_sensordot_1.xmin - tilt_ports,
                                                    rect_plunger_sensordot_1.xmin, 
                                                    rect_plunger_sensordot_1.xmin, 
                                                    rect_plunger_sensordot_1.xmin - tilt_ports], 
                                                    ypts = [rect_plunger_sensordot_1.ymax, 
                                                     rect_plunger_sensordot_1.ymax, 
                                                     rect_plunger_sensordot_1.ymin,
                                                     rect_plunger_sensordot_1.ymax], layer = layer_dot[layer])
    
    sensordot_2_port = D << pg.polygon_ports(xpts = [rect_plunger_sensordot_2.xmax + tilt_ports,
                                                    rect_plunger_sensordot_2.xmax, 
                                                    rect_plunger_sensordot_2.xmax, 
                                                    rect_plunger_sensordot_2.xmax + tilt_ports], 
                                                    ypts = [rect_plunger_sensordot_2.ymin, 
                                                     rect_plunger_sensordot_2.ymin, 
                                                     rect_plunger_sensordot_2.ymax,
                                                     rect_plunger_sensordot_2.ymin], layer = layer_dot[layer])
    
    dot_1_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_1.xmin,
                                                    rect_plunger_dot_1.xmax, 
                                                    rect_plunger_dot_1.xmax, 
                                                    rect_plunger_dot_1.xmin], 
                                                    ypts = [rect_plunger_dot_1.ymin, 
                                                     rect_plunger_dot_1.ymin, 
                                                     rect_plunger_dot_1.ymin -  tilt_ports,
                                                     rect_plunger_dot_1.ymin], layer = layer_dot[layer])
    
    dot_2_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_2.xmin,
                                                    rect_plunger_dot_2.xmin, 
                                                    rect_plunger_dot_2.xmax, 
                                                    rect_plunger_dot_2.xmin], 
                                                    ypts = [rect_plunger_dot_2.ymin, 
                                                     rect_plunger_dot_2.ymin - small_tilt, 
                                                     rect_plunger_dot_2.ymin,
                                                     rect_plunger_dot_2.ymin], layer = layer_dot[layer])
    
    dot_3_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_3.xmin,
                                                    rect_plunger_dot_3.xmin, 
                                                    rect_plunger_dot_3.xmax, 
                                                    rect_plunger_dot_3.xmin], 
                                                    ypts = [rect_plunger_dot_3.ymin, 
                                                     rect_plunger_dot_3.ymin - tilt_ports, 
                                                     rect_plunger_dot_3.ymin,
                                                     rect_plunger_dot_3.ymin], layer = layer_dot[layer])
    
    dot_4_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_4.xmin,
                                                    rect_plunger_dot_4.xmax, 
                                                    rect_plunger_dot_4.xmax, 
                                                    rect_plunger_dot_4.xmin], 
                                                    ypts = [rect_plunger_dot_4.ymax, 
                                                     rect_plunger_dot_4.ymax + tilt_ports, 
                                                     rect_plunger_dot_4.ymax,
                                                     rect_plunger_dot_4.ymax], layer = layer_dot[layer])
    
    dot_5_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_5.xmin,
                                                    rect_plunger_dot_5.xmax, 
                                                    rect_plunger_dot_5.xmax, 
                                                    rect_plunger_dot_5.xmin], 
                                                    ypts = [rect_plunger_dot_5.ymax, 
                                                     rect_plunger_dot_5.ymax + small_tilt, 
                                                     rect_plunger_dot_5.ymax,
                                                     rect_plunger_dot_5.ymax], layer = layer_dot[layer])
    
    dot_6_port = D << pg.polygon_ports(xpts = [rect_plunger_dot_6.xmin,
                                                    rect_plunger_dot_6.xmax, 
                                                    rect_plunger_dot_6.xmin, 
                                                    rect_plunger_dot_6.xmin], 
                                                    ypts = [rect_plunger_dot_6.ymax + tilt_ports, 
                                                     rect_plunger_dot_6.ymax, 
                                                     rect_plunger_dot_6.ymax,
                                                     rect_plunger_dot_6.ymax + tilt_ports], layer = layer_dot[layer])

    spacing_between_barrier_gate_and_lollipop = 30*scale_factor #nm


    # Barrier gates inner dots
    
    width_barrier_rect_plunger_inner_dot = width_rect_plunger_inner_dot*scale_factor #percentage_width_wrt_plunger*width_rect_plunger_inner_dot #nm
    length_barrier_rect_plunger_inner_dot = length_rect_plunger_inner_dot-radius_dot -spacing_between_barrier_gate_and_lollipop #nm
    length_barrier_sensordot = length_barrier_rect_plunger_inner_dot - 30
    width_barrier_sensordot = width_rect_plunger_inner_dot

    length_barrier_cross = length_barrier_rect_plunger_inner_dot - 40*scale_factor
    width_barrier_cross = width_rect_plunger_inner_dot - 5*scale_factor

    radius_arc = width_barrier_cross + 50

    length_barrier_cross_adds = length_barrier_rect_plunger_inner_dot + 100

    barrier_cross_1 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_1.move(destination = (dot1.center[0] + radius_dot + 0.5*(dot2.xmin - dot1.xmax), dot1.center[1]))

    cross_barrier_1_port = D << pg.polygon_ports(xpts = [barrier_cross_1.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                          barrier_cross_1.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_1.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_1.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                                  ypts = [barrier_cross_1.ymin, 
                                                          barrier_cross_1.ymin, 
                                                          barrier_cross_1.ymin - length_barrier_cross_adds - tilt_ports,
                                                          barrier_cross_1.ymin - length_barrier_cross_adds],
                                                  layer = layer_barrier[layer])

    barrier_cross_2 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_2.move(destination = (dot2.center[0] + radius_dot + 0.5*(dot2.xmin - dot1.xmax), dot1.center[1]))

    cross_barrier_2_port = D << pg.polygon_ports(xpts = [barrier_cross_2.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                          barrier_cross_2.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_2.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_2.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                                  ypts = [barrier_cross_2.ymin, 
                                                          barrier_cross_2.ymin, 
                                                          barrier_cross_2.ymin - length_barrier_cross_adds,
                                                          barrier_cross_2.ymin - length_barrier_cross_adds - tilt_ports],
                                                  layer = layer_barrier[layer])

    barrier_cross_3 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_3.move(destination = (dot1.center[0], dot1.center[1] + radius_dot + 0.5*(dot4.ymin - dot1.ymax)))

    barrier_cross_4 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_4.move(destination = (dot2.center[0], dot2.center[1] + radius_dot + 0.5*(dot4.ymin - dot1.ymax)))

    cross_barrier_4_port = D << pg.polygon_ports(xpts = [barrier_cross_4.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                          barrier_cross_4.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_4.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_4.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                                  ypts = [barrier_cross_4.ymax + 2*length_barrier_cross_adds, 
                                                          barrier_cross_4.ymax + 2*length_barrier_cross_adds + tilt_ports, 
                                                          barrier_cross_4.ymax,
                                                          barrier_cross_4.ymax],
                                                  layer = layer_barrier[layer])

    barrier_cross_5 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_5.move(destination = (dot3.center[0], dot3.center[1] + radius_dot + 0.5*(dot4.ymin - dot1.ymax)))

    barrier_cross_6 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_6.move(destination = (dot4.center[0] + radius_dot + 0.5*(dot2.xmin - dot1.xmax), dot4.center[1]))

    cross_barrier_6_port = D << pg.polygon_ports(xpts = [barrier_cross_6.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                          barrier_cross_6.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_6.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_6.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                                  ypts = [barrier_cross_6.ymax + length_barrier_cross_adds, 
                                                          barrier_cross_6.ymax + length_barrier_cross_adds + tilt_ports, 
                                                          barrier_cross_6.ymax,
                                                          barrier_cross_6.ymax],
                                                  layer = layer_barrier[layer])

    barrier_cross_7 = D << pg.cross(length = length_barrier_cross, width = width_barrier_cross, layer = layer_barrier[layer])
    barrier_cross_7.move(destination = (dot5.center[0] + radius_dot + 0.5*(dot2.xmin - dot1.xmax), dot4.center[1]))

    cross_barrier_7_port = D << pg.polygon_ports(xpts = [barrier_cross_7.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                          barrier_cross_7.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_7.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                          barrier_cross_7.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                                  ypts = [barrier_cross_7.ymax + length_barrier_cross_adds + tilt_ports, 
                                                          barrier_cross_7.ymax + length_barrier_cross_adds, 
                                                          barrier_cross_7.ymax,
                                                          barrier_cross_7.ymax],
                                                  layer = layer_barrier[layer])

    cross_3_port = D << pg.polygon_ports(xpts = [  barrier_cross_3.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                    barrier_cross_3.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                    barrier_cross_3.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                    barrier_cross_3.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                            ypts = [barrier_cross_3.ymax, 
                                                    barrier_cross_3.ymax + tilt_ports, 
                                                    barrier_cross_3.ymax,
                                                    barrier_cross_3.ymax], layer = layer_barrier[layer])
    
    cross_5_port = D << pg.polygon_ports(xpts = [  barrier_cross_5.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross,
                                                    barrier_cross_5.xmax - 0.5*length_barrier_cross + 0.5*width_barrier_cross, 
                                                    barrier_cross_5.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross, 
                                                    barrier_cross_5.xmin + 0.5*length_barrier_cross - 0.5*width_barrier_cross], 
                                            ypts = [barrier_cross_5.ymin, 
                                                    barrier_cross_5.ymin, 
                                                    barrier_cross_5.ymin - tilt_ports,
                                                    barrier_cross_5.ymin], layer = layer_barrier[layer])


    barrier_rect_plunger_sensordot_1 = D << pg.rectangle(size = (length_barrier_sensordot + 100, width_barrier_sensordot), layer = layer_barrier[layer])
    barrier_rect_plunger_sensordot_2 = D << pg.rectangle(size = (length_barrier_sensordot, width_barrier_sensordot), layer = layer_barrier[layer])

    barrier_rect_plunger_sensordot_3 = D << pg.rectangle(size = (length_barrier_sensordot, width_barrier_sensordot), layer = layer_barrier[layer])
    barrier_rect_plunger_sensordot_4 = D << pg.rectangle(size = (length_barrier_sensordot, width_barrier_sensordot), layer = layer_barrier[layer])

    barrier_rect_plunger_sensordot_5 = D << pg.rectangle(size = (length_barrier_sensordot, width_barrier_sensordot), layer = layer_barrier[layer])
    barrier_rect_plunger_sensordot_6 = D << pg.rectangle(size = (length_barrier_sensordot + 100, width_barrier_sensordot), layer = layer_barrier[layer])



    barrier_rect_plunger_sensordot_1.move(destination = (sensordot1.center[0] - 0.5*length_barrier_sensordot, sensordot1.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_1.rotate(angle = 270, center = sensordot1.center)
    

    barrier_rect_plunger_sensordot_2.move(destination = (sensordot2.center[0] - 0.5*length_barrier_sensordot, sensordot2.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_2.rotate(angle = 0, center = sensordot2.center)

    barrier_rect_plunger_sensordot_3.move(destination = (sensordot1.center[0] - 0.5*length_barrier_sensordot, sensordot1.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_3.rotate(angle = 0, center = sensordot1.center)

    barrier_rect_plunger_sensordot_4.move(destination = (sensordot2.center[0] - 0.5*length_barrier_sensordot, sensordot2.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_4.rotate(angle = 180, center = sensordot2.center)

    barrier_rect_plunger_sensordot_5.move(destination = (sensordot1.center[0] - 0.5*length_barrier_sensordot, sensordot1.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_5.rotate(angle = 180, center = sensordot1.center)

    barrier_rect_plunger_sensordot_6.move(destination = (sensordot2.center[0] - 0.5*length_barrier_sensordot, sensordot2.center[1] + radius_sensor_dot - 10))
    barrier_rect_plunger_sensordot_6.rotate(angle = 90, center = sensordot2.center)

    barrier_sensordot_1_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_1.xmin,
                                                                barrier_rect_plunger_sensordot_1.xmax, 
                                                                barrier_rect_plunger_sensordot_1.xmax, 
                                                                barrier_rect_plunger_sensordot_1.xmin], 
                                                    ypts = [    barrier_rect_plunger_sensordot_1.ymin, 
                                                                barrier_rect_plunger_sensordot_1.ymin, 
                                                                barrier_rect_plunger_sensordot_1.ymin - tilt_ports,
                                                                barrier_rect_plunger_sensordot_1.ymin],
                                                    layer = layer_barrier[layer])

    barrier_sensordot_2_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_2.xmax + tilt_ports,
                                                                barrier_rect_plunger_sensordot_2.xmax, 
                                                                barrier_rect_plunger_sensordot_2.xmax, 
                                                                barrier_rect_plunger_sensordot_2.xmax + tilt_ports], 
                                                    ypts = [    barrier_rect_plunger_sensordot_2.ymin, 
                                                                barrier_rect_plunger_sensordot_2.ymin, 
                                                                barrier_rect_plunger_sensordot_2.ymax,
                                                                barrier_rect_plunger_sensordot_2.ymin],
                                                    layer = layer_barrier[layer])

    barrier_sensordot_3_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_3.xmin - tilt_ports,
                                                                barrier_rect_plunger_sensordot_3.xmin, 
                                                                barrier_rect_plunger_sensordot_3.xmin, 
                                                                barrier_rect_plunger_sensordot_3.xmin - tilt_ports], 
                                                    ypts = [    barrier_rect_plunger_sensordot_3.ymin, 
                                                                barrier_rect_plunger_sensordot_3.ymin, 
                                                                barrier_rect_plunger_sensordot_3.ymax,
                                                                barrier_rect_plunger_sensordot_3.ymin],
                                                    layer = layer_barrier[layer])
    
    barrier_sensordot_4_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_4.xmax + tilt_ports,
                                                                barrier_rect_plunger_sensordot_4.xmax, 
                                                                barrier_rect_plunger_sensordot_4.xmax, 
                                                                barrier_rect_plunger_sensordot_4.xmax + tilt_ports], 
                                                    ypts = [    barrier_rect_plunger_sensordot_4.ymax, 
                                                                barrier_rect_plunger_sensordot_4.ymax, 
                                                                barrier_rect_plunger_sensordot_4.ymin,
                                                                barrier_rect_plunger_sensordot_4.ymax],
                                                    layer = layer_barrier[layer])
    
    barrier_sensordot_5_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_5.xmin - tilt_ports,
                                                                barrier_rect_plunger_sensordot_5.xmin, 
                                                                barrier_rect_plunger_sensordot_5.xmin, 
                                                                barrier_rect_plunger_sensordot_5.xmin - tilt_ports], 
                                                    ypts = [    barrier_rect_plunger_sensordot_5.ymax, 
                                                                barrier_rect_plunger_sensordot_5.ymax, 
                                                                barrier_rect_plunger_sensordot_5.ymin,
                                                                barrier_rect_plunger_sensordot_5.ymax],
                                                    layer = layer_barrier[layer])
    
    barrier_sensordot_6_port = D << pg.polygon_ports(xpts = [  barrier_rect_plunger_sensordot_6.xmin,
                                                                barrier_rect_plunger_sensordot_6.xmin, 
                                                                barrier_rect_plunger_sensordot_6.xmax, 
                                                                barrier_rect_plunger_sensordot_6.xmin], 
                                                    ypts = [    barrier_rect_plunger_sensordot_6.ymax, 
                                                                barrier_rect_plunger_sensordot_6.ymax + tilt_ports, 
                                                                barrier_rect_plunger_sensordot_6.ymax,
                                                                barrier_rect_plunger_sensordot_6.ymax],
                                                    layer = layer_barrier[layer])



    length_screening_gates = rect_plunger_dot_3.bbox[1, 0] - rect_plunger_dot_1.bbox[0, 0]

    width_screen_rect_plunger_inner_dot = length_screening_gates #percentage_width_wrt_plunger*width_rect_plunger_inner_dot #nm
    length_screen_rect_plunger_inner_dot = length_rect_plunger_inner_dot*0.5 #nm




    # screening gates of array
    screen_rect_plunger_dot_1 = D << pg.rectangle(size = (width_screen_rect_plunger_inner_dot + 20, length_screen_rect_plunger_inner_dot), layer = layer_screening[layer])
    screen_rect_plunger_dot_2 = D << pg.rectangle(size = (width_screen_rect_plunger_inner_dot + 20, length_screen_rect_plunger_inner_dot), layer = layer_screening[layer])

    screen_rect_plunger_gate_1 = D << pg.rectangle(size = (width_barrier_rect_plunger_inner_dot, length_barrier_rect_plunger_inner_dot + 75), layer = layer_screening[layer])
    screen_rect_plunger_gate_2 = D << pg.rectangle(size = (width_barrier_rect_plunger_inner_dot, length_barrier_rect_plunger_inner_dot + 75), layer = layer_screening[layer])

    ramp_sensordot_screening_2 = D << pg.ramp(length = length_barrier_sensordot, width1 = width_barrier_cross, width2 = 2*width_barrier_cross, layer = layer_screening[layer])
    ramp_sensordot_screening_2.center = (0,0)
    ramp_sensordot_screening_2.rotate(-90)
    ramp_sensordot_screening_2.move(destination = (sensordot2.center[0] + radius_sensor_dot + 3*overlap, sensordot2.center[1]))

    ramp_sensordot_screening_1 = D << pg.ramp(length = length_barrier_sensordot, width1 = width_barrier_cross, width2 = 2*width_barrier_cross, layer = layer_screening[layer])
    ramp_sensordot_screening_1.center = (0,0)
    ramp_sensordot_screening_1.rotate(90)
    ramp_sensordot_screening_1.xmax = sensordot1.center[0] - radius_sensor_dot - 3*overlap

    screening_sensordot1_connector = D << pg.polygon_ports(xpts = [ ramp_sensordot_screening_1.xmin,
                                                                    ramp_sensordot_screening_1.xmax ,
                                                                    ramp_sensordot_screening_1.xmax,
                                                                    ramp_sensordot_screening_1.xmin],
                                                            ypts = [ramp_sensordot_screening_1.ymax,
                                                                    ramp_sensordot_screening_1.ymax + tilt_ports,
                                                                    ramp_sensordot_screening_1.ymax,
                                                                    ramp_sensordot_screening_1.ymax],
                                                            layer = layer_screening[layer])

    screening_sensordot2_connector = D << pg.polygon_ports(xpts = [ramp_sensordot_screening_2.xmin,
                                                                    ramp_sensordot_screening_2.xmin ,
                                                                    ramp_sensordot_screening_2.xmax,
                                                                    ramp_sensordot_screening_2.xmin],
                                                            ypts = [ramp_sensordot_screening_2.ymin,
                                                                    ramp_sensordot_screening_2.ymin - tilt_ports,
                                                                    ramp_sensordot_screening_2.ymin,
                                                                    ramp_sensordot_screening_2.ymin],
                                                            layer = layer_screening[layer])


    horizontal_array_screen_gates_1 = Group([screen_rect_plunger_dot_1,screen_rect_plunger_dot_2])
    horizontal_array_screen_gates_1.distribute(direction = 'y', spacing = 2.7*radius_dot + layer_spacing)
    horizontal_array_screen_gates_1.move(destination = (dot1.center[0] - 0.5*width_rect_plunger_inner_dot - 10, dot1.center[1] - 3*radius_dot))


    screen_rect_plunger_gate_1.move(destination = (dot5.center[0] + 0.5*width_rect_plunger_inner_dot/0.5, dot5.center[1] + 3*radius_dot))

    screen_rect_plunger_gate_2.move(destination = (dot2.center[0] - 2*width_barrier_rect_plunger_inner_dot, dot2.center[1] - 3*radius_dot - (length_barrier_rect_plunger_inner_dot + 75)))

    screen_1_port = D << pg.polygon_ports(xpts = [             screen_rect_plunger_gate_2.xmin,
                                                                screen_rect_plunger_gate_2.xmax, 
                                                                screen_rect_plunger_gate_2.xmax, 
                                                                screen_rect_plunger_gate_2.xmin], 
                                                    ypts = [    screen_rect_plunger_gate_2.ymin, 
                                                                screen_rect_plunger_gate_2.ymin, 
                                                                screen_rect_plunger_gate_2.ymin - tilt_ports,
                                                                screen_rect_plunger_gate_2.ymin],
                                                    layer = layer_screening[layer])
    
    screen_2_port = D << pg.polygon_ports(xpts = [             screen_rect_plunger_gate_1.xmin,
                                                                screen_rect_plunger_gate_1.xmin, 
                                                                screen_rect_plunger_gate_1.xmax, 
                                                                screen_rect_plunger_gate_1.xmin], 
                                                    ypts = [    screen_rect_plunger_gate_1.ymax, 
                                                                screen_rect_plunger_gate_1.ymax + tilt_ports, 
                                                                screen_rect_plunger_gate_1.ymax,
                                                                screen_rect_plunger_gate_1.ymax],
                                                    layer = layer_screening[layer])

    #qp(D)


    width_ohmics = width_barrier_rect_plunger_inner_dot
    length_ohmics = length_barrier_rect_plunger_inner_dot - 150

    ohmic_1 = D << pg.rectangle(size = (width_ohmics, length_ohmics), layer = layer_ohmic[layer])
    ohmic_1.move(destination = (sensordot1.center[0] - 0.5*width_barrier_rect_plunger_inner_dot, sensordot1.center[1] + radius_sensor_dot + width_barrier_sensordot + width_ohmics - 20))
    ohmic_2 = D << pg.rectangle(size = (width_ohmics, length_ohmics), layer = layer_ohmic[layer])
    ohmic_2.move(destination = (sensordot1.center[0] - 0.5*width_barrier_rect_plunger_inner_dot, sensordot1.center[1] + radius_sensor_dot + width_barrier_sensordot + width_ohmics - 20)).rotate(angle = -180, center = sensordot1.center)

    ohmic_3 = D << pg.rectangle(size = (width_ohmics, length_ohmics), layer = layer_ohmic[layer])
    ohmic_3.move(destination = (sensordot2.center[0] - 0.5*width_barrier_rect_plunger_inner_dot, sensordot2.center[1] + radius_sensor_dot + width_barrier_sensordot + width_ohmics - 20))
    ohmic_4 = D << pg.rectangle(size = (width_ohmics, length_ohmics), layer = layer_ohmic[layer])
    ohmic_4.move(destination = (sensordot2.center[0] - 0.5*width_barrier_rect_plunger_inner_dot, sensordot2.center[1] + radius_sensor_dot + width_barrier_sensordot + width_ohmics - 20)).rotate(angle = -180, center = sensordot2.center)

    from phidl import set_quickplot_options
    set_quickplot_options(show_ports=None, show_subports=False,  label_aliases=None)

    # Ports for sensor dots

    length_port_rect_plunger_sensordot_1 = width_rect_plunger_sensor_dot
    orientation_port_rect_plunger_sensordot_1 = 0
    port_sensordot_1 = D.add_port(name='SD1', midpoint = [rect_plunger_sensordot_1.xmin, rect_plunger_sensordot_1.ymax - 0.5*width_rect_plunger_sensor_dot],
                                                width = length_port_rect_plunger_sensordot_1, orientation = orientation_port_rect_plunger_sensordot_1)


    length_port_rect_plunger_sensordot_2 = width_rect_plunger_sensor_dot
    orientation_port_rect_plunger_sensordot_2 = 0
    port_sensordot_2 = D.add_port(name='SD2', midpoint = [rect_plunger_sensordot_2.xmax, rect_plunger_sensordot_2.ymax - 0.5*width_rect_plunger_sensor_dot],
                                                width = length_port_rect_plunger_sensordot_2, orientation = orientation_port_rect_plunger_sensordot_2)

    # Ading Ohmic ports

    length_port_ohmic_1 = np.sqrt(width_ohmics**2 + length_ohmics**2)
    orientation_port_ohmic_1 = -45
    port_ohmic_1 = D.add_port(name='O1', midpoint = [ohmic_1.center[0], ohmic_1.center[1]],
                                                width = length_port_ohmic_1, orientation = orientation_port_ohmic_1)

    length_port_ohmic_2 = np.sqrt(width_ohmics**2 + length_ohmics**2)
    orientation_port_ohmic_2 = -135
    port_ohmic_2 = D.add_port(name='O2', midpoint = [ohmic_2.center[0], ohmic_2.center[1]],
                                                width = length_port_ohmic_2, orientation = orientation_port_ohmic_2)

    length_port_ohmic_3 = np.sqrt(width_ohmics**2 + length_ohmics**2)
    orientation_port_ohmic_3 = 45
    port_ohmic_3 = D.add_port(name='O3', midpoint = [ohmic_3.center[0], ohmic_3.center[1]],
                                                width = length_port_ohmic_3, orientation = orientation_port_ohmic_3)

    length_port_ohmic_4 = np.sqrt(width_ohmics**2 + length_ohmics**2)
    orientation_port_ohmic_4 = -45
    port_ohmic_4 = D.add_port(name='O4', midpoint = [ohmic_4.center[0], ohmic_4.center[1]],
                                                width = length_port_ohmic_4, orientation = orientation_port_ohmic_4)

    pu.write_lyp('2x3_inner_structure_props.lyp', layerset = lys)

    D.write_gds(filename = '2x3_inner_structure.gds', # Output GDS file name
            unit = 1e-9,                  # Base unit (1e-6 = microns)
            precision = 1e-9,             # Precision / resolution (1e-9 = nanometers)
            auto_rename = True,           # Automatically rename cells to avoid collisions
            max_cellname_length = 28,     # Max length of cell names
            cellname = 'toplevel'         # Name of output top-level cell
           )



    North_placements = [1, 3, 4, 5, 6, 7, 8, 9, 10]

    South_placements = [3, 4, 5, 6, 7, 8, 10, 11]

    West_placements = [4, 5, 6, 7, 8, 9]

    East_placements = [4, 5, 6, 7, 8, 9]

    N_bond_numbers = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    S_bond_numbers = [3, 4, 5, 6, 7, 8, 9, 10]
    W_bond_numbers = [4, 5, 6, 7, 8, 9]
    E_bond_numbers = [4, 5, 6, 7, 8, 9]

    bond_numbers = [N_bond_numbers, S_bond_numbers, W_bond_numbers, E_bond_numbers]


    list_of_numbers = [North_placements, South_placements, West_placements, East_placements]

    ext_port_N, ext_port_S, ext_port_W, ext_port_E = [], [], [], []

    counter_1 = 0

    for i in list_of_numbers:
        for j in range(len(i)):
            if counter_1 == 0:
                ext_port_N.append('EXT_N' + str(j))
                ext_port_N[j] = D.add_port(name='EXT_N' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]),ymax_ext_box], width = width_external_port, orientation = 270)
            if counter_1 == 1:
                ext_port_S.append('EXT_S' + str(j))
                ext_port_S[j] = D.add_port(name='EXT_S' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]),ymin_ext_box], width = width_external_port, orientation = 90)
            if counter_1 == 2:
                ext_port_W.append('EXT_W' + str(j))
                ext_port_W[j] = D.add_port(name='EXT_W' + str(j), midpoint = [xmin_ext_box, y_midpoint_ext_port_v(i[j])], width = width_external_port, orientation = 0)
            if counter_1 == 3:
                ext_port_E.append('EXT_E' + str(j))
                ext_port_E[j] = D.add_port(name='EXT_E' + str(j), midpoint = [xmax_ext_box, y_midpoint_ext_port_v(i[j])], width = width_external_port, orientation = 180)
        counter_1 += 1

    layer_N = [layer_barrier[layer], layer_dot[layer], layer_barrier[layer], layer_dot[layer], layer_barrier[layer], layer_screening[layer], layer_barrier[layer], layer_dot[layer], layer_barrier[layer]]
    layer_S = [layer_barrier[layer], layer_dot[layer], layer_barrier[layer], layer_screening[layer], layer_dot[layer], layer_barrier[layer], layer_dot[layer], layer_barrier[layer]]
    layer_W = [layer_ohmic[layer], layer_barrier[layer], layer_dot[layer], layer_screening[layer], layer_barrier[layer], layer_ohmic[layer]]
    layer_E = [layer_ohmic[layer], layer_barrier[layer], layer_screening[layer], layer_dot[layer], layer_barrier[layer], layer_ohmic[layer]]

    North_objects_to_route = [cross_3_port.ports['1'], dot_4_port.ports['1'], cross_barrier_6_port.ports['1'],
                            dot_5_port.ports['1'], cross_barrier_4_port.ports['1'], screen_2_port.ports['2'],
                            cross_barrier_7_port.ports['1'], dot_6_port.ports['1'], barrier_sensordot_6_port.ports['2']]

    West_objects_to_route = [port_ohmic_2, barrier_sensordot_5_port.ports['3'], sensordot_1_port.ports['3']
                              , screening_sensordot1_connector.ports['1'], barrier_sensordot_3_port.ports['3'], port_ohmic_1]
    
    South_objects_to_route = [barrier_sensordot_1_port.ports['3'], dot_1_port.ports['3'], cross_barrier_1_port.ports['3'],
                             screen_1_port.ports['3'], dot_2_port.ports['2'], cross_barrier_2_port.ports['3'],
                             dot_3_port.ports['2'], cross_5_port.ports['2']]
    
    East_objects_to_route = [port_ohmic_4, barrier_sensordot_4_port.ports['3'], 
                            screening_sensordot2_connector.ports['2'], sensordot_2_port.ports['3'], barrier_sensordot_2_port.ports['3'],
                            port_ohmic_3]


    ext_ports = [ext_port_N, ext_port_S, ext_port_W, ext_port_E]

    route_ext_N = []
    route_ext_S = []
    route_ext_W = []
    route_ext_E = []


    for i in ext_ports:
        for j in range(len(i)):
            if i == ext_port_N:
                route_ext_N.append('route_ext_N' + str(j))
                route_ext_N[j] = D.add_ref(pr.route_quad(North_objects_to_route[j], i[j], layer = layer_N[j]))
            if i == ext_port_S:
                route_ext_S.append('route_ext_S' + str(j))
                route_ext_S[j] = D.add_ref(pr.route_quad(South_objects_to_route[j], i[j], layer = layer_S[j]))
            if i == ext_port_W:
                route_ext_W.append('route_ext_W' + str(j))
                route_ext_W[j] = D.add_ref(pr.route_quad(West_objects_to_route[j], i[j], layer = layer_W[j]))
            if i == ext_port_E:
                route_ext_E.append('route_ext_E' + str(j))
                route_ext_E[j] = D.add_ref(pr.route_quad(East_objects_to_route[j], i[j], layer = layer_E[j]))

    D_without_local_markers = pg.union(D, by_layer = True)

    D_local_markers = Device()

    d_ref6 = D_local_markers.add_ref(D)


    inner_marker_NE = D_local_markers << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_NE.move(destination = (np.cos(45)*length_to_local, np.cos(45)*length_to_local))

    inner_marker_NW = D_local_markers << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_NW.move(destination = (-np.cos(45)*length_to_local, np.cos(45)*length_to_local))

    inner_marker_SE = D_local_markers << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_SE.move(destination = (np.cos(45)*length_to_local, -np.cos(45)*length_to_local))

    inner_marker_SW = D_local_markers << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_SW.move(destination = (-np.cos(45)*length_to_local, -np.cos(45)*length_to_local))

    D = pg.union(D_local_markers, by_layer = True)

    D_with_local_markers = pg.union(D_local_markers, by_layer = True)

    layer_N_pads = [layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_screening_gates_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad]
    layer_S_pads = [layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad]
    layer_W_pads = [layer_ohmic_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_screening_gates_bond_pad, layer_barrier_gates_bond_pad, layer_ohmic_bond_pad]
    layer_E_pads = [layer_ohmic_bond_pad, layer_barrier_gates_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_ohmic_bond_pad]

    layer_pads = [layer_N_pads, layer_S_pads, layer_W_pads, layer_E_pads]

    return D_with_local_markers, D_without_local_markers, D, list_of_numbers, ext_ports, bond_numbers, layer_pads


def build_1x4resonator(layer, layer_dot, layer_barrier, layer_screening, layer_ohmic, scaling):
    D1 = Device()

    # Define constants 

    angle_resolution = 1

    semi_major_resonator = 500/2*scaling
    semi_minor_resonator = 60/2*scaling

    dot_radius = 120/2*scaling

    sensordot_radius = 1.5*dot_radius

    width_connector = dot_radius*0.5

    length_connector = dot_radius*2

    overlap = 4*scaling

    dot_spacing = width_connector - 2*overlap

    large_tilt = 20*scaling

    small_tilt = 10*scaling


    # build the dots and resonator

    resonator = D1 << pg.ellipse(radii = (semi_major_resonator, semi_minor_resonator), angle_resolution = angle_resolution, layer = layer_dot[layer])

    dot_1 = D1 << pg.circle(radius = dot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    dot_2 = D1 <<  pg.circle(radius = dot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    dot_3 = D1 << pg.circle(radius = dot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    dot_4 = D1 <<  pg.circle(radius = dot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    sensordot_1 = D1 << pg.circle(radius = sensordot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    sensordot_2 = D1 <<  pg.circle(radius = sensordot_radius, angle_resolution = angle_resolution, layer = layer_dot[layer])

    # Group and move the dots and resonator

    dot_array_1 = Group([dot_1, dot_2])
    dot_array_1.distribute(direction = 'x', spacing = dot_spacing)
    dot_array_1.move(destination = ((-semi_major_resonator - 2*dot_spacing - 3*dot_radius), 0))

    dot_array_2 = Group([dot_3, dot_4])
    dot_array_2.distribute(direction = 'x', spacing = dot_spacing)
    dot_array_2.move(destination = ((semi_major_resonator + dot_spacing + dot_radius), 0))

    sensordot_array_1 = Group([sensordot_1, sensordot_2])
    sensordot_array_1.distribute(direction = 'x', spacing = 2*semi_major_resonator + 6*dot_spacing + 8*dot_radius)
    sensordot_array_1.move(destination = ((dot_1.center[0] - dot_radius - dot_spacing - sensordot_radius), 0))

    # Build the connectors to the dots and resonator

    resonator_connector = D1 << pg.compass(size = (width_connector, length_connector), layer = layer_dot[layer])

    dot_1_connector_1 = D1 << pg.rectangle(size = (width_connector, length_connector), layer = layer_dot[layer])

    dot_2_connector_1 = D1 << pg.rectangle(size = (width_connector, length_connector), layer = layer_dot[layer])

    dot_3_connector_1 = D1 << pg.rectangle(size = (width_connector, length_connector), layer = layer_dot[layer])

    dot_4_connector_1 = D1 << pg.rectangle(size = (width_connector, length_connector), layer = layer_dot[layer])

    sensordot_1_connector = D1 << pg.compass(size = (length_connector, width_connector), layer = layer_dot[layer])

    sensordot_2_connector = D1 << pg.compass(size = (length_connector, width_connector), layer = layer_dot[layer])

    # Group and move the connectors

    resonator_connector.move(destination = (0, -semi_minor_resonator + overlap - 0.5*length_connector))

    dot_connector_array_1 = Group([dot_1_connector_1, dot_2_connector_1])
    dot_connector_array_1.distribute(direction = 'x', spacing = dot_spacing + 2*dot_radius - width_connector)
    dot_connector_array_1.move(destination = (dot_1.center[0] - 0.5*width_connector, dot_1.center[1] + dot_radius - overlap))

    dot_connector_array_2 = Group([dot_3_connector_1, dot_4_connector_1])
    dot_connector_array_2.distribute(direction = 'x', spacing = dot_spacing + 2*dot_radius - width_connector)
    dot_connector_array_2.move(destination = (dot_3.center[0] - 0.5*width_connector, dot_3.center[1] + dot_radius - overlap)).rotate(angle = 180, center = (dot_3.center[0] + dot_radius + 0.5*dot_spacing, 0))

    sensordot_connector_array_1 = Group([sensordot_1_connector, sensordot_2_connector])
    sensordot_connector_array_1.distribute(direction = 'x', spacing = 2*semi_major_resonator + 6*dot_spacing + 8*dot_radius + 4*sensordot_radius - 2*overlap)
    sensordot_connector_array_1.move(destination = (sensordot_1.center[0] - 0.5*length_connector - sensordot_radius + overlap,
                                                    sensordot_1.center[1]))
    
    dot_1_connector = D1 << pg.polygon_ports(xpts = [dot_1_connector_1.xmin,
                                                     dot_1_connector_1.xmax,
                                                     dot_1_connector_1.xmax,
                                                     dot_1_connector_1.xmin],
                                             ypts = [dot_1_connector_1.ymax,
                                                     dot_1_connector_1.ymax + large_tilt,
                                                     dot_1_connector_1.ymax,
                                                     dot_1_connector_1.ymax],
                                                            layer = layer_dot[layer])
    
    dot_2_connector = D1 << pg.polygon_ports(xpts = [dot_2_connector_1.xmin,
                                                     dot_2_connector_1.xmax,
                                                     dot_2_connector_1.xmax,
                                                     dot_2_connector_1.xmin],
                                             ypts = [dot_2_connector_1.ymax,
                                                     dot_2_connector_1.ymax + small_tilt,
                                                     dot_2_connector_1.ymax,
                                                     dot_2_connector_1.ymax],
                                                            layer = layer_dot[layer])
    
    dot_3_connector = D1 << pg.polygon_ports(xpts = [dot_4_connector_1.xmin,
                                                     dot_4_connector_1.xmin,
                                                     dot_4_connector_1.xmax,
                                                     dot_4_connector_1.xmin],
                                             ypts = [dot_4_connector_1.ymin,
                                                     dot_4_connector_1.ymin - small_tilt,
                                                     dot_4_connector_1.ymin,
                                                     dot_4_connector_1.ymin],
                                                            layer = layer_dot[layer])
    
    dot_4_connector = D1 << pg.polygon_ports(xpts = [dot_3_connector_1.xmin,
                                                     dot_3_connector_1.xmin,
                                                     dot_3_connector_1.xmax,
                                                     dot_3_connector_1.xmin],
                                             ypts = [dot_3_connector_1.ymin,
                                                     dot_3_connector_1.ymin - large_tilt,
                                                     dot_3_connector_1.ymin,
                                                     dot_3_connector_1.ymin],
                                                            layer = layer_dot[layer])
    

    # Start building the barriers

    barrier_length = length_connector

    barrier_width = width_connector

    # Long barriers for the sensor dots

    vertical_sensordot_barrier = 2.8*barrier_length

    vertical_sensorbarrier_1 = D1 << pg.rectangle(size = (barrier_width, vertical_sensordot_barrier), layer = layer_barrier[layer])

    vertical_sensorbarrier_1.move(destination = (dot_1.center[0] - dot_radius - dot_spacing + 0.5*dot_spacing - 0.5*barrier_width , -0.5*vertical_sensordot_barrier))

    vertical_sensorbarrier_2 = D1 << pg.rectangle(size = (barrier_width, 3*barrier_length), layer = layer_barrier[layer])
 
    vertical_sensorbarrier_2.move(destination = (dot_4.center[0] + dot_radius + dot_spacing - 0.5*dot_spacing - 0.5*barrier_width , -0.5*vertical_sensordot_barrier))

    vertical_barrier_sensordot_1 = D1 << pg.polygon_ports(xpts = [vertical_sensorbarrier_1.xmin,
                                                          vertical_sensorbarrier_1.xmin,
                                                          vertical_sensorbarrier_1.xmax,
                                                          vertical_sensorbarrier_1.xmax],
                                                  ypts = [vertical_sensorbarrier_1.ymin,
                                                          vertical_sensorbarrier_1.ymax,
                                                          vertical_sensorbarrier_1.ymax + large_tilt,
                                                          vertical_sensorbarrier_1.ymin - large_tilt],
                                                            layer = layer_barrier[layer])
    
    vertical_barrier_sensordot_2 = D1 << pg.polygon_ports(xpts = [vertical_sensorbarrier_2.xmin,
                                                          vertical_sensorbarrier_2.xmin,
                                                          vertical_sensorbarrier_2.xmax,
                                                          vertical_sensorbarrier_2.xmax],
                                                  ypts = [vertical_sensorbarrier_2.ymin - large_tilt,
                                                          vertical_sensorbarrier_2.ymax + large_tilt,
                                                          vertical_sensorbarrier_2.ymax,
                                                          vertical_sensorbarrier_2.ymin],
                                                            layer = layer_barrier[layer])


    # Barriers between the dots and resonator

    vertical_barrier_1_1 = D1 << pg.rectangle(size = (barrier_width, 1.5*barrier_length), layer = layer_barrier[layer])
    vertical_barrier_1_1.move(destination = (-0.5*barrier_width, 0)).rotate(angle = 180).move(destination = (dot_2.center[0] - dot_radius - 0.5*dot_spacing, 0.5*barrier_length))

    vertical_barrier_1 = D1 << pg.polygon_ports(xpts = [vertical_barrier_1_1.xmin,
                                                                  vertical_barrier_1_1.xmax,
                                                                  vertical_barrier_1_1.xmax,
                                                                  vertical_barrier_1_1.xmin],
                                                          ypts = [vertical_barrier_1_1.ymin,
                                                                  vertical_barrier_1_1.ymin,
                                                                  vertical_barrier_1_1.ymin - large_tilt,
                                                                  vertical_barrier_1_1.ymin],
                                                            layer = layer_barrier[layer])

    vertical_barrier_2_1 = D1 << pg.rectangle(size = (barrier_width, 2*barrier_length), layer = layer_barrier[layer])
    vertical_barrier_2_1.move(destination = (-0.5*barrier_width, 0)).rotate(angle = 180).move(destination = (dot_2.center[0] + dot_radius + 0.5*dot_spacing, 0.5*barrier_length))

    vertical_barrier_2 = D1 << pg.polygon_ports(xpts = [vertical_barrier_2_1.xmin,
                                                                  vertical_barrier_2_1.xmax,
                                                                  vertical_barrier_2_1.xmax,
                                                                  vertical_barrier_2_1.xmin],
                                                          ypts = [vertical_barrier_2_1.ymin,
                                                                  vertical_barrier_2_1.ymin,
                                                                  vertical_barrier_2_1.ymin - small_tilt,
                                                                  vertical_barrier_2_1.ymin],
                                                            layer = layer_barrier[layer])


    vertical_barrier_3_1 = D1 << pg.rectangle(size = (barrier_width, 2*barrier_length), layer = layer_barrier[layer])
    vertical_barrier_3_1.move(destination = (-0.5*barrier_width, 0)).rotate(angle = 0).move(destination = (dot_3.center[0] - dot_radius - 0.5*dot_spacing, -0.5*barrier_length))

    vertical_barrier_3 = D1 << pg.polygon_ports(xpts = [vertical_barrier_3_1.xmin,
                                                                  vertical_barrier_3_1.xmin,
                                                                  vertical_barrier_3_1.xmax,
                                                                  vertical_barrier_3_1.xmin],
                                                          ypts = [vertical_barrier_3_1.ymax,
                                                                  vertical_barrier_3_1.ymax + small_tilt,
                                                                  vertical_barrier_3_1.ymax,
                                                                  vertical_barrier_3_1.ymax],
                                                            layer = layer_barrier[layer])

    vertical_barrier_4_1 = D1 << pg.rectangle(size = (barrier_width, 1.5*barrier_length), layer = layer_barrier[layer])
    vertical_barrier_4_1.move(destination = (-0.5*barrier_width, 0)).rotate(angle = 0).move(destination = (dot_3.center[0] + dot_radius + 0.5*dot_spacing, -0.5*barrier_length))

    vertical_barrier_4 = D1 << pg.polygon_ports(xpts = [vertical_barrier_4_1.xmin,
                                                                  vertical_barrier_4_1.xmin,
                                                                  vertical_barrier_4_1.xmax,
                                                                  vertical_barrier_4_1.xmin],
                                                          ypts = [vertical_barrier_4_1.ymax,
                                                                  vertical_barrier_4_1.ymax + large_tilt,
                                                                  vertical_barrier_4_1.ymax,
                                                                  vertical_barrier_4_1.ymax],
                                                            layer = layer_barrier[layer])

    # Taper barrier for the resonator

    barrier_sensor_dot_N_S_length = 0.7*barrier_length

    barrier_resonator_connector = D1 << pg.polygon_ports(xpts = [resonator.center[0] - barrier_length,
                                                                 resonator.center[0] - barrier_length,
                                                                 resonator.center[0] + barrier_length,
                                                                 resonator.center[0] + barrier_length + barrier_width],
                                                         ypts = [resonator.center[0] + semi_minor_resonator - 2*overlap,
                                                                 resonator.center[0] + semi_minor_resonator - 2*overlap + barrier_width,
                                                                 resonator.center[0] + semi_minor_resonator - 2*overlap + barrier_width,
                                                                 resonator.center[0] + semi_minor_resonator - 2*overlap],
                                                         layer = layer_barrier[layer])
    


    taper_barrier_sensordot_1 = D1 << pg.taper(length = barrier_sensor_dot_N_S_length, width1 = barrier_width, width2 = barrier_width, layer = layer_barrier[layer])
    taper_barrier_sensordot_1.move(destination = (sensordot_1.center[0] - 0.5*barrier_sensor_dot_N_S_length, sensordot_1.center[1] + sensordot_radius + 0.5*barrier_width - overlap))
    D1_joined_by_layer = pg.union(D1, by_layer = True)


    taper_barrier_sensordot_2 = D1 << pg.taper(length = barrier_sensor_dot_N_S_length, width1 = barrier_width, width2 = barrier_width, layer = layer_barrier[layer])
    taper_barrier_sensordot_2.move(destination = (sensordot_1.center[0] - 0.5*barrier_sensor_dot_N_S_length, sensordot_1.center[1] - sensordot_radius - 0.5*barrier_width + overlap))
    D1_joined_by_layer = pg.union(D1, by_layer = True)

    taper_barrier_sensordot_3 = D1 << pg.taper(length = barrier_sensor_dot_N_S_length, width1 = barrier_width, width2 = barrier_width, layer = layer_barrier[layer])
    taper_barrier_sensordot_3.move(destination = (sensordot_2.center[0] - 0.5*barrier_sensor_dot_N_S_length, sensordot_2.center[1] + sensordot_radius + 0.5*barrier_width - overlap))
    D1_joined_by_layer = pg.union(D1, by_layer = True)

    taper_barrier_sensordot_4 = D1 << pg.taper(length = barrier_sensor_dot_N_S_length, width1 = barrier_width, width2 = barrier_width, layer = layer_barrier[layer])
    taper_barrier_sensordot_4.move(destination = (sensordot_2.center[0] - 0.5*barrier_sensor_dot_N_S_length, sensordot_2.center[1] - sensordot_radius - 0.5*barrier_width + overlap))

    # Build the screening objects

    ramp_screening_resonator = D1 << pg.ramp(length = 0.5*semi_major_resonator, width1 = 2.6*semi_minor_resonator, width2 = 1.4*semi_minor_resonator, layer = layer_screening[layer])
    ramp_screening_resonator.move(destination = (-0.25*semi_major_resonator, -ramp_screening_resonator.center[1])).rotate(180)
    ramp_screening_resonator.ymax = -semi_minor_resonator - 4*overlap

    ramp_sensordot_screening_2 = D1 << pg.ramp(length = 0.5*barrier_length, width1 = 2*barrier_width, width2 = barrier_width, layer = layer_screening[layer])
    ramp_sensordot_screening_2.move(destination = (-0.2*barrier_length, -2*barrier_width)).rotate(angle = -90)
    ramp_sensordot_screening_2.xmin = sensordot_2_connector.xmin + 3*overlap

    ramp_sensordot_screening_1 = D1 << pg.ramp(length = 0.5*barrier_length, width1 = 2*barrier_width, width2 = barrier_width, layer = layer_screening[layer])
    ramp_sensordot_screening_1.move(destination = (-0.2*barrier_length, -2*barrier_width)).rotate(angle = 90)
    ramp_sensordot_screening_1.xmax = sensordot_1_connector.xmax - 3*overlap

    dot_screening_length = dot_2.xmax - dot_1.xmin - 20*overlap
    dot_screening_width = dot_radius

    dot_screening_1 = D1 << pg.rectangle(size = (dot_screening_length, dot_screening_width), layer = layer_screening[layer])
    dot_screening_1.move(destination = (dot_1.xmin + 10*overlap, dot_1.ymax + 0.75*dot_radius))

    dot_screening_2 = D1 << pg.rectangle(size = (dot_screening_length, dot_screening_width), layer = layer_screening[layer])
    dot_screening_2.move(destination = (dot_3.xmin + 10*overlap, dot_3.ymin - 0.75*dot_radius - dot_screening_width))

    # Build connectors for the screenings

    screening_connector_length = length_connector*0.5

    screening_dotconnector_1 = D1 << pg.compass(size = (width_connector, screening_connector_length), layer = layer_screening[layer])
    screening_dotconnector_2 = D1 << pg.compass(size = (width_connector, screening_connector_length), layer = layer_screening[layer])

    screening_dotconnector_1.move(destination = (vertical_barrier_1.xmin + 0.5*barrier_width, dot_screening_1.ymax + 0.5*screening_connector_length))
    screening_dotconnector_2.move(destination = (vertical_barrier_4.xmin + 0.5*barrier_width, dot_screening_2.ymin - 0.5*screening_connector_length))

    screening_dots_1 = D1 << pg.polygon_ports(xpts = [screening_dotconnector_1.xmin,
                                                                screening_dotconnector_1.xmax ,
                                                                screening_dotconnector_1.xmax,
                                                                screening_dotconnector_1.xmin],
                                                            ypts = [screening_dotconnector_1.ymax,
                                                                screening_dotconnector_1.ymax + large_tilt,
                                                                screening_dotconnector_1.ymax,
                                                                screening_dotconnector_1.ymax],
                                                            layer = layer_screening[layer])
    
    screening_dots_2 = D1 << pg.polygon_ports(xpts = [screening_dotconnector_2.xmin,
                                                                screening_dotconnector_2.xmax ,
                                                                screening_dotconnector_2.xmin,
                                                                screening_dotconnector_2.xmin],
                                                            ypts = [screening_dotconnector_2.ymin,
                                                                screening_dotconnector_2.ymin,
                                                                screening_dotconnector_2.ymin - large_tilt,
                                                                screening_dotconnector_2.ymin],
                                                            layer = layer_screening[layer])

    screening_resonator_connector = D1 << pg.polygon_ports(xpts = [ramp_screening_resonator.xmax,
                                                                ramp_screening_resonator.xmax ,
                                                                ramp_screening_resonator.xmax + 4*large_tilt,
                                                                ramp_screening_resonator.xmax + 4*large_tilt],
                                                            ypts = [ramp_screening_resonator.ymin, ramp_screening_resonator.ymax,
                                                            ramp_screening_resonator.ymax, ramp_screening_resonator.ymin],
                                                            layer = layer_screening[layer])
    
    screening_sensordot1_connector = D1 << pg.polygon_ports(xpts = [ramp_sensordot_screening_1.xmin,
                                                                    ramp_sensordot_screening_1.xmax ,
                                                                    ramp_sensordot_screening_1.xmax,
                                                                    ramp_sensordot_screening_1.xmin],
                                                            ypts = [ramp_sensordot_screening_1.ymin,
                                                                    ramp_sensordot_screening_1.ymin,
                                                                    ramp_sensordot_screening_1.ymin - 1.5*large_tilt,
                                                                    ramp_sensordot_screening_1.ymin],
                                                            layer = layer_screening[layer])

    screening_sensordot2_connector = D1 << pg.polygon_ports(xpts = [ramp_sensordot_screening_2.xmin,
                                                                    ramp_sensordot_screening_2.xmin ,
                                                                    ramp_sensordot_screening_2.xmax,
                                                                    ramp_sensordot_screening_2.xmin],
                                                            ypts = [ramp_sensordot_screening_2.ymax,
                                                                    ramp_sensordot_screening_2.ymax + 1.5*large_tilt,
                                                                    ramp_sensordot_screening_2.ymax,
                                                                    ramp_sensordot_screening_2.ymax],
                                                            layer = layer_screening[layer])
    
    # Build connectors for the barriers

    sensordot_1_barrier_connector_1 = D1 << pg.polygon_ports(xpts = [taper_barrier_sensordot_1.xmin - large_tilt,
                                                                    taper_barrier_sensordot_1.xmin,
                                                                    taper_barrier_sensordot_1.xmin,
                                                                    taper_barrier_sensordot_1.xmin - large_tilt],
                                                            ypts = [taper_barrier_sensordot_1.ymin,
                                                                    taper_barrier_sensordot_1.ymax,
                                                                    taper_barrier_sensordot_1.ymin,
                                                                    taper_barrier_sensordot_1.ymin],
                                                            layer = layer_barrier[layer])
    
    sensordot_1_barrier_connector_2 = D1 << pg.polygon_ports(xpts = [taper_barrier_sensordot_2.xmin - large_tilt,
                                                                    taper_barrier_sensordot_2.xmin,
                                                                    taper_barrier_sensordot_2.xmin,
                                                                    taper_barrier_sensordot_2.xmin - large_tilt],
                                                            ypts = [taper_barrier_sensordot_2.ymax,
                                                                    taper_barrier_sensordot_2.ymax,
                                                                    taper_barrier_sensordot_2.ymin,
                                                                    taper_barrier_sensordot_2.ymax],
                                                            layer = layer_barrier[layer])
    
    sensordot_2_barrier_connector_1 = D1 << pg.polygon_ports(xpts = [taper_barrier_sensordot_3.xmax + large_tilt,
                                                                    taper_barrier_sensordot_3.xmax,
                                                                    taper_barrier_sensordot_3.xmax,
                                                                    taper_barrier_sensordot_3.xmax + large_tilt],
                                                            ypts = [taper_barrier_sensordot_3.ymin,
                                                                    taper_barrier_sensordot_3.ymin,
                                                                    taper_barrier_sensordot_3.ymax,
                                                                    taper_barrier_sensordot_3.ymin],
                                                            layer = layer_barrier[layer])
    
    sensordot_2_barrier_connector_2 = D1 << pg.polygon_ports(xpts = [taper_barrier_sensordot_4.xmax + large_tilt,
                                                                    taper_barrier_sensordot_4.xmax,
                                                                    taper_barrier_sensordot_4.xmax,
                                                                    taper_barrier_sensordot_4.xmax + large_tilt],
                                                            ypts = [taper_barrier_sensordot_4.ymax,
                                                                    taper_barrier_sensordot_4.ymax,
                                                                    taper_barrier_sensordot_4.ymin,
                                                                    taper_barrier_sensordot_4.ymax],
                                                            layer = layer_barrier[layer])
    
    # Build the ohmic contacts

    ohmic_length = width_connector

    ohmic_1 = D1 << pg.polygon_ports(xpts = [sensordot_1.center[0] - ohmic_length,
                                             sensordot_1.center[0] + ohmic_length,
                                             sensordot_1.center[0] + ohmic_length,
                                             sensordot_1.center[0] - ohmic_length],
                                     ypts = [sensordot_1.center[1] + sensordot_radius + width_connector - 3*overlap,
                                             sensordot_1.center[1] + sensordot_radius + width_connector - 3*overlap + ohmic_length,
                                             sensordot_1.center[1] + sensordot_radius + width_connector - 3*overlap,
                                             sensordot_1.center[1] + sensordot_radius + width_connector - 3*overlap],
                                    layer = layer_ohmic[layer])
    
    ohmic_2 = D1 << pg.polygon_ports(xpts = [sensordot_1.center[0] - ohmic_length,
                                             sensordot_1.center[0] + ohmic_length,
                                             sensordot_1.center[0] + ohmic_length,
                                             sensordot_1.center[0] - ohmic_length],
                                     ypts = [sensordot_1.center[1] - sensordot_radius - width_connector + 3*overlap,
                                             sensordot_1.center[1] - sensordot_radius - width_connector + 3*overlap - ohmic_length,
                                             sensordot_1.center[1] - sensordot_radius - width_connector + 3*overlap,
                                             sensordot_1.center[1] - sensordot_radius - width_connector + 3*overlap],
                                    layer = layer_ohmic[layer])
    
    ohmic_3 = D1 << pg.polygon_ports(xpts = [sensordot_2.center[0] - ohmic_length,
                                             sensordot_2.center[0] + ohmic_length,
                                             sensordot_2.center[0] - ohmic_length,
                                             sensordot_2.center[0] - ohmic_length],
                                     ypts = [sensordot_2.center[1] + sensordot_radius + width_connector - 3*overlap,
                                             sensordot_2.center[1] + sensordot_radius + width_connector - 3*overlap,
                                             sensordot_2.center[1] + sensordot_radius + width_connector - 3*overlap + ohmic_length,
                                             sensordot_2.center[1] + sensordot_radius + width_connector - 3*overlap],
                                    layer = layer_ohmic[layer])
    
    ohmic_4 = D1 << pg.polygon_ports(xpts = [sensordot_2.center[0] - ohmic_length,
                                             sensordot_2.center[0] + ohmic_length,
                                             sensordot_2.center[0] - ohmic_length,
                                             sensordot_2.center[0] - ohmic_length],
                                     ypts = [sensordot_2.center[1] - sensordot_radius - width_connector + 3*overlap,
                                             sensordot_2.center[1] - sensordot_radius - width_connector + 3*overlap,
                                             sensordot_2.center[1] - sensordot_radius - width_connector + 3*overlap - ohmic_length,
                                             sensordot_2.center[1] - sensordot_radius - width_connector + 3*overlap],
                                    layer = layer_ohmic[layer])
    
    ohmic_resonator = D1 << pg.polygon_ports(xpts = [resonator.center[0] - ohmic_length,
                                                     resonator.center[0] + ohmic_length,
                                                     resonator.center[0] + ohmic_length,
                                                     resonator.center[0] - ohmic_length],
                                             ypts = [resonator.center[1] + semi_minor_resonator + width_connector - 4*overlap,
                                                     resonator.center[1] + semi_minor_resonator + width_connector - 4*overlap,
                                                     resonator.center[1] + semi_minor_resonator + width_connector - 4*overlap + barrier_width,
                                                     resonator.center[1] + semi_minor_resonator + width_connector - 4*overlap + barrier_width],
                                             layer = layer_ohmic[layer])
    

    # Done building all the internal objects, now go to routing and creating relevant bonding pads

    # Create the bonding pads, and route all relevant paths to them

    pu.write_lyp('1x4_inner_structure_props.lyp', layerset = lys)

    D1.write_gds(filename = '1x4_inner_structure.gds', # Output GDS file name
            unit = 1e-9,                  # Base unit (1e-6 = microns)
            precision = 1e-9,             # Precision / resolution (1e-9 = nanometers)
            auto_rename = True,           # Automatically rename cells to avoid collisions
            max_cellname_length = 28,     # Max length of cell names
            cellname = 'toplevel'         # Name of output top-level cell
           )

    ext_ports = [[], [], [], []]

    N_numbers = [2, 4, 5, 6, 7, 8, 9]           # Setting the specific place where the gates on the device are routed to
    S_numbers = [3, 4, 5, 6, 7, 8, 11]
    W_numbers = [2, 3, 5, 6, 7, 8, 10, 11]
    E_numbers = [2, 3, 5, 6, 7, 8, 10, 11]

    N_bond_numbers = [3, 4, 5, 6, 7, 8, 9]      # Setting the specific places where bonding pads should be
    S_bond_numbers = [3, 4, 5, 6, 7, 8, 9]
    W_bond_numbers = [2, 3, 5, 6, 7, 8, 10, 11]
    E_bond_numbers = [2, 3, 5, 6, 7, 8, 10, 11]

    list_of_numbers = [N_numbers, S_numbers, W_numbers, E_numbers]


    counter_1 = 0

    for i in list_of_numbers:
        for j in range(len(i)):
            if counter_1 == 0:
                ext_ports[counter_1].append('EXT_N' + str(j))
                ext_ports[counter_1][j] = D1.add_port(name='EXT_N' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]),ymax_ext_box], width = width_external_port, orientation = 270)
            if counter_1 == 1:
                ext_ports[counter_1].append('EXT_S' + str(j))
                ext_ports[counter_1][j] = D1.add_port(name='EXT_S' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]),ymin_ext_box], width = width_external_port, orientation = 90)
            if counter_1 == 2:
                ext_ports[counter_1].append('EXT_W' + str(j))
                ext_ports[counter_1][j] = D1.add_port(name='EXT_W' + str(j), midpoint = [xmin_ext_box, y_midpoint_ext_port_v(i[j])], width = width_external_port, orientation = 0)
            if counter_1 == 3:
                ext_ports[counter_1].append('EXT_E' + str(j))
                ext_ports[counter_1][j] = D1.add_port(name='EXT_E' + str(j), midpoint = [xmax_ext_box, y_midpoint_ext_port_v(i[j])], width = width_external_port, orientation = 180)
        counter_1 += 1

    layer_N = [layer_dot[layer], layer_screening[layer], layer_dot[layer], layer_ohmic[layer], layer_barrier[layer], layer_barrier[layer], layer_barrier[layer]]
    layer_S = [layer_barrier[layer], layer_barrier[layer], layer_dot[layer], layer_screening[layer], layer_dot[layer], layer_screening[layer], layer_dot[layer]]
    layer_W = [layer_barrier[layer], layer_ohmic[layer], layer_barrier[layer], layer_screening[layer], layer_dot[layer], layer_barrier[layer], layer_ohmic[layer], layer_barrier[layer]]
    layer_E = [layer_barrier[layer], layer_ohmic[layer], layer_barrier[layer], layer_dot[layer], layer_screening[layer], layer_barrier[layer], layer_ohmic[layer], layer_barrier[layer]]

    North_objects_to_route = [dot_1_connector.ports['1'], screening_dots_1.ports['1'], dot_2_connector.ports['1'],
                            ohmic_resonator.ports['3'], barrier_resonator_connector.ports['3'], vertical_barrier_3.ports['2'],
                            vertical_barrier_4.ports['2']]

    South_objects_to_route = [vertical_barrier_1.ports['3'], vertical_barrier_2.ports['3'], resonator_connector.ports['S']
                              , screening_resonator_connector.ports['4'], dot_3_connector.ports['2'], screening_dots_2.ports['2']
                              , dot_4_connector.ports['2']]
    
    West_objects_to_route = [vertical_barrier_sensordot_1.ports['4'], ohmic_2.ports['1'], sensordot_1_barrier_connector_2.ports['3'],
                             screening_sensordot1_connector.ports['3'], sensordot_1_connector.ports['W'], sensordot_1_barrier_connector_1.ports['1'],
                             ohmic_1.ports['1'], vertical_barrier_sensordot_1.ports['2']]
    
    East_objects_to_route = [vertical_barrier_sensordot_2.ports['4'], ohmic_4.ports['2'], sensordot_2_barrier_connector_2.ports['3'], 
                            sensordot_2_connector.ports['E'], screening_sensordot2_connector.ports['2'], sensordot_2_barrier_connector_1.ports['3'],
                            ohmic_3.ports['2'], vertical_barrier_sensordot_2.ports['2']]

    route_ext_N = []
    route_ext_S = []
    route_ext_W = []
    route_ext_E = []

    counter_1 = 0

    for i in ext_ports:
        for j in range(len(i)):
            if counter_1 == 0:
                route_ext_N.append('route_ext_N' + str(j))
                route_ext_N[j] = D1.add_ref(pr.route_quad(North_objects_to_route[j], i[j], layer = layer_N[j]))
            if counter_1 == 1:
                route_ext_S.append('route_ext_S' + str(j))
                route_ext_S[j] = D1.add_ref(pr.route_quad(South_objects_to_route[j], i[j], layer = layer_S[j]))
            if counter_1 == 2:
                route_ext_W.append('route_ext_W' + str(j))
                route_ext_W[j] = D1.add_ref(pr.route_quad(West_objects_to_route[j], i[j], layer = layer_W[j]))
            if counter_1 == 3:
                route_ext_E.append('route_ext_E' + str(j))
                route_ext_E[j] = D1.add_ref(pr.route_quad(East_objects_to_route[j], i[j], layer = layer_E[j]))
        counter_1 += 1
    D_without_local_markers = pg.union(D1, by_layer = True)

    # Adding the local markers

    inner_marker_NE = D1 << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_NE.move(destination = (np.cos(45)*length_to_local, np.cos(45)*length_to_local))

    inner_marker_NW = D1 << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_NW.move(destination = (-np.cos(45)*length_to_local, np.cos(45)*length_to_local))

    inner_marker_SE = D1 << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_SE.move(destination = (np.cos(45)*length_to_local, -np.cos(45)*length_to_local))

    inner_marker_SW = D1 << inner_marker_generator(width_local_marker, length_local_marker, layer_marker)
    inner_marker_SW.move(destination = (-np.cos(45)*length_to_local, -np.cos(45)*length_to_local))

    D_with_local_markers = pg.union(D1, by_layer = True)

    bond_numbers = [N_bond_numbers, S_bond_numbers, W_bond_numbers, E_bond_numbers]

    # Create a list with the relevant layers for the bonding pads, based on the layers of the internal objects

    layer_N_pads = [layer_dot_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad, layer_ohmic_bond_pad, layer_barrier_gates_bond_pad, layer_barrier_gates_bond_pad, layer_barrier_gates_bond_pad]
    layer_S_pads = [layer_barrier_gates_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad]
    layer_W_pads = [layer_barrier_gates_bond_pad, layer_ohmic_bond_pad, layer_barrier_gates_bond_pad, layer_screening_gates_bond_pad, layer_dot_bond_pad, layer_barrier_gates_bond_pad, layer_ohmic_bond_pad, layer_barrier_gates_bond_pad]
    layer_E_pads = [layer_barrier_gates_bond_pad, layer_ohmic_bond_pad, layer_barrier_gates_bond_pad, layer_dot_bond_pad, layer_screening_gates_bond_pad, layer_barrier_gates_bond_pad, layer_ohmic_bond_pad, layer_barrier_gates_bond_pad]

    layer_pads = [layer_N_pads, layer_S_pads, layer_W_pads, layer_E_pads]

    return D_with_local_markers, D_without_local_markers, D1, list_of_numbers, ext_ports, bond_numbers, layer_pads

def build_outer_structure(D1, list_of_numbers, bond_numbers, ext_ports, prod_pads, layer_pads):

    xmax_bonding_box = D1.center[0] + width_bonding_box/2
    xmin_bonding_box = D1.center[0] - width_bonding_box/2

    ymax_bonding_box = D1.center[1] + height_bonding_box/2
    ymin_bonding_box = D1.center[1] - height_bonding_box/2

    spacing_h_bonding_pad = (width_bonding_box - N_ports_per_edge*width_bonding_pad)/(N_ports_per_edge+1)
    spacing_v_bonding_pad = (height_bonding_box - N_ports_per_edge*width_bonding_pad)/(N_ports_per_edge+1)

    def x_min_bonding_pad_h(N):
        x0 = spacing_h_bonding_pad*N+width_bonding_pad*(N-1) + xmin_bonding_box
        return x0
    def y_min_bonding_pad_v(N):
        y0 = spacing_v_bonding_pad*N+width_bonding_pad*(N-1) + ymin_bonding_box
        return y0

    bond_pads_N, bond_pads_S, bond_pads_W, bond_pads_E, bond_pads_port_N, bond_pads_port_S, bond_pads_port_W, bond_pads_port_E = [], [], [], [], [], [], [], []

    mid_port_N, mid_port_S, mid_port_W, mid_port_E = [], [], [], []

    counter_1 = 0

    for i in list_of_numbers:
        for j in range(len(i)):
            if counter_1 == 0:
                bond_pads_N.append("bonding_pad_N" + str(i[j]))
                bond_pads_port_N.append("bond_port_pad_N" + str(i[j]))
                bond_pads_N[j] = D1 << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[0][j]).rotate(90)
                bond_pads_N[j].ymin = ymax_bonding_box
                bond_pads_N[j].xmin = x_min_bonding_pad_h(int(bond_numbers[0][j]))
                bond_pads_port_N[j] = D1.add_port(name='B_N' + str(i[j]), midpoint = [bond_pads_N[j].center[0], bond_pads_N[j].ymin], width = width_port_bonding_pad, orientation = 270)
                D1.add_ref(pr.route_quad(bond_pads_port_N[j], ext_ports[0][j], layer = layer_pads[0][j]))
                mid_port_N.append('mid_port_N' + str(j))
                mid_port_N[j] = D1.add_port(name='MID_N' + str(j), midpoint = [x_midpoint_ext_port_h(i[j]) + 70*(5.5 - i[j]), ymax_ext_box - closing_distance], width = width_middle_port, orientation = 270)
                D1.add_ref(pr.route_quad(ext_ports[0][j], mid_port_N[j], layer = layer_pads[0][j]))
            if counter_1 == 1:
                bond_pads_S.append("bonding_pad_S" + str(i[j]))
                bond_pads_port_S.append("bond_port_pad_S" + str(i[j]))
                bond_pads_S[j] = D1 << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[1][j]).rotate(90)
                bond_pads_S[j].ymax = ymin_bonding_box
                bond_pads_S[j].xmin = x_min_bonding_pad_h(bond_numbers[1][j])
                bond_pads_port_S[j] = D1.add_port(name='B_S' + str(i[j]), midpoint = [bond_pads_S[j].center[0], bond_pads_S[j].ymax], width = width_port_bonding_pad, orientation = 90)
                D1.add_ref(pr.route_quad(bond_pads_port_S[j], ext_ports[1][j], layer = layer_pads[1][j]))
                mid_port_S.append('mid_port_S' + str(j))
                mid_port_S[j] = D1.add_port(name='MID_S' + str(j), midpoint = [x_midpoint_ext_port_h(i[j])  + 70*(5.5 - i[j]), ymin_ext_box + closing_distance], width = width_middle_port, orientation = 90)
                D1.add_ref(pr.route_quad(ext_ports[1][j], mid_port_S[j], layer = layer_pads[1][j]))
            if counter_1 == 2:
                bond_pads_W.append("bonding_pad_W" + str(i[j]))
                bond_pads_port_W.append("bond_port_pad_W" + str(i[j]))
                bond_pads_W[j] = D1 << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[2][j])
                bond_pads_W[j].xmax = xmin_bonding_box
                bond_pads_W[j].ymin = y_min_bonding_pad_v(bond_numbers[2][j])
                bond_pads_port_W[j] = D1.add_port(name='B_W' + str(i[j]), midpoint = [bond_pads_W[j].xmax, bond_pads_W[j].center[1]], width = width_port_bonding_pad, orientation = 0)
                D1.add_ref(pr.route_quad(bond_pads_port_W[j], ext_ports[2][j], layer = layer_pads[2][j]))
                mid_port_W.append('mid_port_W' + str(j))
                mid_port_W[j] = D1.add_port(name='MID_W' + str(j), midpoint = [xmin_ext_box + closing_distance, y_midpoint_ext_port_v(i[j])+ 70*(5.5 - i[j])], width = width_middle_port, orientation = 0)
                D1.add_ref(pr.route_quad(ext_ports[2][j], mid_port_W[j], layer = layer_pads[2][j]))
            if counter_1 == 3:
                bond_pads_E.append("bonding_pad_E" + str(i[j]))
                bond_pads_port_E.append("bond_port_pad_E" + str(i[j]))
                bond_pads_E[j] = D1 << pg.rectangle(size = (length_bonding_pad, width_bonding_pad), layer = layer_pads[3][j])
                bond_pads_E[j].xmin = xmax_bonding_box
                bond_pads_E[j].ymin = y_min_bonding_pad_v(bond_numbers[3][j])
                bond_pads_port_E[j] = D1.add_port(name='B_E' + str(i[j]), midpoint = [bond_pads_E[j].xmin, bond_pads_E[j].center[1]], width = width_port_bonding_pad, orientation = 180)
                D1.add_ref(pr.route_quad(bond_pads_port_E[j], ext_ports[3][j], layer = layer_pads[3][j]))
                mid_port_E.append('mid_port_E' + str(j))
                mid_port_E[j] = D1.add_port(name='MID_E' + str(j), midpoint = [xmax_ext_box - closing_distance, y_midpoint_ext_port_v(i[j]) + 70*(5.5 - i[j])], width = width_middle_port, orientation = 180)
                D1.add_ref(pr.route_quad(ext_ports[3][j], mid_port_E[j], layer = layer_pads[3][j]))
        counter_1 += 1

    counter_1 = 0

    for i in prod_pads:
        for j in range(len(i)):
            if counter_1 == 0:
                protection_pads_N = protection_pads(D1, (bond_pads_N[i[j][0]].xmax + protection_extra - (bond_pads_N[i[j][1]].xmin - protection_extra)), length_bonding_pad,
                    destination = (bond_pads_N[i[j][1]].xmin - protection_extra, bond_pads_N[i[j][1]].ymin + protection_extra))
            if counter_1 == 1:
                protection_pads_S = protection_pads(D1, (bond_pads_S[i[j][0]].xmax + protection_extra - (bond_pads_S[i[j][1]].xmin - protection_extra)), length_bonding_pad,
                    destination = (bond_pads_S[i[j][1]].xmin - protection_extra, bond_pads_S[i[j][1]].ymin - protection_extra))
            if counter_1 == 2:
                protection_pads_W = protection_pads(D1, length_bonding_pad, (bond_pads_W[i[j][0]].ymax + protection_extra - (bond_pads_W[i[j][1]].ymin - protection_extra)),
                    destination = (bond_pads_W[i[j][1]].xmin - protection_extra, bond_pads_W[i[j][1]].ymin - protection_extra))
            if counter_1 == 3:
                protection_pad_E = protection_pads(D1, length_bonding_pad, (bond_pads_E[i[j][0]].ymax + protection_extra - (bond_pads_E[i[j][1]].ymin - protection_extra)),
                    destination = (bond_pads_E[i[j][1]].xmin + protection_extra, bond_pads_E[i[j][1]].ymin - protection_extra))
        counter_1 += 1

    D1_joined_by_layer = pg.union(D1, by_layer = True)

    return D1_joined_by_layer

