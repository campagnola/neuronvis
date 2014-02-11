#!/usr/bin/python
__author__ = 'pbmanis'
"""
hocRender : provide visual rendering for morphology and other attributes
as stored in a "hoc" file.
Usage:

h.loadfile(filename) # standard hoc load
# potentially, you would decorate the membrane with biophysical mechanisms here:
decorate(h)

pg.mkQApp()
pg.dbg()
render = hr.hocRender(h) # where h is the NEURON hoc object (from neuron import h)
render.draw_model(modes=['blob'])
render.getSectionLists(Colors.keys()) # if the sections are named...
render.paint_sections_by_density(self.modelPars.calyxColors, self.modelPars.mechNames['CaPCalyx'])
render.show()

2/3/2014
Portions of this code were taken from neuronvisio (http://neuronvisio.org), specifically, to parse
the hoc file connection structure (specifically: getSectionInfo, and parts of drawModel).

"""


import os, sys, pickle
import pyqtgraph as pg
import numpy as np

# define all commands here.
commands = {
    'sec-type': "Sections colored by type",
    'vm': "Animation of per-section membrane voltage over time.",
    'graph': "Simple wireframe rendering.",
    'cylinders': "Simple cylinder rendering.",
    }

# Handle command line arguments.
# might introduce real command line parsing later..
##########################################################

def error():
    print("Usage: hocRender input_file command")
    print("  input_file may be either *.hoc defining section properties or")
    print("  *.p containing simulation results.\n")
    print("Available commands:")
    for cmd,desc in commands.items():
        print("  "+cmd+":")
        print("\n".join(["    "+line for line in desc.split("\n")]))
    sys.exit(-1)

if len(sys.argv) < 2 or not os.path.isfile(sys.argv[1]):
    error()
input_file = sys.argv[1]

if len(sys.argv) > 2:
    command = sys.argv[2]
    if command not in commands:
        error()
else:
    command = 'sec-type'


# read input file(s)
if input_file.endswith('.p'):
    from sim_result import SimulationResult
    sim_data = SimulationResult(input_file)
    hoc_file = sim_data.hoc_file
elif input_file.endswith('.hoc'):
    sim_data = None
    hoc_file = input_file
else:
    error()

# import here so we can parse commands more quickly 
# (and without neuron garbage)
from hoc_reader import HocReader
from hoc_viewer import HocViewer
hoc = HocReader(hoc_file)
view = HocViewer(hoc)


# Handle commands
##########################################################

section_colors = {
    'axon': 'r', 
    'heminode': 'g', 
    'stalk':'y', 
    'branch': 'b', 
    'neck': 'brown',
    'swelling': 'magenta', 
    'tip': 'powderblue', 
    'parentaxon': 'orange', 
    'synapse': 'k'}

if command == 'sec-type':
    # Color sections by type.
    hoc.read_hoc_section_lists(section_colors.keys())
    surf = view.draw_surface()
    surf.set_group_colors(section_colors, alpha=0.2)
elif command == 'graph':
    g = view.draw_graph()
    hoc.read_hoc_section_lists(section_colors.keys())
    g.set_group_colors(section_colors)
elif command == 'cylinders':
    g = view.draw_cylinders()
    hoc.read_hoc_section_lists(section_colors.keys())
    g.set_group_colors(section_colors)
elif command == 'vm':
    
    # Render animation of membrane voltage
    if sim_data is None:
        raise Exception('Cannot render Vm: no simulation output specified.')

    surf = view.draw_surface()
    start = 375
    stop = 550
    index = start

    def vm_to_color(v):
        """
        Convert an array of Vm to array of representative colors
        """
        color = np.empty((v.shape[0], 4), dtype=float)
        v_min = -80 # mV
        v_range = 100. # mV range in scale
        v = (v - v_min) / v_range
        color[:,0] = v     # R
        color[:,1] = 1.5*abs(v-0.5) # G
        color[:,2] = 1.-v # B
        color[:,3] = 0.1+0.8*v # alpha
        return color

    def set_index(index):
        """
        Set the currently-displayed time index.
        """
        v = sim_data.data['Vm'][:,index]
        color = vm_to_color(v)
        
        # note that we assume sections are ordered the same in the HocReader 
        # as they are in the results data, but really we should use 
        # sim_data.section_map to ensure this is True.
        surf.set_section_colors(color)
        

    def update():
        global index, start, stop, sim_data, surf

        set_index(index)
        
        index += 1
        if index >= stop:
            index = start

    def record(file_name):
        """
        Record a video from *start* to *stop* with the current view
        configuration.
        """
        timer.stop()
        view.begin_video(file_name)
        try:
            for i in range(start, stop):
                set_index(i)
                pg.Qt.QtGui.QApplication.processEvents()
                view.save_frame()
                print("%d / %d" % (i, stop))
        finally:
            view.save_video()

    timer = pg.QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(10.)


if sys.flags.interactive == 0:
    import pyqtgraph as pg
    pg.Qt.QtGui.QApplication.exec_()
