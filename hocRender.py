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


# define all commands here.
commands = {
    'sec-type': "Sections colored by type",
    'vm': "Animation of per-section membrane voltage over time.",
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
    sim_data = pickle.load(input_file)
    try:
        hoc_file = sim_data['hoc_file']
    except KeyError:
        raise Exception(".p file must have 'hoc_file' key!")
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

if command == 'sec-type':
    # Color sections by type.
    section_colors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}
    hoc.read_hoc_section_lists(section_colors.keys())
    surf = view.draw_surface()
    surf.set_group_colors(section_colors, alpha=0.2)
    
elif command == 'vm':
    # Render animation of membrane voltage
    if sim_data is None:
        raise Exception('Cannot render Vm: no simulation output specified.')

    maxcycles = 1
    cycles = 0
    videoOpen = False
    index = 0
    vfile = None
    vdata['data'] = vdata['data'][:,375:550]
    color = np.zeros((vdata['data'].shape[0], 4), dtype=float)

    def update():
        global index, color, render, vdata, videoOpen, vfile, winsize, cycles, maxcycles
        # color goes from 0 to 255
        # vdata range goes from -80 to +20
        # so scale so that -80 to +20

        v0 = -80 # mV
        vr = 100. # mV range in scale
        v = (vdata['data'][:,index] - v0) / vr
        color[:,0] = v     # R
        color[:,1] = 1.5*abs(v-0.5) # G
        color[:,2] = 1.-v # B
        color[:,3] = 0.1+0.8*v # alpha
        render.set_section_colors(color)
        index = (index + 1) % vdata['data'].shape[1]
        if index == 0:
            cycles += 1
            if cycles >= maxcycles:
                timer.stop()
                vfile.release()
                return
        img0 = pg.imageToArray(render.w.readQImage())
        if not videoOpen:
            import cv
            import cv2
            winsize = img0.shape[0:2]
            vfile = cv2.VideoWriter()
            vfile.open(filename='~/Desktop/Python/PyNeuronLibrary/CalyxModel/video.avi',
                        fourcc=cv.CV_FOURCC('M', 'P', '4', 'V'), fps=25, frameSize = winsize,
                        isColor=False)
            videoOpen = True
        if vfile.isOpened():
            print 'writing image'
            vfile.write(img0)


    timer = pg.QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(10.)


if sys.flags.interactive == 0:
    import pyqtgraph as pg
    pg.Qt.QtGui.QApplication.exec_()



# Functions to be used from interactive prompt:
##########################################################

def render_stack(start, stop):
    global index, render, timer
    timer.stop()

    for index in range(start, stop):
        update()
        #success0, img0 = capwin.read()
        img0 = pg.imageToArray(render.w.readQImage())
        # print 'write frame: %d' % index
        # print '    size: ', img0.shape
        # print '    max img0: ', np.max(img0[0:2])
        vfile.write(img0)
        #  render.w.readQImage().save('render-%03d.png' % (index-start))
#render_stack(0, vdata['data'].shape[1])
#    print dir(vfile)

