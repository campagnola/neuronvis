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
render.paintSectionsByDensity(self.modelPars.calyxColors, self.modelPars.mechNames['CaPCalyx'])
render.show()

2/3/2014
Portions of this code were taken from neuronvisio (http://neuronvisio.org), specifically, to parse
the hoc file connection structure (specifically: getSectionInfo, and parts of drawModel).

"""

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import scipy.ndimage
import neuron as h
import neuron as h
from neuron import *
import re
import pickle
import cv
import cv2

class hocRender():
    def __init__(self, h):
        self.h = h
        self.app = pg.mkQApp()
        self.w = gl.GLViewWidget()
        self.w.resize(720,720)
        self.w.show()
        self.w.setWindowTitle('hocRender')
        self.w.setCameraPosition(distance=300.)

        self.g = gl.GLGridItem()
        self.g.scale(2,2,1)
        self.w.addItem(self.g)
        self.sec2coords = {}
        # we need to get the section
        self.cyl2sec = {}
        self.nsec = 0
        for sec in self.h.allsec():
            self.nsec += 1

        # Needed to update the value of a cyl bound to a section
        self.sec2cyl = {}
        self.seg2id = {}
        self.sec2coords = {}
        self.connections = []
        self.n3dpoints_per_sec = {}
        self.Sections = {}
        self.Mechanisms = {} # mechanisms are stored by "axon" section number

        self.Colors = { # colormap
            'b': np.array([0,0,255,255])/255.,
            'blue': np.array([0,0,255,255])/255.,
            'g': np.array([0,255,0,255])/255.,
            'green': np.array([0,255,0,255])/255.,
            'r': np.array([255,0,0,255])/255.,
            'red': np.array([255,0,0,255])/255.,
            'cyan': np.array([0,255,255,255])/255.,
            'c': np.array([0,255,255,255])/255.,
            'm': np.array([255,0,255,255])/255.,
            'magenta': np.array([255,0,255,255])/255.,
            'y': np.array([255,255,0,255])/255.,
            'yellow': np.array([255,255,0,255])/255.,
            'k': np.array([0,0,0,255])/255.,
            'black': np.array([0,0,0,255])/255.,
            'w': np.array([255,255,255,255])/255.,
            'white': np.array([255,255,255,255])/255.,
            'd': np.array([150,150,150,255])/255.,
            'dark': np.array([150,150,150,255])/255.,
            'l': np.array([200,200,200,255])/255.,
            'light': np.array([200,200,200,255])/255.,
            's': np.array([100,100,150,255])/255.,
            'powderblue': np.array([176,230,230,255])/255.,
            'brown': np.array([180,25,25,255])/255.,
            'orange': np.array([255,180,0,255])/255.,
            'pink': np.array([255,190,206,255])/255.,
        }


    def show(self):
        QtGui.QApplication.instance().exec_()


    def get_mechanisms(self, section):
        """
        Get a list of all of the mechanisms inserted into a given section
        Input: section (hoc object)
        Returns: list of mechs:
        Side-effects: None
        """
        mechs = []
        for seg in section:
            for mech in seg:
                mechs.append(mech.name())
        mechs = set(mechs) # Excluding the repeating ones
        return mechs


    def get_density(self, section, mechanism):
        """
        Get density mechanism that may be found the section.
        mechanism is a list ['name', 'gbarname']. This is needed because
        some mechanisms do not adhere to any convention and may have a different
        kind of 'gbarname' than 'gbar<name>_mechname'
        returns the average of the conductance density, as that may range across different
        values in a section (e.g., can vary by segments)
        Input: section (hoc object)
                mechanism mechanism is a list ['name', 'gbarname'].
        Output:
            mean conductance inserted into the section across segments
        Side-effects:
            None
        """
        gmech = []
        for seg in section:
            try:
                x =  eval('seg.%s' % mechanism[0])
                mecbar = '%s_%s' % (mechanism[1], mechanism[0])
                if mecbar in dir(x):
                    gmech.append(eval('x.%s' % mechanism[1]))
                else:
                    print 'hocRender:get_density did not find the mechanism in dir x'
            except:
                print 'hocRender:get_density failed to evaluate the mechanisms... '

#        print gmech
        if len(gmech) == 0:
            gmech = 0.
        return np.mean(gmech)


    def get_sec_info(self, section):
        """
        Get the info of the given section
        modified from: neuronvisio
        """
        info = "<b>Section Name:</b> %s<br/>" %section.name()
        info += "<b>Length [um]:</b> %f<br/>" % section.L
        info += "<b>Diameter [um]:</b> %f<br/>" % section.diam
        info += "<b>Membrane Capacitance:</b> %f<br/>" % section.cm
        info += "<b>Axial Resistance :</b> %f<br/>" % section.Ra
        info += "<b>Number of Segments:</b> %f<br/>" % section.nseg
        mechs = []
        for seg in section:
            for mech in seg:
                mechs.append(mech.name())
        mechs = set(mechs) # Excluding the repeating ones

        mech_info = "<b>Mechanisms in the section</b><ul>"
        for mech_name in mechs:
            s = "<li> %s </li>" % mech_name
            mech_info += s
        mech_info += "</ul>"
        info += mech_info
        return info


    def getSectionLists(self, names):
        """
        Get the list of sections that have names corresponding to the list in names.
        Note that this looks in "axon" again. May break generality
        Side effects (modifies):
            self.Sections
            self.Mechanisms
        returns: Nothing.
        """
        self.Sections = {key: None for key in names}
        self.Mechanisms = {key: None for key in names}
        src = re.compile('axon\[(\d*)\]')
        for name in self.Sections: # for each part
            x=eval("h.%s" % name) # find the associated variable
            axno = []
            for sec in x: # for each section in the list with this name (cell part)
                g = src.match(sec.name())
                axno.append(int(g.groups()[0])) # get the axon # and populate the dictionary
                mech = self.get_mechanisms(sec)
                if self.Mechanisms[name] == None: # fill first
                    self.Mechanisms[name] = [(axno, mech)]
                else:
                    self.Mechanisms[name].append([(axno, mech)])
            self.Sections[name]  = list(axno) # and populate it with the axon #'s of each assigned section


    def draw_model(self, modes=["cylinder"]):
        """
        modified from:neuronvisio
        Draw the model.
        Inputs: modes - a list of potential representations, including:
            cylinder, blob, lines, mesh, and sphere
        Outputs: None (display)
        Side effects:
            modifies cyl2sec, sec2cyl, vertexes and edges.
        """
        self.h.define_shape()

        # map segments (lines) to the section that contains them
        self.segment_to_section = {}
        
        vertexes = []
        connections = []
        for sec in self.h.allsec():
            print sec.name()
            m = re.match(r'axon\[(\d+)\]', sec.name())
            if m is None:
                continue
                #raise Exception('Could not determine section ID from name: %s' % sec.name())
            secid = int(m.groups()[0])
            
            x_sec, y_sec, z_sec, d_sec = self.retrieve_coordinate(sec)
            self.sec2coords[sec.name()] = [x_sec, y_sec, z_sec]
            # Store the section. later.
            radius = sec.diam/2.
            sec_coords_bound = ((x_sec.min(), x_sec.max()),
                                (y_sec.min() - radius,
                                 y_sec.max() + radius),
                                (z_sec.min() - radius,
                                 z_sec.max() + radius))
            self.cyl2sec[sec_coords_bound] = sec
            self.sec2cyl[sec] = sec_coords_bound

            for i,xi in enumerate(x_sec):
                vertexes.append((x_sec[i], y_sec[i], z_sec[i], d_sec[i], secid))
                indx_geom_seg = len(vertexes) - 1
                if len(vertexes) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])


        self.edges  = np.array(connections)
        self.vertexes = np.array(vertexes)

        for mode in modes:
            if mode == 'blob':
                self.drawBlob()
            elif mode == 'volume':
                self.drawVolume()
            else:
                self.drawMeshes(mode)

        
    def makeVolumeData(self):
        """
        Using the current state of vertexes, edges, generates a scalar field
        Returns:
            scalar field, idfield (for mapping) and tranform
        """
        res = 0.4 # resolution of scalar field in microns
        maxdia = 10. # maximum diameter (defines shape of kernel)
        kernel_size = int(maxdia/res) + 1 # width of kernel
        
        # read vertex data
        verlocs = self.vertexes[:, :3]
        x = verlocs[:,0]
        y = verlocs[:,1]
        z = verlocs[:,2]
        d = self.vertexes[:,3]
        sec_id = self.vertexes[:,4]
        
        lines = self.edges
        # decide on dimensions of scalar field
        mx = verlocs[...,:-3].max(axis=0)  # note: skip last 3 due to junk in hoc files
        mn = verlocs.min(axis=0)

        xd = (np.max(x) - np.min(x))
        yd = (np.max(y) - np.min(y))
        zd = (np.max(z[:-3]) - np.min(z))
        nx = int(xd/res) + kernel_size
        ny = int(yd/res) + kernel_size
        nz = int(zd/res) + kernel_size
        
        # prepare blank scalar field for drawing
        scfield = np.zeros((nx, ny, nz), dtype=np.float32)
        scfield[:] = -1000
        
        # array for holding IDs of sections that contribute to each area
        idfield = np.empty((nx, ny, nz), dtype=int)
        idfield[:] = -1
        
        # map vertex locations to voxels
        verlocs -= np.array([[np.min(x), np.min(y), np.min(z)]])
        verlocs *= 1./res

        # Define kernel used to draw scalar field along dendrites
        def cone(i,j,k):
            # value decreases linearly with distance from center of kernel.
            w = kernel_size / 2
            return w - ((i-w)**2 + (j-w)**2 + (k-w)**2)**0.5
        kernel = res * np.fromfunction(cone, (kernel_size,)*3)
        kernel -= kernel.max()

        def array_intersection(arr1, arr2, pos):
            """
            Return slices used to access the overlapping area between two 
            arrays that are offset such that the origin of *arr2* is a *pos* 
            relative to *arr1*.            
            """
            s1 = [0]*3
            s2 = [0]*3
            t1 = [0]*3
            t2 = [0]*3
            pos = map(int, pos)
            for axis in range(3):
                s1[axis] = max(0, -pos[axis])
                s2[axis] = min(arr2.shape[axis], arr1.shape[axis]-pos[axis])
                t1[axis] = max(0, pos[axis])
                t2[axis] = min(arr1.shape[axis], pos[axis]+arr2.shape[axis])
            slice1 = (slice(t1[0],t2[0]), slice(t1[1],t2[1]), slice(t1[2],t2[2]))
            slice2 = (slice(s1[0],s2[0]), slice(s1[1],s2[1]), slice(s1[2],s2[2]))
            return slice1, slice2            

        maplocs = verlocs
        maplocs[:,0] = np.clip(maplocs[:,0], 0, scfield.shape[0]-1)
        maplocs[:,1] = np.clip(maplocs[:,1], 0, scfield.shape[1]-1)
        maplocs[:,2] = np.clip(maplocs[:,2], 0, scfield.shape[2]-1)
        for c in range(lines.shape[0] - 3):
            i = lines[c, 0]
            j = lines[c, 1]
            p1 = maplocs[i].copy()
            p2 = maplocs[j].copy()
            diff = p2-p1
            axis = np.argmax(np.abs(diff))
            dia = d[i]
            nvoxels = abs(int(diff[axis]))+1
            for k in range(nvoxels):
                kern = kernel + (dia/2.0)
                sl1, sl2 = array_intersection(scfield, kern, p1) # find the overlapping area between the field and the kernel
                idfield[sl1] = np.where(scfield[sl1] > kern[sl2], idfield[sl1], sec_id[i])
                scfield[sl1] = np.where(scfield[sl1] > kern[sl2], scfield[sl1], kern[sl2])
                #stamp_array(scfield, kern, p1)
                #stamp_array(idfield, kern, p1)
                dia += (d[j]-d[i]) / nvoxels
                p1 += diff / nvoxels
                
        # return transform relating volume data to original vertex data
        transform = pg.Transform3D()
        w = res * kernel_size / 2 # offset introduced due to kernel
        transform.translate(*(mn-w))
        transform.scale(res, res, res)
        transform.translate(1, 1, 1)
        return scfield, idfield, transform


    def drawVolume(self):
        scfield, idfield, transform = self.makeVolumeData()
        nfdata = np.empty(scfield.shape + (4,), dtype=np.ubyte)
        nfdata[...,0] = 255 #scfield*50
        nfdata[...,1] = 255# scfield*50
        nfdata[...,2] = 255# scfield*50
        nfdata[...,3] = np.clip(scfield*150, 0, 255)
        v = gl.GLVolumeItem(nfdata)
        v.setTransform(transform)
        self.w.addItem(v)


    def drawBlob(self):
        scfield, idfield, transform = self.makeVolumeData()
        scfield = scipy.ndimage.gaussian_filter(scfield, (0.5, 0.5, 0.5))
        #pg.image(scfield)
        verts, faces = pg.isosurface(scfield, level=0.0)
        vertexColors = np.empty((verts.shape[0], 4), dtype=float)
        md = gl.MeshData(vertexes=verts, faces=faces)
        
        # match vertexes to section IDs
        #vox_locations = pg.transformCoordinates(transform, verts, transpose=True).astype(int)
        vox_locations = verts.astype(int)
        # get sction IDs for each vertex
        self.mesh_sec_ids = idfield[vox_locations[:,0], vox_locations[:,1], vox_locations[:,2]] 
        
        mesh = gl.GLMeshItem(meshdata=md, smooth=True, shader='balloon')
        mesh.setTransform(transform)
        mesh.setGLOptions('additive')
        self.mesh = mesh
        self.w.addItem(mesh)


    def show_section(self, sec_id):
        """
        Set the color of sections in the list sec_id
        """
        colors = np.empty((len(self.mesh_sec_ids), 4))
        colors[:] = 0.3
        colors[sec_id] = 1
        self.set_section_colors(colors)
    

    def set_section_colors(self, sec_colors):
        """
        Set the colors of multiple sections
        """
        md = self.mesh.opts['meshdata']
        colors = sec_colors[self.mesh_sec_ids]
        md.setVertexColors(colors)
        self.mesh.meshDataChanged()


    def drawMeshes(self,  mode):
        """
        Draw remaining mesh figures
        mode here can be line, sphere or cylinder

        """
        verlocs = self.vertexes[:, :3]
        x = verlocs[:,0]
        y = verlocs[:,1]
        z = verlocs[:,2]
        d = self.vertexes[:,3]
        dmax = np.max(d)
        lines = np.vstack(self.edges)
        for c in range(len(lines)-3):
            i = lines[c, 0]
            j = lines[c, 1]
            pts = np.vstack([[x[i], x[j]],[y[i], y[j]],[z[i],z[j]]]).transpose()
            if mode == "line":
                plt = gl.GLLinePlotItem(pos=pts, width =(d[i]+d[j])/(2.), color=pg.glColor((int(255.*d[i]/dmax), 128)), connected=True)
                self.w.addItem(plt)
            elif mode == 'sphere':
                md = gl.MeshData.sphere(rows=10, cols=20, radius=d[i]/2.0) # , length=d(i))
                colors = np.ones((md.faceCount(), 4), dtype=float)
                colors[::2,0] = 0
                colors[:,1] = np.linspace(0, 1, colors.shape[0])
                md.setFaceColors(colors)
                m5 = gl.GLMeshItem(meshdata=md, smooth=False, drawEdges=False)
                m5.translate(x[i],y[i],z[i])
                self.w.addItem(m5)
            elif mode == "cylinder":
                cyllen = np.sqrt((x[j]-x[i])**2.0 + (y[j]-y[i])**2.0 + (z[j]-z[i])**2.0)
                md = gl.MeshData.cylinder(rows=1, cols=8, radius=[d[i]/2., d[j]/2.], length=cyllen)
                colors = np.ones((md.faceCount(), 4), dtype=float)
                colors[::2,0] = 0
                colors[:,1] = np.linspace(0, 1, colors.shape[0])
                #colors[:,1] = (int(255.*d[i]/dmax), 128)
                md.setFaceColors(colors)
                m5 = gl.GLMeshItem(meshdata=md, smooth=True, drawEdges=False)
                p2 = pg.Vector(x[j], y[j], z[j])
                p1 = pg.Vector(x[i], y[i], z[i])
                r = pg.Vector(0,0,1)
                axis = pg.QtGui.QVector3D.crossProduct(r, p2-p1)
                ang = r.angle(p2-p1)
                m5.rotate(ang, axis.x(), axis.y(), axis.z())
                m5.translate(x[i],y[i],z[i]+cyllen/2.0) # move into position
                self.w.addItem(m5)


    def retrieve_coordinate(self, sec):
        """Retrieve the coordinates of the section avoiding duplicates"""

        sec.push()
        x, y, z, d = [],[],[],[]

        tot_points = 0
        connect_next = False
        for i in range(int(self.h.n3d())):
            present = False
            x_i = self.h.x3d(i)
            y_i = self.h.y3d(i)
            z_i = self.h.z3d(i)
            d_i = self.h.diam3d(i)
            # Avoiding duplicates in the sec
            if x_i in x:
                ind = len(x) - 1 - x[::-1].index(x_i) # Getting the index of last value
                if y_i == y[ind]:
                    if z_i == z[ind]:
                        present = True

            if not present:
                k =(x_i, y_i, z_i)
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
                d.append(d_i)
        self.h.pop_section()
        #adding num 3d points per section
        self.n3dpoints_per_sec[sec.name()] = len(d)
        return (np.array(x),np.array(y),np.array(z),np.array(d))


    def paintSectionsByType(self, sectionColors, excludeSections = []):
        """
        Color the sections in the reconstruction according to their
        structural type, based on the mapping dictionary in sectionColors
        Inputs: sectionColors, a dictionary of section names and their associated colors
            excludeSections: the names of section types that should not be colored (such as "axon")
        Side-effects: none.
        """
        color = np.zeros((self.nsec, 4), dtype=float)
        for stype in self.Sections:
            if stype in excludeSections: # skip excluded sections and "parent" sections
                continue
            for sno in self.Sections[stype]:
                color[sno, :] = self.Colors[sectionColors[stype]] # map colors
                color[sno, 3] = 0.2 # alpha
        render.set_section_colors(color)
        self.w.setWindowTitle('hocRender: by Section Type')


    def paintSectionsByDensity(self, sectionColors, mechanism, excludeSections = []):
        """
        Color the sections in the reconstruction by the density of the selected mechanism
        in the section
        Inputs: sectionColors, a dictionary of section names and desired colors
                mechanism :should be from modelPars.mechnames['mech'], and is a list with
                    the mech name [0], and the conductance density variable name [1] -
                    for example: ['na', 'gnabar']
                excludeSections: a list of sections that should not be painted.
        """

        color = np.zeros((self.nsec, 4), dtype=float)
        gmax = 0.
        for stype in self.Sections:
            if stype in excludeSections: # everyone is in the synapse and axon, so skip
                continue
            for secno in self.Sections[stype]: # check out sections of a given type
###
### This routine directly references axon, and it should actually reference
### the section as indicated in the sections list (which SHOULD be the same).
### however, this will break the generality of function, as it requires that the primary
### list of sections be called "axon"
###
                ml = self.get_mechanisms(eval('h.axon[%d]' % (secno))) # this is too specific
                if mechanism[0] in ml:
                    gsec = self.get_density(eval('h.axon[%d]' % (secno)), mechanism)
                else:
                    gsec = 0.
                if gsec > gmax:
                    gmax = gsec
                color[secno, :] = self.Colors[sectionColors[stype]] # map colors
                color[secno, 3] = gsec # use alpha, but rescale next
        if gmax > 0:
            color[:,3] = 0.05 + 0.95*color[:,3]/gmax # set alpha for all sections
        else:
            color[:,3] = 0.05
        self.set_section_colors(color)
        self.w.setWindowTitle('hocRender: %s' % (mechanism[0]))


if __name__ == "__main__":
    global index, color, render, vdata, videoOpen, vfile, winsize, cycles, maxcycles
    app = pg.mkQApp()
    #pg.dbg()
    #h.load_file(1, 'Calyx-68cvt2.hoc')
    vfile = None
    h.load_file(1, "Calyx-S53Acvt3.hoc")
    render = hocRender(h)
    render.draw_model(modes=['blob'])
    sectionColors={'axon': 'r', 'heminode': 'g', 'stalk':'y', 'branch': 'b', 'neck': 'brown',
            'swelling': 'magenta', 'tip': 'powderblue', 'parentaxon': 'orange', 'synapse': 'k'}
    render.getSectionLists(sectionColors.keys())
    render.paintSectionsByType(sectionColors)
    #render.paintSectionsByDensity()
    QtGui.QApplication.instance().exec_()


    # get the time course data at every section...
    vdata = pickle.load(open('Canonical/Normal_swellings_14.02.04-18.24.23.p'))

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
    render_stack(0, vdata['data'].shape[1])
#    print dir(vfile)


    QtGui.QApplication.instance().exec_()