__author__ = 'pbmanis'

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import scipy.ndimage
import neuron as h
import neuron as h
from neuron import *
import re

class hocRender():
    def __init__(self, h):
        self.h = h
        self.app = pg.mkQApp()
        self.w = gl.GLViewWidget()
        self.w.show()
        self.w.setWindowTitle('hoc Render')
        self.w.setCameraPosition(distance=200.)

        self.g = gl.GLGridItem()
        self.g.scale(2,2,1)
        self.w.addItem(self.g)
        self.sec2coords = {}
        # we need to get the section
        self.cyl2sec = {}

        # Needed to update the value of a cyl bound to a section
        self.sec2cyl = {}
        self.seg2id = {}
        self.sec2coords = {}
        self.connections = []
        self.n3dpoints_per_sec = {}


    def show(self):
        QtGui.QApplication.instance().exec_()


    def makeStick(self, seclist):
        for i, sec in enumerate(seclist): # self.CalyxStruct[input]['axon']: # set a new value

            if i == 0:
                info = self.get_sec_info(sec)
                print info
                # for nseg in sec:
                #     print sec()
                #     print dir(nseg)
                #     print nseg.diam
                #     print nseg.x
                #     print nseg.sec
                #     for sn in nseg.sec:
                #         print dir(sn)
                #         print sn.x
                #         print sn.diam


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

    def draw_model(self, modes=["cylinder"]):
        """
        modified from:neuronvisio
        Draw the model.
        Params:
        controls - the main gui obj."""

        # Draw the new one
        #self.h.define_shape()
        num_sections = 0
        
        # map segments (lines) to the section that contains them
        self.segment_to_section = {}
        
        vertexes = []
        connections = []
        for sec in self.h.allsec():
            m = re.match(r'axon\[(\d+)\]', sec.name())
            if m is None:
                raise Exception('could not determine section ID from name: %s' % sec.name())
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
                self.drawVolume(x, y, z, d, lines)
            else:
                self.drawMeshes(x, y, z, d, lines, mode)

        
    def makeVolumeData(self):
        
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

    def drawVolume(self, x, y, z, d, lines):
        scfield, idfield, transform = self.makeVolumeData(x, y, z, d, lines)
    
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
        colors = np.empty((len(self.mesh_sec_ids), 4))
        colors[:] = 0.3
        colors[sec_id] = 1
        self.set_section_colors(colors)
    
    def set_section_colors(self, sec_colors):
        md = self.mesh.opts['meshdata']
        colors = sec_colors[self.mesh_sec_ids]
        md.setVertexColors(colors)
        self.mesh.meshDataChanged()

    def drawMeshes(self, x, y, z, d, lines, mode):
        dmax = np.max(d)
        #wmax = 20.
        for c in range(len(lines)-3):
            i = lines[c, 0]
            j = lines[c, 1]
            # if i < 3:
            #     print 'xyzd: %6.1f %6.1f %6.1f %6.2f' % (x[i], y[i], z[i], d[i])
            pts = np.vstack([[x[i], x[j]],[y[i], y[j]],[z[i],z[j]]]).transpose()
            if mode == "line":
                plt = gl.GLLinePlotItem(pos=pts, width =(d[i]+d[j])/(2.), color=pg.glColor((int(255.*d[i]/dmax), 128)), connected=True)
                self.w.addItem(plt)
            elif mode == 'sphere':
                print 'I am a sphere'
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
#                m5.rotate(rz*180./np.pi, (x[j]-x[i]), (y[j]-y[i]), (z[j]-z[i]), local=False)
                m5.translate(x[i],y[i],z[i]+cyllen/2.0) # move into position

                self.w.addItem(m5)
        #self.draw_mayavi(x, y, z, d, self.edges)

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




if __name__ == "__main__":
    app = pg.mkQApp()
    pg.dbg()
    #h.load_file(1, 'Calyx-68cvt2.hoc')
    h.load_file(1, "Calyx-S53Acvt3.hoc")
    render = hocRender(h)
    render.draw_model(modes=['blob'])
    #render.show()
    import user

    import pickle
    vdata = pickle.load(open('data.p'))
    color = np.zeros((vdata['data'].shape[0], 4), dtype=float)
    index = 0
    #vdata['data'] = vdata['data'][:,400:550]
    
    def update():
        global index, color, render, vdata
        v = (vdata['data'][:,index] + 50) * 0.02
        color[:,0] = v
        color[:,1] = v+0.3
        color[:,2] = v-0.3
        color[:,3] = 0.2
        render.set_section_colors(color)
        
        index = (index + 1) % vdata['data'].shape[1]
    
    timer = pg.QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(0)
    
    def render_stack(start, stop):
        global index, render, timer
        timer.stop()
        for index in range(start, stop):
            update()
            render.w.readQImage().save('render-%03d.png' % (index-start))
            