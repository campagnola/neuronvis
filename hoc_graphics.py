import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
import scipy.ndimage

Colors = { # colormap
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




class HocGraphic:
    """
    Methods common to all Hoc graphical representation classes (HocVolume, 
    HocSurface, etc.)
    
    """
    
    def set_section_colors(self, colors):
        """
        Recolor the graphic by section using the *colors* array. This method
        must be reimplemented by HocGraphic subclasses. The order of elements
        in the array must match the order of sections defined in 
        HocReader.sections.
        """
        raise NotImplementedError()
    
    def set_group_colors(self, colors, default_color=(0,0,0,0), alpha=None):
        """
        Color the sections in the reconstruction according to their
        group name.
        Inputs: 
            colors: a dictionary of section group names and their associated colors
            default_color: color to use for any sections that are not included in the
                           groups listed in *colors*.
            alpha: If specified, this overrides the alpha value for all group colors.
        Side-effects: none.
        """
        sec_colors = np.zeros((len(self.h.sections), 4), dtype=float)
        sec_colors[:] = default_color
        
        for group_name, color in colors.items():
            for sec_name in self.h.get_section_group(group_name):
                if isinstance(color, basestring):
                    color = Colors[color]
                index = self.h.sec_index[sec_name]
                sec_colors[index] = color
                if alpha is not None:
                    sec_colors[index, 3] = alpha
        self.set_section_colors(sec_colors)

    def paint_sections_by_density(self, sectionColors, mechanism, excludeSections = []):
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
                color[secno, :] = Colors[sectionColors[stype]] # map colors
                color[secno, 3] = gsec # use alpha, but rescale next
        if gmax > 0:
            color[:,3] = 0.05 + 0.95*color[:,3]/gmax # set alpha for all sections
        else:
            color[:,3] = 0.05
        self.set_section_colors(color)
        #self.w.setWindowTitle('hocRender: %s' % (mechanism[0]))


class HocVolume(gl.GLVolumeItem, HocGraphic):
    """
    Subclass of GLVolumeItem that draws a volume representation of the geometry
    specified by a HocReader.
    
    Input:
        h: HocReader instance
    """
    def __init__(self, h):
        self.h = h
        scfield, idfield, transform = self.h.make_volume_data()
        nfdata = np.empty(scfield.shape + (4,), dtype=np.ubyte)
        nfdata[...,0] = 255 #scfield*50
        nfdata[...,1] = 255# scfield*50
        nfdata[...,2] = 255# scfield*50
        nfdata[...,3] = np.clip(scfield*150, 0, 255)
        super(HocVolume, self).__init__(nfdata)
        self.setTransform(transform)


class HocSurface(gl.GLMeshItem, HocGraphic):
    """
    Subclass of GLMeshItem that draws a surface representation of the geometry
    specified by a HocReader. The surface is generated as an isosurface of
    the scalar field generated by measuring the minimum distance from all
    cell membranes across the volume of the cell.
    
    Input:
        h: HocReader instance
    """
    def __init__(self, h):
        self.h = h
        scfield, idfield, transform = self.h.make_volume_data()
        #scfield = scipy.ndimage.gaussian_filter(scfield, (0.5, 0.5, 0.5))
        #pg.image(scfield)
        verts, faces = pg.isosurface(scfield, level=0.0)
        self.verts = verts
        self.faces = faces
        vertexColors = np.empty((verts.shape[0], 4), dtype=float)
        md = gl.MeshData(vertexes=verts, faces=faces)
        
        # match vertexes to section IDs
        vox_locations = verts.astype(int)
        # get sction IDs for each vertex
        self.vertex_sec_ids = idfield[vox_locations[:,0], vox_locations[:,1], vox_locations[:,2]] 
        
        super(HocSurface, self).__init__(meshdata=md, smooth=True, shader='balloon')
        self.setTransform(transform)
        self.setGLOptions('additive')

    def show_section(self, sec_id, color=(1, 1, 1, 1), bg_color=(1, 1, 1, 0)):
        """
        Set the color of section named *sec_id* to *color*.
        All other sections are colored with *bg_color*.
        """
        colors = np.empty((len(self.vertex_sec_ids), 4))
        colors[:] = bg_color
        colors[sec_id] = color
        self.set_section_colors(colors)
    
    def set_section_colors(self, sec_colors):
        """
        Set the colors of multiple sections.
        
        Input:
            colors: (N,4) float array of (r,g,b,a) colors, in the order that
                    sections are defined by the HocReader.
        """
        colors = sec_colors[self.vertex_sec_ids]
        print colors
        self.opts['meshdata'].setVertexColors(colors)
        self.meshDataChanged()



class HocGraph(gl.GLLinePlotItem, HocGraphic):
    """
    Subclass of GLLinePlotItem that draws a line representation of the geometry
    specified by a HocReader.
    
    Input:
        h: HocReader instance
    """
    def __init__(self, h):
        self.h = h
        
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
            plt = gl.GLLinePlotItem(pos=pts, width =(d[i]+d[j])/(2.), color=pg.glColor((int(255.*d[i]/dmax), 128)), connected=True)
            self.w.addItem(plt)


        
    #def drawMeshes(self,  mode):
        #"""
        #Draw remaining mesh figures
        #mode here can be line, sphere or cylinder

        #"""
        #verlocs = self.vertexes[:, :3]
        #x = verlocs[:,0]
        #y = verlocs[:,1]
        #z = verlocs[:,2]
        #d = self.vertexes[:,3]
        #dmax = np.max(d)
        #lines = np.vstack(self.edges)
        #for c in range(len(lines)-3):
            #i = lines[c, 0]
            #j = lines[c, 1]
            #pts = np.vstack([[x[i], x[j]],[y[i], y[j]],[z[i],z[j]]]).transpose()
            #if mode == "line":
                #plt = gl.GLLinePlotItem(pos=pts, width =(d[i]+d[j])/(2.), color=pg.glColor((int(255.*d[i]/dmax), 128)), connected=True)
                #self.w.addItem(plt)
            #elif mode == 'sphere':
                #md = gl.MeshData.sphere(rows=10, cols=20, radius=d[i]/2.0) # , length=d(i))
                #colors = np.ones((md.faceCount(), 4), dtype=float)
                #colors[::2,0] = 0
                #colors[:,1] = np.linspace(0, 1, colors.shape[0])
                #md.setFaceColors(colors)
                #m5 = gl.GLMeshItem(meshdata=md, smooth=False, drawEdges=False)
                #m5.translate(x[i],y[i],z[i])
                #self.w.addItem(m5)
            #elif mode == "cylinder":
                #cyllen = np.sqrt((x[j]-x[i])**2.0 + (y[j]-y[i])**2.0 + (z[j]-z[i])**2.0)
                #md = gl.MeshData.cylinder(rows=1, cols=8, radius=[d[i]/2., d[j]/2.], length=cyllen)
                #colors = np.ones((md.faceCount(), 4), dtype=float)
                #colors[::2,0] = 0
                #colors[:,1] = np.linspace(0, 1, colors.shape[0])
                ##colors[:,1] = (int(255.*d[i]/dmax), 128)
                #md.setFaceColors(colors)
                #m5 = gl.GLMeshItem(meshdata=md, smooth=True, drawEdges=False)
                #p2 = pg.Vector(x[j], y[j], z[j])
                #p1 = pg.Vector(x[i], y[i], z[i])
                #r = pg.Vector(0,0,1)
                #axis = pg.QtGui.QVector3D.crossProduct(r, p2-p1)
                #ang = r.angle(p2-p1)
                #m5.rotate(ang, axis.x(), axis.y(), axis.z())
                #m5.translate(x[i],y[i],z[i]+cyllen/2.0) # move into position
                #self.w.addItem(m5)

