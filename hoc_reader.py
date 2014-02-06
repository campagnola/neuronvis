import neuron


class HocReader:
    """
    Provides useful methods for reading hoc structures.
    
    Input:
        hoc: a hoc object or a "xxx.hoc" file name.
    """
    def __init__(self, hoc):
        if isinstance(hoc, basestring):
            hoc = neuron.load_file(1, hoc)
            
        self.h = hoc

        # we need to get the section
        self.nsec = len(list(self.h.allsec()))

        # geometry containers
        self.edges = None
        self.vertexes = None
        
        self.sections = {}   # all sections in the hoc  {sec_name: hoc Section}
        self.mechanisms = {} # {sec_name: [mechanism, ...]}
        self.sec_groups = {} # {sec_group_name: [Section, ...]}
        
        # populate self.sections and self.mechanisms
        self._read_section_info()
        

    def get_section(self, sec_name):
        """
        Return the hoc Section object with the given name.
        """
        return self.sections[sec_name]        

    def get_mechanisms(self, section):
        """
        Get a set of all of the mechanisms inserted into a given section
        Input: section (hoc object)
        Returns: set of mechanism names
        Side-effects: None
        """
        return self.mechanisms[section]

    def get_density(self, section, mechanism, variable):
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
                x =  eval('seg.%s' % mechanism)
                mecbar = '%s_%s' % (variable, mechanism)
                if mecbar in dir(x):
                    gmech.append(eval('x.%s' % variable))
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

    def _read_section_info(self):
        # Collect list of all sections and their mechanism names.
        self.sections = {}
        self.mechanisms = {}
        for sec in self.h:
            self.sections[sec.name] = sec
            mechs = set()
            for seg in sec:
                for mech in seg:
                    mechs.add(mech.name())
            self.mechanisms[sec.name] = mechs

    def add_section_group(self, name, sections):
        """
        Declare a grouping of sections (or section names). Sections may be
        grouped by any arbitrary criteria (cell, anatomical type, etc).
        
        Input:
            name: string name of the section group
            sections: list of section names or hoc Section objects.
        
        """
        assert name not in self.sec_groups
        group = []
        for sec in sections:
            group.append(self.get_section(sec))
        self.sec_groups[name] = group

    def get_section_group(self, name):
        """
        Return the list of hoc Section objects in the group *name*.
        """
        return self.sec_groups[name]
        
    def read_hoc_section_list(self, names):
        """
        Add a new section groups from the hoc variables indicated in *names*.
        
        Input:
            names: list of variable names. Each name must refer to a list of 
                   Sections in hoc. If a dict is supplied instead, then it
                   maps {hoc_list_name: section_group_name}.
        Side effects (modifies):
           calls add_section_group
        returns: Nothing.
        """
        # if a list is supplied, then the names of groups to create are 
        # exactly the same as the names of hoc lists.
        if not isinstance(names, dict):
            names = {name:name for name in names}
            
        for hoc_name, group_name in names.items():
            var = getattr(self.h, hoc_name)
            self.add_section_group(group_name, list(var))

    def get_geometry(self):
        """
        modified from:neuronvisio
        Generate structures that describe the geometry of the sections.
        Inputs: None
        Outputs: vertexes: record array containing {pos: (x,y,z), dia, sec_id}
                           for each segment.
                 edges:    array of pairs indicating the indexes of connected
                           vertexes.
        Side effects:
            modifies vertexes and edges.
        """
        
        # return cached geometry if this method has already run.
        if self.vertexes is not None:
            return self.vertexes, self.edges
        
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
            # Store the section. later.
            radius = sec.diam/2.
            sec_coords_bound = ((x_sec.min(), x_sec.max()),
                                (y_sec.min() - radius,
                                 y_sec.max() + radius),
                                (z_sec.min() - radius,
                                 z_sec.max() + radius))

            for i,xi in enumerate(x_sec):
                vertexes.append((x_sec[i], y_sec[i], z_sec[i], d_sec[i], secid))
                indx_geom_seg = len(vertexes) - 1
                if len(vertexes) > 1 and i > 0:
                    connections.append([indx_geom_seg, indx_geom_seg-1])


        self.edges  = np.array(connections)
        self.vertexes = np.array(vertexes)

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
        return (np.array(x),np.array(y),np.array(z),np.array(d))


    def make_volume_data(self):
        """
        Using the current state of vertexes, edges, generates a scalar field
        useful for building isosurface or volumetric renderings.
        Returns:
            * 3D scalar field indicating distance from nearest membrane,
            * 3D field indicating section IDs of nearest membrane,
            * QTransform that maps from 3D array indexes to original vertex 
                coordinates.
        """
        res = 0.4 # resolution of scalar field in microns
        maxdia = 10. # maximum diameter (defines shape of kernel)
        kernel_size = int(maxdia/res) + 1 # width of kernel
        
        vertexes, lines = self.get_geometry()
        
        # read vertex data
        verlocs = vertexes[:, :3]
        x = verlocs[:,0]
        y = verlocs[:,1]
        z = verlocs[:,2]
        d = vertexes[:,3]
        sec_id = vertexes[:,4]
        
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
