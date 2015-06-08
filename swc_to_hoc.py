# -*- coding: utf8 -*-
import numpy as np


class SWC(object):
    """Encapsulates a morphology tree as defined by the SWC standard.
    
    Parameters
    ----------
    filename : str or None
        The name of an swc file to load
    types : dict or None
        A dictionary mapping {type_id: type_name} that describes the type IDs
        in the swc data (second column).
    data : ndarray or None
        Optionally, a data array may be provided instead of an swc file. This
        is used internally.
    """
    def __init__(self, filename=None, types=None, data=None):
        self._dtype = [
            ('id', int), 
            ('type', int), 
            ('x', float), 
            ('y', float), 
            ('z', float), 
            ('r', float), 
            ('parent', int)
        ]
        
        self._id_lookup = None
        self._sections = None

        self.sectypes = {
            0: 'undefined',
            1: 'soma',
            2: 'axon',
            3: 'basal_dendrite',
            4: 'apical_dendrite',
        }
        if types is not None:
            self.sectypes.update(types)
        
        if data is not None:
            self.data = data
        elif filename is not None:
            self.load(filename)
        else:
            raise TypeError("Must initialize with filename or data array.")
        
        if np.any(self.data['parent'] >= self.data['id']):
            raise TypeError("SWC data is not sorted in topological order.")
    
    def load(self, filename):
        self.filename = filename
        self.data = np.loadtxt(filename, dtype=self._dtype)
        
    def copy(self):
        return SWC(self.data.copy())

    @property
    def lookup(self):
        """Return a dict that maps *id* to *index* in the data array.
        """
        if self._id_lookup is None:
            self._id_lookup = {}
            for i, rec in enumerate(self.data):
                self._id_lookup[rec['id']] = i
        return self._id_lookup

    def __getitem__(self, id):
        """Return record for node *id*.
        """
        return self.data[self.lookup[id]]

    def reparent(self, id):
        """Rearrange tree to make *id* the new root parent.
        """
        
    @property
    def sections(self):
        """Return lists of IDs grouped by topological section.
        
        The first item in each list connects to the last item in a previous
        list.
        """
        if self._sections is None:
            sections = []
            sec = []
            
            # find all nodes with nore than 1 child
            branchpts = set()
            endpoints = set(self.data['id'])
            endpoints.add(-1)
            seen = set()
            for r in self.data:
                p = r['parent']
                if p in seen:
                    branchpts.add(p)
                else:
                    seen.add(p)
                    endpoints.remove(p)
            
            # build lists of unbranched node chains
            for r in self.data:
                sec.append(r['id'])
                if r['id'] in branchpts or r['id'] in endpoints:
                    sections.append(sec)
                    sec = []
            
            self._sections = sections
            
        return self._sections
        
    def connect(self, parent_id, swc):
        """Combine this tree with another by attaching the root of *swc* as a 
        child of *parent_id*.
        """
        data = swc.data.copy()
        shift = self.data['id'].max() + 1 - data['id'].min()
        data['id'] += shift
        rootmask = data['parent'] == -1
        data['parent'] += shift
        data['parent'][rootmask] = parent_id
        
        self.data = np.concatenate([self.data, data])
        self.sort()
        
    def set_type(self, typ):
        self.data['type'] = typ
        
    def write_hoc(self, filename, types=None):
        """Write data to a HOC file.
        
        Each node type is written to a separate section list.
        """
        hoc = []
        
        sectypes = self.sectypes.copy()
        for t in np.unique(self.data['type']):
            if t not in sectypes:
                sectypes[t] = 'type_%d' % t

        # create section lists
        for t in sectypes.values():
            hoc.extend(['objref %s' % t,
                        '%s = new SectionList()' % t])
        hoc.append('')
            
        # create sections
        sects = self.sections
        hoc.append('create sections[%d]' % len(sects))
        sec_ids = {}
        for i, sec in enumerate(sects):
            # remember hoc index for this section
            endpt = self[sec[-1]]['id']
            sec_id = len(sec_ids)
            sec_ids[endpt] = sec_id
            
            # add section to list
            hoc.append('access sections[%d]' % sec_id)
            typ = self[sec[0]]['type']
            hoc.append('%s.append()' % sectypes[typ])
            
            # connect section to parent
            p = self[sec[0]]['parent']
            if p != -1:
                hoc.append('connect sections[%d](0), sections[%d](1)' % (sec_id, sec_ids[p]))

            # set up geometry for this section
            hoc.append('sections[%d] {' % sec_id)
            for seg in sec:
                rec = self[seg]
                hoc.append('  pt3dadd(%f, %f, %f, %f)' % (rec['x'], rec['y'], rec['z'], rec['r']))
            hoc.append('}')
            
            hoc.append('')
        
        open(filename, 'w').write('\n'.join(hoc))

    def sort(self):
        """Sort the tree in topological order.
        """
        recs = list(self.data)
        recs.sort(cmp=lambda a, b: 1 if b['id'] in self.path(a['id']) else -1)
        self.data = np.array(recs, dtype=self._dtype)
        
        self._id_lookup = None
        self._sections = None
        
    def path(self, node):
        path = [node]
        while True:
            node = self[node]['parent']
            if node < 0:
                return path
            path.append(node)

    def scale(self, x, y, z, r):
        self.data['x'] *= x
        self.data['y'] *= y
        self.data['z'] *= z
        self.data['r'] *= r
        
    def branch(self, id):
        """Return the branch beginning at *id*.
        """
    
    def topology(self):
        """Print the tree topology.
        """
        path = []
        indent = ''
        secparents = [self[s[0]]['parent'] for s in self.sections]
        
        for i, sec in enumerate(self.sections):
            p = secparents[i]
            if p != -1:
                ind = path.index(p)
                path = path[:ind+1]
                indent = indent[:(ind+1) * 3]
            path.append(self[sec[-1]]['id'])

            # look ahead to see whether subsequent sections are children
            if p in secparents[i+1:]:
                this_indent = indent[:-2] + u"├─ "
                indent =      indent[:-2] + u"│  │  "
            else:
                this_indent = indent[:-2] + u"└─ "
                indent =      indent[:-2] + u"   │  "
                
                
            typ = self.sectypes[self[sec[0]]['type']]
            if len(sec) > 10:
                secstr = "%s,...%s" % (str(tuple(sec[:3]))[:-1], str(tuple(sec[-3:]))[1:])
            else:
                secstr = str(tuple(sec))
            print "%ssections[%d] type=%s parent=%d %s" % (this_indent, i, typ, p, secstr)


if __name__ == '__main__':
    soma = SWC('data/cellbody.swc', types={1:'soma', 2:'axon', 3:'dendrite'})
    soma.set_type(1)
    axon = SWC('data/axonnonscaled.swc')
    axon.set_type(2)
    dend = SWC('data/dendnonscaled.swc')
    dend.set_type(3)
    
    soma.connect(10, dend)
    soma.connect(20, axon)
    
    soma.scale(0.11, 0.11, 0.06, 0.11)
    
    soma.write_hoc('test.hoc')
    soma.topology()
