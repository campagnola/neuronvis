import os, pickle

class SimulationResult(object):
    """
    Class used to formalize structure and storage of results from NEURON 
    simulations.
    """
    version = 1
    
    def __init__(self, file_name=None):
        if file_name is not None:
            self.load(file_name)

    def load(self, file_name):
        """
        Load results from a .p file.        
        """
        self.file_name = file_name
        results = pickle.load(open(file_name, 'rb'))
        
        typ = results.get('_result_type', None)
        if typ is None:
            raise Exception('Unknown result type (no _result_type key)')
        elif typ != type(self).__name__:
            raise Exception('Unknown result type "%s"' % typ)
        
        ver = results['_result_version']
        if ver > self.version:
            print("Warning: results file is from newer version of %s: %d" % (type(self).__name__, ver))
        
        self.hoc_file = results['hoc_file']
        self.data = results['data']
        self.time = results['time']
        self.section_map = results['section_map']
        
    def get_hoc_file(self):
        """
        Return the full path to the hoc file for this result file. 
        hoc files with relative path names should be interpreted as being 
        relative to the path of the result file.
        """
        if self.hoc_file.startswith('/'):
            return self.hoc_file
        
        path = os.path.dirname(self.file_name)
        return os.path.join(path, self.hoc_file)
        
    
    def save(self, file_name):
        """
        Save simulation results to file. 
        This instance must have the following attributes in order to save:
        
            data: {'measurement': array(N,time), ...}
                  where 'measurement' might be 'Vm', etc..
            time: array of time values in simulation
            hoc_file: the hoc file that defines the geometry of the simulation
            section_map: list of section names in the same order that sections 
                         appear in the data arrays.
        """
        self.file_name = file_name
        
        path = os.path.dirname(self.file_name)
        hoc = os.path.relpath(self.hoc_file, path)
        
        fh = open(file_name, 'wb')
        pickle.dump({
            'hoc_file': self.hoc_file,
            'data': self.data,
            'time': self.time,
            'section_map': self.section_map,
            '_result_version': self.version,
            '_result_type': type(self).__name__,
            },
            fh)
