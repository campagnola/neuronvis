import pyqtgraph as pg
import pyqtgraph.opengl as gl

from hoc_reader import HocReader
from hoc_graphics import *

class HocViewer(gl.GLViewWidget):
    """
    Subclass of GLViewWidget that displays data from HocReader.
    This is a convenience class implementing boilerplate display code.
    
    Input:
        h: HocReader instance or "xxxx.hoc" file name
    """
    def __init__(self, h):
        if isinstance(h, basestring):
            h = HocReader(h)
        self.h = h
        pg.mkQApp()  # make sure there is a QApplication before instantiating any QWidgets.
        super(HocRender, self).__init__()
        self.resize(720,720)
        self.show()
        self.setWindowTitle('hocRender')
        self.setCameraPosition(distance=300.)

        self.g = gl.GLGridItem()
        self.g.scale(2,2,1)
        self.addItem(self.g)
        
        self.graphics = []

    def show(self):
        """
        Start the Qt event loop.
        """
        QtGui.QApplication.exec_()

    def draw_volume(self):
        """
        Add a HocVolume graphic to this view.
        
        Returns:  HocVolume instance
        """
        g = HocVolume(self.h)
        self.graphics.append(g)
        self.addItem(g)
        return g

    def draw_surface(self):
        """
        Add a HocSurface graphic to this view.
        
        Returns:  HocSurface instance
        """
        g = HocSurface(self.h)
        self.graphics.append(g)
        self.addItem(g)
        return g

    def draw_volume(self):
        """
        Add a HocGraph graphic to this view.
        
        Returns:  HocGraph instance
        """
        g = HocGraph(self.h)
        self.graphics.append(g)
        self.addItem(g)
        return g

