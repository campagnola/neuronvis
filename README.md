neuronvis
=========

Routines in Python for visualization of morphology files (".hoc" format) from NEURON. Part of pipeline for construction of biophysically and morphologically detailed models.

This requires: pyqtgraph (www.pyqtgraph.org)

Options are to color by section type (e.g., morphological region), by total conductance
of a channel (klt, for example) across sections, or to make a movie of vm as a function
of time and position over the structure, based on results from a prior simulation.

Note: there are limits regarding the structure of the hoc file that can be read. Some preprocessing of the files is recommened to keep the import clean. 

Original from Luke Campagnola's repository of the same name.
Some code based on neuron_visio project from http://michelemattioni.me/neuronvisio/ (but note: we use pyqtgraph, which requires PyQt and OpenGL, rather than the mayavi and traits tools).


