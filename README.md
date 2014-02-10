neuronvis
=========

Routines in Python for visualization of morphology files (".hoc" format) from NEURON.
Requires pyqtgraph (www.pyqtgraph.org)

Options are to color by section type (e.g., morphological region), by total conductance
of a channel (klt, for example) across sections, or to make a movie of vm as a function
of time and position over the structure, based on results from a prior simulation.

Note: there are limits regarding the structure of the hoc file that can be read. Some preprocessing of the files is recommened to keep things clean

Original from Luke Campagnola's repository of the same name.
Some code based on neuron_visio project from http://michelemattioni.me/neuronvisio/.


