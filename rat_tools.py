import glob

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

        

class RatTools(QtGui.QWidget):
    def __init__(self, parent):
        super(RatTools, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()    
        
        #NB: ###########NEED TO MOVE THESE TO INDIVIDUALLY CALLABLE OR AUTO-CALLABLE BUTTONS
        #self.animal.rhd_to_tsf()
        #self.animal.tsf_to_lfp()
        #self.animal.lfp_compress()
        #self.animal.bin_to_npy()
        #self.animal.rhd_digital_save()
