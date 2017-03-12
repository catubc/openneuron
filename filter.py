
import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class Filter(QtGui.QWidget):
    def __init__(self, parent):
        super(Filter, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        row_index = 0
        
        self.low_cutoff = QLineEdit('0.1');                #parent.start_time = self.start_time
        self.low_cutoff .setMaximumWidth(50)
        self.lowcutoff_lbl = QLabel('Low cutoff (Hz)', self)
        self.lowcutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.lowcutoff_lbl, row_index, 0)
        layout.addWidget(self.low_cutoff, row_index, 1)



        self.high_cutoff = QLineEdit('110');                #parent.start_time = self.start_time
        self.high_cutoff .setMaximumWidth(50)
        self.high_cutoff_lbl = QLabel('High cutoff (Hz)', self)
        self.high_cutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.high_cutoff_lbl, row_index, 2)
        layout.addWidget(self.high_cutoff, row_index, 3); row_index+=1



        #MAKE BUTTONS             
        self.button1 = QPushButton('Butterworth - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.fltr_ephys)
        layout.addWidget(self.button1, 5, 0)

        self.button1 = QPushButton('Wavelet - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.fltr_ephys)
        layout.addWidget(self.button1, 5, 2)

        self.button1 = QPushButton('Butterworth - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.fltr_ephys)
        layout.addWidget(self.button1, 6, 0)

        self.button1 = QPushButton('Cheby - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.fltr_ephys)
        layout.addWidget(self.button1, 6, 2)


        self.setLayout(layout)
   
    def fltr_ephys(self):
        
        self.tsf_filename = QtGui.QFileDialog.getOpenFileName(self, "Select tsf file (*.tsf)", self.root_dir,"TSF (*.tsf)")
        #self.button_set_sua_file_lbl.setText(self.parent.sua_file.replace(os.path.dirname(self.parent.sua_file),''))
        #self.parent.setWindowTitle(self.parent.sua_file)

        filter_ephys(self)
        
