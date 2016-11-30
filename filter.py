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
        
        parent.start_time = QLineEdit('0.5');                #parent.start_time = self.start_time
        parent.start_time.setMaximumWidth(50)
        start_time_lbl = QLabel('low_pass (hz):', self)
        start_time_lbl.setMaximumWidth(100)

        parent.end_time = QLineEdit('6.0');                  #parent.end_time = self.end_time
        parent.end_time.setMaximumWidth(50)
        end_time_lbl = QLabel('high_pass (hz):', self)
        end_time_lbl.setMaximumWidth(100)
      
        layout.addWidget(start_time_lbl, 0,0)
        layout.addWidget(parent.start_time, 0,1)
        layout.addWidget(end_time_lbl,0,2)
        layout.addWidget(parent.end_time,0,3)

        #MAKE BUTTONS             
        self.button1 = QPushButton('Butterworth - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 5, 0)

        self.button1 = QPushButton('Wavelet - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 5, 2)

        self.button1 = QPushButton('Butterworth - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 6, 0)

        self.button1 = QPushButton('Cheby - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 6, 2)


        self.setLayout(layout)
   
    def filter_imaging(self):
        print "...pop up window to batch select invidiual files..."

