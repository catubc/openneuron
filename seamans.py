
import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class Seamans(QtGui.QWidget):
    def __init__(self, parent):
        super(Seamans, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        
        row_index = 0
        #*******************************************************************************************
        #*******************************************************************************************
        
        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ncs to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ncs_convert)
        layout.addWidget(self.button1, row_index, 0)
        
        self.parent.make_hp = QLineEdit('False');                #parent.start_time = self.start_time
        self.parent.make_hp.setMaximumWidth(50)
        parent.make_hp_lbl = QLabel('Make High Pass .tsf', self)
        self.parent.make_hp.setMaximumWidth(100)
        layout.addWidget(self.parent.make_hp_lbl, row_index,1)
        layout.addWidget(self.parent.make_hp, row_index,2)
        
        
        row_index+=1

        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ntt to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ntt_convert)
        layout.addWidget(self.button1, row_index, 0);row_index+=1



        self.subsample_tsf = QPushButton('Subsample .tsf file')
        self.subsample_tsf.setMaximumWidth(250)
        self.subsample_tsf.clicked.connect(self.tsf_subs)
        layout.addWidget(self.subsample_tsf, row_index, 0)
        
        
        #self.button1 = QPushButton('Concatenate .lfp.zip files')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.multi_lfp_zip)
        #layout.addWidget(self.button1, 6, 0)


        #self.button1 = QPushButton('Compress .tsf files')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.comp_lfp)
        #layout.addWidget(self.button1, 7, 0)
        
        
        #self.button1 = QPushButton('Make Track Wide HighPass .tsf')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.concatenate_tsf)
        #layout.addWidget(self.button1, 7, 0)

        

        self.setLayout(layout)

        

    def ncs_convert(self):
        self.parent.root_dir = '/media/cat/12TB/in_vivo/jeremy/' 
        ncs_to_tsf(self, QtGui.QFileDialog.getOpenFileNames(self.parent, 'Load Experiment', self.parent.root_dir, "*.ncs"))

    def ntt_convert(self):
        self.parent.root_dir = '/media/cat/12TB/in_vivo/jeremy/' 
        ntt_to_tsf(self, QtGui.QFileDialog.getOpenFileNames(self.parent, 'Load Experiment', self.parent.root_dir, "*.ntt"))

    def tsf_subs(self):
        
        self.root_dir = '/media/cat/12TB/in_vivo/jeremy/' 
        
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)") 

        tsf_subsample(self)
        
        
        
        
        
