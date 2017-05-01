import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *


#LFP ANALYSIS TOOLBOX
class IntanTools(QtGui.QWidget):
    def __init__(self, parent):
        super(IntanTools, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        
        row_index=0
        
        
        #*************************************************************
        #*********************** CONVERT TO TSF   ********************
        #*************************************************************
        self.button_ephys_to_tsf = QPushButton('.rhd -> .tsf (+ dig chs)')
        self.button_ephys_to_tsf.setMaximumWidth(250)
        self.button_ephys_to_tsf.clicked.connect(self.rhdtotsf)
        layout.addWidget(self.button_ephys_to_tsf, row_index, 0); row_index+=1

        self.button_tsf_to_digchs = QPushButton('.tsf -> extract digital chs')
        self.button_tsf_to_digchs.setMaximumWidth(250)
        self.button_tsf_to_digchs.clicked.connect(self.tsf_extract_digchannels)
        layout.addWidget(self.button_tsf_to_digchs, row_index, 0); row_index+=1


        #*************************************************************
        #********************* CONVERT DIGITAL CHS *******************
        #*************************************************************
        self.button_digital_convert = QPushButton('Save external channels only')
        self.button_digital_convert.setMaximumWidth(250)
        self.button_digital_convert.clicked.connect(self.rhd_digital_convert)
        layout.addWidget(self.button_digital_convert, row_index, 0); row_index+=1
        
        
        self.setLayout(layout)


    def rhdtotsf(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileNames(self, "RHD (*.rhd)", self.root_dir,"RHD (*.rhd)") 
        rhd_to_tsf(self.selected_recording)   #Send as list in case we revert back to multiple recordings at once need to rmemeber this


    def tsf_extract_digchannels(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileNames(self, "tsf (*.tsf)", self.root_dir,"tsf (*.tsf)")[0]
        tsf_to_digchs(self)   #Send as list in case we revert back to multiple recordings at once need to rmemeber this
        
        
    def rhd_digital_convert(self):
        rhd_digital_save(self)
        
        
        
        
