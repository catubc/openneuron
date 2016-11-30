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
        self.root_dir = '/media/cat/12TB/in_vivo/tim/'
        
        row_index=0
        
        #*************************************************************
        #*********************** SELECT RECORDING ********************
        #*************************************************************
        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording2)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.selected_recording  = '/media/cat/12TB/in_vivo/tim/cat'
        self.select_recording_lbl = QLabel(self.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1); row_index+=1
        
        
        #*************************************************************
        #*********************** CONVERT TO TSF   ********************
        #*************************************************************
        self.button_rhd_to_tsf = QPushButton('Convert .rhd to .tsf')
        self.button_rhd_to_tsf.setMaximumWidth(250)
        self.button_rhd_to_tsf.clicked.connect(self.rhdtotsf)
        layout.addWidget(self.button_rhd_to_tsf, row_index, 0); row_index+=1


        #*************************************************************
        #********************* CONVERT DIGITAL CHS *******************
        #*************************************************************
        self.button_digital_convert = QPushButton('Save external channels only')
        self.button_digital_convert.setMaximumWidth(250)
        self.button_digital_convert.clicked.connect(self.rhd_digital_convert)
        layout.addWidget(self.button_digital_convert, row_index, 0); row_index+=1
        
        
        self.setLayout(layout)

    def slct_recording2(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileNames(self, "RHD (*.rhd)", self.root_dir,"RHD (*.rhd)") 
        #path_name, file_name = os.path.split(self.selected_recording)
        self.select_recording_lbl.setText(self.selected_recording[0])
        

    def rhdtotsf(self):
        rhd_to_tsf(self.selected_recording)   #Send as list in case we revert back to multiple recordings at once need to rmemeber this
    

    def rhd_digital_convert(self):
        rhd_digital_save(self.selected_recording)
        
