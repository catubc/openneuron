import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *


#LFP ANALYSIS TOOLBOX
class Track(QtGui.QWidget):
    def __init__(self, parent):
        super(Track, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        self.root_dir = '/media/cat/'

        row_index= 0
        #*****************************************************************
        parent.start_time = QLineEdit('0.5');                #parent.start_time = self.start_time
        parent.start_time.setMaximumWidth(50)
        start_time_lbl = QLabel('low_pass (hz):', self)
        start_time_lbl.setMaximumWidth(100)

        parent.end_time = QLineEdit('6.0');                  #parent.end_time = self.end_time
        parent.end_time.setMaximumWidth(50)
        end_time_lbl = QLabel('high_pass (hz):', self)
        end_time_lbl.setMaximumWidth(100)
      
        layout.addWidget(start_time_lbl, row_index,0)
        layout.addWidget(parent.start_time, row_index,1)
        layout.addWidget(end_time_lbl,row_index,2)
        layout.addWidget(parent.end_time,row_index,3); row_index+=1


        self.button_tsf_to_lfp = QPushButton('Convert .raw .tsf to low-pass (1khz) .tsf')
        self.button_tsf_to_lfp.setMaximumWidth(350)
        self.button_tsf_to_lfp.clicked.connect(self.tsftolfp)
        layout.addWidget(self.button_tsf_to_lfp, row_index, 0); row_index+=1
        

        self.button_lfpzip_to_lptsf = QPushButton('Convert .lpf.zip to _lp.tsf')
        self.button_lfpzip_to_lptsf.setMaximumWidth(250)
        self.button_lfpzip_to_lptsf.clicked.connect(self.lfpzip_to_lptsf)
        layout.addWidget(self.button_lfpzip_to_lptsf, row_index, 0); row_index+=1


        self.button1 = QPushButton('Concatenate multiple .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_tsf)
        layout.addWidget(self.button1, row_index, 0); row_index+=1


        self.button1 = QPushButton('Concatenate .lfp.zip files')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_lfp_zip)
        layout.addWidget(self.button1, row_index, 0); row_index+=1


        self.button1 = QPushButton('Time compress (1Khz) .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.comp_lfp)
        layout.addWidget(self.button1, row_index, 0)
        
        self.parent.compress_factor = QLineEdit('50');                  #parent.end_time = self.end_time
        self.parent.compress_factor.setMaximumWidth(50)
        self.parent.compress_factor_lbl = QLabel('Compress Factor (multiples of 25):', self)
        layout.addWidget(self.parent.compress_factor_lbl, row_index, 1)
        layout.addWidget(self.parent.compress_factor, row_index, 2)
        
        self.setLayout(layout)

    
    def tsftolfp(self):
        tsf_files = QtGui.QFileDialog.getOpenFileNames(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf_to_lfp(tsf_files)
    

    def lfpzip_to_lptsf(self):
        lfpzip_file = QtGui.QFileDialog.getOpenFileName(self, "LFP ZIP (*.lfp.zip)", self.root_dir,"*.lfp.zip (*.lfp.zip)")
        lfp_to_lptsf(lfpzip_file)
        

    def multi_tsf(self):
        print "... selecting multiple recordings ..."
        out_files = QtGui.QFileDialog.getOpenFileNames(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        self.tsf_files = out_files
        
        concatenate_tsf(self) 


    def multi_lfp_zip(self):
        print "... selecting multiple .lfp.zip recording directories ..."
        dialog = FileDialog()   #**** SELECT MULTIPLE DIRECTORIES, NOT INDIVIDUAL FIELS
        dialog.exec_()
        
        self.parent.animal.tsf_files = dialog.out_files
        concatenate_lfp_zip(self) 
            

    def comp_lfp(self):
        print "... selecting single lfp recording ..."
        self.parent.animal.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        compress_lfp(self) 

