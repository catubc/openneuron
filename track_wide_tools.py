import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *


#LFP ANALYSIS TOOLBOX
class TrackWideTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TrackWideTools, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        self.root_dir = '/media/cat/All.Data.3TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/'

        row_index= 0
        #*****************************************************************
        self.low_cutoff = QLineEdit('0.1');                #parent.start_time = self.start_time
        self.low_cutoff.setMaximumWidth(50)
        low_cutoff_lbl = QLabel('low_pass (hz):', self)
        low_cutoff_lbl.setMaximumWidth(100)
      
        layout.addWidget(low_cutoff_lbl, row_index,0)
        layout.addWidget(self.low_cutoff, row_index,1)
        
        self.high_cutoff = QLineEdit('110.0');                  #parent.end_time = self.end_time
        self.high_cutoff.setMaximumWidth(50)
        high_cutoff_lbl = QLabel('high_pass (hz):', self)
        high_cutoff_lbl.setMaximumWidth(100)

        layout.addWidget(high_cutoff_lbl,row_index,2)
        layout.addWidget(self.high_cutoff,row_index,3)

        #% of electrodes
        self.n_electrodes = QLineEdit('1.0');                
        self.n_electrodes.setMaximumWidth(50)
        self.n_electrodes_lbl = QLabel('%Electrodes', self)
        self.n_electrodes_lbl.setMaximumWidth(100)
        layout.addWidget(self.n_electrodes_lbl, row_index, 4)
        layout.addWidget(self.n_electrodes, row_index, 5); row_index+=1
        
        
        self.button_tsf_to_lfp = QPushButton('Convert .tsf to low-pass @ 1Khz (apply filter)')
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
        
        self.compress_factor = QLineEdit('50');                  #parent.end_time = self.end_time
        self.compress_factor.setMaximumWidth(50)
        self.compress_factor_lbl = QLabel('Compress Factor (multiples of 25):', self)
        layout.addWidget(self.compress_factor_lbl, row_index, 1)
        layout.addWidget(self.compress_factor, row_index, 2)
        
        self.setLayout(layout)
    
    
    def tsftolfp(self):
        self.tsf_files = QtGui.QFileDialog.getOpenFileNames(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf_to_lfp(self)
    
    
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
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        compress_lfp(self) 
    