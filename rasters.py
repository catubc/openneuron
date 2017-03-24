
import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class Rasters(QtGui.QWidget):
    def __init__(self, parent):
        super(Rasters, self).__init__(parent)
        self.parent = parent
        

        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        #self.parent.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/'
        #self.parent.root_dir = '/media/cat/2TB/in_vivo/nick/ptc21_tr5c/tsf_files/'
        #self.parent.root_dir = '/media/cat/All.Data.3TB/in_vivo/nick/ptc21/tr5c/'

        #***************************** Cat DATA *************************************
        #ptc21 tr5c 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/recordings/61-tr5c-blankscreen/61-tr5c-blankscreen_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/recordings/61-tr5c-blankscreen/61-tr5c-blankscreen_alltrack_lfp_50compressed.ptcs'


        #***************************** MOUSE DATA *************************************
        #2017-02_03 VISUAL
        self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_hp_butter_alltrack.ptcs'
        self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed_4.0threshold_3clusters.ptcs'

        #2017-01_31 BARREL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_lfp_100hz_alltrack_50compressed.ptcs'

        #2017_01_30 AUDITORY - 2 clusters only !?
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_30_auditory_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170130_164612_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_30_auditory_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170130_164612_lfp_250hz_alltrack_50compressed.ptcs'

        #2016_07_26 AUDITORY   - Multiople clusters
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_lfp_250hz_alltrack_50compressed.ptcs'


        #********************************* CHRONIC V1 ************************************
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack.tsf'
        
        self.selected_sort_sua = self.parent.sua_file



        layout = QtGui.QGridLayout()

        row_index=0

        #**************************************************************************************
        #*********************************** LOAD FILES  **************************************
        #**************************************************************************************
        self.preprocess_lbl = QLabel('LOAD FILES', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1


        #Set root dir; button, then label
        self.button_set_sua_file = QPushButton('Single Unit File')
        self.button_set_sua_file.setMaximumWidth(200)
        self.button_set_sua_file.clicked.connect(self.set_sua_file)
        layout.addWidget(self.button_set_sua_file, row_index, 0)
        
        #self.parent.button_set_sua_file  = os.getcwd()
        self.button_set_sua_file_lbl = QLabel(os.path.split(self.parent.sua_file)[1], self)
        layout.addWidget(self.button_set_sua_file_lbl, row_index, 1); row_index+=1

        #Set LFP event file
        self.button_set_lfp_event_file = QPushButton('LPF Event File')
        self.button_set_lfp_event_file.setMaximumWidth(200)
        self.button_set_lfp_event_file.clicked.connect(self.set_lfp_event_file)
        layout.addWidget(self.button_set_lfp_event_file, row_index, 0)
        
        #self.parent.set_lfp_event_file = os.getcwd()
        self.button_set_lfp_event_file_lbl = QLabel(os.path.split(self.parent.lfp_event_file)[1], self)
        layout.addWidget(self.button_set_lfp_event_file_lbl, row_index, 1); row_index+=1

        #***************************************************************************
        self.preprocess_lbl = QLabel('RASTERS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        
        self.button_rasters = QPushButton('Plot Rasters')
        self.button_rasters.setMaximumWidth(200)
        self.button_rasters.clicked.connect(self.plt_rasters)
        layout.addWidget(self.button_rasters, row_index, 0)
        
        #***********************************************************************************************************
        #***********************************************************************************************************

        self.setLayout(layout)

    
    #**************************** LOAD FILE ROUTINES *******************************
    def set_sua_file(self):
        self.parent.sua_file = QtGui.QFileDialog.getOpenFileName(self, "Select SUA file (*.ptcs)", self.parent.root_dir,"PTCS (*.ptcs)")
        self.button_set_sua_file_lbl.setText(self.parent.sua_file.replace(os.path.dirname(self.parent.sua_file),''))
        #self.parent.setWindowTitle(self.parent.sua_file)
        self.selected_sort_sua = self.parent.sua_file
        

    def set_lfp_event_file(self):
        self.parent.lfp_event_file = QtGui.QFileDialog.getOpenFileName(self, "Select LFP event file (*.ptcs)", self.parent.sua_file,"PTCS (*.ptcs)")
        self.button_set_lfp_event_file_lbl.setText(self.parent.lfp_event_file.replace(os.path.dirname(self.parent.lfp_event_file),''))
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))
    
    
    def set_recording(self):
        temp_dir = '/media/cat/2TB/in_vivo/nick/ptc21/'

        self.parent.recording_dir = QtGui.QFileDialog.getExistingDirectory(self, "Select recording", temp_dir)
        
        self.rec_path = os.path.dirname(self.parent.recording_dir)
        self.rec_name = self.parent.recording_dir.replace(self.rec_path,'')
        self.button_set_recording_lbl.setText(self.rec_name)
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))

    def plt_rasters(self):
        
        Plot_rasters(self)
