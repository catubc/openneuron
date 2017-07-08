import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class LFP(QtGui.QWidget):
    def __init__(self, parent):
        print "... initializing LFP tools..."
        super(LFP, self).__init__(parent)
        layout = QtGui.QGridLayout()
        
        self.parent = parent
        
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        #self.parent.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering/'
        #self.selected_recording = self.parent.root_dir
        #self.selected_sort_sua = self.parent.root_dir
        #self.selected_sort_lfp = self.parent.root_dir
        
        #CAT VISUAL
        #self.selected_recording = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/recordings/61-tr5c-blankscreen/61-tr5c-blankscreen_alltrack_lfp.tsf'
        #self.selected_sort_sua = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/recordings/61-tr5c-blankscreen/61-tr5c-blankscreen_alltrack.ptcs'
        #self.selected_sort_lfp = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/recordings/61-tr5c-blankscreen/61-tr5c-blankscreen_alltrack_lfp_50compressed.ptcs'
        
        #MOUSE AUDITORY
        self.selected_recording = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_lfp_250hz_alltrack.tsf'
        self.selected_sort_sua = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_hp_butter_alltrack.ptcs'
        self.selected_sort_lfp = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_lfp_250hz_alltrack_50compressed.ptcs'

        #MOUSE VISUAL
        #self.selected_recording = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack.tsf'
        #self.selected_sort_sua = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_hp_butter_alltrack.ptcs'
        #self.selected_sort_lfp = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed_4.0threshold_3clusters.ptcs'

        #MOUSE BARREL
        #self.selected_recording = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_lfp_100hz_alltrack.tsf'
        #self.selected_sort_sua = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_hp_butter_alltrack.ptcs'
        #self.selected_sort_lfp = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_lfp_100hz_alltrack_50compressed.ptcs'

        #self.selected_recording = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_lfp_250hz_alltrack.tsf'
        #self.selected_sort_sua = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_hp_butter_alltrack.ptcs'
        #self.selected_sort_lfp = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_lfp_250hz_alltrack_50compressed.ptcs'


        #********************************* CHRONIC V1 ************************************
        #self.selected_sort_sua = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_hp_butter_alltrack.ptcs'
        #self.selected_sort_lfp = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack_notch_50compressed.ptcs'
        #self.selected_recording = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack.tsf'
        
       



        row_index = 0
        #***********************************************************************************
        #******************************* SELECT RECORDING / SORT ***************************
        #***********************************************************************************
        self.preprocess_lbl = QLabel('LOAD FILES', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        
        #Select recording
        self.button_select_recording = QPushButton('Recording File')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        self.select_recording_lbl = QLabel(os.path.split(self.selected_recording)[1], self)
        layout.addWidget(self.select_recording_lbl, row_index,1); row_index+=1

        #Select sort
        self.button_select_lfpsort = QPushButton('LFP Sort File')
        self.button_select_lfpsort.setMaximumWidth(200)
        self.button_select_lfpsort.clicked.connect(self.slct_sort_lfp)
        layout.addWidget(self.button_select_lfpsort, row_index, 0)

        self.select_sort_lfp_lbl = QLabel(os.path.split(self.selected_sort_lfp)[1], self)
        layout.addWidget(self.select_sort_lfp_lbl, row_index,1); row_index+=1

        #Select sort
        self.button_select_suasort = QPushButton('SUA Sort File')
        self.button_select_suasort.setMaximumWidth(200)
        self.button_select_suasort.clicked.connect(self.slct_sort_sua)
        layout.addWidget(self.button_select_suasort, row_index, 0)
                
        self.selected_sort_sua_lbl = QLabel(os.path.split(self.selected_sort_sua)[1], self)
        layout.addWidget(self.selected_sort_sua_lbl, row_index,1); row_index+=1
        
        #***********************************************************************************
        #************************** VIEW SPECGRAM AND RASTERS ******************************
        #***********************************************************************************
        self.preprocess_lbl = QLabel('SPECGRAMS AND RASTERS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1


        self.specgram_ch = QLineEdit('63');                #parent.start_time = self.start_time
        self.specgram_ch.setMaximumWidth(50)
        specgram_ch_lbl = QLabel('specgram_ch:', self)
        specgram_ch_lbl.setMaximumWidth(100)
        layout.addWidget(specgram_ch_lbl, row_index,0)
        layout.addWidget(self.specgram_ch, row_index,1)

        self.time_start = QLineEdit('0');                #parent.start_time = self.start_time
        self.time_start.setMaximumWidth(50)
        time_start_lbl = QLabel('t_start:', self)
        time_start_lbl.setMaximumWidth(100)
        layout.addWidget(time_start_lbl, row_index,2)
        layout.addWidget(self.time_start, row_index,3)
        
        self.time_end = QLineEdit('10');                #parent.start_time = self.start_time
        self.time_end.setMaximumWidth(50)
        time_end_lbl = QLabel('t_end:', self)
        time_end_lbl.setMaximumWidth(100)
        layout.addWidget(time_end_lbl, row_index,4)
        layout.addWidget(self.time_end, row_index,5)
              
        self.specgram_db_clip = QLineEdit('-60');                #parent.start_time = self.start_time
        self.specgram_db_clip.setMaximumWidth(50)
        specgram_db_clip_lbl = QLabel('db_clip:', self)
        specgram_db_clip_lbl.setMaximumWidth(100)
        layout.addWidget(specgram_db_clip_lbl, row_index,6)
        layout.addWidget(self.specgram_db_clip, row_index,7); row_index+=1  

        #****************************************************************************************************************
        #************************************************** FUNCTIONS ***************************************************
        #****************************************************************************************************************

        self.button10 = QPushButton('View Spectrum')
        self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.button10, row_index, 0)
        self.button10.clicked.connect(self.view_Spectrum); row_index+=1  
        
        self.button10 = QPushButton('View Spegram')
        self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.button10, row_index, 0)
        self.button10.clicked.connect(self.view_Specgram); row_index+=1  
        
        
        self.button3 = QPushButton('Specgram + LFP Rasters')
        self.button3.setMaximumWidth(200)
        layout.addWidget(self.button3, row_index, 0)
        self.button3.clicked.connect(self.view_lfp_raster); row_index+=1
 
        
        self.tfr_spec = QPushButton('TFR Spegram + Rasters')
        self.tfr_spec.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.tfr_spec, row_index, 0)
        self.tfr_spec.clicked.connect(self.view_Specgram_tfr); row_index +=1
        
        
        self.all_rasters = QPushButton('SUA and LFP Rasters')
        layout.addWidget(self.all_rasters, row_index, 0)
        self.all_rasters.clicked.connect(self.view_allrasters)    
        
        
        
        self.setLayout(layout)


    def slct_recording(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.selected_recording,"TSF (*.tsf)") 
        path_name, file_name = os.path.split(self.selected_recording)
        self.select_recording_lbl.setText(file_name)
        self.recName = self.selected_recording
        #self.tsf = Tsf_file(self.selected_recording)
        #self.tsf.read_ec_traces()

  
    def slct_sort_lfp(self):
        self.selected_sort_lfp =  QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.selected_recording,"PTCS (*.ptcs)")
        path_name, file_name = os.path.split(self.selected_sort_lfp)
        self.select_sort_lfp_lbl.setText(file_name)

        
    def slct_sort_sua(self):
        self.selected_sort_sua =  QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.selected_recording,"PTCS (*.ptcs)")
        path_name, file_name = os.path.split(self.selected_sort_sua)
        self.selected_sort_sua_lbl.setText(file_name)
        
        
    #***************************************************************************************************
    #************************************** DEFINE FUNCTIONS *******************************************
    #***************************************************************************************************

    def view_Spectrum(self):
        #Specgram_syncindex(self, 0) #Use regular specgram
        Spectrum(self) #Use regular specgram
        

    def view_Specgram(self):
        #Specgram_syncindex(self, 0) #Use regular specgram
        Specgram_syncindex(self) #Use regular specgram
        
        
    def view_lfp_raster(self):
        #self.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        self.ptcsName = self.selected_sort_lfp
        plot_rasters(self)

        
    def view_Specgram_tfr(self):
        Specgram_syncindex_tfr(self)     #Use TFR specgram
            

    def view_allrasters(self):
        plot_all_rasters(self)


