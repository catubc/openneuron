import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class TracesTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TracesTools, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        #self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        self.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering/'
        #self.parent.root_dir = '/media/cat/All.Data.3TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/'
        #self.parent.root_dir = '/media/cat/All.Data.3TB/in_vivo/nick/ptc21/tr5c/'
        
#***************************** Cat DATA *************************************
        ##ptc21 tr5c 
        self.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/synch_sort/61-tr5c-blankscreen_alltrack.ptcs'
        self.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/synch_sort/61-tr5c-blankscreen_alltrack_lfp_50compressed.ptcs'
        self.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/synch_sort/61-tr5c-blankscreen_alltrack_lfp.tsf'


        #***************************** MOUSE DATA *************************************
        #2017-02_03 VISUAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed_4.0threshold_3clusters.ptcs'
        ##self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack.tsf'

        #2017-01_31 BARREL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_lfp_250hz_alltrack_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_31_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170131_164034_lfp_250hz_alltrack.tsf'

        #2017_01_30 AUDITORY - 2 clusters only !?
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_30_auditory_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170130_164612_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_30_auditory_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170130_164612_lfp_250hz_alltrack_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_30_auditory_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170130_164612_lfp_250hz_alltrack.tsf'

        #2017-01_26 BARREL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_lfp_100hz_alltrack_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2017_01_26_barrel_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170126_153637_lfp_100hz_alltrack.tsf'
        

        #************************************************ OLDER DATA **********************************
    
        #2016_08_31 VSD VISUAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_08_31_vsd_visual/sort_alltrack/track1_spontaneous_1_160831_213746_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_08_31_vsd_visual/sort_alltrack/track1_spontaneous_1_160831_213746_lfp_250hz_alltrack_lowcut0.1_highcut110.0_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2016_08_31_vsd_visual/sort_alltrack/track1_spontaneous_1_160831_213746_lfp_250hz_alltrack.tsf'



        #************************************************ SUBCORTICAL LFP CLUSTERS **********************************
        #2016_07_26 AUDITORY   - Multiople clusters - *********** SUBCORTICAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_lfp_250hz_alltrack_50compressed.ptcs'#

        #2016_07_15 VISUAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd_visual/sort_alltrack_spontaneous/track1_150Hz_1st_spontaneous_10iso_160715_181445_hp_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd_visual/sort_alltrack_spontaneous/track1_150Hz_1st_spontaneous_10iso_160715_181445_lfp_250hz_alltrack_50compressed.ptcs'


        ##********************************* CHRONIC V1 ************************************
        #ALL EPOCHS
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack.tsf'
        
        #EPOCH #0
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_hp_butter_alltrack_epoch_0.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack_notch_50compressed_epoch_0.ptcs2.npz'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack.tsf'

        #EPOCH #1
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_hp_butter_alltrack_epoch_1.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack_notch_50compressed_epoch_1.ptcs2.npz'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/Mouse_42/sorted_11_12_iso/spontaneous_42_iso_170311_103957_lfp_250hz_alltrack.tsf'

        #*************************************************************************************************************************************************************************
        #*************************************************************************************************************************************************************************
        #*************************************************************************************************************************************************************************
        
        row_index = 0

        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.selected_recording  = self.root_dir
        self.select_recording_lbl = QLabel(self.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1)
        
        #Rec length
        self.rec_length = QLabel('0', self)
        self.rec_length.setMaximumWidth(100)
        self.rec_length_lbl = QLabel('Rec Length:', self)
        self.rec_length_lbl.setMaximumWidth(100)
        layout.addWidget(self.rec_length_lbl, row_index, 3)
        layout.addWidget(self.rec_length, row_index, 4); row_index+=1

        #Select recording
        self.button_select_sort = QPushButton('Select Event File')
        self.button_select_sort.setMaximumWidth(200)
        self.button_select_sort.clicked.connect(self.slct_event_file)
        layout.addWidget(self.button_select_sort, row_index, 0)
        
        self.selected_sort  = self.root_dir
        self.select_sort_lbl = QLabel(self.selected_sort, self)
        layout.addWidget(self.select_sort_lbl, row_index,1)
    

        #% of electrodes
        self.n_electrodes = QLineEdit('1.0');                #parent.start_time = self.start_time
        self.n_electrodes.setMaximumWidth(50)
        self.n_electrodes_lbl = QLabel('%Electrodes', self)
        self.n_electrodes_lbl.setMaximumWidth(100)
        layout.addWidget(self.n_electrodes_lbl, row_index, 5)
        layout.addWidget(self.n_electrodes, row_index, 6)
        
        
        self.low_cutoff = QLineEdit('0.1');                #parent.start_time = self.start_time
        self.low_cutoff.setMaximumWidth(50)
        self.low_cutoff_lbl = QLabel('Lowcut', self)
        self.low_cutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.low_cutoff_lbl, row_index, 7)
        layout.addWidget(self.low_cutoff, row_index, 8)


        self.high_cutoff = QLineEdit('110');                #parent.start_time = self.start_time
        self.high_cutoff.setMaximumWidth(50)
        self.high_cutoff_lbl = QLabel('Highcut', self)
        self.high_cutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.high_cutoff_lbl, row_index,9)
        layout.addWidget(self.high_cutoff, row_index, 10); row_index+=1
        
        #View traces
        self.button_view_traces = QPushButton('View Traces')
        self.button_view_traces.setMaximumWidth(200)
        self.button_view_traces.clicked.connect(self.vw_traces)
        layout.addWidget(self.button_view_traces, row_index, 0)

        self.start_time = QLineEdit('0');                #parent.start_time = self.start_time
        self.start_time.setMaximumWidth(50)
        self.start_time_lbl = QLabel('Start (sec):', self)
        self.start_time_lbl.setMaximumWidth(100)
        layout.addWidget(self.start_time_lbl, row_index,1)
        layout.addWidget(self.start_time, row_index, 2)

        self.end_time = QLineEdit('60');                  #parent.end_time = self.end_time
        self.end_time.setMaximumWidth(50)
        self.end_time_lbl = QLabel('End (sec):', self)
        self.end_time_lbl.setMaximumWidth(100)
        layout.addWidget(self.end_time_lbl, row_index, 3)
        layout.addWidget(self.end_time, row_index,4)
        
        self.probe_penentration = QLineEdit('1.0');                  #parent.end_time = self.end_time
        self.probe_penentration.setMaximumWidth(50)
        self.probe_penentration_lbl = QLabel('% Probe in Tissue', self)
        self.probe_penentration_lbl.setMaximumWidth(150)
        layout.addWidget(self.probe_penentration_lbl, row_index, 5)
        layout.addWidget(self.probe_penentration, row_index,6)
        
        self.voltage_scale = QLineEdit('10.0');                  #parent.end_time = self.end_time
        self.voltage_scale.setMaximumWidth(50)
        self.voltage_scale_lbl = QLabel('Voltage Scaling', self)
        self.voltage_scale_lbl.setMaximumWidth(150)
        layout.addWidget(self.voltage_scale_lbl, row_index, 7)
        layout.addWidget(self.voltage_scale, row_index,8)
        
        self.setLayout(layout)

    def slct_event_file(self):
        self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, ".ptcs or .txt or .npy)", self.selected_recording,"PTCS, TXT, NPY (*.ptcs *.txt *.npy)")
        path_name, file_name = os.path.split(self.selected_sort)
        self.select_sort_lbl.setText(file_name)

        #Reload session list and reset session box
        #self.comboBox_selected_unit.clear()
        if '.ptcs' in self.selected_sort: 
            self.Sort = Ptcs(self.selected_sort)
            n_units = len(self.Sort.units)
            #self.selected_unit_spikes_lbl.setText(str(len(self.Sort.units[0])))
            
        elif '.txt' in self.selected_sort: 
            n_units = 1
            self.triggers = np.loadtxt(self.selected_sort)
            #self.selected_unit_spikes_lbl.setText(str(len(self.triggers)))
        
        elif '.npy' in self.selected_sort: 
            n_units = 1
            self.camera_pulses = np.load(self.selected_sort)
            
            parse_camera_pulses(self); print "...# of triggers: ", len(self.triggers)
           
            print self.triggers/25000.
            print self.triggers_length/25000.
            


    def slct_recording(self):
        #self.selected_recording =  QtGui.QFileDialog.getOpenFileName(self, 'Load File', self.selected_recording)
        self.selected_recording =  QtGui.QFileDialog.getOpenFileName(self, ".tsf", self.selected_recording,"*.tsf")

        path_name, file_name = os.path.split(self.selected_recording)

        self.select_recording_lbl.setText(file_name)

        #self.tsf = Tsf_file(self.selected_recording)
        #self.tsf.read_ec_traces()
        
        #print "...len rec: ", self.tsf.n_vd_samples/float(self.tsf.SampleFrequency)
        #self.rec_length.setText(str(self.tsf.n_vd_samples/float(self.tsf.SampleFrequency)))
   
   
    def vw_traces(self):
        view_traces(self)
