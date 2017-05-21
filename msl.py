
import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class MSL(QtGui.QWidget):
    def __init__(self, parent):
        super(MSL, self).__init__(parent)
        self.parent = parent
        
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        #self.parent.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/'
        #self.parent.root_dir = '/media/cat/2TB/in_vivo/nick/ptc21_tr5c/tsf_files/'
        #self.parent.root_dir = '/media/cat/All.Data.3TB/in_vivo/nick/ptc21/tr5c/'

        '''
        #************************************************************************************************************************************************************
        #************************************************************************************************************************************************************
        #*************************************************** Cat DATA ************************************************************************
        #************************************************************************************************************************************************************
        #************************************************************************************************************************************************************
        '''
        
        #**************************** PTC17 ***********************
        #ptc17 tr1 - 1202um electrode
        self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr1/sort_all_track/01-tr1-driftgrating_alltrack_sync.ptcs2.npz'
        self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr1/sort_all_track/01-tr1-driftgrating_lp_alltrack_desynch_zeroed_50compressed.ptcs'
        self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr1/sort_all_track/01-tr1-driftgrating_lp_alltrack_desynch_zeroed.tsf'
        
        #ptc17 tr2 - 1202um electrode
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr2b/sort_all_track/44-tr2b-blankscreen_alltrack_sync.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr2b/sort_all_track/44-tr2b-blankscreen_lp_alltrack_desynch_zeroed_0.7si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc17/tr2b/sort_all_track/44-tr2b-blankscreen_lp_alltrack_desynch_zeroed_0.7si.tsf'
        
        #**************************** PTC18 ***********************
        ##ptc18 tr1 - 1325um electrode
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr1/sort_alltrack/01-tr1-driftbar_longbar_alltrack.ptcs' #IS THIS SYNC PERIOD ONLY?!??!?!
        #MUA
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr1/sort_alltrack/01-tr1-driftbar_longbar_alltrack_MUA_5x.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr1/sort_alltrack/01-tr1-driftbar_longbar_lp_alltrack_desynch_zeroed_0.5si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr1/sort_alltrack/01-tr1-driftbar_longbar_lp_alltrack_desynch_zeroed_0.5si_50compressed.tsf'

        ##ptc18 tr2 - 1325um electrode
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr2c/sort_alltrack/47-tr2c-csd_alltrack.ptcs' #IS THIS SYNC PERIOD ONLY !?!??!?!
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr2c/sort_alltrack/47-tr2c-csd_lp_alltrack_notch_desynch_zeroed_0.7si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc18/tr2c/sort_alltrack/47-tr2c-csd_lp_alltrack_notch_desynch_zeroed_0.7si.tsf'

        
        #***************** PTC20 *************************
        #ptc20 tr1 - 1325um electrode
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr1/sort_alltrack/01-tr1-csd_alltrack_sync.ptcs2.npz'  #HAS THIS BEEN DESCYNCHED!??!?!
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr1/sort_alltrack/01-tr1-csd_lp_alltrack_notch_desynch_zeroed_0.4si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr1/sort_alltrack/01-tr1-csd_lp_alltrack_notch_desynch_zeroed_0.4si.tsf'
    
        #ptc20 tr2 - NOT GOOD TRACK; 1325um electrode **********************************
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr2/sort_alltrack/34-tr2-blankscreen_alltrack.ptcs'    #HAS THIS BEEN DESCYNCHED!??!?!
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr2/sort_alltrack/34-tr2-blankscreen_lp_alltrack_notch_desynch_zeroed_0.5si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr2/sort_alltrack/34-tr2-blankscreen_lp_alltrack_notch_desynch_zeroed_0.5si_50compressed.tsf'
    
        #ptc20 tr3  - 1700um electrode
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr3/sort_alltrack/58-tr3-spontaneous_alltrack.ptcs'    #HAS THIS BEEN DESCYNCHED!??!?!
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr3/sort_alltrack/58-tr3-spontaneous_lp_alltrack_notch_desynch_zeroed_0.5si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc20/tr3/sort_alltrack/58-tr3-spontaneous_lp_alltrack_notch_desynch_zeroed_0.5si_50compressed.tsf'


        #*************** PTC21 *******************
        #ptc21 tr1 - 1755um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr1/sort_alltrack/01-tr1-csd_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr1/sort_alltrack/01-tr1-csd_lp_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr1/sort_alltrack/01-tr1-csd_lp_alltrack_notch_50compressed.tsf'

        #ptc21 tr2 - 1755um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr2/sort_alltrack/13-tr2-spontaneous_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr2/sort_alltrack/13-tr2-spontaneous_lp_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr2/sort_alltrack/13-tr2-spontaneous_lp_alltrack_notch_50compressed.tsf'

        #ptc21 tr5c  - 1755um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_lp_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_lp_alltrack_notch_50compressed.tsf'

        #ptc21 tr6b - 1755um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr6b/sort_alltrack/80-t6b-csd_binocular_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr6b/sort_alltrack/80-t6b-csd_binocular_lp_alltrack_notch_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr6b/sort_alltrack/80-t6b-csd_binocular_lp_alltrack_notch_50compressed.tsf'
    

        #****************** PTC22 **********************
        #tr1 - 1202um electrode 


        #tr2 - 1202um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr2/sort_alltrack/23-tr2-csd_alltrack.ptcs'
        #SYNC SUA
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr2/sort_alltrack/23-tr2-csd_alltrack_0.7si_sync.ptcs2.npz'
        #SYNC MUA:
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr2/sort_alltrack/23-tr2-csd_alltrack_MUA_5x_0.7si_sync.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr2/sort_alltrack/23-tr2-csd_lp_alltrack_notch_desynch_zeroed_0.6si_50compressed.ptcs'
        #elf.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr2/sort_alltrack/23-tr2-csd_lp_alltrack_notch_desynch_zeroed_0.6si_50compressed.tsf'
        
                
        #tr4b - 1202um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_alltrack_0.8si_sync.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_lp_alltrack_notch_desynch_zeroed_0.7si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_lp_alltrack_notch_desynch_zeroed_0.7si.tsf'
    

        #tr5b - 1202um electrode 
        #self.parent.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_alltrack_0.8si_sync.ptcs2.npz'
        #self.parent.lfp_event_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_lp_alltrack_notch_desynch_zeroed_0.7si_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc22/tr4b/sort_alltrack/40-tr4b-csd_lp_alltrack_notch_desynch_zeroed_0.7si.tsf'
    
    
        '''
        #************************************************************************************************************************************************************
        #************************************************************************************************************************************************************
        #***************************************************************** MOUSE DATA ***********************************************************
        #************************************************************************************************************************************************************
        #************************************************************************************************************************************************************
        '''
        
        #2017-02_03 VISUAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed_4.0threshold_3clusters.ptcs'
        ###self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2017_02_03_visual_ephys_ophys/sort_alltrack_spontaneous/track_1_spontaneous_1_170203_172405_lfp_250hz_alltrack_50compressed.ptcs'
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
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_26_vsd_auditory/sort_alltrack2/track2_spontaneous_1iso_160726_215426_lfp_250hz_alltrack_50compressed.tsf'


        #2016_07_20 BARREL  
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_20_vsd_barrel/sort_track2/track2_iso75_spontaneous_160720_220333_hp_butter_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_20_vsd_barrel/sort_track2/track2_iso75_spontaneous_160720_220333_lfp_250hz_alltrack_50compressed.ptcs' 
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_20_vsd_barrel/sort_track2/track2_iso75_spontaneous_160720_220333_lfp_250hz_alltrack_50compressed.tsf'


        #2016_07_15 VISUAL
        #self.parent.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd_visual/sort_alltrack_spontaneous/track1_150Hz_1st_spontaneous_10iso_160715_181445_hp_alltrack.ptcs'
        #self.parent.lfp_event_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd_visual/sort_alltrack_spontaneous/track1_150Hz_1st_spontaneous_10iso_160715_181445_lfp_250hz_alltrack_50compressed.ptcs'
        #self.parent.lfp_tsf_file = '/media/cat/12TB/in_vivo/tim/cat/2016_07_15_vsd_visual/sort_alltrack_spontaneous/track1_150Hz_1st_spontaneous_10iso_160715_181445_lfp_250hz_alltrack_50compressed.tsf'

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

        

        layout = QtGui.QGridLayout()

        row_index=0

        #**************************************************************************************
        #*********************************** LOAD FILES  **************************************
        #**************************************************************************************
        self.preprocess_lbl = QLabel('MSL ANALYSIS', self)
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
        layout.addWidget(self.button_set_sua_file_lbl, row_index, row_index, 1, 9); row_index+=1

        #Set LFP event file
        self.button_set_lfp_event_file = QPushButton('LPF Event File')
        self.button_set_lfp_event_file.setMaximumWidth(200)
        self.button_set_lfp_event_file.clicked.connect(self.set_lfp_event_file)
        layout.addWidget(self.button_set_lfp_event_file, row_index, 0)
        self.button_set_lfp_event_file_lbl = QLabel(os.path.split(self.parent.lfp_event_file)[1], self)
        layout.addWidget(self.button_set_lfp_event_file_lbl, row_index, row_index, 1, 9); row_index+=1

        #Set LFP event file
        self.button_set_lfp_tsf = QPushButton('LPF TSF (uncomp.) ')
        self.button_set_lfp_tsf.setMaximumWidth(200)
        self.button_set_lfp_tsf.clicked.connect(self.set_lfp_tsf_file)
        layout.addWidget(self.button_set_lfp_tsf, row_index, 0)
        self.button_set_lfp_tsf_lbl = QLabel(os.path.split(self.parent.lfp_tsf_file)[1], self)
        layout.addWidget(self.button_set_lfp_tsf_lbl, row_index, row_index, 1, 9); row_index+=1

        #**************************************************************************************
        #*********************************** SET LFP PARAMETERS *******************************
        #**************************************************************************************
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1

        self.preprocess_lbl = QLabel('PARAMETERS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        
        parent.lfp_cluster = QLineEdit('0');                   
        parent.lfp_cluster.setMaximumWidth(50)
        lfp_cluster_lbl = QLabel('LFP Cluster #:', self)
        layout.addWidget(lfp_cluster_lbl, row_index,0)
        layout.addWidget(parent.lfp_cluster, row_index,1)

        parent.time_chunks = QLineEdit('1');               
        parent.time_chunks.setMaximumWidth(50)
        time_chunks_lbl = QLabel('# of time_chunks:', self)
        layout.addWidget(time_chunks_lbl,row_index,2)
        layout.addWidget(parent.time_chunks,row_index,3)

        self.cmap = QLineEdit('viridis');               
        self.cmap.setMaximumWidth(250)
        cmap_lbl = QLabel('Colour Map', self)
        layout.addWidget(cmap_lbl,row_index,4)
        layout.addWidget(self.cmap,row_index,5)
                
        parent.lock_window_start = QLineEdit('-100');               
        parent.lock_window_start.setMaximumWidth(50)
        lock_window_start_lbl = QLabel('start window (ms):', self)
        layout.addWidget(lock_window_start_lbl,row_index,6)
        layout.addWidget(parent.lock_window_start,row_index,7)

        parent.lock_window_end = QLineEdit('100');               
        parent.lock_window_end.setMaximumWidth(50)
        lock_window_end_lbl = QLabel('end window (ms):', self)
        layout.addWidget(lock_window_end_lbl,row_index,8)
        layout.addWidget(parent.lock_window_end,row_index,9); row_index+=1


        self.min_spikes = QLineEdit('1');               
        self.min_spikes.setMaximumWidth(50)
        self.min_spikes_lbl = QLabel('min_spikes:', self)
        layout.addWidget(self.min_spikes_lbl,row_index,0)
        layout.addWidget(self.min_spikes,row_index,1)

        self.sigma_width = QLineEdit('20');               
        self.sigma_width.setMaximumWidth(50)
        self.sigma_width_lbl = QLabel('Sigma Width:', self)
        layout.addWidget(self.sigma_width_lbl,row_index,2)
        layout.addWidget(self.sigma_width,row_index,3)
        
        self.min_fire_rate = QLineEdit('0.001');               
        self.min_fire_rate.setMaximumWidth(50)
        self.min_fire_rate_lbl = QLabel('min_rate:', self)
        layout.addWidget(self.min_fire_rate_lbl,row_index,4)
        layout.addWidget(self.min_fire_rate,row_index,5)
        
       
        self.starting_cell = QLineEdit('0');               
        self.starting_cell.setMaximumWidth(50)
        self.starting_cell_lbl = QLabel('start cell:', self)
        layout.addWidget(self.starting_cell_lbl,row_index,6)
        layout.addWidget(self.starting_cell,row_index,7)
        
        
        self.ending_cell = QLineEdit('10');               
        self.ending_cell.setMaximumWidth(50)
        self.ending_cell_lbl = QLabel('end cell:', self)
        layout.addWidget(self.ending_cell_lbl,row_index,8)
        layout.addWidget(self.ending_cell,row_index,9); row_index+=1
        
        

        #**************************************************************************************
        #*********************************** MSL PARAMETERS *******************************
        #**************************************************************************************
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1

        self.preprocess_lbl = QLabel('MSL PLOTS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        
        
        self.sliding_window_length = QLineEdit('60');               
        self.sliding_window_length.setMaximumWidth(50)
        self.sliding_window_length_lbl = QLabel('Window (mins)', self)
        layout.addWidget(self.sliding_window_length_lbl,row_index,0)
        layout.addWidget(self.sliding_window_length,row_index,1)
        

        self.sliding_window_step = QLineEdit('1');               
        self.sliding_window_step.setMaximumWidth(50)
        self.sliding_window_step_lbl = QLabel('Step size (mins)', self)
        layout.addWidget(self.sliding_window_step_lbl,row_index,2)
        layout.addWidget(self.sliding_window_step,row_index,3)

        self.multiple_units = QLineEdit('0, 1, 2');               
        self.multiple_units.setMaximumWidth(50)
        self.multiple_units_lbl = QLabel('Selected Units', self)
        layout.addWidget(self.multiple_units_lbl,row_index,4)
        layout.addWidget(self.multiple_units,row_index,5); row_index+=1



        #************* PREPROCESSING: Rasters and Histograms ****************************
        self.button_lfp_rasters = QPushButton('Cell LFP Rasters')
        self.button_lfp_rasters.setMaximumWidth(200)
        self.button_lfp_rasters.clicked.connect(self.compute_lfp_rasters)
        layout.addWidget(self.button_lfp_rasters, row_index, 0)
        
        
        self.button_lfp_histograms = QPushButton('Cell LFP Histograms')
        self.button_lfp_histograms.setMaximumWidth(200)
        self.button_lfp_histograms.clicked.connect(self.compute_lfp_histograms)
        layout.addWidget(self.button_lfp_histograms, row_index, 1)
        
        
        self.button_peth = QPushButton('PETH & Rasters')
        self.button_peth.setMaximumWidth(200)
        self.button_peth.clicked.connect(self.view_peth)
        layout.addWidget(self.button_peth, row_index, 2)
        
        
        self.button_spikerates = QPushButton('Epoch Spike-Rates')
        self.button_spikerates.setMaximumWidth(200)
        self.button_spikerates.clicked.connect(self.msl_spikerates)
        layout.addWidget(self.button_spikerates, row_index, 3); row_index+=1
        
        
        self.button_msl_depth = QPushButton('Compute MSL Depth')
        self.button_msl_depth.setMaximumWidth(200)
        self.button_msl_depth.clicked.connect(self.compute_msl_depth)
        layout.addWidget(self.button_msl_depth, row_index, 0)

        self.button_msl = QPushButton('Compute MSL Chunks')
        self.button_msl.setMaximumWidth(200)
        self.button_msl.clicked.connect(self.compute_msl_chunks)
        layout.addWidget(self.button_msl, row_index, 1)
        
        self.button_MSL_drift = QPushButton('Single Cell MSL Drift')
        self.button_MSL_drift.setMaximumWidth(200)
        self.button_MSL_drift.clicked.connect(self.view_msl_drift)
        layout.addWidget(self.button_MSL_drift, row_index, 2)
        
        self.button_allcell_drift = QPushButton('All Cell MSL Stats')
        self.button_allcell_drift.setMaximumWidth(200)
        self.button_allcell_drift.clicked.connect(self.view_allcell_drift)
        layout.addWidget(self.button_allcell_drift, row_index, 4)
        
        self.button_drift_trends = QPushButton('MSL Drift Trends')
        self.button_drift_trends.setMaximumWidth(200)
        self.button_drift_trends.clicked.connect(self.view_drift_trends)
        layout.addWidget(self.button_drift_trends, row_index, 5)

        self.button_drift_movie = QPushButton('MSL Drift Movies')
        self.button_drift_movie.setMaximumWidth(200)
        self.button_drift_movie.clicked.connect(self.view_drift_movie)
        layout.addWidget(self.button_drift_movie, row_index, 6); row_index+=1


        
        #************** SINGLE CELL MSL FUNCTIONS ***********

        self.button_msl_continuous = QPushButton('MSL Sliding Time Window')
        self.button_msl_continuous.setMaximumWidth(200)
        self.button_msl_continuous.clicked.connect(self.view_msl_continuous)
        layout.addWidget(self.button_msl_continuous, row_index, 0)
        
        self.button_msl_continuous_single = QPushButton('MSL Sliding Single Cell')
        self.button_msl_continuous_single.setMaximumWidth(200)
        self.button_msl_continuous_single.clicked.connect(self.compute_msl_continuous_single)
        layout.addWidget(self.button_msl_continuous_single, row_index, 1); row_index+=1
        

        self.button_msl_discrete_single = QPushButton('Single Cell MSL - Discrete')
        self.button_msl_discrete_single.setMaximumWidth(200)
        self.button_msl_discrete_single.clicked.connect(self.view_msl_discrete_single)
        layout.addWidget(self.button_msl_discrete_single, row_index,0)

        self.button_msl_continuous_single = QPushButton('MSL Sliding Win - Single Unit')
        self.button_msl_continuous_single.setMaximumWidth(200)
        self.button_msl_continuous_single.clicked.connect(self.view_msl_continuous_single)
        layout.addWidget(self.button_msl_continuous_single, row_index, 1)


        self.button_msl_continuous_multi = QPushButton('MSL Sliding Win - Multi Unit')
        self.button_msl_continuous_multi.setMaximumWidth(200)
        self.button_msl_continuous_multi.clicked.connect(self.view_msl_continuous_multi)
        layout.addWidget(self.button_msl_continuous_multi, row_index, 2)
        

        self.button_msl_single_lfpevent = QPushButton('MSL Single Unit - Single Event')     
        self.button_msl_single_lfpevent.setMaximumWidth(200)
        self.button_msl_single_lfpevent.clicked.connect(self.view_msl_single_lfpevent)
        layout.addWidget(self.button_msl_single_lfpevent, row_index, 3); 
        
        
                
        self.max_error = QLineEdit('10');               
        self.max_error.setMaximumWidth(50)
        self.max_error_lbl = QLabel('Max Error (ms)', self)
        layout.addWidget(self.max_error_lbl,row_index,4)
        layout.addWidget(self.max_error, row_index,5)
        
        
        
        self.ylim_bottom = QLineEdit('-40');               
        self.ylim_bottom.setMaximumWidth(50)
        self.ylim_bottom_lbl = QLabel('Ylim (bottom)', self)
        layout.addWidget(self.ylim_bottom_lbl,row_index,6)
        layout.addWidget(self.ylim_bottom, row_index,7)
        
        
        self.ylim_top = QLineEdit('20');               
        self.ylim_top.setMaximumWidth(50)
        self.ylim_top_lbl = QLabel('Ylim (top)', self)
        layout.addWidget(self.ylim_top_lbl,row_index,8)
        layout.addWidget(self.ylim_top, row_index,9); row_index+=1


        self.button_msl_state_space = QPushButton('MSL - State Space')     #WHAT DOES THIS DO AGAIN? 
        self.button_msl_state_space.setMaximumWidth(200)
        self.button_msl_state_space.clicked.connect(self.view_msl_state_space)
        layout.addWidget(self.button_msl_state_space, row_index, 0)
        
        
        self.button_window_pval = QPushButton('MSL - Windowed P-Val Plots')     #WHAT DOES THIS DO AGAIN? 
        self.button_window_pval.setMaximumWidth(200)
        self.button_window_pval.clicked.connect(self.view_window_pval)
        layout.addWidget(self.button_window_pval, row_index, 2); row_index+=1

        
        
        #**************************************************************************************
        #*********************************** SPIKING STATISTICS *******************************
        #**************************************************************************************
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1

        self.preprocess_lbl = QLabel('SPIKING STATISTICS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        

        self.button_percentages = QPushButton('% Lock: Depth and F-Rate')
        self.button_percentages.setMaximumWidth(200)
        self.button_percentages.clicked.connect(self.sua_lockpercentage)
        layout.addWidget(self.button_percentages, row_index, 0)
        
        self.start_window = QLineEdit('-100');               
        self.start_window.setMaximumWidth(50)
        self.start_window_lbl = QLabel('Start Window (ms)', self)
        layout.addWidget(self.start_window_lbl,row_index,1)
        layout.addWidget(self.start_window, row_index,2)

        self.end_window = QLineEdit('100');               
        self.end_window.setMaximumWidth(50)
        self.end_window_lbl = QLabel('End Window (ms)', self)
        layout.addWidget(self.end_window_lbl,row_index,3)
        layout.addWidget(self.end_window,row_index,4); row_index+=1
              
        
        self.button1 = QPushButton('Compute All P-Values')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.compute_msl_Pvals)
        layout.addWidget(self.button1, row_index, 0)


        self.view_pvalue = QPushButton('View P-Value')
        self.view_pvalue.setMaximumWidth(200)
        self.view_pvalue.clicked.connect(self.view_pval)
        layout.addWidget(self.view_pvalue, row_index, 1)

                
        self.vmin_value = QLineEdit('0.01');               
        self.vmin_value.setMaximumWidth(50)
        self.vmin_value_lbl = QLabel('Vmin Value', self)
        layout.addWidget(self.vmin_value_lbl,row_index,2)
        layout.addWidget(self.vmin_value,row_index,3)

        self.vmax_value = QLineEdit('1.0');               
        self.vmax_value.setMaximumWidth(50)
        self.vmax_value_lbl = QLabel('Vmax Value', self)
        layout.addWidget(self.vmax_value_lbl,row_index,4)
        layout.addWidget(self.vmax_value,row_index,5); row_index+=1

        
        
        self.specgram_ch = QLineEdit('9');               
        self.specgram_ch.setMaximumWidth(50)
        self.specgram_ch_lbl = QLabel('Specgram Ch', self)
        layout.addWidget(self.specgram_ch_lbl,row_index,8)
        layout.addWidget(self.specgram_ch,row_index,9); 


        #**************************************************************************************
        #***************************** CSD ANALYSIS **************************************
        #**************************************************************************************
        
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1
        self.preprocess_lbl = QLabel('CSD ANALYSIS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        #Set recording name
        self.button_set_recording = QPushButton('CSD Recording')
        self.button_set_recording.setMaximumWidth(200)
        self.button_set_recording.clicked.connect(self.set_recording)
        layout.addWidget(self.button_set_recording, row_index, 0)
        
        self.parent.button_set_recording = ''
        self.button_set_recording_lbl = QLabel(self.parent.button_set_recording, self)
        layout.addWidget(self.button_set_recording_lbl, row_index,1); row_index+=1

        self.button_csd_rasters = QPushButton('CSD: Stimulus + Rasters')
        self.button_csd_rasters.setMaximumWidth(200)
        self.button_csd_rasters.clicked.connect(self.view_csd_rasters)
        layout.addWidget(self.button_csd_rasters, row_index, 0)

        self.button_csd_histogram = QPushButton('CSD: SUA Histograms')
        self.button_csd_histogram.setMaximumWidth(200)
        self.button_csd_histogram.clicked.connect(self.view_csd_histogram)
        layout.addWidget(self.button_csd_histogram, row_index, 1); row_index+=1


        #**************************************************************************************
        #***************************** NAT SCENE ANALYSIS **************************************
        #**************************************************************************************
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1
        self.preprocess_lbl = QLabel('NAT SCENE ANALYSIS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.button_set_natscene_rec = QPushButton('Nat Scene Recording')
        self.button_set_natscene_rec.setMaximumWidth(200)
        self.button_set_natscene_rec.clicked.connect(self.set_nat_scene_rec)
        layout.addWidget(self.button_set_natscene_rec, row_index, 0)

        self.parent.button_set_natscene_rec = ''
        self.button_set_natscene_rec_lbl = QLabel(self.parent.button_set_natscene_rec, self)
        layout.addWidget(self.button_set_natscene_rec_lbl, row_index,1); row_index+=1
        
        self.button_natscene_rasters = QPushButton('Nat Scene: Stimulus + Rasters')
        self.button_natscene_rasters.setMaximumWidth(200)
        self.button_natscene_rasters.clicked.connect(self.view_natscene_rasters)
        layout.addWidget(self.button_natscene_rasters, row_index, 0)

        self.button_natscene_seqisi = QPushButton('Sequential ISI Distributions')
        self.button_natscene_seqisi.setMaximumWidth(200)
        self.button_natscene_seqisi.clicked.connect(self.view_natscene_seqisi)
        layout.addWidget(self.button_natscene_seqisi, row_index, 1)
        
        self.button_natscene_pairisi = QPushButton('Pairwise ISI Distributions')
        self.button_natscene_pairisi.setMaximumWidth(200)
        self.button_natscene_pairisi.clicked.connect(self.view_natscene_pairisi)
        layout.addWidget(self.button_natscene_pairisi, row_index, 2)

        
        self.bin_width = QLineEdit('0.001');               
        self.bin_width.setMaximumWidth(50)
        self.bin_width_lbl = QLabel('Bin Width (sec)', self)
        layout.addWidget(self.bin_width_lbl,row_index,4)
        layout.addWidget(self.bin_width,row_index,5)
        
        self.seqisi_start = QLineEdit('0');               
        self.seqisi_start.setMaximumWidth(50)
        self.seqisi_start_lbl = QLabel('Time Start (sec)', self)
        layout.addWidget(self.seqisi_start_lbl,row_index,6)
        layout.addWidget(self.seqisi_start,row_index,7)
        
        self.seqisi_end = QLineEdit('5.5');               
        self.seqisi_end.setMaximumWidth(50)
        self.seqisi_end_lbl = QLabel('Time End (sec)', self)
        layout.addWidget(self.seqisi_end_lbl,row_index,8)
        layout.addWidget(self.seqisi_end,row_index,9); row_index+=1


        #**************************************************************************************
        #***************************** Pairise spiking analysis **************************************
        #**************************************************************************************
        layout.addWidget(QLabel('', self), row_index,0); row_index+=1
        self.preprocess_lbl = QLabel('DISTRIBUTION ANALYSIS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.button_triplet_sequences = QPushButton('Triplet Sequences')
        self.button_triplet_sequences.setMaximumWidth(200)
        self.button_triplet_sequences.clicked.connect(self.view_triplet_sequences)
        layout.addWidget(self.button_triplet_sequences, row_index, 0); row_index+=1

        self.button_nspikes_perlfp = QPushButton('Nspikes per LFP event')
        self.button_nspikes_perlfp.setMaximumWidth(200)
        self.button_nspikes_perlfp.clicked.connect(self.view_nspikes_histograms)
        layout.addWidget(self.button_nspikes_perlfp, row_index, 0); row_index+=1

        self.button_isi_perlfp = QPushButton('ISI from LFP events')
        self.button_isi_perlfp.setMaximumWidth(200)
        self.button_isi_perlfp.clicked.connect(self.view_isi_histograms)
        layout.addWidget(self.button_isi_perlfp, row_index, 0); row_index+=1


        self.setLayout(layout)

    
    #**************************** LOAD FILE ROUTINES *******************************
    def set_sua_file(self):
        #self.parent.sua_file = QtGui.QFileDialog.getOpenFileName(self, "*.ptcs (*.ptcs)", self.parent.root_dir,"ptcs (*.ptcs; *.npz)")
        self.parent.sua_file = QtGui.QFileDialog.getOpenFileName(self, "*.ptcs (*.ptcs)", self.parent.root_dir, "*.ptcs *.npz")
        self.button_set_sua_file_lbl.setText(self.parent.sua_file.replace(os.path.dirname(self.parent.sua_file),''))
        #self.parent.setWindowTitle(self.parent.sua_file)


    def set_lfp_event_file(self):
        self.parent.lfp_event_file = QtGui.QFileDialog.getOpenFileName(self, "ptcs (*.ptcs)", self.parent.sua_file,"*.ptcs *.npz")
        self.button_set_lfp_event_file_lbl.setText(self.parent.lfp_event_file.replace(os.path.dirname(self.parent.lfp_event_file),''))
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))
    
    def set_lfp_tsf_file(self):
        self.parent.lfp_tsf_file  = QtGui.QFileDialog.getOpenFileName(self, "tsf (*.tsf)", self.parent.sua_file,"*.tsf")
        self.button_set_lfp_tsf_lbl.setText(self.parent.lfp_tsf_file.replace(os.path.dirname(self.parent.lfp_tsf_file),''))
    
    
    
    def set_recording(self):
        temp_dir = '/media/cat/2TB/in_vivo/nick/ptc21/'

        self.parent.recording_dir = QtGui.QFileDialog.getExistingDirectory(self, "Select recording", temp_dir)
        
        self.rec_path = os.path.dirname(self.parent.recording_dir)
        self.rec_name = self.parent.recording_dir.replace(self.rec_path,'')
        self.button_set_recording_lbl.setText(self.rec_name)
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))


    def set_nat_scene_rec(self):
        temp_dir = self.parent.root_dir

        self.parent.recording_dir = QtGui.QFileDialog.getExistingDirectory(self, "Select recording", temp_dir)
        
        self.rec_path = os.path.dirname(self.parent.recording_dir)
        self.rec_name = self.parent.recording_dir.replace(self.rec_path,'')
        self.button_set_natscene_rec_lbl.setText(self.rec_name)
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))
        

    #***************************** COMPUTE MSL ************************************
    def compute_lfp_rasters(self):          #This should be the first function called, it computes each cell rasters for each LFP event 
        Compute_LFP_rasters(self)
        
    
    def compute_lfp_histograms(self):
        Compute_LFP_histograms(self)
    
    
    #************************** VIEWING FUCNTIONS?!
    
    def view_peth(self):
        peth_scatter_plots(self)
    
    def msl_spikerates(self):
        compute_msl_spikerates(self)

        
    def compute_msl_chunks(self, main_widget):
        Compute_MSL_chunks(self)                    #MSL CHUNKS!!!
    

    def compute_msl_depth(self, main_widget):
        Compute_MSL_depth(self)                     #MSL DEPTH

    def view_msl_single_lfpevent(self):
        compute_msl_single_lfpevent(self)

    
    def view_msl_continuous(self):                      #Compute MSL SLIDING WINDOW - ALL CELLS  ***************************
        compute_msl_continuous(self)


    def compute_msl_continuous_single(self):
        compute_msl_continuous_single_cell(self)        #Compute MSL SLIDING WINDOW - single cell
        


    #****************************** STATE SPACE ANALYSIS **********************************
    
    def view_msl_state_space(self):
        compute_msl_state_space(self)
    
    def view_window_pval(self):
        compute_window_pval(self)
    
    #***************************** PLOT MSL ************************************
    
    def view_msl_discrete_single(self):
        plt_msl_discrete_single(self)               #Plot MSL discrete for a single cell
        
    
    def view_msl_continuous_single(self):           #Plot MSL drift for a single cell
        plot_msl_continuous_single(self)
    
    
    def view_msl_continuous_multi(self):            #Plot MSL drift for multiple cells
        plot_msl_continuous_multi_unit(self)
    
    
    def view_msl_drift(self):
        cell_msl_drift(self)

        
    def view_allcell_drift(self):
        all_cell_msl_stats(self)       
        
        
    def view_drift_trends(self):
        drift_trends(self)
    
    
    def view_drift_movie(self):
        drift_movies(self)


    #*************************************************************************
    #***************************** STATISTICS FUNCTIONS **********************
    #*************************************************************************
    def sua_lockpercentage(self):
        sua_lock_percentage(self)

            
    def compute_msl_Pvals(self):
        compute_msl_pvals(self)                 #Compute P-values from 2-value KS test


    def view_pval(self):
        view_msl_pval(self)                     #View P-Values



    #*************************************************************************
    #********************* CSD ANALYSYS AND OTHER STIMULI ********************
    #*************************************************************************
    def view_csd_rasters(self):
        compute_csd_rasters(self)


    def view_csd_histogram(self):
        compute_csd_histogram(self)


    def view_triplet_sequences(self):
        compute_triplet_sequences(self)


    def view_nspikes_histograms(self):
        compute_nspikes_histograms(self)


    def view_isi_histograms(self):
        compute_isi_histograms(self)


    #********************** NAT SCENE ANALYSIS ************************
    def view_natscene_rasters(self):
        compute_natscene_rasters(self)


    def view_natscene_rasters(self):
        compute_natscene_rasters(self)


    def view_natscene_seqisi(self):
        compute_natscene_seqisi(self)


    def view_natscene_pairisi(self):
        compute_natscene_pairisi(self)
