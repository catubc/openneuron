import glob, os, sys

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

sys.path.append('/home/cat/code/')
import TSF.TSF as TSF
import PTCS.PTCS as PTCS

from analysis import *


#LFP ANALYSIS TOOLBOX
class TrackWideTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TrackWideTools, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        #self.root_dir = '/media/cat/All.Data.3TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/'
        self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        #self.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering'

        row_index= 0
        #*****************************************************************

        self.preprocess_lbl = QLabel('CONVERT /CONCATENATE', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1


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

        self.button_tsf_to_dat = QPushButton('Concatenate multiple .tsf to .dat')
        self.button_tsf_to_dat.setMaximumWidth(250)
        self.button_tsf_to_dat.clicked.connect(self.multi_tsf_dat)
        layout.addWidget(self.button_tsf_to_dat, row_index, 0); row_index+=1
        
        self.button1 = QPushButton('Concatenate _lp.tsf files')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_lfp_zip)
        layout.addWidget(self.button1, row_index, 0); row_index+=1

        
        self.lfpcsd_button = QPushButton('Convert LFP->CSD')
        self.lfpcsd_button.setMaximumWidth(250)
        self.lfpcsd_button.clicked.connect(self.lfp_to_csd)
        layout.addWidget(self.lfpcsd_button, row_index, 0); row_index+=1
        
        #***********************************************************************************
        #***********************************************************************************
        #***********************************************************************************

        layout.addWidget(QLabel('', self), row_index,0); row_index+=1

        self.preprocess_lbl = QLabel('TIME COMPRESSION', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.button12 = QPushButton('Convert to LFP (XX Khz-> 1Khz)')
        self.button12.setMaximumWidth(250)
        self.button12.clicked.connect(self.subsample_data)
        layout.addWidget(self.button12, row_index, 0)
        
        #self.subsample_factor = QLineEdit('25');                  #parent.end_time = self.end_time
        #self.subsample_factor.setMaximumWidth(50)
        #self.subsample_factor_lbl = QLabel('Compress Factor (multiples of 25): ', self)
        #layout.addWidget(self.subsample_factor_lbl, row_index, 1)
        #layout.addWidget(self.subsample_factor, row_index, 2); 
        row_index+=1

        
        self.button1 = QPushButton('Time compress (1Khz) .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.comp_lfp)
        layout.addWidget(self.button1, row_index, 0)
        
        self.compress_factor = QLineEdit('50');                  #parent.end_time = self.end_time
        self.compress_factor.setMaximumWidth(50)
        self.compress_factor_lbl = QLabel('Compress Factor (multiples of 25): ', self)
        layout.addWidget(self.compress_factor_lbl, row_index, 1)
        layout.addWidget(self.compress_factor, row_index, 2); row_index+=1
        
        self.button_subsample_tsf = QPushButton('Clip + Subsample .tsf')
        self.button_subsample_tsf.setMaximumWidth(250)
        self.button_subsample_tsf.clicked.connect(self.clip_subsample_tsf)
        layout.addWidget(self.button_subsample_tsf, row_index, 0)

        self.top_channel = QLineEdit('0');                  #parent.end_time = self.end_time
        self.top_channel.setMaximumWidth(50)
        self.top_channel_lbl = QLabel('Top Ch:', self)
        layout.addWidget(self.top_channel_lbl, row_index, 1)
        layout.addWidget(self.top_channel, row_index, 2)
        
        self.bottom_channel = QLineEdit('63');                  #parent.end_time = self.end_time
        self.bottom_channel.setMaximumWidth(50)
        self.bottom_channel_lbl = QLabel('Top Ch:', self)
        layout.addWidget(self.bottom_channel_lbl, row_index, 3)
        layout.addWidget(self.bottom_channel, row_index, 4)
        
        self.total_channels = QLineEdit('10');                  #parent.end_time = self.end_time
        self.total_channels.setMaximumWidth(50)
        self.total_channels_lbl = QLabel('Total Chs:', self)
        layout.addWidget(self.total_channels_lbl, row_index, 5)
        layout.addWidget(self.total_channels, row_index, 6); row_index+=1
        
                
        self.zero_desynch_states = QPushButton('Zero Desynch Periods')
        self.zero_desynch_states.setMaximumWidth(250)
        self.zero_desynch_states.clicked.connect(self.zero_desynch)
        layout.addWidget(self.zero_desynch_states, row_index, 0)
        
        
        self.specgram_ch = QLineEdit('9');                  #parent.end_time = self.end_time
        self.specgram_ch.setMaximumWidth(50)
        self.specgram_ch_lbl = QLabel('Spec Channel: ', self)
        layout.addWidget(self.specgram_ch_lbl, row_index, 1)
        layout.addWidget(self.specgram_ch, row_index, 2)
            
        self.sync_limit = QLineEdit('0.4');                  #parent.end_time = self.end_time
        self.sync_limit.setMaximumWidth(50)
        self.sync_limit_lbl = QLabel('Synchrony Index Limit: ', self)
        layout.addWidget(self.sync_limit_lbl, row_index, 3)
        layout.addWidget(self.sync_limit, row_index, 4); row_index+=1    
        
        
        
        self.make_ptcs2 = QPushButton('Make Synch .ptcs2')
        self.make_ptcs2.setMaximumWidth(250)
        self.make_ptcs2.clicked.connect(self.do_sync_ptcs2)
        layout.addWidget(self.make_ptcs2, row_index, 0)

        
        #self.print_footer = QPushButton('Print Footer .tsf')
        #self.print_footer.setMaximumWidth(250)
        #self.print_footer.clicked.connect(self.view_footer)
        #layout.addWidget(self.print_footer, row_index, 0)
         
        #self.print_header= QPushButton('Print Header .tsf')
        #self.print_header.setMaximumWidth(250)
        #self.print_header.clicked.connect(self.view_header)
        #layout.addWidget(self.print_header, row_index, 1)
                
        #***************************************
        
        
        self.setLayout(layout)
    
    
    def tsftolfp(self):
        self.tsf_files = QtGui.QFileDialog.getOpenFileNames(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf_to_lfp(self)
    
    
    def lfpzip_to_lptsf(self):
        lfpzip_file = QtGui.QFileDialog.getOpenFileName(self, "lfp.zip (*lfp.zip, *.txt)", self.root_dir,"*lfp.zip or .txt (*lfp.zip; *.txt)")
        lfp_to_lptsf(lfpzip_file)
        

    def multi_tsf(self):
        print "... selecting multiple recordings ..."
        out_files = QtGui.QFileDialog.getOpenFileNames(self, "*", self.root_dir, "*")
        self.tsf_files = out_files
        
        print out_files[0]
        
        self.file_list = out_files[0]
        if '.txt' in out_files[0]:
            with open(out_files[0]) as f: 
                self.tsf_files = [line.rstrip('\n') for line in f]
        
        print "... concatenate multiple recordings ..."
        concatenate_tsf(self) 
    
    def multi_tsf_dat(self):
        print "... selecting multiple recordings ..."
        out_files = QtGui.QFileDialog.getOpenFileNames(self, "*", self.root_dir, "*")
        self.tsf_files = out_files
        
        print out_files[0]
        
        if '.txt' in out_files[0]:
            with open(out_files[0]) as f: 
                self.tsf_files = [line.rstrip('\n') for line in f]
        
        print "... concatenate multiple recordings ..."
        concatenate_tsf_to_dat(self)
        
    def multi_lfp_zip(self):
        print "... selecting multiple _lp.tsf recording directories ..."
        #dialog = FileDialog()   #**** SELECT MULTIPLE DIRECTORIES, NOT INDIVIDUAL FIELS
        #dialog.exec_()

        in_files = QtGui.QFileDialog.getOpenFileNames(self, "*", self.root_dir, "*")

        self.tsf_files = in_files
        concatenate_lfp_zip(self) 


    def subsample_save(self):
        print "... selecting single tsf recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        
        subsample_channels_tsf(self) 
    
    def lfp_to_csd(self):
        print "... selecting single tsf recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        
        convert_lfp_to_csd(self)

    def clip_subsample_tsf(self):
        print "... selecting single tsf recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        subsample_channels_tsf(self)
        
        
    def subsample_data(self):
        print "...subsampling data (converting 25Khz->1Khz LFP data)..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        subsample_lfp(self) 

        
    def comp_lfp(self):
        print "... selecting single lfp recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        compress_lfp(self) 
    
    
    def zero_desynch(self):

        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")

        zero_out_desynch_periods(self)

    def do_sync_ptcs2(self):

        print "...select .ptcs file..."
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.root_dir,"PTCS (*.ptcs)")

        print "...select .tsf file..."
        self.lfp_tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")

        make_sync_ptcs2(self)
    
    
    #def view_footer(self):
        
        #print "... selecting single lfp recording ..."
        #self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        #tsf = TSF.TSF(self.tsf_file)
        
        #tsf.read_footer()
    
    #def view_header(self):
        
        #print "... selecting single lfp recording ..."
        #self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        #tsf = TSF.TSF(self.tsf_file)
        
        #tsf.print_header()
    
