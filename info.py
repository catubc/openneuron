import glob, os, sys

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

sys.path.append('/home/cat/code/')
import TSF.TSF as TSF
import PTCS.PTCS as PTCS

from analysis import *


#LFP ANALYSIS TOOLBOX
class Info(QtGui.QWidget):
    def __init__(self, parent):
        super(Info, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        #self.root_dir = '/media/cat/All.Data.3TB/in_vivo/tim/cat/2016_05_27_gcamp/tsf_files/'
        self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        #self.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering'

        row_index= 0

        #**************************************************************************************
        #*********************************** TSF **************************************
        #**************************************************************************************

        self.preprocess_lbl = QLabel('TSF FILE INFO', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.print_footer = QPushButton('Print Footer .tsf')
        self.print_footer.setMaximumWidth(250)
        self.print_footer.clicked.connect(self.view_footer)
        layout.addWidget(self.print_footer, row_index, 0)
         
        self.print_header= QPushButton('Print Header .tsf')
        self.print_header.setMaximumWidth(250)
        self.print_header.clicked.connect(self.view_header)
        layout.addWidget(self.print_header, row_index, 1); row_index+=1
                
        #***************************************
        
        self.preprocess_lbl = QLabel('PTCS FILE INFO', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.print_header_ptcs= QPushButton('Print Header .ptcs')
        self.print_header_ptcs.setMaximumWidth(250)
        self.print_header_ptcs.clicked.connect(self.view_header_ptcs)
        layout.addWidget(self.print_header_ptcs, row_index, 0); row_index+=1
                
        
        self.setLayout(layout)
    
    
    def tsftolfp(self):
        self.tsf_files = QtGui.QFileDialog.getOpenFileNames(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf_to_lfp(self)
    
    
    def lfpzip_to_lptsf(self):
        lfpzip_file = QtGui.QFileDialog.getOpenFileName(self, "LFP ZIP (*.lfp.zip)", self.root_dir,"*.lfp.zip (*.lfp.zip)")
        lfp_to_lptsf(lfpzip_file)
        
    
    def multi_tsf(self):
        print "... selecting multiple recordings ..."
        out_files = QtGui.QFileDialog.getOpenFileNames(self, "*", self.root_dir, "*")
        self.tsf_files = out_files
        
        print out_files[0]
        
        if '.txt' in out_files[0]:
            with open(out_files[0]) as f: 
                self.tsf_files = [line.rstrip('\n') for line in f]
                    
        print "... concatenate multiple recordings ..."
        concatenate_tsf(self) 
    
    
    def multi_lfp_zip(self):
        print "... selecting multiple .lfp.zip recording directories ..."
        #dialog = FileDialog()   #**** SELECT MULTIPLE DIRECTORIES, NOT INDIVIDUAL FIELS
        #dialog.exec_()

        in_files = QtGui.QFileDialog.getOpenFileNames(self, "*", self.root_dir, "*")

        self.tsf_files = in_files
        concatenate_lfp_zip(self) 
            
    
    def comp_lfp(self):
        print "... selecting single lfp recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        compress_lfp(self) 
    
    
    
    def view_footer(self):
        
        print "... selecting single lfp recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf = TSF.TSF(self.tsf_file)
        
        tsf.read_footer()
    
    def view_header(self):
        
        print "... selecting single lfp recording ..."
        self.tsf_file = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.root_dir,"TSF (*.tsf)")
        tsf = TSF.TSF(self.tsf_file)
        
        tsf.print_header()
    
    def view_header_ptcs(self):
        
        print "... selecting single lfp recording ..."
        self.ptcs_file = QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.root_dir,"PTCS (*.ptcs)")
        sort = PTCS.PTCS(self.ptcs_file)
        
        sort.print_header()
    
