import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from analysis import *

class CountMatrix(QtGui.QWidget):
    def __init__(self, parent):
        super(CountMatrix, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()

        self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        self.sua_file = '/media/cat/12TB/in_vivo/tim/cat/2017_04_26_barrel_ephys/track1/sort_alltrack_track1/track1_10ms_10sec_600sec_1_170426_152509_hp_butter_alltrack_goodUnits.ptcs'
        self.lfp_ptcs = '/media/cat/12TB/in_vivo/tim/cat/2017_04_26_barrel_ephys/track1/sort_alltrack_track1/track1_10ms_10sec_600sec_1_170426_152509_lfp_250hz_alltrack_50compressed_reduced_10chs_new_2units.ptcs'
        self.lfp_file = '/media/cat/12TB/in_vivo/tim/cat/2017_04_26_barrel_ephys/track1/sort_alltrack_track1/track1_10ms_10sec_600sec_1_170426_152509_hp_butter_alltrack_goodUnits_cell_rasters_lfp0.npy'
        
        #self.root_dir = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/'
        #self.sua_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_alltrack.ptcs'
        #self.lfp_ptcs = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_lp_alltrack_notch_50compressed_4Units.ptcs'
        #self.lfp_file = '/media/cat/8TB/in_vivo/nick/lfp_clustering/ptc21/tr5c/sort_alltrack/55-tr5c-csd_alltrack_cell_rasters_lfp0.npy'
        
        
        row_index=0
        #self.title = QLineEdit('');                #parent.start_time = self.start_time
        #self.title.setMaximumWidth(0)
        #self.title_lbl = QLabel('COUNT MATRIX TOOLS',self);row_index+=1
        ##self.title_len_lbl = QLabel(' Rec len: ' + str(int(parent.animal.rec_length))+ " (sec)", self)

        self.preprocess_lbl = QLabel('COUNT MATRIX TOOLS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        
        
        #Set root dir; button, then label
        self.button_set_sua_file = QPushButton('Single Unit File')
        self.button_set_sua_file.setMaximumWidth(200)
        self.button_set_sua_file.clicked.connect(self.set_sua_file)
        layout.addWidget(self.button_set_sua_file, row_index, 0)
        self.button_set_sua_file_lbl = QLabel(os.path.split(self.sua_file)[1], self)
        layout.addWidget(self.button_set_sua_file_lbl, row_index, row_index, 1, 9); row_index+=1
        
        
        #Set root dir; button, then label
        self.button_set_lfp_file = QPushButton('LFP Sort')
        self.button_set_lfp_file.setMaximumWidth(200)
        self.button_set_lfp_file.clicked.connect(self.set_lfp_ptcs)
        layout.addWidget(self.button_set_lfp_file, row_index, 0)
        self.button_set_lfp_file_lbl = QLabel(os.path.split(self.lfp_ptcs)[1], self)
        layout.addWidget(self.button_set_lfp_file_lbl, row_index, row_index, 1, 9); row_index+=1
        
        
        #Set root dir; button, then label
        self.button_set_raster_file = QPushButton('LEC Raster File')
        self.button_set_raster_file.setMaximumWidth(200)
        self.button_set_raster_file.clicked.connect(self.set_lfp_file)
        layout.addWidget(self.button_set_raster_file, row_index, 0)
        self.button_set_raster_file_lbl = QLabel(os.path.split(self.lfp_file)[1], self)
        layout.addWidget(self.button_set_raster_file_lbl, row_index, row_index, 1, 9); row_index+=1
        
        
        self.interpolation = QLineEdit('none');                   #parent.block_save = self.block_save
        self.interpolation.setMaximumWidth(50)
        interpolation_lbl = QLabel('interp:', self.interpolation)
        interpolation_lbl.setMaximumWidth(100)
        layout.addWidget(interpolation_lbl, row_index,3)
        layout.addWidget(self.interpolation, row_index,4)    

        self.zoom = QLineEdit('300');                   #parent.block_save = self.block_save
        self.zoom.setMaximumWidth(50)
        zoom_lbl = QLabel('zoom (ms):', self.zoom)
        zoom_lbl.setMaximumWidth(100)
        layout.addWidget(zoom_lbl, row_index,6)
        layout.addWidget(self.zoom, row_index,7)
        
        self.bin_len = QLineEdit('10');                   #parent.block_save = self.block_save
        self.bin_len.setMaximumWidth(50)
        bin_len_lbl = QLabel('Bin Length:', self.zoom)
        bin_len_lbl.setMaximumWidth(100)
        layout.addWidget(bin_len_lbl, row_index,8)
        layout.addWidget(self.bin_len, row_index,9);row_index+=1
        

        #MAKE BUTTONS
        self.button4 = QPushButton('LEC Count Matrix')
        self.button4.setMaximumWidth(200)
        layout.addWidget(self.button4, row_index, 0)
        self.button4.clicked.connect(self.view_count_matrix)
        
        self.start_cell = QLineEdit('0');                   #parent.start_cell = self.start_cell
        self.start_cell.setMaximumWidth(50)
        layout.addWidget(self.start_cell, row_index,2)

        start_cell_lbl = QLabel('1st_cell:', self)
        start_cell_lbl.setMaximumWidth(60)
        layout.addWidget(start_cell_lbl, row_index,1)

        self.end_cell = QLineEdit('1');                     #parent.end_cell = self.end_cell
        self.end_cell.setMaximumWidth(50)
        layout.addWidget(self.end_cell,row_index,4)

        end_cell_lbl = QLabel('2nd_cell:', self)
        end_cell_lbl.setMaximumWidth(60)
        layout.addWidget(end_cell_lbl,row_index,3)

        self.other_cell = QLineEdit('2');                     #parent.end_cell = self.end_cell
        self.other_cell.setMaximumWidth(50)
        layout.addWidget(self.other_cell,row_index,6)
        other_cell_lbl = QLabel('3rd_cell:', self)
        other_cell_lbl.setMaximumWidth(60)
        layout.addWidget(other_cell_lbl,row_index,5); 
        
        self.lfp_selected = QLineEdit('0');                     #parent.end_cell = self.end_cell
        self.lfp_selected.setMaximumWidth(50)
        layout.addWidget(self.lfp_selected,row_index,8)
        lfp_selected_lbl = QLabel('LEC #:', self)
        lfp_selected_lbl.setMaximumWidth(60)
        layout.addWidget(lfp_selected_lbl,row_index,7); row_index+=1

        self.button4 = QPushButton('AllSpike Count Matrix')
        self.button4.setMaximumWidth(200)
        layout.addWidget(self.button4, row_index, 0)
        self.button4.clicked.connect(self.view_count_matrix_allspikes); row_index+=1
        
        self.button41 = QPushButton('View CMatrix Chunks')
        self.button41.setMaximumWidth(200)
        layout.addWidget(self.button41, row_index, 0)
        self.button41.clicked.connect(self.view_count_matrix_chunks)
        
        self.time_chunks = QLineEdit('4');               
        self.time_chunks.setMaximumWidth(50)
        time_chunks_lbl = QLabel('# of time_chunks:', self)
        layout.addWidget(time_chunks_lbl, row_index,1)
        layout.addWidget(self.time_chunks, row_index,2); row_index+=1
        

        self.button5 = QPushButton('Compute CMatrix All')
        self.button5.setMaximumWidth(200)
        layout.addWidget(self.button5, row_index, 0)
        self.button5.clicked.connect(self.compute_all_cell_count_matrix); row_index+=1
        
        self.button6 = QPushButton('CMatrix to .png')
        self.button6.setMaximumWidth(200)
        layout.addWidget(self.button6, row_index, 0)
        self.button6.clicked.connect(self.view_all_count_matrix)

        self.isi_binning = QLineEdit('10');                   #parent.block_save = self.block_save
        self.isi_binning.setMaximumWidth(50)
        isi_binning_lbl = QLabel('isi_binning (ms):', self.isi_binning)
        isi_binning_lbl.setMaximumWidth(110)
        layout.addWidget(isi_binning_lbl, row_index,1)
        layout.addWidget(self.isi_binning, row_index,2)        

        self.cmatrix_nspikes = QLineEdit('1000');                   #parent.block_save = self.block_save
        self.cmatrix_nspikes.setMaximumWidth(50)
        cmatrix_nspikes_lbl = QLabel('min_#spks:', self.cmatrix_nspikes)
        
        cmatrix_nspikes_lbl.setMaximumWidth(100)
        layout.addWidget(cmatrix_nspikes_lbl, row_index,3)
        layout.addWidget(self.cmatrix_nspikes, row_index,4)
        

        self.setLayout(layout)

    def set_sua_file(self):
        self.sua_file = QtGui.QFileDialog.getOpenFileName(self, "*.ptcs (*.ptcs)", self.root_dir, "*.ptcs")
        self.button_set_sua_file_lbl.setText(self.sua_file.replace(os.path.dirname(self.sua_file),''))
        self.root_dir = os.path.split(self.sua_file)[0]

    def set_lfp_ptcs(self):
        self.lfp_ptcs = QtGui.QFileDialog.getOpenFileName(self, "*.ptcs (*.ptcs)", self.root_dir, "*.ptcs")
        self.button_set_lfp_file_lbl.setText(self.sua_file.replace(os.path.dirname(self.sua_file),''))
        self.root_dir = os.path.split(self.sua_file)[0]

    def set_lfp_file(self):
        self.lfp_file = QtGui.QFileDialog.getOpenFileName(self, "*.npy (*.npy)", self.root_dir, "*.npy")
        self.button_set_lfp_file_lbl.setText(self.lfp_file.replace(os.path.dirname(self.lfp_file),''))
        self.root_dir = os.path.split(self.sua_file)[0]
        

    def view_count_matrix(self):
        cell_count_matrix_LEC(self)

    def view_count_matrix_allspikes(self):
        cell_count_matrix_allspikes(self)
        

    def view_count_matrix_chunks(self):
        cell_count_matrix_chunks(self)
        
    def compute_all_cell_count_matrix(self):
        all_cell_count_matrix(self)

    def view_all_count_matrix(self):
        view_all_cell_count_matrix(self)
        
