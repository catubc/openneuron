import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class CountMatrix(QtGui.QWidget):
    def __init__(self, parent):
        super(CountMatrix, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()
        
        self.title = QLineEdit('');                #parent.start_time = self.start_time
        self.title.setMaximumWidth(0)
        self.title_lbl = QLabel('COUNT MATRIX TOOLS',self)
        self.title_len_lbl = QLabel(' Rec len: ' + str(int(parent.animal.rec_length))+ " (sec)", self)
        
        parent.interpolation = QLineEdit('none');                   #parent.block_save = self.block_save
        parent.interpolation.setMaximumWidth(50)
        interpolation_lbl = QLabel('interp:', parent.interpolation)
        interpolation_lbl.setMaximumWidth(100)

        parent.zoom = QLineEdit('300');                   #parent.block_save = self.block_save
        parent.zoom.setMaximumWidth(50)
        zoom_lbl = QLabel('zoom (ms):', parent.zoom)
        zoom_lbl.setMaximumWidth(100)

        parent.start_cell = QLineEdit('0');                   #parent.start_cell = self.start_cell
        parent.start_cell.setMaximumWidth(50)
        start_cell_lbl = QLabel('1st_cell:', self)
        start_cell_lbl.setMaximumWidth(60)
        
        parent.end_cell = QLineEdit('1');                     #parent.end_cell = self.end_cell
        parent.end_cell.setMaximumWidth(50)
        end_cell_lbl = QLabel('2nd_cell:', self)
        end_cell_lbl.setMaximumWidth(60)

        parent.other_cell = QLineEdit('2');                     #parent.end_cell = self.end_cell
        parent.other_cell.setMaximumWidth(50)
        other_cell_lbl = QLabel('3rd_cell:', self)
        other_cell_lbl.setMaximumWidth(60)

        parent.time_chunks = QLineEdit('4');               
        parent.time_chunks.setMaximumWidth(50)
        time_chunks_lbl = QLabel('# of time_chunks:', self)
        
        parent.isi_binning = QLineEdit('10');                   #parent.block_save = self.block_save
        parent.isi_binning.setMaximumWidth(50)
        isi_binning_lbl = QLabel('isi_binning (ms):', parent.isi_binning)
        isi_binning_lbl.setMaximumWidth(110)

        parent.cmatrix_nspikes = QLineEdit('1000');                   #parent.block_save = self.block_save
        parent.cmatrix_nspikes.setMaximumWidth(50)
        cmatrix_nspikes_lbl = QLabel('min_#spks:', parent.cmatrix_nspikes)
        cmatrix_nspikes_lbl.setMaximumWidth(100)



        #ADD TO LAYOUT
        layout.addWidget(self.title_lbl, 0,0)
        layout.addWidget(self.title, 0,1)
        layout.addWidget(self.title_len_lbl, 0,2)
        layout.addWidget(interpolation_lbl, 0,3)
        layout.addWidget(parent.interpolation, 0,4)    
        layout.addWidget(zoom_lbl, 0,6)
        layout.addWidget(parent.zoom, 0,7)
                        
        layout.addWidget(start_cell_lbl, 2,1)
        layout.addWidget(parent.start_cell, 2,2)
        layout.addWidget(end_cell_lbl,2,3)
        layout.addWidget(parent.end_cell,2,4)
        layout.addWidget(other_cell_lbl,2,5)
        layout.addWidget(parent.other_cell,2,6)

        layout.addWidget(time_chunks_lbl,3,1)
        layout.addWidget(parent.time_chunks,3,2)
        
        layout.addWidget(isi_binning_lbl, 7,1)
        layout.addWidget(parent.isi_binning, 7,2)        
        layout.addWidget(cmatrix_nspikes_lbl, 7,3)
        layout.addWidget(parent.cmatrix_nspikes, 7,4)

                 
        
        #MAKE BUTTONS
        self.button4 = QPushButton('View CMatrix Cell')
        self.button4.setMaximumWidth(200)
        layout.addWidget(self.button4, 2, 0)
        self.button4.clicked.connect(self.view_count_matrix)

        self.button41 = QPushButton('View CMatrix Chunks')
        self.button41.setMaximumWidth(200)
        layout.addWidget(self.button41, 3, 0)
        self.button41.clicked.connect(self.view_count_matrix_chunks)

        self.button5 = QPushButton('Compute CMatrix All')
        self.button5.setMaximumWidth(200)
        layout.addWidget(self.button5, 6, 0)
        self.button5.clicked.connect(self.compute_all_cell_count_matrix)

        self.button6 = QPushButton('CMatrix to .png')
        self.button6.setMaximumWidth(200)
        layout.addWidget(self.button6, 7, 0)
        self.button6.clicked.connect(self.view_all_count_matrix)


        self.setLayout(layout)

    def view_count_matrix(self):
        cell_count_matrix(self.parent)

    def view_count_matrix_chunks(self):
        cell_count_matrix_chunks(self.parent)
        
    def compute_all_cell_count_matrix(self):
        all_cell_count_matrix(self.parent)

    def view_all_count_matrix(self):
        view_all_cell_count_matrix(self.parent)
        
