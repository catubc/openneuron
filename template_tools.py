import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *


class TemplateTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TemplateTools, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        
        row_index = 0
        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.selected_recording  = '/media/cat/8TB/in_vivo/nick'
        self.select_recording_lbl = QLabel(self.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1); row_index+=1

        #Select recording
        self.button_select_sort = QPushButton('Select Sort')
        self.button_select_sort.setMaximumWidth(200)
        self.button_select_sort.clicked.connect(self.slct_sort)
        layout.addWidget(self.button_select_sort, row_index, 0)
        
        self.parent.selected_sort  = '/media/cat/12TB/in_vivo/'
        self.select_sort_lbl = QLabel(self.parent.selected_sort, self)
        layout.addWidget(self.select_sort_lbl, row_index,1)
        
        #Color picker
        self.comboBox3 = QtGui.QComboBox(self)
        self.comboBox3.addItem("blue")
        self.comboBox3.addItem("red")
        self.comboBox3.addItem("green")
        self.comboBox3.addItem("magenta")
        self.comboBox3.addItem("cyan")
        self.comboBox3.addItem("brown")
        self.comboBox3.addItem("pink")
        self.comboBox3.addItem("orange")
        
        layout.addWidget(self.comboBox3, row_index,6)
        self.comboBox3.activated[str].connect(self.slct_colour); self.selected_colour="blue"

        #Prefilter templates
        self.low_cutoff = QLineEdit('0.0');                #parent.start_time = self.start_time
        self.low_cutoff.setMaximumWidth(50)
        self.low_cutoff_lbl = QLabel('Lowcut Filter', self)
        self.low_cutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.low_cutoff_lbl, row_index,7)
        layout.addWidget(self.low_cutoff, row_index, 8); row_index+=1
        
        
        
        #View traces
        self.button_view_templates = QPushButton('View Templates')
        self.button_view_templates.setMaximumWidth(200)
        self.button_view_templates.clicked.connect(self.vw_templates)
        layout.addWidget(self.button_view_templates, row_index, 0)

        self.selected_unit = QLineEdit('0');                #parent.start_time = self.start_time
        self.selected_unit.setMaximumWidth(50)
        self.selected_unit_lbl = QLabel('Unit:', self)
        self.selected_unit_lbl.setMaximumWidth(100)
        layout.addWidget(self.selected_unit_lbl, row_index,1)
        layout.addWidget(self.selected_unit, row_index, 2)

        self.n_electrodes = QLineEdit('1.0');                #parent.start_time = self.start_time
        self.n_electrodes.setMaximumWidth(50)
        self.n_electrodes_lbl = QLabel('%Electrodes', self)
        self.n_electrodes_lbl.setMaximumWidth(100)
        layout.addWidget(self.n_electrodes_lbl, row_index,3)
        layout.addWidget(self.n_electrodes, row_index, 4)
        
        self.n_sample_pts = QLineEdit('20');                #parent.start_time = self.start_time
        self.n_sample_pts.setMaximumWidth(50)
        self.n_sample_pts_lbl = QLabel('#Sample Pts:', self)
        self.n_sample_pts_lbl.setMaximumWidth(100)
        layout.addWidget(self.n_sample_pts_lbl, row_index,5)
        layout.addWidget(self.n_sample_pts, row_index, 6)
        
        self.voltage_scale = QLineEdit('10.0');                  #parent.end_time = self.end_time
        self.voltage_scale.setMaximumWidth(50)
        self.voltage_scale_lbl = QLabel('Voltage Scaling', self)
        self.voltage_scale_lbl.setMaximumWidth(150)
        layout.addWidget(self.voltage_scale_lbl, row_index, 7)
        layout.addWidget(self.voltage_scale, row_index,8); row_index+=1
        
        self.button_view_csd = QPushButton('View CSD')
        self.button_view_csd.setMaximumWidth(200)
        self.button_view_csd.clicked.connect(self.vw_csd)
        layout.addWidget(self.button_view_csd, row_index, 0)
        
        
        self.start_ch = QLineEdit('0');               
        self.start_ch.setMaximumWidth(50)
        self.start_ch_lbl = QLabel('Starting ch:', self)
        self.start_ch_lbl.setMaximumWidth(100)
        layout.addWidget(self.start_ch_lbl, row_index,1)
        layout.addWidget(self.start_ch, row_index, 2)
        
        
        self.end_ch = QLineEdit('32');                
        self.end_ch.setMaximumWidth(50)
        self.end_ch_lbl = QLabel('Ending ch:', self)
        self.end_ch_lbl.setMaximumWidth(100)
        layout.addWidget(self.end_ch_lbl, row_index,3)
        layout.addWidget(self.end_ch, row_index, 4); row_index+=1
        
        
        self.button_view_all_csd = QPushButton('View all CSD')
        self.button_view_all_csd.setMaximumWidth(200)
        self.button_view_all_csd.clicked.connect(self.vw_all_csd)
        layout.addWidget(self.button_view_all_csd, row_index, 0)
        
        
        
        self.snr_value = QLineEdit('3.0');                
        self.snr_value.setMaximumWidth(50)
        self.snr_value_lbl = QLabel('SNR factor:', self)
        self.snr_value_lbl.setMaximumWidth(100)
        layout.addWidget(self.snr_value_lbl, row_index,3)
        layout.addWidget(self.snr_value, row_index, 4); row_index+=1
        
        self.setLayout(layout)

    def slct_recording(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.selected_recording,"TSF (*.tsf)") 
        path_name, file_name = os.path.split(self.selected_recording)
        self.select_recording_lbl.setText(file_name)

        self.tsf = Tsf_file(self.selected_recording)
        self.tsf.read_ec_traces()

  
    def slct_sort(self):
        self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.selected_recording,"PTCS (*.ptcs)")
        path_name, file_name = os.path.split(self.selected_sort)
        self.select_sort_lbl.setText(file_name)
   
    def vw_templates(self):
        view_templates(self)
    
    def vw_csd(self):
        view_csd(self)

    def vw_all_csd(self):
        view_all_csd(self)


    def slct_colour(self, text):
        self.selected_colour = text