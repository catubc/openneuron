import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *


class TracesTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TracesTools, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat'
        
        row_index = 0

        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.parent.selected_recording  = self.parent.root_dir
        self.select_recording_lbl = QLabel(self.parent.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1)
        
        #Rec length
        self.parent.rec_length = QLabel('0', self)
        self.parent.rec_length.setMaximumWidth(100)
        self.parent.rec_length_lbl = QLabel('Rec Length:', self)
        self.parent.rec_length_lbl.setMaximumWidth(100)
        layout.addWidget(self.parent.rec_length_lbl, row_index, 3)
        layout.addWidget(self.parent.rec_length, row_index, 4); row_index+=1

        #Select recording
        self.button_select_sort = QPushButton('Select Event File')
        self.button_select_sort.setMaximumWidth(200)
        self.button_select_sort.clicked.connect(self.slct_event_file)
        layout.addWidget(self.button_select_sort, row_index, 0)
        
        self.selected_sort  = self.parent.root_dir
        self.select_sort_lbl = QLabel(self.selected_sort, self)
        layout.addWidget(self.select_sort_lbl, row_index,1)
    


        #% of electrodes
        self.n_electrodes = QLineEdit('1.0');                #parent.start_time = self.start_time
        self.n_electrodes.setMaximumWidth(50)
        self.n_electrodes_lbl = QLabel('%Electrodes', self)
        self.n_electrodes_lbl.setMaximumWidth(100)
        layout.addWidget(self.n_electrodes_lbl, row_index,5)
        layout.addWidget(self.n_electrodes, row_index, 6)
        
        self.low_cutoff = QLineEdit('10.0');                #parent.start_time = self.start_time
        self.low_cutoff.setMaximumWidth(50)
        self.low_cutoff_lbl = QLabel('Lowcut Filter', self)
        self.low_cutoff_lbl.setMaximumWidth(100)
        layout.addWidget(self.low_cutoff_lbl, row_index,7)
        layout.addWidget(self.low_cutoff, row_index, 8); row_index+=1
                
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
        self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, ".ptcs or .txt or .npy)", self.parent.selected_recording,"PTCS, TXT, NPY (*.ptcs *.txt *.npy)")
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
        self.parent.selected_recording =  QtGui.QFileDialog.getOpenFileName(self, 'Load File', self.parent.selected_recording)
        path_name, file_name = os.path.split(self.parent.selected_recording)

        self.select_recording_lbl.setText(file_name)

        self.tsf = Tsf_file(self.parent.selected_recording)
        self.tsf.read_ec_traces()
        print "...len rec: ", self.tsf.n_vd_samples/float(self.tsf.SampleFrequency)
        self.parent.rec_length.setText(str(self.tsf.n_vd_samples/float(self.tsf.SampleFrequency)))
   
   
    def vw_traces(self):
        view_traces(self)
