from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *


class EventTriggeredEphys(QtGui.QWidget):
    def __init__(self, parent):
        super(EventTriggeredEphys, self).__init__(parent)
        self.parent = parent
        layout = QtGui.QGridLayout()    


        #Default: july 11, 2016 experiment
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/' 
        #self.parent.root_dir = '/media/cat/500GB/in_vivo/tim/cat/'
        self.parent.n_sec = 3
        
        row_index = 0   #Keep track of button/box row
        
        #**************************************************************************************
        #******************************** SELECT ANIMAL & SESSION *****************************
        #**************************************************************************************
        self.vid_analysis_lbl = QLabel('EVENT TRIGGERED EPHYS', self)
        self.vid_analysis_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.vid_analysis_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.vid_analysis_lbl, row_index, 0); row_index+=1
        
        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording1)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.selected_recording  = self.parent.root_dir
        self.select_recording_lbl = QLabel(self.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1); row_index+=1

        #Select recording
        self.button_select_sort = QPushButton('Select Event File')
        self.button_select_sort.setMaximumWidth(200)
        self.button_select_sort.clicked.connect(self.slct_event_file)
        layout.addWidget(self.button_select_sort, row_index, 0)
        
        self.selected_sort  = self.parent.root_dir
        self.select_sort_lbl = QLabel(self.selected_sort, self)
        layout.addWidget(self.select_sort_lbl, row_index,1)
    

        selected_unit_lbl = QLabel('Unit #:', self)
        layout.addWidget(selected_unit_lbl, row_index, 5)
        
        #self.filter_list = ['nofilter', 'butterworth', 'chebyshev']
        self.comboBox_selected_unit = QtGui.QComboBox(self)
        self.comboBox_selected_unit.addItem('')
        layout.addWidget(self.comboBox_selected_unit, row_index, 6)
        self.comboBox_selected_unit.activated[str].connect(self.slct_unit); self.selected_filter = "" #Set default

        self.selected_unit_spikes_lbl = QLabel('', self)
        layout.addWidget(self.selected_unit_spikes_lbl, row_index, 7)


        for k in range(2): layout.addWidget(QLabel(' '*40, self), row_index,0); row_index+=1

        #**************************************************************************************
        #***************************** PRE-PROCESSING HEADER **********************************
        #**************************************************************************************

        self.preprocess_lbl = QLabel('PRE-PROCESSING', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1
        

        self.button2 = QPushButton('Filter ==>')
        self.button2.setMaximumWidth(200)
        self.button2.clicked.connect(self.fltr_npy)
        layout.addWidget(self.button2, row_index, 0)
        

        self.filter_list = ['butterworth', 'chebyshev']
        self.comboBox_filter = QtGui.QComboBox(self)
        for filter_ in self.filter_list[1:]: self.comboBox_filter.addItem(filter_)
        layout.addWidget(self.comboBox_filter, row_index, 1)
        self.comboBox_filter.activated[str].connect(self.slct_filter); self.selected_filter = "butterworth" #Set default


        self.lowcut = QLineEdit('0.1')
        self.lowcut.setMaximumWidth(50)
        lowcut_lbl = QLabel('Low Cutoff (HZ):', self)
        layout.addWidget(lowcut_lbl, row_index, 2)
        layout.addWidget(self.lowcut, row_index, 3)

       
        self.highcut = QLineEdit('6.0')
        self.highcut.setMaximumWidth(50)
        highcut_lbl = QLabel('High Cutoff (HZ):', self)
        layout.addWidget(highcut_lbl, row_index, 4)
        layout.addWidget(self.highcut, row_index, 5); row_index+=1

        
        #**************************************************************************************
        #******************************** COMPUTE AVERAGES ***********************************
        #**************************************************************************************
        for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
        row_index+=1
                
        self.preprocess_lbl = QLabel('EVENT TRIGGERED AVERAGES', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        row_index+=1
        
        self.button_lfp_template = QPushButton('LFP Template')
        self.button_lfp_template.setMaximumWidth(200)
        self.button_lfp_template.clicked.connect(self.lfp_avg)
        layout.addWidget(self.button_lfp_template, row_index, 0)
        


        self.excluded_trials = QLineEdit('');                #parent.start_time = self.start_time
        self.excluded_trials.setMaximumWidth(100)
        self.excluded_trials_lbl = QLabel('Excluded Trials:', self)
        self.excluded_trials_lbl.setMaximumWidth(100)
        layout.addWidget(self.excluded_trials_lbl, row_index,3)
        layout.addWidget(self.excluded_trials, row_index, 4)
        
                
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
        
                
        self.button31 = QPushButton('CSD')
        self.button31.setMaximumWidth(200)
        self.button31.clicked.connect(self.csd_avg)
        layout.addWidget(self.button31, row_index, 0)
        
        self.snr_value = QLineEdit('3.0');                
        self.snr_value.setMaximumWidth(50)
        self.snr_value_lbl = QLabel('SNR factor:', self)
        self.snr_value_lbl.setMaximumWidth(100)
        layout.addWidget(self.snr_value_lbl, row_index,3)
        layout.addWidget(self.snr_value, row_index, 4); row_index+=1
        

        self.button31 = QPushButton('Single Unit Histogram')
        self.button31.setMaximumWidth(200)
        self.button31.clicked.connect(self.static_stm_mouse_lever)
        layout.addWidget(self.button31, row_index, 0); row_index+=1



        #self.block_save = QLineEdit('10')
        #self.block_save.setMaximumWidth(50)
        #self.block_save_lbl = QLabel('Block Ave & Mid-Mask Pixels:', self)
        #layout.addWidget(self.block_save_lbl, row_index,4)
        #layout.addWidget(self.block_save, row_index,5)
        
        #self.midline_mask = QLineEdit('5')
        #self.midline_mask.setMaximumWidth(50)
        #layout.addWidget(self.midline_mask, row_index,6); row_index+=1
        
        
        
        #*************************************************************************
        row_index+=1
        for o in range(4):
            for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
            row_index+=1
        
        
        self.setLayout(layout)


    def slct_recording1(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, ".tsf", self.selected_recording,"(*.tsf)") 
        path_name, file_name = os.path.split(self.selected_recording)
        self.select_recording_lbl.setText(file_name)



    def slct_event_file(self):
        self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, ".ptcs or .txt or .npy)", self.selected_recording,"PTCS, TXT, NPY (*.ptcs *.txt *.npy)")
        path_name, file_name = os.path.split(self.selected_sort)
        self.select_sort_lbl.setText(file_name)

        #Reload session list and reset session box
        self.comboBox_selected_unit.clear()
        if '.ptcs' in self.selected_sort: 
            self.Sort = Ptcs(self.selected_sort)
            n_units = len(self.Sort.units)
            self.selected_unit_spikes_lbl.setText(str(len(self.Sort.units[0])))
            
        elif '.txt' in self.selected_sort: 
            n_units = 1
            self.triggers = np.loadtxt(self.selected_sort)
            self.selected_unit_spikes_lbl.setText(str(len(self.triggers)))
        
        elif '.npy' in self.selected_sort: 
            n_units = 1
            self.camera_pulses = np.load(self.selected_sort)
            
            parse_camera_pulses(self); print "...# of triggers: ", len(self.triggers)
            
            self.selected_unit_spikes_lbl.setText(str(len(self.triggers)))
            
            print self.triggers/25000.
            print self.triggers_length/25000.
            
        for unit in range(n_units):
            self.comboBox_selected_unit.addItem(str(unit))
        self.slct_unit('0')


    def slct_unit(self, text):
        self.selected_unit = text

        if '.ptcs' in self.selected_sort: self.selected_unit_spikes_lbl.setText(str(len(self.Sort.units[int(text)])))
        else: self.selected_unit_spikes_lbl.setText(str(len(self.triggers)))


    def lfp_avg(self):
        compute_lfp_triggered_template(self)


    def csd_avg(self):
        compute_csd_event_triggered(self)

    def select_dff_filter(self,text):
        self.selected_dff_filter=text
        
            
    def select_dff_method(self, text):
        self.dff_method = text
            

    def static_stm_mouse_lever(self):
        view_static_stm_events(self)


    def video_stm_mouse_lever(self):
        view_video_stm(self)
        
        
    def bn_to_npy(self):
        convert_bin_to_npy(self)
            
    def fltr_npy(self):
        filter_for_event_trigger_analysis(self)
        
        
    def slct_filter(self, text):
        self.selected_filter = text
        print self.selected_filter


    def vw_activity(self):
        view_spontaneous_activity(self)
        

    def select_dim_red(self, text):
        self.selected_dim_red = text
        print text


    def st_space(self):
        compute_dim_reduction(self)


    def plt_distribution(self):
        plot_3D_distribution_new(self)
        
    
    def vw_ave_points(self):
        view_ave_points(self)

    def vw_all_points(self):
        view_all_points(self)

    def list_points(self):
        print self.parent.glwindow.glWidget.points_selected

    def clr_points(self):
        self.parent.glwindow.glWidget.points_selected = []

    def vid_points(self):
        print "...making vid starting at frame: ", self.parent.glwindow.glWidget.points_selected[0]
        video_points(self)

