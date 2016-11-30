from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *



class EventTriggeredImaging(QtGui.QWidget):
    def __init__(self, parent):
        super(EventTriggeredImaging, self).__init__(parent)
        self.parent = parent
        layout = QtGui.QGridLayout()    


        #Default: july 11, 2016 experiment
        #self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/' 
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        self.parent.n_sec = 3
        
        row_index = 0   #Keep track of button/box row
        
        #**************************************************************************************
        #******************************** SELECT ANIMAL & SESSION *****************************
        #**************************************************************************************
        self.vid_analysis_lbl = QLabel('EVENT TRIGGERED IMAGING', self)
        self.vid_analysis_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.vid_analysis_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.vid_analysis_lbl, row_index, 0); row_index+=1
        
        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording_event)
        layout.addWidget(self.button_select_recording, row_index, 0)
        
        #self.parent.selected_recording  = os.getcwd()
        self.selected_recording  = self.parent.root_dir
        self.select_recording_lbl = QLabel(self.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1)

        self.button_bin_to_npy = QPushButton('Convert .bin -> .npy')
        self.button_bin_to_npy.setMaximumWidth(200)
        self.button_bin_to_npy.clicked.connect(self.bn_to_npy)
        layout.addWidget(self.button_bin_to_npy, row_index, 5)

        self.n_pixels = QLineEdit('128')
        self.n_pixels.setMaximumWidth(50)
        n_pixels_lbl = QLabel('No. Pixels:', self)
        layout.addWidget(n_pixels_lbl, row_index, 6)
        layout.addWidget(self.n_pixels, row_index, 7); row_index+=1
        
        
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
        self.button2.clicked.connect(self.fltr_npy_event)
        layout.addWidget(self.button2, row_index, 0)
        

        self.filter_list = ['nofilter', 'butterworth', 'chebyshev']
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
        #************************************ COMPUTE DFF *************************************
        #**************************************************************************************
        for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
        row_index+=1
                
        self.preprocess_lbl = QLabel('DFF-COMPUTATION', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.button3 = QPushButton('Compute DFF Pipe ====>')
        self.button3.setMaximumWidth(200)
        self.button3.clicked.connect(self.dff_compute)
        layout.addWidget(self.button3, row_index, 0)

        self.comboBox_select_dff_filter = QtGui.QComboBox(self)
        for filter_ in self.filter_list: self.comboBox_select_dff_filter.addItem(filter_)
        layout.addWidget(self.comboBox_select_dff_filter, row_index,1)
        self.comboBox_select_dff_filter.activated[str].connect(self.select_dff_filter); self.selected_dff_filter = "nofilter" #Set default

        dff_list = ['globalAverage', 'slidingWindow']
        self.comboBox_select_dff_method = QtGui.QComboBox(self)
        for dff_ in dff_list: self.comboBox_select_dff_method.addItem(dff_)
        layout.addWidget(self.comboBox_select_dff_method, row_index,2)
        self.comboBox_select_dff_method.activated[str].connect(self.select_dff_method); self.dff_method = "globalAverage"
        
        
        self.n_sec_window = QLineEdit('3.0')
        self.n_sec_window.setMaximumWidth(50)
        self.n_sec_window_lbl = QLabel('#Sec Window:', self)
        layout.addWidget(self.n_sec_window_lbl, row_index,4)
        layout.addWidget(self.n_sec_window, row_index,5)
        
        self.control_lbl = QLabel('Control', self)
        layout.addWidget(self.control_lbl, row_index, 6)
        
        self.comboBox_select_control = QtGui.QComboBox(self)
        control_txt = ["no", "yes"]
        for ctxt in control_txt: self.comboBox_select_control.addItem(ctxt)
        layout.addWidget(self.comboBox_select_control, row_index,7)
        self.comboBox_select_control.activated[str].connect(self.select_control); self.selected_control = "no" #Set default
       
        row_index+=1
        
        self.button31 = QPushButton('Static STM')
        self.button31.setMaximumWidth(200)
        self.button31.clicked.connect(self.static_stm_event_trigger)
        layout.addWidget(self.button31, row_index, 0)
        
        self.button32 = QPushButton('Video STM')
        self.button32.setMaximumWidth(200)
        self.button32.clicked.connect(self.video_stm_mouse_lever)
        layout.addWidget(self.button32, row_index, 1)        
        

        self.window_start = QLineEdit('-3.0')
        self.window_start.setMaximumWidth(50)
        layout.addWidget(self.window_start, row_index, 2)
        
        self.window_end = QLineEdit('+3.0')
        self.window_end.setMaximumWidth(50)
        layout.addWidget(self.window_end, row_index, 3)
        

        self.block_save = QLineEdit('10')
        self.block_save.setMaximumWidth(50)
        self.block_save_lbl = QLabel('Block Ave:', self)
        layout.addWidget(self.block_save_lbl, row_index,4)
        layout.addWidget(self.block_save, row_index,5)

        self.midline_mask = QLineEdit('5')
        self.midline_mask.setMaximumWidth(50)
        self.midline_mask_lbl = QLabel('Mid-Mask Pixels:', self)
        layout.addWidget(self.midline_mask_lbl, row_index,6)
        layout.addWidget(self.midline_mask, row_index,7)
        
        
        self.vmax_default = QLineEdit('0.0')
        self.vmax_default.setMaximumWidth(50)
        layout.addWidget(self.vmax_default, row_index,8)
        self.vmin_default = QLineEdit('0.0')
        self.vmin_default.setMaximumWidth(50)
        layout.addWidget(self.vmin_default, row_index,9); row_index+=1
        
        
        #**************************************************************************************
        #************************************ VIEW ACTIVITY ***********************************
        #**************************************************************************************


        self.button_view_activity = QPushButton('View Spontaneous Activity')
        self.button_view_activity.setMaximumWidth(200)
        self.button_view_activity.clicked.connect(self.vw_activity)
        layout.addWidget(self.button_view_activity, row_index, 0)

        self.starting_frame = QLineEdit('500')
        self.starting_frame.setMaximumWidth(100)
        self.starting_frame_lbl = QLabel('Starting Frame:', self)
        layout.addWidget(self.starting_frame_lbl, row_index, 1)
        layout.addWidget(self.starting_frame, row_index, 2)
        
        self.number_frame = QLineEdit('1000')
        self.number_frame.setMaximumWidth(100)
        self.number_frame_lbl = QLabel('No. of Frames:', self)
        layout.addWidget(self.number_frame_lbl, row_index, 3)
        layout.addWidget(self.number_frame, row_index, 4); row_index+=1
        
        
        for k in range(2): 
            layout.addWidget(QLabel(' '*40, self), row_index,0); row_index+=1
        
        
        #**************************************************************************************
        #**************************** DIMENSIONALITY REDUCTION ********************************
        #**************************************************************************************
        #Select animal
        self.select_lbl = QLabel('STATE SPACE ANALYSIS', self)
        self.select_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold))
        self.select_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.select_lbl, row_index,0); row_index+=1
        
        self.button_state_space = QPushButton('State Space Distributions')
        self.button_state_space.setMaximumWidth(200)
        self.button_state_space.clicked.connect(self.st_space)
        layout.addWidget(self.button_state_space, row_index, 0)

        self.dim_red_list = ['PCA', 'MDS', 'tSNE', 'tSNE_Barnes_Hut']
        self.comboBox_dim_red = QtGui.QComboBox(self)
        for dim_red in self.dim_red_list: 
            self.comboBox_dim_red.addItem(dim_red)
        layout.addWidget(self.comboBox_dim_red, row_index,1)
        self.comboBox_dim_red.activated[str].connect(self.select_dim_red); self.selected_dim_red = self.dim_red_list[0] #Set default
        

        self.button_plot_distribution = QPushButton('Plot Distribution')
        self.button_plot_distribution.setMaximumWidth(200)
        self.button_plot_distribution.clicked.connect(self.plt_distribution)
        layout.addWidget(self.button_plot_distribution, row_index, 2)
        
        self.scaling_factor = QLineEdit('1.0')
        self.scaling_factor.setMaximumWidth(100)
        self.scaling_factor_lbl = QLabel('Scaling Factor:', self)
        layout.addWidget(self.scaling_factor_lbl, row_index, 3)
        layout.addWidget(self.scaling_factor, row_index, 4); row_index+=1
        
        
        #**************************************************************************************
        #******************************* PICKING POINTS  **************************************
        #**************************************************************************************
                
        self.button_list_points = QPushButton('List Selected Points')
        self.button_list_points.setMaximumWidth(200)
        self.button_list_points.clicked.connect(self.list_points)
        layout.addWidget(self.button_list_points, row_index, 0)
        
        self.button_clear_points = QPushButton('Clear Selected Points')
        self.button_clear_points.setMaximumWidth(200)
        self.button_clear_points.clicked.connect(self.clr_points)
        layout.addWidget(self.button_clear_points, row_index, 1)

        self.button_view_ave_points = QPushButton('View Ave of Selected Points')
        self.button_view_ave_points.setMaximumWidth(250)
        self.button_view_ave_points.clicked.connect(self.vw_ave_points)
        layout.addWidget(self.button_view_ave_points, row_index, 2)


        self.button_view_all_points = QPushButton('View All Selected Points')
        self.button_view_all_points.setMaximumWidth(250)
        self.button_view_all_points.clicked.connect(self.vw_all_points)
        layout.addWidget(self.button_view_all_points, row_index, 3); row_index+=1
                
        
        #**************************************************************************************
        #******************************* PICKING POINTS  **************************************
        #**************************************************************************************
        
        self.button_video_points = QPushButton('Make Video From 1st Point')
        self.button_video_points.setMaximumWidth(250)
        self.button_video_points.clicked.connect(self.vid_points)
        layout.addWidget(self.button_video_points, row_index, 0)


        
        #*************************************************************************
        row_index+=1
        for o in range(4):
            for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
            row_index+=1
        
        
        self.setLayout(layout)


    def slct_recording_event(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, ".bin, *.tif, .npy", self.selected_recording,"(*.bin *.tif *.npy)") 
        path_name, file_name = os.path.split(self.selected_recording)
        self.select_recording_lbl.setText(file_name)


    def slct_event_file(self):
        self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, ".ptcs or .txt)", self.selected_recording,"PTCS, TXT (*.ptcs *.txt)")
        path_name, file_name = os.path.split(self.selected_sort)
        self.select_sort_lbl.setText(file_name)

        #Reload session list and reset session box
        self.comboBox_selected_unit.clear()
        if '.ptcs' in self.selected_sort: 
            self.Sort = Ptcs(self.selected_sort)
            n_units = len(self.Sort.units)
            self.selected_unit_spikes_lbl.setText(str(len(self.Sort.units[0])))
        else: 
            n_units = 1
            self.text_Sort = np.loadtxt(self.selected_sort)
            self.selected_unit_spikes_lbl.setText(str(len(self.text_Sort)))

        for unit in range(n_units):
            self.comboBox_selected_unit.addItem(str(unit))
        self.slct_unit('0')


    def slct_unit(self, text):
        self.selected_unit = text

        if '.ptcs' in self.selected_sort: self.selected_unit_spikes_lbl.setText(str(len(self.Sort.units[int(text)])))
        else: self.selected_unit_spikes_lbl.setText(str(len(self.text_Sort)))


    def dff_compute(self):
        compute_dff_events(self)

    def static_stm_event_trigger(self):
        view_static_stm_events(self)


    def select_dff_filter(self,text):
        self.selected_dff_filter=text
            
    def select_dff_method(self, text):
        self.dff_method = text
    
    def select_control(self, text):
        self.selected_control = text


    def video_stm_mouse_lever(self):
        pass
        view_video_stm(self)
        
        
    def bn_to_npy(self):
        convert_bin_to_npy(self)
            
    def fltr_npy_event(self):
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
        plot_3D_distribution_new(self)      #Plot 3D points a
        
    
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
        

