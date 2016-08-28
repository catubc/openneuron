import sip
sip.setapi('QString', 2) #Sets the qt string to native python strings so can be read without weird stuff

import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *


from mouse import *
from mouse_lever import *
from cat import *
from rat import *

from analysis import *
from mouse_lever_analysis import *

np.set_printoptions(suppress=True)      #Supress scientific notation printing


class Load(QtGui.QWidget):
    def __init__(self, parent):
        super(Load, self).__init__(parent)

        layout = QtGui.QGridLayout()
        self.parent=parent
        
        parent.animal_name = QLineEdit(parent.animal_name_text) #'2016_07_11_vsd'
        animal_name_text_lbl = QLabel('Exp name:', parent.animal_name)

        parent.rec_name = QLineEdit(parent.rec_name_text) #'track1_150Hz_iso1.0_spontaneous.rhd'
        rec_name_text_lbl = QLabel('Rec name:', parent.rec_name)

        layout.addWidget(animal_name_text_lbl, 0,0)
        layout.addWidget(parent.animal_name, 0,1)
        layout.addWidget(rec_name_text_lbl,1,0)
        layout.addWidget(parent.rec_name,1,1)
                
        #MAKE BUTTONS             
        self.button1 = QPushButton('place holder')
        self.button1.setMaximumWidth(200)
        #self.button1.clicked.connect(self.view_stm)
        layout.addWidget(self.button1, 5, 0)

        self.setLayout(layout)


    def load_mouse(self, main_widget):
        
        print("....loading experiment ...")
        main_widget.root_dir = '/media/cat/12TB/in_vivo/tim/cat/' 
        
        main_widget.animal = Mouse(QtGui.QFileDialog.getExistingDirectory(main_widget, 'Load Experiment', main_widget.root_dir).replace(main_widget.root_dir,''), main_widget.root_dir)
        main_widget.animal_name_text=main_widget.animal.name
        
        self.parent.setWindowTitle(main_widget.animal.name)
        
        
        
        #OTHER OPTIONS FOR LOADING
        #1) QFileDialog.getExistingDirectory(...)
        #2) QFileDialog.getOpenFileName(...)
        #3) QFileDialog.getOpenFileNames(...)
        #4) QFileDialog.getSaveFileName(...)
        #self.default_parameters() 

    def load_mouse_lever(self, main_widget):
        
        print("....loading experiment ...")
        main_widget.root_dir = '/media/cat/12TB/in_vivo/tim/yuki/' 
        main_widget.n_sec = 3
        
        main_widget.animal = Mouse_lever(QtGui.QFileDialog.getExistingDirectory(main_widget, 'Load Experiment', main_widget.root_dir).replace(main_widget.root_dir,''), main_widget.root_dir, main_widget.n_sec)
        
        main_widget.animal_name_text=main_widget.animal.name
        
        self.parent.setWindowTitle(main_widget.animal.name)
        


    def load_cat(self, main_widget):
        
        print("....loading experiment ...")
        
        main_widget.root_dir = '/media/cat/8TB/in_vivo/nick/' 

        main_widget.animal = Cat(QtGui.QFileDialog.getExistingDirectory(main_widget, 'Load Experiment', main_widget.root_dir).replace(main_widget.root_dir,''), main_widget.root_dir)
        main_widget.animal_name_text=main_widget.animal.name

        self.parent.setWindowTitle(main_widget.animal.name)


    def load_rat(self, main_widget):
        
        print("....loading experiment ...")
        
        main_widget.root_dir = '/media/cat/8TB/in_vivo/seamans/' 

        main_widget.animal = Rat(QtGui.QFileDialog.getExistingDirectory(main_widget, 'Load Experiment', main_widget.root_dir).replace(main_widget.root_dir,''), main_widget.root_dir)
        main_widget.animal_name_text=main_widget.animal.name

        self.parent.setWindowTitle(main_widget.animal.name)
        

    def select_recording(self, main_widget):
        print "... selecting recording ..."
        print main_widget.animal_name.text()

        main_widget.animal.recName = QtGui.QFileDialog.getOpenFileName(self, 'Load File', main_widget.root_dir+main_widget.animal.name+'/rhd_files/')
        
        #Simplified loading
        if True:
            self.parent.recName = main_widget.animal.recName
        else:
            main_widget.rec_name_text = main_widget.animal.recName.replace(main_widget.animal.home_dir+main_widget.animal.name+'/rhd_files/', '')
            
            self.parent.setWindowTitle(main_widget.animal.name+'/'+main_widget.rec_name_text)

            self.parent.animal.load_tsf_header(main_widget.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','_hp.tsf'))
            self.parent.animal.rec_length = self.parent.animal.tsf.n_vd_samples/float(self.parent.animal.tsf.SampleFrequency)

        #self.default_parameters()



class EventTriggeredAnalysis(QtGui.QWidget):
    def __init__(self, parent):
        super(EventTriggeredAnalysis, self).__init__(parent)
        self.parent = parent
        layout = QtGui.QGridLayout()    


        #Default: july 11, 2016 experiment
        #self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/' 
        self.parent.root_dir = '/media/cat/500GB/in_vivo/tim/cat/'
        self.parent.n_sec = 3
        
        row_index = 0   #Keep track of button/box row
        
        #**************************************************************************************
        #******************************** SELECT ANIMAL & SESSION *****************************
        #**************************************************************************************
        self.vid_analysis_lbl = QLabel('EVENT TRIGGERED ANALYSIS', self)
        self.vid_analysis_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.vid_analysis_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.vid_analysis_lbl, row_index, 0); row_index+=1
        
        #Select recording
        self.button_select_recording = QPushButton('Select Recording')
        self.button_select_recording.setMaximumWidth(200)
        self.button_select_recording.clicked.connect(self.slct_recording)
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


        layout.addWidget(QLabel(' '*40, self), row_index,0); row_index+=1

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
        
        
        self.n_sec_window = QLineEdit('3')
        self.n_sec_window.setMaximumWidth(50)
        self.n_sec_window_lbl = QLabel('#Sec Window:', self)
        layout.addWidget(self.n_sec_window_lbl, row_index,4)
        layout.addWidget(self.n_sec_window, row_index,5)        
        
        row_index+=1
        
        self.button31 = QPushButton('Static STM')
        self.button31.setMaximumWidth(200)
        self.button31.clicked.connect(self.static_stm_mouse_lever)
        layout.addWidget(self.button31, row_index, 0)
        
        self.button32 = QPushButton('Video STM')
        self.button32.setMaximumWidth(200)
        self.button32.clicked.connect(self.video_stm_mouse_lever)
        layout.addWidget(self.button32, row_index, 1)        
        
        self.block_save = QLineEdit('10')
        self.block_save.setMaximumWidth(50)
        self.block_save_lbl = QLabel('Block Ave & Mid-Mask Pixels:', self)
        layout.addWidget(self.block_save_lbl, row_index,4)
        layout.addWidget(self.block_save, row_index,5)
        
        self.midline_mask = QLineEdit('5')
        self.midline_mask.setMaximumWidth(50)
        layout.addWidget(self.midline_mask, row_index,6); row_index+=1
        
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


    def slct_recording(self):
        self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, ".bin, .npy", self.selected_recording,"(*.bin *.tif)") 
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
        


class MouseLeverTools(QtGui.QWidget):
    def __init__(self, parent):
        super(MouseLeverTools, self).__init__(parent)
        self.parent = parent

        layout = QtGui.QGridLayout()    

        #Mouse IA1 as default experiment for Greg experiments
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/yuki/' 
        self.parent.n_sec = 3
        self.parent.exp_type = 'mouse_lever'
        
        row_index = 0   #Keep track of button/box row


        #**************************************************************************************
        #******************************** SELECT ANIMAL & SESSION *****************************
        #**************************************************************************************
        #Select animal
        self.select_lbl = QLabel('Select Experiment----->', self)
        self.select_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        layout.addWidget(self.select_lbl, row_index,1)

        self.comboBox_select_animal = QtGui.QComboBox(self)
        file_names = sorted(glob.glob(self.parent.root_dir+"*"))
        for file_name in file_names:
            self.comboBox_select_animal.addItem(file_name.replace(self.parent.root_dir,''))
        layout.addWidget(self.comboBox_select_animal, row_index,2)
        self.comboBox_select_animal.activated[str].connect(self.select_animal); self.selected_animal = file_names[0].replace(self.parent.root_dir,'')

        #Select session
        self.select_lbl = QLabel('Select Session----->', self)
        self.select_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        layout.addWidget(self.select_lbl, row_index,3)

        self.comboBox_select_session = QtGui.QComboBox(self)
        file_names = sorted(glob.glob(self.parent.root_dir + self.selected_animal + "/tif_files/*"))
        for file_name in file_names:
            self.comboBox_select_session.addItem(file_name.replace(self.parent.root_dir+self.selected_animal+"/tif_files/",''))
        layout.addWidget(self.comboBox_select_session, row_index,4)
        self.comboBox_select_session.activated[str].connect(self.select_session); self.selected_session = file_names[0].replace(self.parent.root_dir + self.selected_animal + "/tif_files/",'')
        
        #Load Mouse from default values
        self.parent.animal_name_text=self.selected_animal.replace(self.parent.root_dir,'')
        self.parent.animal = Mouse_lever(self.parent.animal_name_text, self.parent.root_dir, self.parent.n_sec)
        self.parent.setWindowTitle(self.parent.animal.name)

        #Find movie if available:
        self.movie_available_lbl = QLabel('', self)
        self.movie_available_lbl.setFont(QtGui.QFont("", 12) )
        if os.path.exists(self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v')==True: 
            self.movie_available_lbl.setText("Movie exists") 
        else: 
            self.movie_available_lbl.setText("No movie")
        layout.addWidget(self.movie_available_lbl, row_index,5)
        
        row_index+=1

    

        #**************************************************************************************
        #********************************** PRE-PROCESSING  ************************************
        #**************************************************************************************

        self.preprocess_lbl = QLabel('PRE-PROCESSING', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        #Convert .tif to .npy
        self.button1 = QPushButton('Convert .tif -> .npy')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.cvrt_tif_npy)
        layout.addWidget(self.button1, row_index, 0)
        
        #Align sessions
        self.button1 = QPushButton('Align Sessions')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.algn_images)
        layout.addWidget(self.button1, row_index, 1)
        
        #Other Functions
        self.button1 = QPushButton('Other Functions')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.algn_images)
        layout.addWidget(self.button1, row_index, 2); row_index +=1
        
        
        #**************************************************************************************
        #************************************ FILTERING ***************************************
        #**************************************************************************************
        
        self.button2 = QPushButton('Pre-Filter Images')
        self.button2.setMaximumWidth(200)
        self.button2.clicked.connect(self.fltr_mouse_lever)
        layout.addWidget(self.button2, row_index, 0)

        self.filter_list = ['nofilter', 'butterworth', 'chebyshev']
        self.comboBox_filter = QtGui.QComboBox(self)
        for filter_ in self.filter_list[1:]: self.comboBox_filter.addItem(filter_)
        layout.addWidget(self.comboBox_filter, row_index,1)
        self.comboBox_filter.activated[str].connect(self.select_filter); self.selected_filter = "butterworth" #Set default

        parent.filter_low = QLineEdit('0.1')
        parent.filter_low.setMaximumWidth(50)
        filter_low_lbl = QLabel('Low Cutoff (HZ):', self)
        layout.addWidget(filter_low_lbl, row_index,2)
        layout.addWidget(parent.filter_low, row_index,3)
       
        parent.filter_high = QLineEdit('6.0')
        parent.filter_high.setMaximumWidth(50)
        filter_high_lbl = QLabel('High Cutoff (HZ):', self)
        layout.addWidget(filter_high_lbl, row_index,4)
        layout.addWidget(parent.filter_high, row_index,5); row_index+=1
        
        
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
        self.button3.clicked.connect(self.dff_mouse_lever)
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
        
        #Select reward code; self.selected_session contains full file name including session name
        self.code_list = ['02', '04', '07']
        self.locs_44threshold = []; self.code_44threshold = []
        temp_file = self.parent.root_dir+self.selected_animal + '/tif_files/'+self.selected_session+'/'+self.selected_session
        print temp_file
        if os.path.exists(temp_file + '_locs44threshold.npy')==True: 
            self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
            self.code_44threshold = np.load(temp_file+'_code44threshold.npy')
        
        self.selected_code = '02'; self.n_codes = np.count_nonzero(self.code_44threshold == self.selected_code)         
        
        self.comboBox_select_code = QtGui.QComboBox(self)
        for code_ in self.code_list: self.comboBox_select_code.addItem(code_)
        layout.addWidget(self.comboBox_select_code, row_index, 3)
        self.comboBox_select_code.activated[str].connect(self.select_reward_code); 
        
        self.n_codes_lbl = QLabel(str(self.n_codes), self)
        layout.addWidget(self.n_codes_lbl, row_index, 4); row_index+=1
        
        
        self.button31 = QPushButton('Static STM')
        self.button31.setMaximumWidth(200)
        self.button31.clicked.connect(self.static_stm_mouse_lever)
        layout.addWidget(self.button31, row_index, 0)
        
        self.button32 = QPushButton('Video STM')
        self.button32.setMaximumWidth(200)
        self.button32.clicked.connect(self.video_stm_mouse_lever)
        layout.addWidget(self.button32, row_index, 1)        
        
        trial_lbl = QLabel('Select Trial:', self)
        layout.addWidget(trial_lbl, row_index,2)
                
        self.comboBox_select_trial = QtGui.QComboBox(self)
        n_trials_in = self.load_stm_name()
        self.selected_trial = ''
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #Redundant method
        layout.addWidget(self.comboBox_select_trial, row_index,3)
        self.comboBox_select_trial.activated[str].connect(self.select_trial)
        
        self.block_save = QLineEdit('10')
        self.block_save.setMaximumWidth(50)
        self.block_save_lbl = QLabel('Block Ave & Mid-Mask Pixels:', self)
        layout.addWidget(self.block_save_lbl, row_index,4)
        layout.addWidget(self.block_save, row_index,5)
        
        self.midline_mask = QLineEdit('5')
        self.midline_mask.setMaximumWidth(50)
        layout.addWidget(self.midline_mask, row_index,6)
               
        row_index+=1
        
        
        #**************************************************************************************
        #*********************************** VIDEO TOOLS **************************************
        #**************************************************************************************
        
        for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
        row_index+=1

        self.vid_analysis_lbl = QLabel('VIDEO ANALYSIS', self)
        self.vid_analysis_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.vid_analysis_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.vid_analysis_lbl, row_index, 0); row_index+=1
        
        #Convert videos; 
        self.button4 = QPushButton('Convert: .m4v => .npy')
        self.button4.setMaximumWidth(200)
        self.button4.clicked.connect(self.conv_video)
        layout.addWidget(self.button4, row_index, 0)
        
        #self.comboBox_select_video = QtGui.QComboBox(self)
        #layout.addWidget(self.comboBox_select_video, row_index,1)
        #self.comboBox_select_video.activated[str].connect(self.select); self.choice2 = ""

        self.button41 = QPushButton('Find start/end of vid')
        self.button41.setMaximumWidth(200)
        self.button41.clicked.connect(self.fnd_start_end)
        layout.addWidget(self.button41, row_index, 1)

        self.button42 = QPushButton('Plot blue_light ROI')
        self.button42.setMaximumWidth(200)
        self.button42.clicked.connect(self.bl_light_roi)
        layout.addWidget(self.button42, row_index, 2)

        self.blue_light_std = QLineEdit('1.0')
        self.blue_light_std.setMaximumWidth(50)
        self.blue_light_std_lbl = QLabel('std of blue light):', self)
        layout.addWidget(self.blue_light_std_lbl, row_index, 3)
        layout.addWidget(self.blue_light_std, row_index, 4); row_index+=1
        

        self.button431 = QPushButton('Movies: Single-Trial + [Ca]')
        self.button431.setMaximumWidth(200)
        self.button431.clicked.connect(self.evt_movies_ca)
        layout.addWidget(self.button431, row_index, 0)

        self.button43 = QPushButton('Movies: Multi-Trial')
        self.button43.setMaximumWidth(200)
        self.button43.clicked.connect(self.evt_movies)
        layout.addWidget(self.button43, row_index, 1)

        self.n_trials_movies = QLineEdit('12')
        self.n_trials_movies.setMaximumWidth(50)
        self.n_trials_movies_lbl = QLabel('# of trials:', self)
        layout.addWidget(self.n_trials_movies_lbl, row_index, 2)
        layout.addWidget(self.n_trials_movies, row_index, 3); row_index+=1

        

        #**************************************************************************************
        #********************************** STROKE TOOLS **************************************
        #**************************************************************************************
        for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
        row_index+=1
        
        self.post_processing_lbl = QLabel('STROKE TOOLS', self)
        self.post_processing_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.post_processing_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.post_processing_lbl, row_index, 0); row_index+=1
        
        #Load mouse from disk
        self.button5 = QPushButton('Load Mouse')
        self.button5.setMaximumWidth(200)
        self.button5.clicked.connect(self.load_m_lever)
        layout.addWidget(self.button5, row_index, 0); row_index+=1


        #Separate sessions into epoch periods
        self.button6 = QPushButton('Set Epochs (days)')
        self.button6.setMaximumWidth(200)
        self.button6.clicked.connect(self.load_m_lever)
        layout.addWidget(self.button6, row_index, 0)

        pre_stroke_days = 21;  post_stroke_days = 14;  post_post_stroke_days = 42
        
        parent.pre_stroke_days = QLineEdit('21');                   
        parent.pre_stroke_days.setMaximumWidth(50)
        pre_stroke_days_lbl = QLabel('# days pre_stroke:', self)
        layout.addWidget(pre_stroke_days_lbl, row_index,1)
        layout.addWidget(parent.pre_stroke_days, row_index,2)
        
        parent.post_stroke_days = QLineEdit('14');                     
        parent.post_stroke_days.setMaximumWidth(50)
        post_stroke_days_lbl = QLabel('# days post_stroke:', self)
        layout.addWidget(post_stroke_days_lbl, row_index,3)
        layout.addWidget(parent.post_stroke_days, row_index,4)

        parent.post_post_stroke_days = QLineEdit('42');               
        parent.post_post_stroke_days.setMaximumWidth(50)
        post_post_stroke_days_lbl = QLabel('# days pp_stroke:', self)
        layout.addWidget(post_post_stroke_days_lbl, row_index,5)
        layout.addWidget(parent.post_post_stroke_days, row_index,6)   ; row_index+=1     
        
        
        #Dim reduction: options
        self.button7 = QPushButton('Dimension Reduction:')
        self.button7.setMaximumWidth(200)
        self.button7.clicked.connect(self.dim_red_mouse_lever)
        layout.addWidget(self.button7, row_index, 0)
        self.comboBox3 = QtGui.QComboBox(self)
        self.comboBox3.addItem("PCA")
        self.comboBox3.addItem("MDS")
        self.comboBox3.addItem("tSNE")
        self.comboBox3.addItem("tSNE Barnes-Hutt")
        
        layout.addWidget(self.comboBox3, row_index,1); row_index+=1
        self.comboBox3.activated[str].connect(self.style_choice3); self.choice3="PCA"
        

        #KMeans: options
        self.button8 = QPushButton('Kmeans:')
        self.button8.setMaximumWidth(200)
        self.button8.clicked.connect(self.kmeans_mouse_lever)
        layout.addWidget(self.button8, row_index, 0)

        parent.kmeans_clusters = QLineEdit('16');                   
        parent.kmeans_clusters.setMaximumWidth(50)
        kmeans_clusters_lbl = QLabel('# of clusters:', self)
        layout.addWidget(kmeans_clusters_lbl, row_index,1)
        layout.addWidget(parent.kmeans_clusters, row_index,2)
        
        
        #View clusters in dim reduction space
        self.button9 = QPushButton('View Clusters - 3D')
        self.button9.setMaximumWidth(200)
        self.button9.clicked.connect(self.clusters_plot_mouse_lever)
        layout.addWidget(self.button9, row_index, 3); row_index+=1
        

        #Select a cluster from 2D trace groupings
        self.button10 = QPushButton('Select Cluster - 2D')
        self.button10.setMaximumWidth(200)
        self.button10.clicked.connect(self.select_cluster_mouse_lever)
        layout.addWidget(self.button10, row_index, 0)

        #Plot Cluster Motiff
        self.button11 = QPushButton('Plot Cluster Motiff')
        self.button11.setMaximumWidth(200)
        self.button11.clicked.connect(self.plot_motiff_mouse_lever)
        layout.addWidget(self.button11, row_index, 1); row_index+=1
        
              
        #Movies - single trial
        self.button12 = QPushButton('Movies - Single Trial')
        self.button12.setMaximumWidth(200)
        self.button12.clicked.connect(self.movies_single_mouse_lever)
        layout.addWidget(self.button12, row_index, 0); row_index+=1
        
              
        #1D plots
        self.button13 = QPushButton('1D - Plots')
        self.button13.setMaximumWidth(200)
        self.button13.clicked.connect(self.plot_1D_mouse_lever)
        layout.addWidget(self.button13, row_index, 0)        
        
        #Pixel plots
        self.button14 = QPushButton('Pixel - Plots')
        self.button14.setMaximumWidth(200)
        self.button14.clicked.connect(self.plot_pixel_mouse_lever)
        layout.addWidget(self.button14, row_index, 1)   
        
        #Make Movies 
        self.button15 = QPushButton('Make Movies')
        self.button15.setMaximumWidth(200)
        self.button15.clicked.connect(self.make_movies_mouse_lever)
        layout.addWidget(self.button15, row_index, 2)   ; row_index+=1
        
        
        self.setLayout(layout)


    def cvrt_tif_npy(self):
        print "...convert .tif to .npy... CODE BELOW (BUT NOT IMPLEMENTED)..."

        if (os.path.exists(self.tif_file[:-4] +'.npy')==False) and (os.path.exists(self.tif_file[:-4] +'_aligned.npy')==False):
            print "...read: ", self.tif_file
            images_raw = tiff.imread(self.tif_file)

            print "... saving .npy"
            np.save(self.tif_file[:-4], images_raw)

    def fltr_mouse_lever(self):
        filter_data(self)
    
    def movie_available(self, text):
        
        self.comboBox_movie_available = QtGui.QComboBox(self)
        if os.path.exists(self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v')==True: 
            self.comboBox_movie_available.addItem("Movie exists")
        else: 
            self.comboBox_movie_available.addItem("No movie")
        layout.addWidget(self.comboBox_movie_available, row_index,5)
        self.comboBox_select_session.activated[str].connect(self.movie_available)
        


    #SELECT ANIMAL
    def select_animal(self, text):
        print "...animal: ", text
        self.selected_animal=text
        self.parent.animal_name_text = text
        
        #Reload animal object
        self.parent.animal = Mouse_lever(self.parent.animal_name_text, self.parent.root_dir, self.parent.n_sec)
        self.parent.setWindowTitle(self.parent.animal.name)

        #Reload session list and reset session box
        self.comboBox_select_session.clear()
        file_names = sorted(glob.glob(self.parent.root_dir+self.parent.animal.name+"/tif_files/*"))
        for file_name in sorted(file_names):
            self.comboBox_select_session.addItem(file_name.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",''))
        self.selected_session = file_names[0].replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
        self.select_session(self.selected_session)
       
        #Reset code functions
        self.comboBox_select_code.setCurrentIndex(0)
        self.locs_44threshold = []; self.code_44threshold = []
        temp_file = self.selected_session+'/'+self.selected_session.replace(self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/','')
        if os.path.exists(temp_file + '_locs44threshold.npy')==True: 
            self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
            self.code_44threshold = np.load(temp_file+'_code44threshold.npy')

        self.select_reward_code('02')   #Reset reward code to '02' for new animal load
    
        print self.selected_animal
        print self.selected_session
        print self.selected_code

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
     
        #Find movie availabliltiy
        if os.path.exists(self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v')==True: 
            self.movie_available_lbl.setText("Movie exists")
        else: 
            self.movie_available_lbl.setText("No movie")
        
    #SELECT SESSION
    def select_session(self, text):
        print "...session: ", text
        self.selected_session=text
        self.parent.animal.recName = text; self.parent.rec_name_text= text
        
        print "... self.selected_animal: ", self.selected_animal
        print "... self.selected_session: ", self.selected_session
        
        #CALL CODE COUNTING FUNCTIONS
        self.locs_44threshold = []; self.code_44threshold = []
        temp_file = self.parent.root_dir+self.parent.animal.name+"/tif_files/" +self.selected_session+'/'+self.selected_session
        if os.path.exists(temp_file + '_locs44threshold.npy')==True: 
            self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
            self.code_44threshold = np.load(temp_file+'_code44threshold.npy')

        print self.code_44threshold
        
        self.comboBox_select_code.setCurrentIndex(0)
        self.selected_code = '02'; self.n_codes = np.count_nonzero(self.code_44threshold == self.selected_code)         
        self.n_codes_lbl.setText(str(self.n_codes))
        print "\n\n"

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
            
        #Find movie availabliltiy
        if os.path.exists(self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v')==True: 
            self.movie_available_lbl.setText("Movie exists")
        else: 
            self.movie_available_lbl.setText("No movie")
            
    #SELECT REWARD CODE
    def select_reward_code(self, text):
        print "...code selected: ", text 
        self.selected_code = text
        
        #CALL CODE COUNTING FUNCTIONS
        self.locs_44threshold = []; self.code_44threshold = []
        temp_file = self.parent.root_dir + self.parent.animal.name + "/tif_files/"+ self.selected_session+'/'+self.selected_session
        print "...temp_file: ", temp_file
        
        if os.path.exists(temp_file + '_locs44threshold.npy')==True: 
            self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
            self.code_44threshold = np.load(temp_file+'_code44threshold.npy')

    
        print self.code_44threshold     
        self.n_codes = np.count_nonzero(self.code_44threshold == self.selected_code)         
        self.n_codes_lbl.setText(str(self.n_codes))
        
        print "..reseting reward code, # codes: ", self.n_codes
        print "\n\n"

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
        

    def select_trial(self, text):
        print "...selecting trial: ", text
        self.selected_trial = text
        

    def load_stm_name(self):        

        if self.selected_dff_filter == 'nofilter':
            self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
                str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
        else:
            self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
                str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

        filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'

        #print "...resetting stm_name to: ", filename

        if os.path.exists(filename)==True: 
            data = np.load(filename,  mmap_mode='r+')
            print data.shape
            return data.shape[0]
        
        return 0
        

    def select_dff_filter(self,text):
        self.selected_dff_filter=text
        
        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
            
        
    def select_filter(self, text):
        self.selected_filter = text

    def select_dff_method(self, text):
        self.dff_method = text
        
        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
            

    def style_choice2(self, text):
        self.choice2 = text
        
    def style_choice3(self, text):
        self.choice3 = text

    def static_stm_mouse_lever(self):
        view_static_stm(self)


    def video_stm_mouse_lever(self):
        view_video_stm(self)
        

    def algn_images(self):
        
        #load specific session .npy images and align to 1000 frame of 1st session
        print "... alligning session ... NOT CURRENTLY IMPLEMENTED..."
        
        
    def conv_video(self):
        convert_video(self)

    def fnd_start_end(self):
        find_start_end(self)

    def bl_light_roi(self):
        plot_blue_light_roi(self)

    def evt_movies(self):
        event_triggered_movies(self)

    def evt_movies_ca(self):
        event_triggered_movies_Ca(self)

    def kmeans_mouse_lever(self):
        print "...kmeans..."
        
    def clusters_plot_mouse_lever(self):
        print "... visualizing clustered distributions in 3D post dim-reduction..."
        
        plot_pyqt(app, dim_red_data, mouse.cluster_labels)

    def select_cluster_mouse_lever(self):
        print "... selecting cluster post dim-reduction..."
 
        plot_traces_pyqt(app, mouse)

    def plot_motiff_mouse_lever(self):
        print "... plotting mouse lever motiff ..."
        
        plot_selected_cluster_DFF(mouse)

    def plot_1D_mouse_lever(self):
        print "... plotting 1D ..."
        
        plot_1D(mouse, generic_mask_indexes)   

    def plot_pixel_mouse_lever(self):
        
        print "... plotting 1D ..."
        
        pixel_plot(mouse, generic_mask_indexes)

    def make_movies_mouse_lever(self):
        print "... make movies ..."
        make_movies(mouse, stacks, v_max, v_min)


    def movies_single_mouse_lever(self):
        print "... make movies signle trial ..."
        make_movies_singletrials(mouse, data_chunks)


    def dim_red_mouse_lever(self):
        print "... dim reduction..."
        dim_reduction (self.parent.animal, text) 


    def dff_mouse_lever(self):
        compute_dff_mouse_lever(self)

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant

    def convert_files(self):
        self.parent.animal.preprocess_mouse_lever()


    def load_m_lever(self):
        self.parent.animal.load()

        
    def chunk_mouse_lever(self):
        #Call chunking algorithm
        self.parent.animal.chunk_data(pre_stroke_days, post_stroke_days, post_post_stroke_days)

    ##THIS DOESN"T WORK PROPERLY
    #def paintEvent(self, event):
        #painter = QtGui.QPainter()
        #painter.begin(self)
        #self.hourColor = QtGui.QColor(255, 255, 0)
        #painter.drawLine(10, 150, 200, 150)
        
        
class CatTools(QtGui.QWidget):
    def __init__(self, parent):
        super(CatTools, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()    
        
        #NB: ###########NEED TO MOVE THESE TO INDIVIDUALLY CALLABLE OR AUTO-CALLABLE BUTTONS
        #self.animal.rhd_to_tsf()
        #self.animal.tsf_to_lfp()
        #self.animal.lfp_compress()
        #self.animal.bin_to_npy()
        #self.animal.rhd_digital_save()


class RatTools(QtGui.QWidget):
    def __init__(self, parent):
        super(RatTools, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()    
        
        #NB: ###########NEED TO MOVE THESE TO INDIVIDUALLY CALLABLE OR AUTO-CALLABLE BUTTONS
        #self.animal.rhd_to_tsf()
        #self.animal.tsf_to_lfp()
        #self.animal.lfp_compress()
        #self.animal.bin_to_npy()
        #self.animal.rhd_digital_save()



class MSL(QtGui.QWidget):
    def __init__(self, parent):
        super(MSL, self).__init__(parent)
        self.parent = parent
        
        self.parent.root_dir = '/media/cat/12TB/in_vivo/tim/cat/'
        
        layout = QtGui.QGridLayout()

        row_index=0

        #**************************************************************************************
        #*********************************** LOAD FILES  **************************************
        #**************************************************************************************
        #Set root dir; button, then label
        self.button_set_sua_file = QPushButton('Single Unit File')
        self.button_set_sua_file.setMaximumWidth(200)
        self.button_set_sua_file.clicked.connect(self.set_sua_file)
        layout.addWidget(self.button_set_sua_file, row_index, 0)
        
        self.parent.button_set_sua_file  = os.getcwd()
        self.button_set_sua_file_lbl = QLabel(self.parent.button_set_sua_file, self)
        layout.addWidget(self.button_set_sua_file_lbl, row_index,1); row_index+=1

        #Set recording name; button, then label
        self.button_set_lfp_event_file = QPushButton('LPF Event File')
        self.button_set_lfp_event_file.setMaximumWidth(200)
        self.button_set_lfp_event_file.clicked.connect(self.set_lfp_event_file)
        layout.addWidget(self.button_set_lfp_event_file, row_index, 0)
        
        self.parent.set_lfp_event_file = ''
        self.button_set_lfp_event_file_lbl = QLabel(self.parent.set_lfp_event_file, self)
        layout.addWidget(self.button_set_lfp_event_file_lbl, row_index,1); row_index+=1
        

        #**************************************************************************************
        #*********************************** SET LFP #  **************************************
        #**************************************************************************************
        parent.start_lfp = QLineEdit('0');                   
        parent.start_lfp.setMaximumWidth(50)
        start_lfp_lbl = QLabel('start_lfp_cluster:', self)
        layout.addWidget(start_lfp_lbl, row_index,0)
        layout.addWidget(parent.start_lfp, row_index,1)
               
        parent.end_lfp = QLineEdit('1');                     
        parent.end_lfp.setMaximumWidth(50)
        end_lfp_lbl = QLabel('end_lfp_cluster:', self)
        layout.addWidget(end_lfp_lbl,row_index,2)
        layout.addWidget(parent.end_lfp,row_index,3)


        parent.time_chunks = QLineEdit('1');               
        parent.time_chunks.setMaximumWidth(50)
        time_chunks_lbl = QLabel('# of time_chunks:', self)
        layout.addWidget(time_chunks_lbl,row_index,4)
        layout.addWidget(parent.time_chunks,row_index,5)
                
        parent.lock_window = QLineEdit('200');               
        parent.lock_window.setMaximumWidth(50)
        lock_window_lbl = QLabel('lock_window (ms):', self)
        layout.addWidget(lock_window_lbl,row_index,6)
        layout.addWidget(parent.lock_window,row_index,7)


        parent.n_spikes = QLineEdit('250');               
        parent.n_spikes.setMaximumWidth(50)
        n_spikes_lbl = QLabel('min_spikes:', self)
        layout.addWidget(n_spikes_lbl,row_index,8)
        layout.addWidget(parent.n_spikes,row_index,9)


        #MAKE BUTTONS             
        self.button1 = QPushButton('Compute MSL')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.view_msl)
        layout.addWidget(self.button1, 5, 0)

        self.button1 = QPushButton('View P-vals')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.view_msl_Pvals)
        layout.addWidget(self.button1, 6, 0)
        
        self.setLayout(layout)


    def set_sua_file(self):
        self.parent.sua_file = QtGui.QFileDialog.getOpenFileName(self.parent, "Select SUA file (*.ptcs)", self.parent.root_dir,"PTCS (*.ptcs)")
        self.button_set_sua_file_lbl.setText(self.parent.sua_file)
        #self.parent.setWindowTitle(self.parent.sua_file)

    def set_lfp_event_file(self):
        self.parent.lfp_event_file = QtGui.QFileDialog.getOpenFileName(self.parent, "Select LFP event file (*.ptcs)", self.parent.sua_file,"PTCS (*.ptcs)")
        self.button_set_lfp_event_file_lbl.setText(self.parent.lfp_event_file.replace(self.parent.root_dir,''))
        #self.parent.setWindowTitle(self.parent.button_set_sua_file_lbl.replace(self.parent.root_dir,''))
        
    def view_msl(self, main_widget):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        
        msl_plots(self)

            
    def view_msl_Pvals(self):
        compute_msl_pvals(self.parent)


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
        

#class EventTriggered(QtGui.QWidget):
    #def __init__(self, parent):
        #super(EventTriggered, self).__init__(parent)
        ##layout = QtGui.QFormLayout()
        #layout = QtGui.QGridLayout()
        #self.parent = parent


        #row_index = 0
        ##*****************************************        
        #self.vid_analysis_lbl = QLabel('EVENT TRIGGERED ANALYSIS', self)
        #self.vid_analysis_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        #self.vid_analysis_lbl.setStyleSheet('color: blue')
        #layout.addWidget(self.vid_analysis_lbl, row_index, 0); row_index+=1
        
        ##Select recording
        #self.button_select_recording = QPushButton('Select Recording')
        #self.button_select_recording.setMaximumWidth(200)
        #self.button_select_recording.clicked.connect(self.slct_recording)
        #layout.addWidget(self.button_select_recording, row_index, 0)
        
        ##self.parent.selected_recording  = os.getcwd()
        #self.selected_recording  = '/media/cat/8TB/in_vivo/nick'
        #self.select_recording_lbl = QLabel(self.selected_recording, self)
        #layout.addWidget(self.select_recording_lbl, row_index,1); row_index+=1

        ##Select recording
        #self.button_select_sort = QPushButton('Select Event File')
        #self.button_select_sort.setMaximumWidth(200)
        #self.button_select_sort.clicked.connect(self.slct_sort)
        #layout.addWidget(self.button_select_sort, row_index, 0)
        
        #self.parent.selected_sort  = '/media/cat/12TB/in_vivo/'
        #self.select_sort_lbl = QLabel(self.parent.selected_sort, self)
        #layout.addWidget(self.select_sort_lbl, row_index,1); row_index+=1
    
    
        #self.button2 = QPushButton('Filter ==>')
        #self.button2.setMaximumWidth(200)
        #self.button2.clicked.connect(self.fltr_npy)
        #layout.addWidget(self.button2, row_index, 1)
        

        #self.filter_list = ['nofilter', 'butterworth', 'chebyshev']
        #self.comboBox_filter = QtGui.QComboBox(self)
        #for filter_ in self.filter_list[1:]: self.comboBox_filter.addItem(filter_)
        #layout.addWidget(self.comboBox_filter, row_index, 2)
        #self.comboBox_filter.activated[str].connect(self.slct_filter); self.selected_filter = "butterworth" #Set default


        #parent.filter_low = QLineEdit('0.1')
        #parent.filter_low.setMaximumWidth(50)
        #filter_low_lbl = QLabel('Low Cutoff (HZ):', self)
        #layout.addWidget(filter_low_lbl, row_index, 3)
        #layout.addWidget(parent.filter_low, row_index, 4)

       
        #parent.filter_high = QLineEdit('6.0')
        #parent.filter_high.setMaximumWidth(50)
        #filter_high_lbl = QLabel('High Cutoff (HZ):', self)
        #layout.addWidget(filter_high_lbl, row_index, 5)
        #layout.addWidget(parent.filter_high, row_index, 6); row_index+=1
        
        
        
        ##Prefilter templates
        #self.low_cutoff = QLineEdit('0.1');                #parent.start_time = self.start_time
        #self.low_cutoff.setMaximumWidth(50)
        #self.low_cutoff_lbl = QLabel('Lowcut', self)
        #self.low_cutoff_lbl.setMaximumWidth(100)
        #layout.addWidget(self.low_cutoff_lbl, row_index,0)
        #layout.addWidget(self.low_cutoff, row_index, 1)

        #self.low_cutoff = QLineEdit('6.0');                #parent.start_time = self.start_time
        #self.low_cutoff.setMaximumWidth(50)
        #self.low_cutoff_lbl = QLabel('Highcut', self)
        #self.low_cutoff_lbl.setMaximumWidth(100)
        #layout.addWidget(self.low_cutoff_lbl, row_index,2)
        #layout.addWidget(self.low_cutoff, row_index, 3); row_index+=1

        ##*************************************************************************************
        ##******************************* MAKE STM ********************************************
        ##*************************************************************************************

        
        #parent.start_time = QLineEdit('-0.2');                #parent.start_time = self.start_time
        #parent.start_time.setMaximumWidth(50)
        #start_time_lbl = QLabel('start_time:', self)
        #start_time_lbl.setMaximumWidth(80)
        #layout.addWidget(start_time_lbl, row_index,0)
        #layout.addWidget(parent.start_time, row_index,1)

        #parent.end_time = QLineEdit('+0.2');                  #parent.end_time = self.end_time
        #parent.end_time.setMaximumWidth(50)
        #end_time_lbl = QLabel('end_time:', self)
        #end_time_lbl.setMaximumWidth(80)
        #layout.addWidget(end_time_lbl,row_index,2)
        #layout.addWidget(parent.end_time,row_index,3); row_index+=1
        
        #parent.start_cell = QLineEdit('0');                   #parent.start_cell = self.start_cell
        #parent.start_cell.setMaximumWidth(50)
        #start_cell_lbl = QLabel('start_cell:', self)
        #start_cell_lbl.setMaximumWidth(80)
        #layout.addWidget(start_cell_lbl, row_index,0)
        #layout.addWidget(parent.start_cell, row_index,1)
        
        #parent.end_cell = QLineEdit('1');                     #parent.end_cell = self.end_cell
        #parent.end_cell.setMaximumWidth(50)
        #end_cell_lbl = QLabel('end_cell:', self)
        #end_cell_lbl.setMaximumWidth(80)
        #layout.addWidget(end_cell_lbl,row_index,2)
        #layout.addWidget(parent.end_cell,row_index,3); row_index+=1
        
        
        
        ##MAKE BUTTONS             
        #self.button1 = QPushButton('STM')
        #self.button1.setMaximumWidth(200)
        #self.button1.clicked.connect(self.view_stm)
        #layout.addWidget(self.button1, row_index, 0)

        #self.button2 = QPushButton('STM Movies')
        #self.button2.setMaximumWidth(200)
        #layout.addWidget(self.button2, row_index, 1)
        #self.button2.clicked.connect(self.view_stm_movie)
                
        #self.button10 = QPushButton('Compute STM')
        #self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        #layout.addWidget(self.button10, row_index, 2)
        #self.button10.clicked.connect(self.compute_cell_stm)        

        #self.button3 = QPushButton('Cell Rasters')
        #self.button3.setMaximumWidth(200)
        #layout.addWidget(self.button3, row_index, 3)
        #self.button3.clicked.connect(self.view_cell_raster)


        #self.setLayout(layout)



    #def slct_recording(self):
        #self.selected_recording = QtGui.QFileDialog.getOpenFileName(self, "TSF (*.tsf)", self.selected_recording,"TSF (*.tsf)") 
        #path_name, file_name = os.path.split(self.selected_recording)
        #self.select_recording_lbl.setText(file_name)

        #self.tsf = Tsf_file(self.selected_recording)
        #self.tsf.read_ec_traces()

  
    #def slct_sort(self):
        #self.selected_sort =  QtGui.QFileDialog.getOpenFileName(self, "PTCS (*.ptcs)", self.selected_recording,"PTCS (*.ptcs)")
        #path_name, file_name = os.path.split(self.selected_sort)
        #self.select_sort_lbl.setText(file_name)
   
   
        
    #def fltr_npy(self):
        #filter_single_file(self)   
   
   
    #def view_stm(self):
        #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        #sta_maps(self.parent)
        
    #def view_stm_movie(self):
        #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        #sta_movies(self.parent)

    #def view_cell_raster(self):
        #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        #plot_rasters(self.parent)

    #def compute_cell_stm(self):
        #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        #compute_sta(self.parent) #NEEDS TO BE FIXED FROM (self.animal, self.animal.ptcsName) TO CURRENT FORM


class FileDialog(QtGui.QFileDialog):
    def __init__(self, *args, **kwargs):
        super(FileDialog, self).__init__(*args, **kwargs)
        self.setOption(QtGui.QFileDialog.DontUseNativeDialog, True)
        self.setFileMode(QtGui.QFileDialog.ExistingFiles)
        #self.setFileMode(QtGui.QFileDialog.DirectoryOnly)

        self.out_files=[]

    def accept(self):
        self.out_files = self.selectedFiles()   #not sure why can't index into this directly, but seems ok to use "out_files"
        self.deleteLater()
        #super(FileDialog, self).accept()
        
        
#LFP ANALYSIS TOOLBOX
class Track(QtGui.QWidget):
    def __init__(self, parent):
        super(Track, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        
        row_index= 0
        parent.start_time = QLineEdit('0.5');                #parent.start_time = self.start_time
        parent.start_time.setMaximumWidth(50)
        start_time_lbl = QLabel('low_pass (hz):', self)
        start_time_lbl.setMaximumWidth(100)

        parent.end_time = QLineEdit('6.0');                  #parent.end_time = self.end_time
        parent.end_time.setMaximumWidth(50)
        end_time_lbl = QLabel('high_pass (hz):', self)
        end_time_lbl.setMaximumWidth(100)
      
        layout.addWidget(start_time_lbl, row_index,0)
        layout.addWidget(parent.start_time, row_index,1)
        layout.addWidget(end_time_lbl,row_index,2)
        layout.addWidget(parent.end_time,row_index,3); row_index+=1

        self.button_tsf_to_lfp = QPushButton('Convert .tsf to low-pass (1khz) .tsf')
        self.button_tsf_to_lfp.setMaximumWidth(250)
        self.button_tsf_to_lfp.clicked.connect(self.tsftolfp)
        layout.addWidget(self.button_tsf_to_lfp, row_index, 0); row_index+=1

        self.button1 = QPushButton('Concatenate multiple .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_tsf)
        layout.addWidget(self.button1, row_index, 0); row_index+=1

        self.button1 = QPushButton('Concatenate .lfp.zip files')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_lfp_zip)
        layout.addWidget(self.button1, row_index, 0); row_index+=1


        self.button1 = QPushButton('Time compress (1Khz) .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.comp_lfp)
        layout.addWidget(self.button1, row_index, 0)
        
        self.parent.compress_factor = QLineEdit('50');                  #parent.end_time = self.end_time
        self.parent.compress_factor.setMaximumWidth(50)
        self.parent.compress_factor_lbl = QLabel('Compress Factor (multiples of 25):', self)
        layout.addWidget(self.parent.compress_factor_lbl, row_index, 1)
        layout.addWidget(self.parent.compress_factor, row_index, 2)
        
        self.setLayout(layout)

    
    def tsftolfp(self):
        tsf_files = QtGui.QFileDialog.getOpenFileNames(self.parent, "TSF (*.tsf)", self.parent.root_dir,"TSF (*.tsf)")
        tsf_to_lfp(tsf_files)
    

    def multi_tsf(self):
        print "... selecting multiple recordings ..."
        out_files = QtGui.QFileDialog.getOpenFileNames(self.parent, "TSF (*.tsf)", self.parent.root_dir,"TSF (*.tsf)")
        self.parent.animal.tsf_files = out_files
        
        concatenate_tsf(self) 


    def multi_lfp_zip(self):
        print "... selecting multiple .lfp.zip recording directories ..."
        dialog = FileDialog()   #**** SELECT MULTIPLE DIRECTORIES, NOT INDIVIDUAL FIELS
        dialog.exec_()
        
        self.parent.animal.tsf_files = dialog.out_files
        concatenate_lfp_zip(self) 
            

    def comp_lfp(self):
        print "... selecting single lfp recording ..."
        self.parent.animal.tsf_file = QtGui.QFileDialog.getOpenFileName(self.parent, "TSF (*.tsf)", self.parent.root_dir,"TSF (*.tsf)")
        compress_lfp(self) 


#LFP ANALYSIS TOOLBOX
class IntanTools(QtGui.QWidget):
    def __init__(self, parent):
        super(IntanTools, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        
        row_index=0
        
        self.button_rhd_to_tsf = QPushButton('Convert .rhd to .tsf')
        self.button_rhd_to_tsf.setMaximumWidth(250)
        self.button_rhd_to_tsf.clicked.connect(self.rhdtotsf)
        layout.addWidget(self.button_rhd_to_tsf, row_index, 0); row_index+=1
        
        self.setLayout(layout)
        

    def rhdtotsf(self):
        rhd_files = QtGui.QFileDialog.getOpenFileNames(self.parent, "RHD (*.rhd)", self.parent.root_dir,"RHD (*.rhd)")
        rhd_to_tsf(rhd_files) 


#LFP ANALYSIS TOOLBOX
class Seamans(QtGui.QWidget):
    def __init__(self, parent):
        super(Seamans, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        
        row_index = 0
        #*******************************************************************************************
        #*******************************************************************************************
        
        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ncs to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ncs_convert)
        layout.addWidget(self.button1, row_index, 0)
        
        self.parent.make_hp = QLineEdit('False');                #parent.start_time = self.start_time
        self.parent.make_hp.setMaximumWidth(50)
        parent.make_hp_lbl = QLabel('Make High Pass .tsf', self)
        self.parent.make_hp.setMaximumWidth(100)
        layout.addWidget(self.parent.make_hp_lbl, row_index,1)
        layout.addWidget(self.parent.make_hp, row_index,2)
        
        
        row_index+=1

        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ntt to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ntt_convert)
        layout.addWidget(self.button1, row_index, 0)

        #self.button1 = QPushButton('Concatenate .lfp.zip files')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.multi_lfp_zip)
        #layout.addWidget(self.button1, 6, 0)


        #self.button1 = QPushButton('Compress .tsf files')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.comp_lfp)
        #layout.addWidget(self.button1, 7, 0)
        
        
        #self.button1 = QPushButton('Make Track Wide HighPass .tsf')
        #self.button1.setMaximumWidth(250)
        #self.button1.clicked.connect(self.concatenate_tsf)
        #layout.addWidget(self.button1, 7, 0)

        self.setLayout(layout)

        

    def ncs_convert(self):
        self.parent.root_dir = '/media/cat/12TB/in_vivo/jeremy/' 
        ncs_to_tsf(self, QtGui.QFileDialog.getOpenFileNames(self.parent, 'Load Experiment', self.parent.root_dir, "*.ncs"))

    def ntt_convert(self):
        self.parent.root_dir = '/media/cat/12TB/in_vivo/jeremy/' 
        ntt_to_tsf(self, QtGui.QFileDialog.getOpenFileNames(self.parent, 'Load Experiment', self.parent.root_dir, "*.ntt"))


        
#LFP ANALYSIS TOOLBOX
class Filter(QtGui.QWidget):
    def __init__(self, parent):
        super(Filter, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        
        parent.start_time = QLineEdit('0.5');                #parent.start_time = self.start_time
        parent.start_time.setMaximumWidth(50)
        start_time_lbl = QLabel('low_pass (hz):', self)
        start_time_lbl.setMaximumWidth(100)

        parent.end_time = QLineEdit('6.0');                  #parent.end_time = self.end_time
        parent.end_time.setMaximumWidth(50)
        end_time_lbl = QLabel('high_pass (hz):', self)
        end_time_lbl.setMaximumWidth(100)
      
        layout.addWidget(start_time_lbl, 0,0)
        layout.addWidget(parent.start_time, 0,1)
        layout.addWidget(end_time_lbl,0,2)
        layout.addWidget(parent.end_time,0,3)

        #MAKE BUTTONS             
        self.button1 = QPushButton('Butterworth - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 5, 0)

        self.button1 = QPushButton('Wavelet - Ephys (.tsf)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 5, 2)

        self.button1 = QPushButton('Butterworth - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 6, 0)

        self.button1 = QPushButton('Cheby - Imaging (.npy)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_imaging)
        layout.addWidget(self.button1, 6, 2)


        self.setLayout(layout)
   
    def filter_imaging(self):
        print "...pop up window to batch select invidiual files..."



#LFP ANALYSIS TOOLBOX
class TemplatesTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TemplatesTools, self).__init__(parent)
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

#LFP ANALYSIS TOOLBOX
class TraceTools(QtGui.QWidget):
    def __init__(self, parent):
        super(TraceTools, self).__init__(parent)
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
        self.parent.selected_recording  = '/media/cat/8TB/in_vivo/nick'
        self.select_recording_lbl = QLabel(self.parent.selected_recording, self)
        layout.addWidget(self.select_recording_lbl, row_index,1)
        
        #Rec length
        self.parent.rec_length = QLabel('0', self)
        self.parent.rec_length.setMaximumWidth(100)
        self.parent.rec_length_lbl = QLabel('Rec Length:', self)
        self.parent.rec_length_lbl.setMaximumWidth(100)
        layout.addWidget(self.parent.rec_length_lbl, row_index, 3)
        layout.addWidget(self.parent.rec_length, row_index, 4)


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


#LFP ANALYSIS TOOLBOX
class LFP(QtGui.QWidget):
    def __init__(self, parent):
        super(LFP, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        
        row_index = 0
        parent.specgram_ch = QLineEdit('5');                #parent.start_time = self.start_time
        parent.specgram_ch.setMaximumWidth(50)
        specgram_ch_lbl = QLabel('specgram_ch:', self)
        specgram_ch_lbl.setMaximumWidth(100)
        layout.addWidget(specgram_ch_lbl, row_index,0)
        layout.addWidget(parent.specgram_ch, row_index,1)
      
        
        parent.specgram_db_clip = QLineEdit('-40');                #parent.start_time = self.start_time
        parent.specgram_db_clip.setMaximumWidth(50)
        specgram_db_clip_lbl = QLabel('specgram_db_clip:', self)
        specgram_db_clip_lbl.setMaximumWidth(100)
        layout.addWidget(specgram_db_clip_lbl, row_index,2)
        layout.addWidget(parent.specgram_db_clip, row_index,3); row_index+=1  
                      
        self.button3 = QPushButton('LFP Rasters')
        self.button3.setMaximumWidth(200)
        layout.addWidget(self.button3, row_index, 0)
        self.button3.clicked.connect(self.view_lfp_raster); row_index+=1

        self.button10 = QPushButton('View Spegram')
        self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.button10, row_index, 0)
        self.button10.clicked.connect(self.view_Specgram)     
        
        self.setLayout(layout)

   
    def view_stm(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        sta_maps(self.parent)
        
    def view_stm_movie(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        sta_movies(self.parent)

    def view_lfp_raster(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        plot_rasters(self.parent)

    def compute_lfp_stm(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        compute_sta(self.parent) #NEEDS TO BE FIXED FROM (self.animal, self.animal.ptcsName) TO CURRENT FORM

    def view_Specgram(self):
        Specgram_syncindex(self)
            
            
class Window(QtGui.QMainWindow):

    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        #self.setWindowTitle("OpenNeuron")
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))


        #Set widget to show up with viewbox
        toolMenu = QtGui.QMenuBar()
        toolMenu.setNativeMenuBar(False) # <--Sets the menu with the widget; otherwise shows up as global (i.e. at top desktop screen)
        self.setMenuBar(toolMenu)

        #***** TEXT PARAMETERS FIELDS ******
        if True: 
            #Mouse July 11 as default experiment
            self.root_dir = '/media/cat/12TB/in_vivo/tim/cat/' 
            self.animal_name_text='2016_07_11_vsd' 
            self.rec_name_text='track1_150Hz_iso1.0_spontaneous.rhd'
            
            self.animal = Mouse(self.animal_name_text, self.root_dir)
            self.animal.recName =self.root_dir+self.animal_name_text+'/rhd_files/'+ self.rec_name_text
            self.setWindowTitle(self.animal.name+'/'+self.animal.recName.replace(self.animal.home_dir+self.animal.name+'/rhd_files/', ''))
            self.animal.load_tsf_header(self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','_hp.tsf'))
            self.exp_type = 'mouse'

        else: 
            #Cat ptc20 as default experiment
            self.root_dir = '/media/cat/8TB/in_vivo/nick/' 
            self.animal_name_text='ptc20' 
            self.rec_name_text='71-tr3-blankscreen/71-tr3-blankscreen.tsf'

            self.animal = Cat(self.animal_name_text, self.root_dir)
            self.animal.recName =self.root_dir+self.animal_name_text+'/'+ self.rec_name_text
            self.setWindowTitle(self.animal.name+'/'+self.animal.recName.replace(self.root_dir+self.animal_name_text+'/', ''))
            self.animal.load_tsf_header(self.animal.recName)
            self.exp_type = 'cat'

        #Load default experiment

        self.animal.rec_length = self.animal.tsf.n_vd_samples/float(self.animal.tsf.SampleFrequency)


        #Menu Item Lists
        self.make_menu()

        #LOAD CENTRAL WIDGET
        self.central_widget = QtGui.QStackedWidget()
        self.setCentralWidget(self.central_widget)
        
        #SET DEFAULT WIDGET TO PROCESS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        
        self.show()

    def make_menu(self):
        
        #FILE MENUS
        loadMouse = QtGui.QAction("&Load Mouse", self)
        loadMouse.setStatusTip('Load Mouse')
        loadMouse.triggered.connect(self.ld_mouse)

        loadMouseLever = QtGui.QAction("&Load Mouse-Lever", self)
        loadMouseLever.setStatusTip('Load Mouse-Lever')
        loadMouseLever.triggered.connect(self.ld_mouse_lever)

        loadCat = QtGui.QAction("&Load Cat", self)
        loadCat.setStatusTip('Load Cat')
        loadCat.triggered.connect(self.ld_cat)
                
        loadRat = QtGui.QAction("&Load Rat", self)
        loadRat.setStatusTip('Load Rat')
        loadRat.triggered.connect(self.ld_rat)
                
        loadRecording = QtGui.QAction("&Select Recording", self)
        loadRecording.setStatusTip('Select Recording')
        loadRecording.triggered.connect(self.ld_rec)


        #PROCESSING MENUS
        trackTools = QtGui.QAction("&Track Wide Tools", self)
        trackTools.setStatusTip('Track Wide Tools')
        trackTools.triggered.connect(self.track_tools)
        
        intanTools = QtGui.QAction("&Intan Conversion Tools", self)
        intanTools.setStatusTip('Intan Conversion Tools')
        intanTools.triggered.connect(self.intan_tools)


        preprocessExperiment = QtGui.QAction("&Process Data", self)
        preprocessExperiment.setStatusTip('Process Data')
        preprocessExperiment.triggered.connect(self.prepreprocess_data)


        seamansData = QtGui.QAction("&Seamans' Lab Data", self)
        seamansData.setStatusTip('Seamans Lab Data')
        seamansData.triggered.connect(self.seamans_data)


        filterData = QtGui.QAction("&Filter Data", self)
        filterData.setStatusTip('Filter Data')
        filterData.triggered.connect(self.fltr_data)


        #IMAGING TOOLS MENUS
        mouseTools = QtGui.QAction("&Event Triggered Imaging", self)
        mouseTools.setStatusTip('Event Triggered Imaging')
        mouseTools.triggered.connect(self.mouse_tools)

        mouseLeverTools = QtGui.QAction("&Mouse-Lever", self)
        mouseLeverTools.setStatusTip('Mouse-Lever')
        mouseLeverTools.triggered.connect(self.mouse_lever_tools)

        catTools = QtGui.QAction("&Cat", self)
        catTools.setStatusTip('Cat')
        catTools.triggered.connect(self.cat_tools)

        ratTools = QtGui.QAction("&Rat", self)
        ratTools.setStatusTip('Rat')
        ratTools.triggered.connect(self.rat_tools)


        #EPHYS TOOLS MENUS
        View_Traces = QtGui.QAction("&View Traces", self)
        View_Traces.setStatusTip('View Traces')
        View_Traces.triggered.connect(self.view_rtraces)
        
        View_Templates = QtGui.QAction("&View Templates", self)
        View_Templates.setStatusTip('View Templates')
        View_Templates.triggered.connect(self.view_templts)
        
        #Event_Triggered_Maps = QtGui.QAction("&Cell STM", self)
        #Event_Triggered_Maps.setStatusTip('Cell Analysis')
        #Event_Triggered_Maps.triggered.connect(self.event_triggered_analysis)

        LFP_Analysis = QtGui.QAction("&LFP Tools", self)
        LFP_Analysis.setStatusTip('LFP Analysis')
        LFP_Analysis.triggered.connect(self.lfp_analysis)

        MSL_Analysis = QtGui.QAction("&MSL Tools", self)
        MSL_Analysis.setStatusTip('MSL')
        MSL_Analysis.triggered.connect(self.msl_analysis)

        Count_Matrix = QtGui.QAction("&Count Matrix", self)
        Count_Matrix.setStatusTip('Count Matrix')
        Count_Matrix.triggered.connect(self.view_count_matrix)


        exitApplication = QtGui.QAction("&Exit Application", self)
        exitApplication.setStatusTip('Exit')
        exitApplication.triggered.connect(self.close_application)
        
        #MAKE MENUS
        mainMenu = self.menuBar()
        
        fileMenu = mainMenu.addMenu('&Load')
        fileMenu.addAction(loadMouse)
        fileMenu.addAction(loadMouseLever)
        fileMenu.addAction(loadCat)
        fileMenu.addAction(loadRat)
        fileMenu.addAction(loadRecording)
        fileMenu.addAction(exitApplication)

        fileMenu = mainMenu.addMenu('Pre-Process')
        fileMenu.addAction(intanTools)
        fileMenu.addAction(trackTools)
        fileMenu.addAction(seamansData)
        fileMenu.addAction(filterData)
        #fileMenu.addAction(preprocessExperiment)

        fileMenu = mainMenu.addMenu('Imaging Analysis')
        #fileMenu.addAction(Event_Triggered_Maps)
        fileMenu.addAction(mouseTools)
        fileMenu.addAction(mouseLeverTools)
        
        #fileMenu.addAction(catTools)
        #fileMenu.addAction(ratTools)



        fileMenu = mainMenu.addMenu('Ephys Analysis')
        fileMenu.addAction(View_Traces)
        fileMenu.addAction(View_Templates)

        fileMenu.addAction(LFP_Analysis)
        fileMenu.addAction(MSL_Analysis)
        fileMenu.addAction(Count_Matrix)

    #************* LOAD FILE MENUS *****************
    def ld_mouse(self):

        self.load_widget.load_mouse(self)   #Pass main widget to subwidgets as it contains needed parameters.

        self.exp_type = 'mouse'

        #RESTART Process widget with updated info; SEEMS THERE IS A BETTER WAY TO DO THIS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        self.central_widget.setCurrentWidget(self.load_widget)


    def ld_mouse_lever(self):

        self.load_widget.load_mouse_lever(self)   #Pass main widget to subwidgets as it contains needed parameters.

        self.exp_type = 'mouse_lever'

        #RESTART Process widget with updated info; SEEMS THERE IS A BETTER WAY TO DO THIS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        self.central_widget.setCurrentWidget(self.load_widget)

        
    def ld_cat(self):

        self.load_widget.load_cat(self)   #Pass main widget to subwidgets as it contains needed parameters.
        self.exp_type = 'cat'

        #RESTART Process widget with updated info; SEEMS THERE IS A BETTER WAY TO DO THIS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        self.central_widget.setCurrentWidget(self.load_widget)


    def ld_rat(self):

        self.load_widget.load_rat(self)   #Pass main widget to subwidgets as it contains needed parameters.
        self.exp_type = 'rat'

        #RESTART Process widget with updated info; SEEMS THERE IS A BETTER WAY TO DO THIS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        self.central_widget.setCurrentWidget(self.load_widget)
        
        
    def ld_rec(self):
        self.load_widget.select_recording(self)   #Pass main widget to subwidgets as it contains needed parameters.
 
        #RESTART Process widget with updated info; SEEMS THERE SHOULD BETTER WAY TO DO THIS
        self.load_widget = Load(self)
        self.central_widget.addWidget(self.load_widget)
        self.central_widget.setCurrentWidget(self.load_widget)
                
    def close_application(self):
        print("whooaaaa so custom!!!")
        
        sys.exit()


    #********** EXP TOOLS MENUS *******************
    def mouse_tools(self):
        mouse_widget = EventTriggeredAnalysis(self)
        self.central_widget.addWidget(mouse_widget)  
        self.central_widget.setCurrentWidget(mouse_widget)

    def mouse_lever_tools(self):
        mouse_lever_widget = MouseLeverTools(self)
        self.central_widget.addWidget(mouse_lever_widget)  
        self.central_widget.setCurrentWidget(mouse_lever_widget)
    
    def cat_tools(self):
        cat_widget = CatTools(self)
        self.central_widget.addWidget(cat_widget)  
        self.central_widget.setCurrentWidget(cat_widget)

    def rat_tools(self):
        rat_widget = RatTools(self)
        self.central_widget.addWidget(rat_widget)  
        self.central_widget.setCurrentWidget(rat_widget)




    #********** ANALYSIS MENUS *******************
    
    def view_rtraces(self):
        traces_widget = TraceTools(self)
        self.central_widget.addWidget(traces_widget)  
        self.central_widget.setCurrentWidget(traces_widget)
    
    def view_templts(self):
        templates_widget = TemplatesTools(self)
        self.central_widget.addWidget(templates_widget)  
        self.central_widget.setCurrentWidget(templates_widget)
        
    
    #def event_triggered_analysis(self):
        #event_widget = EventTriggered(self)
        #self.central_widget.addWidget(event_widget)  
        #self.central_widget.setCurrentWidget(event_widget)

    def lfp_analysis(self):
        lfp_widget = LFP(self)
        self.central_widget.addWidget(lfp_widget)  
        self.central_widget.setCurrentWidget(lfp_widget)
    
    def msl_analysis(self):
        msl_widget = MSL(self)
        self.central_widget.addWidget(msl_widget)  
        self.central_widget.setCurrentWidget(msl_widget)

    def view_count_matrix(self):
        count_matrix_widget = CountMatrix(self)
        self.central_widget.addWidget(count_matrix_widget)  
        self.central_widget.setCurrentWidget(count_matrix_widget)



    #************ PROCESSING MENUS ***************
    def track_tools(self):
        track_widget = Track(self)
        self.central_widget.addWidget(track_widget)      
        self.central_widget.setCurrentWidget(track_widget)

    def intan_tools(self):
        intan_widget = IntanTools(self)
        self.central_widget.addWidget(intan_widget)      
        self.central_widget.setCurrentWidget(intan_widget)
        


    def seamans_data(self):
        seamans_widget = Seamans(self)
        self.central_widget.addWidget(seamans_widget)      
        self.central_widget.setCurrentWidget(seamans_widget)


    def prepreprocess_data(self):  #NOT CURENTLY USED...
        print "....processing experiment: ", self.fileName
        
        #MUST BE EXPERIMENT SPECIFIC
        

    def fltr_data(self):
        filter_widget = Filter(self)
        self.central_widget.addWidget(filter_widget)  
        self.central_widget.setCurrentWidget(filter_widget)
        
        

def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())


run()
