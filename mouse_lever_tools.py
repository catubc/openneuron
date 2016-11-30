import glob

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *



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
        #************************* PRE-PROCESSING BATCH / ALL SESSIONS ************************
        #**************************************************************************************

        self.preprocess_lbl = QLabel('PREPROCESS ALL SESSIONS', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1


        self.button_preprocess_mlever = QPushButton('Preprocess All Sessions')
        self.button_preprocess_mlever.setMaximumWidth(200)
        self.button_preprocess_mlever.clicked.connect(self.preprocess_mlever)
        layout.addWidget(self.button_preprocess_mlever, row_index, 0)
        
        
        self.button_filter_mlever = QPushButton('Filter All Sessions')
        self.button_filter_mlever.setMaximumWidth(200)
        self.button_filter_mlever.clicked.connect(self.filter_mlever)
        layout.addWidget(self.button_filter_mlever, row_index, 1); row_index+=1
        
        #**************************************************************************************
        #************************* PREPROCESSING SINGLE SESSION  ******************************
        #**************************************************************************************

        
        for k in range(6): layout.addWidget(QLabel(' '*40, self), row_index,k)
        row_index+=1
        
        self.preprocess_lbl = QLabel('PREPROCESS SINGLE SESSION', self)
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
        
        
        #************************************ FILTERING ***************************************
        
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
                
        self.preprocess_lbl = QLabel('[Ca] EVENT TRIGGERED', self)
        self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        self.preprocess_lbl.setStyleSheet('color: blue')
        layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        self.button3 = QPushButton('Compute DFF Pipe ====>')
        self.button3.setMaximumWidth(200)
        self.button3.clicked.connect(self.dff_mouse_lever)
        layout.addWidget(self.button3, row_index, 0)

        #Filter dropdown box
        self.comboBox_select_dff_filter = QtGui.QComboBox(self)
        for filter_ in self.filter_list: self.comboBox_select_dff_filter.addItem(filter_)
        layout.addWidget(self.comboBox_select_dff_filter, row_index,1)
        self.comboBox_select_dff_filter.activated[str].connect(self.select_dff_filter); self.selected_dff_filter = "nofilter" #Set default

        #DFF method dropdown box
        dff_list = ['globalAverage', 'slidingWindow']
        self.comboBox_select_dff_method = QtGui.QComboBox(self)
        for dff_ in dff_list: self.comboBox_select_dff_method.addItem(dff_)
        layout.addWidget(self.comboBox_select_dff_method, row_index,2)
        self.comboBox_select_dff_method.activated[str].connect(self.select_dff_method); self.dff_method = "globalAverage"
        
        #Reward code dropdown box
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
        
        # Number of trials for reward code
        self.n_codes_lbl = QLabel(str(self.n_codes), self)
        layout.addWidget(self.n_codes_lbl, row_index, 4)
        
        
        self.n_sec_window = QLineEdit('3')
        self.n_sec_window.setMaximumWidth(50)
        self.n_sec_window_lbl = QLabel('Window (sec):', self)
        layout.addWidget(self.n_sec_window_lbl, row_index,6)
        layout.addWidget(self.n_sec_window, row_index,7); row_index+=1
        
        
        #**************************************************************************************
        #************************************ VIEW DFF *************************************
        #**************************************************************************************

        
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
        
        self.block_save = QLineEdit('1')
        self.block_save.setMaximumWidth(50)
        self.block_save_lbl = QLabel('Block Ave:', self)
        layout.addWidget(self.block_save_lbl, row_index,4)
        layout.addWidget(self.block_save, row_index,5)
        
        self.midline_mask = QLineEdit('1')
        self.midline_mask.setMaximumWidth(50)
        self.midline_mask_lbl = QLabel('MidlineMask:', self)
        layout.addWidget(self.midline_mask_lbl, row_index,6)
        layout.addWidget(self.midline_mask, row_index,7)

        self.stm_start_time = QLineEdit('-.3')
        self.stm_start_time.setMaximumWidth(50)
        self.stm_start_time_lbl = QLabel('Start time:', self)
        layout.addWidget(self.stm_start_time_lbl, row_index,8)
        layout.addWidget(self.stm_start_time, row_index,9)
        
        self.stm_end_time = QLineEdit('+.3')
        self.stm_end_time.setMaximumWidth(50)
        self.stm_end_time_lbl = QLabel('End time:', self)
        layout.addWidget(self.stm_end_time_lbl, row_index,10)
        layout.addWidget(self.stm_end_time, row_index,11); row_index+=1


        #************************************ FILTERING ***************************************

        #self.stm_activity = QPushButton('Activity Triggered STM')
        #self.stm_activity.setMaximumWidth(200)
        #self.stm_activity.clicked.connect(self.view_stm_activity)
        #layout.addWidget(self.stm_activity, row_index, 0)
        
        
        self.mask_start_frame = QLineEdit('0')
        self.mask_start_frame.setMaximumWidth(50)
        self.mask_start_frame_lbl = QLabel('Mask start frame:', self)
        layout.addWidget(self.mask_start_frame_lbl, row_index,1)
        layout.addWidget(self.mask_start_frame, row_index,2)
        
        self.mask_end_frame = QLineEdit('6')
        self.mask_end_frame.setMaximumWidth(50)
        self.mask_end_frame_lbl = QLabel('End time:', self)
        layout.addWidget(self.mask_end_frame_lbl, row_index,3)
        layout.addWidget(self.mask_end_frame, row_index,4)        
        
        self.mask_width = QLineEdit('12')
        self.mask_width.setMaximumWidth(50)
        self.mask_width_lbl = QLabel('Mask width (pixels):', self)
        layout.addWidget(self.mask_width_lbl, row_index,5)
        layout.addWidget(self.mask_width, row_index,6)
        
        self.mask_percentile = QLineEdit('99.99')
        self.mask_percentile.setMaximumWidth(50)
        self.mask_percentile_lbl = QLabel('Mask Percentile:', self)
        layout.addWidget(self.mask_percentile_lbl, row_index,7)
        layout.addWidget(self.mask_percentile, row_index,8)
        
        self.mask_power = QLineEdit('1.2')
        self.mask_power.setMaximumWidth(50)
        self.mask_power_lbl = QLabel('Mask Power:', self)
        layout.addWidget(self.mask_power_lbl, row_index,9)
        layout.addWidget(self.mask_power, row_index,10); row_index+=1
        
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

        
        for k in range(2): 
            layout.addWidget(QLabel(' '*40, self), row_index,0)
            row_index+=1


        ##**************************************************************************************
        ##***************************** COMBINED [Ca] & VIDEO **********************************
        ##**************************************************************************************
        
                
        #self.preprocess_lbl = QLabel('COMBINED [Ca] & VIDEO', self)
        #self.preprocess_lbl.setFont(QtGui.QFont("Times", 12, QtGui.QFont.Bold) )
        #self.preprocess_lbl.setStyleSheet('color: blue')
        #layout.addWidget(self.preprocess_lbl, row_index, 0); row_index+=1

        
        #self.button_ca_stm = QPushButton('Static [Ca] + STM')
        #self.button_ca_stm.setMaximumWidth(200)
        #self.button_ca_stm.clicked.connect(self.static_stm_ca_mouse_lever)
        #layout.addWidget(self.button_ca_stm, row_index, 0)
        
        #self.button_ca_stm_video = QPushButton('Video [Ca] + STM')
        #self.button_ca_stm_video.setMaximumWidth(200)
        #self.button_ca_stm_video.clicked.connect(self.video_stm_ca_mouse_lever)
        #layout.addWidget(self.button_ca_stm_video, row_index, 1)        
        
        #row_index+=1
                
        
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
        epoch_days_lbl = QLabel('Epoch days: ', self)
        layout.addWidget(epoch_days_lbl, row_index, 0)
        #self.button6 = QPushButton('Set Epochs (days)')
        #self.button6.setMaximumWidth(200)
        #self.button6.clicked.connect(self.load_m_lever)
        #layout.addWidget(self.button6, row_index, 0)

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


    #*****************************************************************************************************
    #************************************* SELECT SESSIONS FUNCTIONS *************************************
    #*****************************************************************************************************


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
        self.tif_file = self.parent.root_dir+self.parent.animal.name+"/tif_files/" +self.selected_session+'/'+self.selected_session+'.tif'
        if os.path.exists(self.tif_file[:-4] + '_locs44threshold.npy')==True: 
            self.locs_44threshold = np.load(self.tif_file[:-4]+'_locs44threshold.npy')
            self.code_44threshold = np.load(self.tif_file[:-4]+'_code44threshold.npy')

        print "... code_44threshold: ", self.code_44threshold
        
        self.comboBox_select_code.setCurrentIndex(0)
        self.selected_code = '02'; self.n_codes = np.count_nonzero(self.code_44threshold == self.selected_code)         
        self.n_codes_lbl.setText(str(self.n_codes))
        #print "\n\n"

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant
            
        #Find movie availability
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
        
    #SELECT TRIAL
    def select_trial(self, text):
        print "...selecting trial: ", text
        self.selected_trial = text
        

    #*****************************************************************************************************
    #************************************* PREPROCESS ALL DATA ******************************************
    #*****************************************************************************************************


    def preprocess_mlever(self):
        #PREPROCESS ALL DATA: Load filenames, sessions, etc.
        self.parent.animal.preprocess_mouse_lever()
        

    def filter_mlever(self):
        #FILTER ALL DATA
        
        file_names = sorted(glob.glob(self.parent.root_dir+self.parent.animal.name+"/tif_files/*"))
        for file_name in sorted(file_names):
            #print file_name.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
       
            self.selected_session = file_name.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
            self.select_session(self.selected_session)        
            
            filter_data(self)

   
    def algn_images(self):
        #ALIGN ALL IMAGES

        file_names = sorted(glob.glob(self.parent.root_dir+self.parent.animal.name+"/tif_files/*"))
        for file_name in sorted(file_names):
            print file_name.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
       
            self.selected_session = file_name.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
            self.select_session(self.selected_session)        
            
            self.animal.sessfilter_data(self)

        self.align_images()                 #Align raw_images to first session frame 1000


        #load specific session .npy images and align to 1000 frame of 1st session
        print "... alligning session ... NOT CURRENTLY IMPLEMENTED..."
        

    def load_m_lever(self):
        self.parent.animal.load()


    #*****************************************************************************************************
    #************************************* PROCESS INDIVIDUAL SESSION ************************************
    #*****************************************************************************************************


    def dff_mouse_lever(self):
        compute_dff_mouse_lever(self)
        
        #self.parent.animal.process_sessions(self)
        

        #Reset trial value
        n_trials_in = self.load_stm_name()
        self.comboBox_select_trial.clear()
        for k in range(n_trials_in):
            self.comboBox_select_trial.addItem(str(k))
            self.selected_trial = '0' #crappy method; redundant


    def cvrt_tif_npy(self):
        print "...convert .tif to .npy..."

        #Check to see if .npy or _aligned.npy file exists
        if (os.path.exists(self.tif_file[:-4] +'.npy')==False) and (os.path.exists(self.tif_file[:-4] +'_aligned.npy')==False):
            print "...read: ", self.tif_file
            images_raw = tiff.imread(self.tif_file)

            print "... saving .npy"
            np.save(self.tif_file[:-4], images_raw)
        else:
            print "... .npy file exists ..."


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
            

    def static_stm_mouse_lever(self):
        view_static_stm(self)


    #def view_stm_activity(self):
    #    make_stm_motion_mask(self)  #NOT USED ANYMORE; WAS DEVELOPED TO TEST OUT MASKING



    def video_stm_mouse_lever(self):
        view_video_stm(self)


    def style_choice2(self, text):
        self.choice2 = text
        
        
    def style_choice3(self, text):
        self.choice3 = text



    #*************************************************************************************************
    #*********************************** VIDEO PROCESSING TOOLS **************************************
    #*************************************************************************************************


    def conv_video(self):
        convert_video(self)


    def fnd_start_end(self):
        find_start_end(self)


    def bl_light_roi(self):
        plot_blue_light_roi(self)


    def evt_movies(self):
        event_triggered_movies_multitrial(self)


    def evt_movies_ca(self):
        event_triggered_movies_single_Ca(self)


    #*************************************************************************************************
    #****************************************** STROKE TOOLS *****************************************
    #*************************************************************************************************

        
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



    def chunk_mouse_lever(self):
        #Call chunking algorithm
        self.parent.animal.chunk_data(pre_stroke_days, post_stroke_days, post_post_stroke_days)

