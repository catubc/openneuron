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
        main_widget.rec_name_text = main_widget.animal.recName.replace(main_widget.animal.home_dir+main_widget.animal.name+'/rhd_files/', '')
        
        self.parent.setWindowTitle(main_widget.animal.name+'/'+main_widget.rec_name_text)

        self.parent.animal.load_tsf_header(main_widget.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','_hp.tsf'))
        self.parent.animal.rec_length = self.parent.animal.tsf.n_vd_samples/float(self.parent.animal.tsf.SampleFrequency)

        #self.default_parameters()



class MouseTools(QtGui.QWidget):
    def __init__(self, parent):
        super(MouseTools, self).__init__(parent)
        self.parent = parent

        self.styleChoice = QtGui.QLabel("Windows Vista", self)

        comboBox = QtGui.QComboBox(self)
        comboBox.addItem("motif")
        comboBox.addItem("Windows")
        comboBox.addItem("cde")
        comboBox.addItem("Plastique")
        comboBox.addItem("Cleanlooks")
        comboBox.addItem("windowsvista")
        comboBox.move(50, 250)
        
        self.styleChoice.move(50,150)
        comboBox.activated[str].connect(self.style_choice)

        self.show()


    def style_choice(self, text):
        print "...text: ", text
        self.styleChoice.setText(text)
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create(text))




class MouseLeverTools(QtGui.QWidget):
    def __init__(self, parent):
        super(MouseLeverTools, self).__init__(parent)
        self.parent = parent
        
        layout = QtGui.QGridLayout()    

        #********* PRE-PROCESSING *************
        
        self.button1 = QPushButton('Pre-Process:')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.pprocess_mouse_lever)
        self.button1_lbl = QLabel('.tif -> .npy; align sessions', self)
        layout.addWidget(self.button1, 0, 0)
        layout.addWidget(self.button1_lbl, 0,1)
        
        self.button1 = QPushButton('Filter (aligned) Images:')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.filter_mouse_lever)
        layout.addWidget(self.button1, 1, 0)

        parent.filter_low = QLineEdit('0.1')
        parent.filter_low.setMaximumWidth(50)
        filter_low_lbl = QLabel('Filter img (low):', self)
        layout.addWidget(filter_low_lbl, 1,1)
        layout.addWidget(parent.filter_low, 1,2)
       
        parent.filter_high = QLineEdit('0.1')
        parent.filter_high.setMaximumWidth(50)
        filter_high_lbl = QLabel('Filter img (high):', self)
        layout.addWidget(filter_high_lbl, 1,3)
        layout.addWidget(parent.filter_high, 1,4)
        
        #DFF Drop down box choices
        self.button2 = QPushButton('Compute DFF:')
        self.button2.setMaximumWidth(200)
        self.button2.clicked.connect(self.dff_mouse_lever)
        layout.addWidget(self.button2, 2, 0)
        self.comboBox1 = QtGui.QComboBox(self)
        self.comboBox1.addItem("Global Average")
        self.comboBox1.addItem("Sliding Window (3sec)")
        self.comboBox1.addItem("User Defined Period")
        
        layout.addWidget(self.comboBox1, 2,1)
        
        self.comboBox1.activated[str].connect(self.style_choice1)
        
        #layout.addItem(QtGui.QSpacerItem(3,4), 0, 4)
        
        #*************POST PROCESSING ****************
        #Load mouse from disk
        self.button1 = QPushButton('Load Mouse')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.load_m_lever)
        layout.addWidget(self.button1, 6, 0)


        #Separate sessions into epoch periods
        self.button1 = QPushButton('Load Epochs (days)')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.load_m_lever)
        layout.addWidget(self.button1, 7, 0)

        pre_stroke_days = 21;  post_stroke_days = 14;  post_post_stroke_days = 42
        
        parent.pre_stroke_days = QLineEdit('21');                   
        parent.pre_stroke_days.setMaximumWidth(50)
        pre_stroke_days_lbl = QLabel('# days pre_stroke:', self)
        layout.addWidget(pre_stroke_days_lbl, 7,1)
        layout.addWidget(parent.pre_stroke_days, 7,2)
        
        parent.post_stroke_days = QLineEdit('14');                     
        parent.post_stroke_days.setMaximumWidth(50)
        post_stroke_days_lbl = QLabel('# days post_stroke:', self)
        layout.addWidget(post_stroke_days_lbl, 7,3)
        layout.addWidget(parent.post_stroke_days, 7,4)

        parent.post_post_stroke_days = QLineEdit('42');               
        parent.post_post_stroke_days.setMaximumWidth(50)
        post_post_stroke_days_lbl = QLabel('# days pp_stroke:', self)
        layout.addWidget(post_post_stroke_days_lbl, 7,5)
        layout.addWidget(parent.post_post_stroke_days, 7,6)        
        
        
        #Dim reduction: options
        self.button3 = QPushButton('Dimension Reduction:')
        self.button3.setMaximumWidth(200)
        self.button3.clicked.connect(self.dim_red_mouse_lever)
        layout.addWidget(self.button3, 8, 0)
        self.comboBox2 = QtGui.QComboBox(self)
        self.comboBox2.addItem("PCA")
        self.comboBox2.addItem("MDS")
        self.comboBox2.addItem("tSNE")
        self.comboBox2.addItem("tSNE Barnes-Hutt")
        
        layout.addWidget(self.comboBox2, 8,1)
        self.comboBox2.activated[str].connect(self.style_choice2)
        

        #KMeans: options
        self.button4 = QPushButton('Kmeans:')
        self.button4.setMaximumWidth(200)
        self.button4.clicked.connect(self.kmeans_mouse_lever)
        layout.addWidget(self.button4, 9, 0)

        parent.kmeans_clusters = QLineEdit('16');                   
        parent.kmeans_clusters.setMaximumWidth(50)
        kmeans_clusters_lbl = QLabel('# of clusters:', self)
        layout.addWidget(kmeans_clusters_lbl, 9,1)
        layout.addWidget(parent.kmeans_clusters, 9,2)
        
        
        #View clusters in dim reduction space
        self.button5 = QPushButton('View Clusters - 3D')
        self.button5.setMaximumWidth(200)
        self.button5.clicked.connect(self.clusters_plot_mouse_lever)
        layout.addWidget(self.button5, 9, 3)
        

        #Select a cluster from 2D trace groupings
        self.button6 = QPushButton('Select Cluster - 2D')
        self.button6.setMaximumWidth(200)
        self.button6.clicked.connect(self.select_cluster_mouse_lever)
        layout.addWidget(self.button6, 10, 0)

        #Plot Cluster Motiff
        self.button7 = QPushButton('Plot Cluster Motiff')
        self.button7.setMaximumWidth(200)
        self.button7.clicked.connect(self.plot_motiff_mouse_lever)
        layout.addWidget(self.button7, 10, 1)
        
              
        #Movies - single trial
        self.button8 = QPushButton('Movies - Single Trial')
        self.button8.setMaximumWidth(200)
        self.button8.clicked.connect(self.movies_single_mouse_lever)
        layout.addWidget(self.button8, 15, 0)
        
              
        #1D plots
        self.button9 = QPushButton('1D - Plots')
        self.button9.setMaximumWidth(200)
        self.button9.clicked.connect(self.plot_1D_mouse_lever)
        layout.addWidget(self.button9, 15, 0)        
        
        #Pixel plots
        self.button10 = QPushButton('Pixel - Plots')
        self.button10.setMaximumWidth(200)
        self.button10.clicked.connect(self.plot_pixel_mouse_lever)
        layout.addWidget(self.button10, 15, 1)   
        
        #Make Movies 
        self.button11 = QPushButton('Make Movies')
        self.button11.setMaximumWidth(200)
        self.button11.clicked.connect(self.make_movies_mouse_lever)
        layout.addWidget(self.button11, 15, 2)   
        
                 
                
        self.setLayout(layout)


    def filter_mouse_lever(self):
        print "...filtering aligned images..."

    def style_choice1(self, text):
        print "...text: ", text
        #self.styleChoice.setText(text)
        #QtGui.QApplication.setStyle(QtGui.QStyleFactory.create(text))

        compute_DFF()   #saves processed data to disk; 
                        #Needs to work with different options
        
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

    def style_choice2(self, text):
        print "...text: ", text
        
        dim_red_data = dim_reduction(TEST)  #Return some arrays; May wish to pack them differently;
        

    def dim_red_mouse_lever(self):
        print "... dim reduction..."
        
        dim_reduction (self.parent.animal, text) 

    def dff_mouse_lever(self):
        print "... dff computation..."
        
        #NB: COMPUTE DFF ON FILTERED DATA IF IT EXISTS; OTHERWISE ALIGNED DATA....

    def pprocess_mouse_lever(self):
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
        
        layout = QtGui.QGridLayout()
     
        parent.start_lfp = QLineEdit('0');                   
        parent.start_lfp.setMaximumWidth(50)
        start_lfp_lbl = QLabel('start_lfp_cluster:', self)
       
        parent.end_lfp = QLineEdit('1');                     
        parent.end_lfp.setMaximumWidth(50)
        end_lfp_lbl = QLabel('end_lfp_cluster:', self)

        parent.time_chunks = QLineEdit('1');               
        parent.time_chunks.setMaximumWidth(50)
        time_chunks_lbl = QLabel('# of time_chunks:', self)
        
        parent.lock_window = QLineEdit('200');               
        parent.lock_window.setMaximumWidth(50)
        lock_window_lbl = QLabel('lock_window (ms):', self)

        parent.n_spikes = QLineEdit('250');               
        parent.n_spikes.setMaximumWidth(50)
        n_spikes_lbl = QLabel('min_spikes:', self)

        #ADD TO LAYOUT
        layout.addWidget(start_lfp_lbl, 0,0)
        layout.addWidget(parent.start_lfp, 0,1)
        layout.addWidget(end_lfp_lbl,0,2)
        layout.addWidget(parent.end_lfp,0,3)
        
        layout.addWidget(time_chunks_lbl,1,0)
        layout.addWidget(parent.time_chunks,1,1)

        layout.addWidget(lock_window_lbl,2,0)
        layout.addWidget(parent.lock_window,2,1)
        
        layout.addWidget(n_spikes_lbl,6,2)
        layout.addWidget(parent.n_spikes,6,3)


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


    def view_msl(self, main_widget):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        msl_plots(self.parent)

            
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
        

class Cell(QtGui.QWidget):
    def __init__(self, parent):
        super(Cell, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        
        self.title = QLineEdit('');                #parent.start_time = self.start_time
        self.title.setMaximumWidth(0)
        self.title_lbl = QLabel('CELL STM TOOLS', self)

        parent.start_time = QLineEdit('-0.2');                #parent.start_time = self.start_time
        parent.start_time.setMaximumWidth(50)
        start_time_lbl = QLabel('start_time:', self)
        start_time_lbl.setMaximumWidth(80)

        parent.end_time = QLineEdit('+0.2');                  #parent.end_time = self.end_time
        parent.end_time.setMaximumWidth(50)
        end_time_lbl = QLabel('end_time:', self)
        end_time_lbl.setMaximumWidth(80)
        
        parent.start_cell = QLineEdit('0');                   #parent.start_cell = self.start_cell
        parent.start_cell.setMaximumWidth(50)
        start_cell_lbl = QLabel('start_cell:', self)
        start_cell_lbl.setMaximumWidth(80)
        
        parent.end_cell = QLineEdit('1');                     #parent.end_cell = self.end_cell
        parent.end_cell.setMaximumWidth(50)
        end_cell_lbl = QLabel('end_cell:', self)
        end_cell_lbl.setMaximumWidth(80)



        #ADD TEXTBOXES TO LAYOUT
        layout.addWidget(self.title_lbl, 0,0)
        layout.addWidget(self.title, 0,1)

        layout.addWidget(start_time_lbl, 1,0)
        layout.addWidget(parent.start_time, 1,1)
        layout.addWidget(end_time_lbl,1,2)
        layout.addWidget(parent.end_time,1,3)
        
        layout.addWidget(start_cell_lbl, 2,0)
        layout.addWidget(parent.start_cell, 2,1)
        layout.addWidget(end_cell_lbl,2,2)
        layout.addWidget(parent.end_cell,2,3)
        
        
        
        #MAKE BUTTONS             
        self.button1 = QPushButton('STM')
        self.button1.setMaximumWidth(200)
        self.button1.clicked.connect(self.view_stm)
        layout.addWidget(self.button1, 4, 0)

        self.button2 = QPushButton('STM Movies')
        self.button2.setMaximumWidth(200)
        layout.addWidget(self.button2, 4, 1)
        self.button2.clicked.connect(self.view_stm_movie)
                
        self.button10 = QPushButton('Compute STM')
        self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.button10, 4, 2)
        self.button10.clicked.connect(self.compute_cell_stm)        

        self.button3 = QPushButton('Cell Rasters')
        self.button3.setMaximumWidth(200)
        layout.addWidget(self.button3, 4, 3)
        self.button3.clicked.connect(self.view_cell_raster)


        self.setLayout(layout)
   
    def view_stm(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
        sta_maps(self.parent)
        
    def view_stm_movie(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        sta_movies(self.parent)

    def view_cell_raster(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        plot_rasters(self.parent)

    def compute_cell_stm(self):
        self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
        compute_sta(self.parent) #NEEDS TO BE FIXED FROM (self.animal, self.animal.ptcsName) TO CURRENT FORM


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
        self.button1 = QPushButton('Concatenate multiple .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_tsf)
        layout.addWidget(self.button1, 5, 0)

        self.button1 = QPushButton('Concatenate .lfp.zip files')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.multi_lfp_zip)
        layout.addWidget(self.button1, 6, 0)


        self.button1 = QPushButton('Time compress .tsf (1Khz)')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.comp_lfp)
        layout.addWidget(self.button1, 7, 0)
        
        self.setLayout(layout)

        

    def multi_tsf(self):

        print "... selecting multiple recordings ..."

        dialog = FileDialog()
        dialog.exec_()
        
        self.parent.animal.tsf_files = dialog.out_files
        
        #animal_analysis.py has this function.
        concatenate_tsf(self) 


    def multi_lfp_zip(self):

        print "... selecting multiple recordings ..."

        dialog = FileDialog()
        dialog.exec_()
        
        self.parent.animal.tsf_files = dialog.out_files
        
        #animal_analysis.py has this function.
        concatenate_lfp_zip(self) 
            

    def comp_lfp(self):

        print "... selecting multiple recordings ..."

        dialog = FileDialog()
        dialog.exec_()
        
        self.parent.animal.tsf_files = dialog.out_files
        
        #animal_analysis.py has this function.
        compress_lfp(self) 


#LFP ANALYSIS TOOLBOX
class Seamans(QtGui.QWidget):
    def __init__(self, parent):
        super(Seamans, self).__init__(parent)
        
        self.parent = parent
        layout = QtGui.QGridLayout()
        
        #parent.start_time = QLineEdit('0.5');                #parent.start_time = self.start_time
        #parent.start_time.setMaximumWidth(50)
        #start_time_lbl = QLabel('low_pass (hz):', self)
        #start_time_lbl.setMaximumWidth(100)

        #parent.end_time = QLineEdit('6.0');                  #parent.end_time = self.end_time
        #parent.end_time.setMaximumWidth(50)
        #end_time_lbl = QLabel('high_pass (hz):', self)
        #end_time_lbl.setMaximumWidth(100)
      
        #layout.addWidget(start_time_lbl, 0,0)
        #layout.addWidget(parent.start_time, 0,1)
        #layout.addWidget(end_time_lbl,0,2)
        #layout.addWidget(parent.end_time,0,3)

        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ncs to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ncs_convert)
        layout.addWidget(self.button1, 5, 0)

        #MAKE BUTTONS             
        self.button1 = QPushButton('Convert .ntt to .tsf')
        self.button1.setMaximumWidth(250)
        self.button1.clicked.connect(self.ntt_convert)
        layout.addWidget(self.button1, 6, 0)

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

        print "... selecting multiple recordings ..."

        self.parent.root_dir = '/media/cat/8TB/in_vivo/seamans/' 


        ncs_to_tsf(self, QtGui.QFileDialog.getOpenFileNames(self.parent, 'Load Experiment', self.parent.root_dir, "*.ncs"))

        #QFileDialog.getOpenFileNames(...)
        #self.QFileDialog.getOpenFileName(self, "Select File", "", "*.ncs")
        
        #animal_analysis.py has this function.
        #ncs_to_tsf(self) 


 

    def ntt_convert(self):

        print "... selecting multiple recordings ..."

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
        
        pass
        

#LFP ANALYSIS TOOLBOX
class LFP(QtGui.QWidget):
    def __init__(self, parent):
        super(LFP, self).__init__(parent)
        #layout = QtGui.QFormLayout()
        layout = QtGui.QGridLayout()

        self.parent = parent
        
        parent.specgram_ch = QLineEdit('55');                #parent.start_time = self.start_time
        parent.specgram_ch.setMaximumWidth(50)
        specgram_ch_lbl = QLabel('specgram_ch:', self)
        specgram_ch_lbl.setMaximumWidth(100)
      
        #ADD TO LAYOUT
        layout.addWidget(specgram_ch_lbl, 0,0)
        layout.addWidget(parent.specgram_ch, 0,1)
                
        self.button3 = QPushButton('LFP Rasters')
        self.button3.setMaximumWidth(200)
        layout.addWidget(self.button3, 6, 0)
        self.button3.clicked.connect(self.view_lfp_raster)

        self.button10 = QPushButton('View Spegram')
        self.button10.setLayoutDirection(QtCore.Qt.RightToLeft)
        layout.addWidget(self.button10, 7, 0)
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
        if False: 
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


        #EXPERIMENT TOOLS MENUS
        mouseTools = QtGui.QAction("&Mouse", self)
        mouseTools.setStatusTip('Mouse')
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

        #PROCESSING MENUS
        trackTools = QtGui.QAction("&Track Tools", self)
        trackTools.setStatusTip('Track Tools')
        trackTools.triggered.connect(self.track_tools)

        preprocessExperiment = QtGui.QAction("&Process Data", self)
        preprocessExperiment.setStatusTip('Process Data')
        preprocessExperiment.triggered.connect(self.prepreprocess_data)


        seamansData = QtGui.QAction("&Seamans' Lab Data", self)
        seamansData.setStatusTip('Seamans Lab Data')
        seamansData.triggered.connect(self.seamans_data)


        filterData = QtGui.QAction("&Filter Data", self)
        filterData.setStatusTip('Filter Data')
        filterData.triggered.connect(self.filter_data)


        #ANALYSIS MENUS
        Cell_Analysis = QtGui.QAction("&Cell STM", self)
        Cell_Analysis.setStatusTip('Cell Analysis')
        Cell_Analysis.triggered.connect(self.cell_analysis)

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

        fileMenu = mainMenu.addMenu('Experiment Tools')
        fileMenu.addAction(mouseTools)
        fileMenu.addAction(mouseLeverTools)
        fileMenu.addAction(catTools)
        fileMenu.addAction(ratTools)


        fileMenu = mainMenu.addMenu('Pre-Process')
        fileMenu.addAction(trackTools)
        fileMenu.addAction(seamansData)
        fileMenu.addAction(filterData)
        #fileMenu.addAction(preprocessExperiment)

        fileMenu = mainMenu.addMenu('Analysis')
        fileMenu.addAction(Cell_Analysis)
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
        mouse_widget = MouseTools(self)
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
    def cell_analysis(self):
        cell_widget = Cell(self)
        self.central_widget.addWidget(cell_widget)  
        self.central_widget.setCurrentWidget(cell_widget)

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


    def seamans_data(self):
        seamans_widget = Seamans(self)
        self.central_widget.addWidget(seamans_widget)      
        self.central_widget.setCurrentWidget(seamans_widget)


    def prepreprocess_data(self):  #NOT CURENTLY USED...
        print "....processing experiment: ", self.fileName
        
        #MUST BE EXPERIMENT SPECIFIC
        

    def filter_data(self):
        filter_widget = Filter(self)
        self.central_widget.addWidget(filter_widget)  
        self.central_widget.setCurrentWidget(filter_widget)
        
        

def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())


run()
