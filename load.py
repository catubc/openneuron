from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

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


