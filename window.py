import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from load import Load
from event_triggered_ephys import EventTriggeredEphys
from event_triggered_imaging import EventTriggeredImaging
from event_triggered_imaging_mcd import EventTriggeredImagingMCD

from mouse_lever_tools import MouseLeverTools
from cat_tools import CatTools
from rat_tools import RatTools
from msl import MSL
from count_matrix import CountMatrix
from file_dialog import FileDialog
from track import Track
from intan_tools import IntanTools
from seamans import Seamans
from template_tools import TemplateTools
from traces_tools import TracesTools
from lfp import LFP
         
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
        self.animal_name_text='' 
        self.rec_name_text=''
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

        if False: 
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

        #self.animal.rec_length = self.animal.tsf.n_vd_samples/float(self.animal.tsf.SampleFrequency)


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
        ophysTools = QtGui.QAction("&Event Triggered Imaging", self)
        ophysTools.setStatusTip('Event Triggered Imaging')
        ophysTools.triggered.connect(self.ophys_tools)

        ophysToolsMCD = QtGui.QAction("&Event Triggered Imaging - MCD", self)
        ophysToolsMCD.setStatusTip('Event Triggered Imaging MCD')
        ophysToolsMCD.triggered.connect(self.ophys_tools_mcd)

        ephysTools = QtGui.QAction("&Event Triggered Ephys", self)
        ephysTools.setStatusTip('Event Triggered Ephys')
        ephysTools.triggered.connect(self.ephys_tools)
        
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
        View_Traces = QtGui.QAction("&Traces Tools", self)
        View_Traces.setStatusTip('Traces Tools')
        View_Traces.triggered.connect(self.view_rtraces)
        
        View_Templates = QtGui.QAction("&Template Tools", self)
        View_Templates.setStatusTip('Template Tools')
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

        fileMenu = mainMenu.addMenu('Event Triggered Analysis')
        #fileMenu.addAction(Event_Triggered_Maps)
        fileMenu.addAction(ophysTools)
        fileMenu.addAction(ophysToolsMCD)
        fileMenu.addAction(ephysTools)
        fileMenu.addAction(mouseLeverTools)
        #fileMenu.addAction(catTools)
        #fileMenu.addAction(ratTools)

        fileMenu = mainMenu.addMenu('Ephys Tools')
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
    def ophys_tools(self):
        ophys_widget = EventTriggeredImaging(self)
        self.central_widget.addWidget(ophys_widget)  
        self.central_widget.setCurrentWidget(ophys_widget)
        
    def ophys_tools_mcd(self):
        ophys_widget = EventTriggeredImagingMCD(self)
        self.central_widget.addWidget(ophys_widget)  
        self.central_widget.setCurrentWidget(ophys_widget)
        
        

    def ephys_tools(self):
        ephys_widget = EventTriggeredEphys(self)
        self.central_widget.addWidget(ephys_widget)  
        self.central_widget.setCurrentWidget(ephys_widget)

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
        traces_widget = TracesTools(self)
        self.central_widget.addWidget(traces_widget)  
        self.central_widget.setCurrentWidget(traces_widget)
    
    def view_templts(self):
        templates_widget = TemplateTools(self)
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
        
