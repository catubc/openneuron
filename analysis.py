#Code for doing simultaneous imaging and ephys analysis

import struct, array, csv, os
import numpy as np
import matplotlib.pyplot as plt
import time
import random
import glob
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib.path import Path
import matplotlib.animation as animation
import scipy.ndimage as ndimage
from scipy.signal import butter, filtfilt, cheby1
from load_intan_rhd_format import *

from numpy import nan

np.set_printoptions(precision=20, threshold=nan, edgeitems=None, suppress=None)

from openglclasses import *     #Custom plotting functions


class Object_empty(object):
    def __init__(self):
        pass
            



class Ptcs(object):
    """Polytrode clustered spikes file neuron record"""
    def __init__(self, file_name):
        
        f = open(file_name, "rb")
        self.sorted_file = file_name
        self.name = file_name
        self.full_path =file_name
        # call the appropriate method:
        self.VER2FUNC = {1: self.readHeader, 2: self.readHeader, 3: self.readHeader}

        self.readHeader(f)
        
        self.nid = []  #Make unique unit id list for loading later.
        
        self.loadData(self.nsamplebytes, f)
        
        f.close()

    def __getstate__(self):
        """Instance methods must be excluded when pickling"""
        d = self.__dict__.copy()
        try: del d['VER2FUNC']
        except KeyError: pass
        return d

    def readHeader(self, f):
        """Read in neuron record of .ptcs file version 3. 'zpos' field was replaced
        by 'sigma' field.
        nid: int64 (signed neuron id, could be -ve, could be non-contiguous with previous)
        ndescrbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment, defaults to 0)
        descr: ndescrbytes of ASCII text
        (padded with null bytes if needed for 8 byte alignment)
        clusterscore: float64
        xpos: float64 (um)
        ypos: float64 (um)
        sigma: float64 (um) (Gaussian spatial sigma)
        nchans: uint64 (num chans in template waveforms)
        chanids: nchans * uint64 (0 based IDs of channels in template waveforms)
        maxchanid: uint64 (0 based ID of max channel in template waveforms)
        nt: uint64 (num timepoints per template waveform channel)
        nwavedatabytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavedata: nwavedatabytes of nsamplebytes sized floats
        (template waveform data, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nwavestdbytes: uint64 (nbytes, keep as multiple of 8 for nice alignment)
        wavestd: nwavestdbytes of nsamplebytes sized floats
        (template waveform standard deviation, laid out as nchans * nt, in uV,
        padded with null bytes if needed for 8 byte alignment)
        nspikes: uint64 (number of spikes in this neuron)
        spike timestamps: nspikes * uint64 (us, should be sorted)
        """

        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.nneurons = int(np.fromfile(f, dtype=np.uint64, count=1)) # nneurons
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes
        self.nsamplebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsamplebytes
        self.samplerate = int(np.fromfile(f, dtype=np.uint64, count=1)) # samplerate
        self.npttypebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # npttypebytes

        self.pttype = f.read(self.npttypebytes).rstrip('\0 ') # pttype

        self.nptchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nptchans
        self.chanpos = np.fromfile(f, dtype=np.float64, count=self.nptchans*2) # chanpos
        self.chanpos.shape = self.nptchans, 2 # reshape into rows of (x, y) coords
        self.nsrcfnamebytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nsrcfnamebytes
        self.srcfname = f.read(self.nsrcfnamebytes).rstrip('\0 ') # srcfname
        # maybe convert this to a proper Python datetime object in the Neuron:
        self.datetime = float(np.fromfile(f, dtype=np.float64, count=1)) # datetime (days)
        self.ndatetimestrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndatetimestrbytes
        self.datetimestr = f.read(self.ndatetimestrbytes).rstrip('\0 ') # datetimestr


    def loadData(self, n_bytes, f):
        #call the appropriate method:
        #self.VER2FUNC = {1: self.read_ver_1, 2:self.read_ver_2, 3:self.read_ver_3}
        self.nsamplebytes = n_bytes
        self.wavedtype = {2: np.float16, 4: np.float32, 8: np.float64}[self.nsamplebytes]

        self.n_units=self.nneurons
        self.units=[None]*self.n_units
        self.uid = [None]*self.n_units  #Unique id for full track sorts
        self.n_sorted_spikes = [None]*self.n_units
        self.ptp=np.zeros((self.n_units), dtype=np.float32)
        self.size = []
        self.maxchan = []
        self.sigma = []
        self.xpos = []
        self.ypos = []
        self.wavedata = []

        for k in range(self.n_units):
            self.readUnit(f)
            self.units[k]= self.spikes

            if 'martin' in self.full_path:
                self.uid[k]= self.nid
            else: #All other sorts are from Nick's SS so should be the same
                self.uid[k]= self.nid-1
               
            #print "SAMPLERATE: ", self.samplerate
            #if ptcs_flag: #Martin's data has wrong flag for saves
            
            #CONVERT UNITS TO TIMESTEPS
            #self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps
            
            
            #else:
            #    self.units[k]=[x*self.samplerate/2/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps

            self.n_sorted_spikes[k] = len(self.units[k])
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
            self.sigma.append(self.zps)
            self.xpos.append(self.xps)
            self.ypos.append(self.yps)
            self.wavedata.append(self.wvdata)

            #self.ptp[k]=max(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) - \
            #            min(self.wavedata[np.where(self.chans==self.maxchanu)[0][0]]) #compute PTP of template;


        f.close()


    def readUnit(self,f):
        self.nid = int(np.fromfile(f, dtype=np.int64, count=1)) # nid
        self.ndescrbytes = int(np.fromfile(f, dtype=np.uint64, count=1)) # ndescrbytes
        self.descr = f.read(self.ndescrbytes).rstrip('\0 ') # descr

        if self.descr:
            try:
                self.descr = eval(self.descr) # might be a dict
            except: pass

        self.clusterscore = float(np.fromfile(f, dtype=np.float64, count=1)) # clusterscore
        self.xps = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.yps = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zps = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um) #Replaced by spatial sigma
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) #NB: Some errors here from older .ptcs formats
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid

        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wvdata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes

        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes) # spike timestamps (us):

        # convert from unsigned to signed int for calculating intervals:
        #self.spikes = np.asarray(self.spikes, dtype=np.float64)

            
    def read_wave(self, f):
        """Read wavedata/wavestd bytes"""
        # nwavedata/nwavestd bytes, padded:
        nbytes = int(np.fromfile(f, dtype=np.uint64, count=1))
        fp = f.tell()
        count = nbytes // self.nsamplebytes # trunc to ignore any pad bytes
        X = np.fromfile(f, dtype=self.wavedtype, count=count) # wavedata/wavestd (uV)
        if nbytes != 0:
            X.shape = self.nchans, self.nt # reshape
        f.seek(fp + nbytes) # skip any pad bytes
        return nbytes, X

    def rstrip(s, strip):
        """What I think str.rstrip should really do"""
        if s.endswith(strip):
            return s[:-len(strip)] # strip it
        else:
            return s

    def read(self):
        self.nid = self.parse_id()
        with open(self.fname, 'rb') as f:
            self.spikes = np.fromfile(f, dtype=np.int64) # spike timestamps (us)
        self.nspikes = len(self.spikes)

class Tsf_file(object):

    def __init__(self, file_name):
        
        self.read_header(file_name)
        
    def read_header(self, file_name):
        
        self.fin = open(file_name, "rb")
        
        self.header = self.fin.read(16)
        self.iformat = struct.unpack('i',self.fin.read(4))[0] 
        self.SampleFrequency = struct.unpack('i',self.fin.read(4))[0] 
        self.n_electrodes = struct.unpack('i',self.fin.read(4))[0] 
        self.n_vd_samples = struct.unpack('i',self.fin.read(4))[0] 
        self.vscale_HP = struct.unpack('f',self.fin.read(4))[0] 

        if self.iformat==1001:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Siteloc = struct.unpack(str(2*self.n_electrodes)+'h', self.fin.read(2*self.n_electrodes*2))
        if self.iformat==1002:
            self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
            for i in range(self.n_electrodes):
                self.Siteloc[i*2] = struct.unpack('h', self.fin.read(2))[0]
                self.Siteloc[i*2+1] = struct.unpack('h', self.fin.read(2))[0]
                self.Readloc[i] = struct.unpack('i', self.fin.read(4))[0]

    def read_ec_traces(self):
        print " ... reading data, #chs: ", self.n_electrodes, " nsamples: ", self.n_vd_samples, " len: ", float(self.n_vd_samples)/float(self.SampleFrequency), " sec."
        self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        self.ec_traces.shape = self.n_electrodes, self.n_vd_samples

        self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 
        
        #print "No. ground truth cell spikes: ", self.n_cell_spikes
        if (self.n_cell_spikes>0):
            if (self.iformat==1001):
                self.vertical_site_spacing = struct.unpack('i',self.fin.read(4))[0] 
                self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 

            self.fake_spike_times =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_assignment =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            self.fake_spike_channels =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
        
        self.fin.close()

    def read_trace(self, channel):
        #Load single channel 

        indent = 16+20+self.n_electrodes*8

        self.fin.seek(indent+channel*2*self.n_vd_samples, os.SEEK_SET)         #Not 100% sure this indent is correct.
        self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_vd_samples)
        self.fin.close()
    
    def save_tsf(self, file_name):
        
        fout = open(file_name, 'wb')
        print "...saving: ",  file_name
        fout.write(self.header)
        fout.write(struct.pack('i', self.iformat))
        fout.write(struct.pack('i', self.SampleFrequency))
        fout.write(struct.pack('i', self.n_electrodes))
        fout.write(struct.pack('i', self.n_vd_samples))
        fout.write(struct.pack('f', self.vscale_HP))
        for i in range (self.n_electrodes):
            fout.write(struct.pack('h', self.Siteloc[i*2]))
            fout.write(struct.pack('h', self.Siteloc[i*2+1]))
            fout.write(struct.pack('i', i+1))                 #CAREFUL, SOME FILES MAY USE ReadLoc values..

        self.ec_traces.tofile(fout)

        fout.write(struct.pack('i', self.n_cell_spikes))

        try:
            self.subsample
        except NameError:
            self.subsample = 1.0

        fout.write(struct.pack('i', self.subsample))

        fout.close()



class Probe(object):      

    def __init__(self):

        print "...loading probe..."

        self.name = "NeuroNexus 64Ch probe"         #Hardwired, but should add options here...
    
        self.load_layout()
    
    def load_layout(self):
        ''' Load intan probe map layout 
        '''

        self.n_electrodes = 64

        #Fixed location array for NeurNexus probe layotus
        self.Siteloc = np.zeros((self.n_electrodes,2), dtype=np.int16) #Read as 1D array
        for i in range (self.n_electrodes):
            self.Siteloc[i][0]=30*(i%2)
            self.Siteloc[i][1]=i*23


        #A64 Omnetics adaptor
        adaptor_map = []
        adaptor_map.append([34,35,62,33,60,54,57,55,10,8,11,5,32,3,30,31])
        adaptor_map.append([64,58,63,56,61,59,52,50,15,13,6,4,9,2,7,1])
        adaptor_map.append([53,51,49,47,45,36,37,38,27,28,29,20,18,16,14,12])
        adaptor_map.append([48,46,44,42,40,39,43,41,24,22,26,25,23,21,19,17])

        adaptor_layout1=[]      #Concatenated rows
        for maps in adaptor_map:
            adaptor_layout1.extend(maps)

        #Intan adapter - if inserted right-side up
        intan_map = []
        intan_map.append(list(reversed([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])))     #NB: need to reverse these arrays:  list(reversed(...))
        intan_map.append(list(reversed([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])))
        intan_map.append(list(reversed([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])))
        intan_map.append(list(reversed([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])))

        intan_layout1=[]
        for maps in intan_map:
            intan_layout1.extend(maps)


        #Intan adapter - if inserted upside-down; no need to reverse
        intan_map = []
        intan_map.append([48,50,52,54,56,58,60,62,0,2,4,6,8,10,12,14])
        intan_map.append([49,51,53,55,57,59,61,63,1,3,5,7,9,11,13,15])
        intan_map.append([47,45,43,41,39,37,35,33,31,29,27,25,23,21,19,17])
        intan_map.append([46,44,42,40,38,36,34,32,30,28,26,24,22,20,18,16])

        intan_layout2=[]
        for maps in intan_map:
            intan_layout2.extend(maps)


        #A1x64 probe layout
        a = [27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,28,29,30,31,32]
        b = [37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,36,35,34,33]
        probe_map = a+b
        probe_map[::2] = a
        probe_map[1::2] = b
        
        self.layout = []
        for i in range(len(probe_map)):
            self.layout.append(intan_layout1[adaptor_layout1.index(probe_map[i])])
 

#DUPLICATE FUNCTION WITH TSF CLASS FUNCTION; May still need it for stand alone functions; but LIKELY OBSOLETE... ERASE!!!!!!!!!!!!
def save_tsf_single(tsf, file_name):
    
    fout = open(file_name, 'wb')
    print file_name
    fout.write(tsf.header)
    fout.write(struct.pack('i', tsf.iformat))
    fout.write(struct.pack('i', tsf.SampleFrequency))
    fout.write(struct.pack('i', tsf.n_electrodes))
    fout.write(struct.pack('i', tsf.n_vd_samples))
    fout.write(struct.pack('f', tsf.vscale_HP))
    for i in range (tsf.n_electrodes):
        fout.write(struct.pack('h', tsf.Siteloc[i*2]))
        fout.write(struct.pack('h', tsf.Siteloc[i*2+1]))
        fout.write(struct.pack('i', i+1))                 #CAREFUL, SOME FILES MAY USE ReadLoc values..

    tsf.ec_traces.tofile(fout)

    fout.write(struct.pack('i', tsf.n_cell_spikes))
    fout.close()

def convert_bin_to_npy(self):
    
    filename = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording
    print filename
    if os.path.exists(filename+'.bin')==False: 
        print "...loading .bin file..."
        data = np.fromfile(file_out+'.bin', dtype=np.int16)
        print "...reshaping array..."
        data = data.reshape((-1, 128, 128))
        print "...saving .npy array..."
        np.save(file_out, data)
    else:
        print "... .npy file already exists..."

def convert_video(self):
    
    print "Loading file: ", self.selected_session
    import cv2

    print self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v'
    
    vid = cv2.VideoCapture(self.parent.animal.home_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.m4v')
    length = vid.get(cv2.cv.CV_CAP_PROP_FRAME_COUNT)
    width  = vid.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH)
    height = vid.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)
    fps    = vid.get(cv2.cv.CV_CAP_PROP_FPS)

    #NB: CAN ALSO INDEX DIRECTLY INTO .WMV FILES:
    #time_length = 30.0
    #fps=25
    #frame_seq = 749
    #frame_no = (frame_seq /(time_length*fps))

    ##The first argument of cap.set(), number 2 defines that parameter for setting the frame selection.
    ##Number 2 defines flag CV_CAP_PROP_POS_FRAMES which is a 0-based index of the frame to be decoded/captured next.
    ##The second argument defines the frame number in range 0.0-1.0
    #cap.set(2,frame_no);

    print length, width, height, fps
    if length==0: print "... no movie file... returning..."; return
        
    time.sleep(.5)
    #Show video
    if True:
        data = []
        ctr=0
        now = time.time()
        while True:
            vid.grab()
            retval, image = vid.retrieve()
            if not retval: break
            #cv2.imshow("Test", image)
            #cv2.waitKey(1)
            image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            
            data.append(image)
            ctr+=1; print ctr

        np.save(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session, data)


def find_start_end(self):
    
    self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_blue_light_frames.npy'
    
    #if os.path.exists(self.blue_light_filename)==True: 
    #    print "...Blue Light Boundaries already found... returning..."
    #    return

    global coords, images_temp, ax, fig, cid
    
    #Re-Compute frames per second relative session.reclength; ensure makes sense ~15Hz
    #First, load movie data
    movie_data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy')
    print movie_data.shape
    
    #Second, load imaging data
    temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session
    
    print "...loading .npy img file..."
    data = np.load(temp_file+'_aligned.npy')
    print data.shape

    #********Check to see if img_rate was ~30.00Hz; otherwise skip
    self.img_rate = np.load(temp_file+'_img_rate.npy') #LOAD IMG_RATE
    self.abstimes = np.load(temp_file+'_abstimes.npy')
    
    print "...movie fps: ", float(movie_data.shape[0])/self.abstimes[-1]
    
    #Select pixels from 2 pics: 30th frame and 300th frame; should be able to see light differences
    images_temp = movie_data[300]
    coords=[]
    self.coords = coords
    fig = plt.figure()
    
    ax=plt.subplot(1,2,1)
    plt.imshow(movie_data[30],cmap=plt.get_cmap('gray') )

    ax = plt.subplot(1,2,2)
    ax.imshow(images_temp, cmap=plt.get_cmap('gray') )
    ax.set_title("Click blue light")
    cid = fig.canvas.mpl_connect('button_press_event', on_click_single_frame)
    plt.show()
    

def plot_blue_light_roi(self):
    
    plotting = True
    movie_filename = self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy'
    movie_data = np.load(movie_filename)
    print movie_data.shape
    
    t = np.arange(0,movie_data.shape[0],1)/15.
    blue_light_roi = []
    for k in range(len(movie_data)):
        blue_light_roi.append(np.sum(movie_data[k, int(self.coords[0][0])-1:int(self.coords[0][0])+1, 
                                int(self.coords[0][1])-1:int(self.coords[0][1])+1]))
    
    blight_ave = np.average(blue_light_roi)
    print "std...", self.blue_light_std.text()
    blight_std = np.std(blue_light_roi)*float(self.blue_light_std.text())
    print blight_std
    
    if plotting: 
        ax = plt.subplot(2,1,1)
        plt.plot(t, blue_light_roi)
        plt.plot([t[0],t[-1]], [blight_ave,blight_ave], color='red')
        plt.plot([t[0],t[-1]], [blight_ave-blight_std,blight_ave-blight_std], color='green')
    
    blue_light_roi = np.int32(blue_light_roi)
    indexes = np.where(blue_light_roi>(blight_ave-blight_std))[0]
    blue_light_roi= blue_light_roi[indexes[0]:indexes[-1]]  #SOME OF THE FRAMES DIP BELOW SO JUST ASSUME OK 
    
    print indexes[0], indexes[-1]
    print "...no of frames w. blue light: ", len(indexes)
    movie_rate = float(len(blue_light_roi))/self.abstimes[-1]
    print "...movie fps: ", movie_rate
    
    if abs(15-movie_rate)>0.02: 
        print "************* movie frame rate incorrect *************"
    else:
        #Save frames for blue light
        np.save(self.blue_light_filename, np.arange(indexes[0],indexes[-1],1))
        vid_rate_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_vid_rate.npy'
        np.savetxt(vid_rate_filename, [movie_rate])

    if plotting: 
        ax = plt.subplot(2,1,2)
        plt.plot(blue_light_roi)
        plt.ylim(bottom=0)
        plt.show()
    
    

def event_triggered_movies_multitrial(self):
    """ Make multiple mouse lever pull trials video"""
    
    #Load imaging data
    temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session    
    self.img_rate = np.load(temp_file+'_img_rate.npy') #LOAD IMG_RATE
    self.abstimes = np.load(temp_file+'_abstimes.npy')

    self.abstimes = np.load(temp_file+'_abstimes.npy')
    self.abspositions = np.load(temp_file+'_abspositions.npy')
    self.abscodes = np.load(temp_file+'_abscodes.npy')
    self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
    self.code_44threshold = np.load(temp_file+'_code44threshold.npy')

    print self.abspositions
    print self.abscodes
    print self.locs_44threshold
    print self.code_44threshold


    #Load original .npy movie data and index only during blue_light_frames
    movie_data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy')
    self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'_blue_light_frames.npy'
    self.movie_data = movie_data[np.load(self.blue_light_filename)]

    print "... movie_data.shape: ", self.movie_data.shape
    
    movie_times = np.linspace(0, self.abstimes[-1], self.movie_data.shape[0])
    print movie_times
    
    #NB: indexes are not just for '04' codes but for value in selected_code 
    indexes_04 = np.where(self.code_44threshold==self.selected_code)
    times_04 = self.locs_44threshold[indexes_04]
    print times_04
    
    #Find movie frame centred on time_index
    movie_04frame_locations = []
    for time_index in times_04: 
        movie_04frame_locations.append(find_nearest(movie_times, time_index))
    
    self.movie_04frame_locations = movie_04frame_locations
    
    print "... frame event triggers: ", self.movie_04frame_locations
    
    #Load original .npy movie data and index only during blue_light_frames
    print "... movie_data.shape: ", self.movie_data.shape
  
    temp_img_rate = 15
    self.movie_stack = []
    for frame in self.movie_04frame_locations:
        self.movie_stack.append(self.movie_data[frame-3*temp_img_rate: frame+3*temp_img_rate])



    make_movies_from_triggers(self)

def make_movies_from_triggers(self):
    
    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    self.movie_stack = self.movie_stack[:int(self.n_trials_movies.text())]
    
    n_cols = int(np.sqrt(len(self.movie_stack))+0.999999); n_rows = n_cols-1    #Assume on less row needed (always the case unless perfect square
    if (n_rows*n_cols)<(len(self.movie_stack)): n_rows = n_cols                 #If not enough rows, add 1

    fig = plt.figure()
    im=[]
    for k in range(len(self.movie_stack)):
        im.append([])
        ax = plt.subplot(n_rows,n_cols,k+1)
        
        ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
        
        im[k] = plt.imshow(self.movie_stack[k][0], cmap=plt.get_cmap('gray'), interpolation='none')
        
    def updatefig(j):
        print j
        plt.suptitle("Frame: "+str(j)+"  " +str(format(float(j)/15-3.,'.2f'))+"sec", fontsize = 20)

        # set the data in the axesimage object
        for k in range(len(self.movie_stack)):
            im[k].set_array(self.movie_stack[k][j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.movie_stack[0])), interval=100, blit=False, repeat=True)

    if True:
        ani.save(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'_'+str(len(self.movie_stack))+'.mp4', writer=writer)
    plt.show()


def event_triggered_movies_single_Ca(self):
    """ Load [Ca] imaging and behavioural camera data and align to selected trial"""

    self.parent.n_sec = float(self.n_sec_window.text())
    #**************************************
    #Read [Ca] data
    #**************************************
    temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session    
    self.img_rate = np.load(temp_file+'_img_rate.npy') #imaging rate
    
    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_stm.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_stm.npy"
    print "...stm_name: ", self.traces_filename
    
    data = np.load(self.traces_filename, mmap_mode='r+')
    print data.shape
    
    #Mask [Ca] stack
    print "...selected trial for stm: ", self.selected_trial
    self.ca_stack = quick_mask(self, data[int(self.selected_trial)])
    self.start_time = -self.parent.n_sec; self.end_time = self.parent.n_sec
    
    
    #**************************************
    #Load behaviour camera data
    #**************************************
    vid_rate_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_vid_rate.npy'
    self.vid_rate = np.loadtxt(vid_rate_filename)

    self.abstimes = np.load(temp_file+'_abstimes.npy')
    self.locs_44threshold = np.load(temp_file+'_locs44threshold.npy')
    self.code_44threshold = np.load(temp_file+'_code44threshold.npy')
    print self.locs_44threshold
    print self.code_44threshold

    print "...self.selected_code: ", self.selected_code
    print "...self.selected_trial: ", self.selected_trial
    
    indexes = np.where(self.code_44threshold==self.selected_code)[0]
    print indexes
    self.selected_locs_44threshold = self.locs_44threshold[indexes][int(self.selected_trial)]
    self.selected_code_44threshold = self.code_44threshold[indexes][int(self.selected_trial)]

    #Load original movie data and index only during blue_light_frames
    movie_data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy')
    self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_blue_light_frames.npy'
    self.movie_data = movie_data[np.load(self.blue_light_filename)]
    
    #Find movie frame corresponding to lever pull trigger
    movie_times = np.linspace(0, self.abstimes[-1], self.movie_data.shape[0])
    self.movie_04frame_locations = find_nearest(movie_times, self.selected_locs_44threshold)
    print "... frame event triggers: ", self.movie_04frame_locations

    #Make movie stack
    self.movie_stack = self.movie_data[self.movie_04frame_locations-self.parent.n_sec*self.vid_rate: self.movie_04frame_locations+self.parent.n_sec*self.vid_rate]

    #Interpolate movie stack to match [Ca] imaging rate
    new_stack = []
    for frame in range(len(self.movie_stack)-1):
        new_stack.append(self.movie_stack[frame])
        new_stack.append((np.int16(self.movie_stack[frame])+np.int16(self.movie_stack[frame+1]))/2.)
        
    new_stack.append(self.movie_stack[-1]);  new_stack.append(self.movie_stack[-1])
    self.movie_stack = np.uint8(new_stack)
    
    print self.ca_stack.shape
    print self.movie_stack.shape

    make_movies_ca(self)


def make_movies_ca(self):
    
    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
  
    fig = plt.figure()
    im = []

    #[Ca] stack
    ax = plt.subplot(2,1,1)
    v_max = np.nanmax(np.ma.abs(self.ca_stack)); v_min = -v_max
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.ca_stack[0], vmin=v_min, vmax = v_max, cmap=plt.get_cmap('jet'), interpolation='none'))

    #Camera stack
    ax = plt.subplot(2,1,2)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.movie_stack[0], cmap=plt.get_cmap('gray'), interpolation='none'))

    def updatefig(j):
        print "...frame: ", j
        plt.suptitle(self.selected_dff_filter+'  ' +self.dff_method + "\nFrame: "+str(j)+"  " +str(format(float(j)/self.img_rate-self.parent.n_sec,'.2f'))+"sec", fontsize = 15)

        # set the data in the axesimage object
        im[0].set_array(self.ca_stack[j])
        im[1].set_array(self.movie_stack[j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.movie_stack)), interval=100, blit=False, repeat=True)

    if True:
        ani.save(self.parent.root_dir+self.parent.animal.name+"/movie_files/"+self.selected_session+'_'+str(len(self.movie_stack))+'_'+str(self.selected_trial)+'trial.mp4', writer=writer)
    plt.show()



def filter_data(self):
    """ Filter _aligned.npy files for lever_pull analysis.  
        NB: mean value of stack is added back into filtered data - so it isn't a purely filtered 
    """
    
    plotting = True
    #self.filter_list = ['No Filter', 'Butterworth', 'Chebyshev']

    rec_name = self.selected_session.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')
    images_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+rec_name+'/'+rec_name+'_aligned.npy'
    
    filter_type = self.selected_filter
    lowcut = float(self.parent.filter_low.text())
    highcut = float(self.parent.filter_high.text())
    fs = self.parent.animal.img_rate
    print "... frame rate: ", fs, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

    #Check to see if data already exists
    if os.path.exists(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'):
        print "...data already filtered..."
        return

    #Load aligned images
    print "... loading aligned imgs..."
    images_aligned = np.load(images_file)
    
    #Save mean of images_aligned if not already done
    if os.path.exists(images_file[:-4]+'_mean.npy')==False: 
        images_aligned_mean = np.mean(images_aligned, axis=0)
        np.save(images_file[:-4]+'_mean', images_aligned_mean)
    else:
        images_aligned_mean = np.load(images_file[:-4]+'_mean.npy')
            
    #Load mask - filter only datapoints inside mask
    n_pixels = len(images_aligned[0])
    generic_coords = np.loadtxt(self.parent.animal.home_dir + self.parent.animal.name+'/genericmask.txt')
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)): generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Filter selection and parameters
    if filter_type == 'butterworth':
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        order = 2
        b, a = butter(order, [low, high], btype='band')
    elif filter_type == 'chebyshev':
        nyq = fs / 2.0
        order = 4
        rp = 0.1
        Wn = [lowcut / nyq, highcut / nyq]
        b, a = cheby1(order, rp, Wn, 'bandpass', analog=False)
    
    
    #Load individual pixel time courses; SWITCH TO UNRAVEL HERE****
    import time
   
    filtered_array = np.zeros(images_aligned.shape, dtype=np.float16)
    now = time.time(); start_time = now
    cutoff=n_pixels
    for p1 in range(n_pixels):
        print "...row: ", p1, " ... time: ", time.time()-now,
        now=time.time(); n_pixels_in=0
        for p2 in range(n_pixels):
            if generic_mask_indexes[p1,p2]==False:
                filtered_array[:,p1,p2] = np.float16(filtfilt(b, a, images_aligned[:,p1,p2])); n_pixels_in+=1   #filter pixel inside mask
        
        print " # pixels filtered: ", n_pixels_in
    print "...total filter time: ", time.time()-start_time

    plt.imshow(filtered_array[1000]); plt.show()        #Check the 1000th frame see what it looks like

    print "... saving filtered data...", filtered_array.shape
    np.save(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz', filtered_array+np.float16(images_aligned_mean))
    print "...DONE..."


def filter_for_event_trigger_analysis(self):
    """ Filter a single file of data. 
        NB: Mean of data saved separately - i.e. does NOT added it back to filtered data.
        NB: Very similar fucntion to filter_data; may wish to combine syntax eventually
    """
        
    main_dir = os.path.dirname(os.path.dirname(self.selected_recording)[:-1])   #Strip file name and 'tif_files' directory 
    
    #Load parameters
    images_file = self.selected_recording[:-4]+'.npy'
    filter_type = self.selected_filter
    lowcut = float(self.lowcut.text())
    highcut = float(self.highcut.text())
    fs = np.loadtxt(main_dir+'/img_rate.txt')     #imaging rate
    print "... frame rate: ", fs, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

    #Check to see if data already exists
    if os.path.exists(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'):
        print "...data already filtered..."; return

    #Load aligned images
    print "... loading images..."
    images = np.load(images_file)
    print images.shape
    
    
    #Save mean of images_aligned if not already done
    if os.path.exists(images_file[:-4]+'_mean.npy')==False: 
        images_mean = np.mean(images, axis=0)
        np.save(images_file[:-4]+'_mean', images_mean)
    else:
        images_mean = np.load(images_file[:-4]+'_mean.npy')


    #Load mask - filter only datapoints inside mask
    n_pixels = len(images[0])
    generic_coords = np.loadtxt(main_dir+'/genericmask.txt')
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)): generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Filter selection and parameters
    if filter_type == 'butterworth':
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        order = 2
        b, a = butter(order, [low, high], btype='band')
    elif filter_type == 'chebyshev':
        nyq = fs / 2.0
        order = 4
        rp = 0.1
        Wn = [lowcut / nyq, highcut / nyq]
        b, a = cheby1(order, rp, Wn, 'bandpass', analog=False)
    
    
    #Load individual pixel time courses; SWITCH TO UNRAVEL HERE****
    import time

    pixels = []
    for p1 in range(n_pixels):
        for p2 in range(n_pixels):
            if generic_mask_indexes[p1,p2]==False:
                pixels.append(images[:,p1,p2])
    
    print "...filtering parallel..."
    start_time = time.time()

    #Use parmap
    import parmap
    filtered_pixels = parmap.map(do_filter, pixels, b, a, processes=30)

    #Use multiprocessing pool:
    #import multiprocessing
    #pool = multiprocessing.Pool() #use all available cores, otherwise specify the number you want as an argument
    #for i in range(len(pixels)):
    #    pool.apply_async(do_filter, args=(pixels[k],b,a))
    #pool.close()
    #pool.join()

    print "...total filter time: ", time.time()-start_time

    print "...reconstructing data ..."
    ctr=0
    filtered_array = np.zeros(images.shape, dtype=np.float16)
    for p1 in range(n_pixels):
        for p2 in range(n_pixels):
            if generic_mask_indexes[p1,p2]==False:
                filtered_array[:,p1,p2] = filtered_pixels[ctr]; ctr+=1

    plt.imshow(filtered_array[1000]); plt.show()        #Check the 1000th frame see what it looks like

    print "... saving filtered data...", filtered_array.shape
    np.save(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz', filtered_array)
    print "...DONE..."

def do_filter(pixel, b, a):
    return np.float16(filtfilt(b, a, pixel))



def filter_single_file(self):
    """ Filter a single file of data. 
        NB: Mean of data saved separately - i.e. does NOT added it back to filtered data.
        NB: Very similar fucntion to filter_data; may wish to combine syntax eventually
    """
        
    plotting = False
    #self.filter_list = ['No Filter', 'Butterworth', 'Chebyshev']

    #Load parameters
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter
    lowcut = float(self.parent.filter_low.text())
    highcut = float(self.parent.filter_high.text())
    fs = np.loadtxt(self.parent.root_dir+self.selected_animal+'/img_rate.txt')
    print "... frame rate: ", fs, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

    #Check to see if data already exists
    if os.path.exists(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'):
        print "...data already filtered..."
        return

    #Load aligned images
    print "... loading aligned imgs..."
    images_aligned = np.load(images_file)
    
    #Save mean of images_aligned if not already done
    if os.path.exists(images_file[:-4]+'_mean.npy')==False: 
        images_aligned_mean = np.mean(images_aligned, axis=0)
        np.save(images_file[:-4]+'_mean', images_aligned_mean)
    else:
        images_aligned_mean = np.load(images_file[:-4]+'_mean.npy')


    #Load mask - filter only datapoints inside mask
    n_pixels = len(images_aligned[0])
    generic_coords = np.loadtxt(self.parent.root_dir+self.selected_animal+'/genericmask.txt')
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)): generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Filter selection and parameters
    if filter_type == 'butterworth':
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq
        order = 2
        b, a = butter(order, [low, high], btype='band')
    elif filter_type == 'chebyshev':
        nyq = fs / 2.0
        order = 4
        rp = 0.1
        Wn = [lowcut / nyq, highcut / nyq]
        b, a = cheby1(order, rp, Wn, 'bandpass', analog=False)
    
    
    #Load individual pixel time courses; SWITCH TO UNRAVEL HERE****
    import time
   
    filtered_array = np.zeros(images_aligned.shape, dtype=np.float16)
    now = time.time(); start_time = now
    cutoff=n_pixels
    for p1 in range(n_pixels):
        print "...row: ", p1, " ... time: ", time.time()-now,
        now=time.time(); n_pixels_in=0
        for p2 in range(n_pixels):
            if generic_mask_indexes[p1,p2]==False:
                filtered_array[:,p1,p2] = np.float16(filtfilt(b, a, images_aligned[:,p1,p2])); n_pixels_in+=1   #filter pixel inside mask
        
        print " # pixels filtered: ", n_pixels_in
    print "...total filter time: ", time.time()-start_time

    plt.imshow(filtered_array[1000]); plt.show()        #Check the 1000th frame see what it looks like

    print "... saving filtered data...", filtered_array.shape
    np.save(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz', filtered_array)
    print "...DONE..."


def view_spontaneous_activity(self):
    
    start_frame = int(self.starting_frame.text())
    n_frames = int(self.number_frame.text())
    print "...start frame: ", start_frame, "   n_frames: ", n_frames

    #Load data from disk
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter; lowcut = float(self.parent.filter_low.text()); highcut = float(self.parent.filter_high.text())
    filtered_file = images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'
        
    data= np.load(filtered_file,  mmap_mode='r+')
    print data.shape

    #Load stack and mean of filtered data
    self.stack = data[start_frame:start_frame+n_frames]
    self.stack_mean = self.stack/np.load(images_file[:-4]+'_mean.npy')

    print self.stack.shape

    make_spontaneous_movies(self)

def view_ave_points(self):

    #Load dim reduced filename and data
    start_frame = int(self.starting_frame.text())
    n_frames = int(self.number_frame.text())
    print "...selected frames: ", start_frame, "   n_frames: ", n_frames
    
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter; lowcut = float(self.parent.filter_low.text()); highcut = float(self.parent.filter_high.text())
    filtered_file = images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'
        
    data= np.load(filtered_file,  mmap_mode='r+')
    print data.shape

    #Select only clustered frames
    selected_frames = np.int16(self.parent.glwindow.glWidget.points_selected)
    self.stack = data[selected_frames]
    self.stack_mean = np.average(self.stack, axis=0)
    plt.imshow(self.stack_mean)
    plt.show()

def view_all_points(self):
    
    #Load dim reduced filename and data
    start_frame = int(self.starting_frame.text())
    n_frames = int(self.number_frame.text())
    print "...selected frames: ", start_frame, "   n_frames: ", n_frames
    
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter; lowcut = float(self.parent.filter_low.text()); highcut = float(self.parent.filter_high.text())
    filtered_file = images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'
        
    data= np.load(filtered_file,  mmap_mode='r+')
    print data.shape

    #Select only clustered frames
    selected_frames = np.int16(self.parent.glwindow.glWidget.points_selected)
    data = data[selected_frames]
    

    n_cols = int(np.sqrt(len(data))+0.999999); n_rows = n_cols-1    #Assume on less row needed (always the case unless perfect square
    if (n_rows*n_cols)<(len(data)): n_rows = n_cols                 #If not enough rows, add 1
    
    v_max = np.nanmax(np.ma.abs(data)); v_min = -v_max; print v_min, v_max
    
    for k in range(len(data)):
        ax = plt.subplot(n_rows,n_cols,k+1)
        ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
              
        #plt.imshow(data[k], vmin=v_min, vmax=v_max, interpolation='none')
        plt.imshow(data[k], interpolation='none')
        plt.title(str(selected_frames[k]),fontsize=10)
    plt.show()
    

def video_points(self):
    
    start_frame = int(self.parent.glwindow.glWidget.points_selected[0])
    n_frames = int(self.number_frame.text())
    
    print "...start frame: ", start_frame, "   n_frames: ", n_frames
    
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter; lowcut = float(self.parent.filter_low.text()); highcut = float(self.parent.filter_high.text())
    filtered_file = images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'
        
    data= np.load(filtered_file,  mmap_mode='r+')
    print data.shape

    #Load stack and mean of filtered data
    self.stack = data[start_frame:start_frame+n_frames]
    self.stack_mean = self.stack/np.load(images_file[:-4]+'_mean.npy')

    print self.stack.shape

    make_spontaneous_movies(self)
    
    
def make_spontaneous_movies(self):
    
    file_movie = self.parent.root_dir+self.selected_animal+"/movie_files/"+self.selected_recording+'.mp4'

    
    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=50, metadata=dict(artist='Me'), bitrate=10000)
  
 
    fig = plt.figure()

    im = [] 
    #raw stack
    ax = plt.subplot(1,2,1)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    v_max1 = np.nanmax(np.ma.abs(self.stack)); v_min1 = -v_max1; print v_min1, v_max1
    im.append(plt.imshow(self.stack[0], vmin = v_min1, vmax=v_max1, cmap=plt.get_cmap('jet'), interpolation='none'))

    #stack + mean frame
    ax = plt.subplot(1,2,2)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    v_max2 = np.nanmax(np.ma.abs(self.stack_mean)); v_min2 = -v_max2; print v_min2, v_max2
    im.append(plt.imshow(self.stack_mean[0], vmin=v_min2, vmax=v_max2, cmap=plt.get_cmap('jet'), interpolation='none'))

    frame_offset = int(self.parent.glwindow.glWidget.points_selected[0])+int(self.starting_frame.text())

    def updatefig(j):
        print j
        plt.suptitle("Frame: "+str(j+frame_offset)+"  " +str(format(float(j)/150,'.2f'))+"sec", fontsize = 15)

        # set the data in the axesimage object
        im[0].set_array(self.stack[j])
        im[1].set_array(self.stack_mean[j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.stack)), interval=100, blit=False, repeat=True)

    if True:
        ani.save(file_movie, writer=writer)
    plt.show()


def compute_dim_reduction(self):
    
    start_frame = int(self.starting_frame.text())
    n_frames = int(self.number_frame.text())
    
    images_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filter_type = self.selected_filter; lowcut = float(self.parent.filter_low.text()); highcut = float(self.parent.filter_high.text())
    self.filtered_file = images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz.npy'
        
    data= np.load(self.filtered_file,  mmap_mode='r+')
    print data.shape

    #Load stack and mean of filtered data
    self.stack = data[start_frame:start_frame+n_frames]
    print self.stack.shape
    self.stack = self.stack.reshape(self.stack.shape[0],-1)
    print self.stack.shape
    
    dim_reduction_stack(self)


def compute_dff_events(self):

    compress_factor = 50.   #Needed to uncompress the LFP compressed sorts

    print "... control computation: ", self.selected_control

    print "\n\n... dff computation event triggers..."

    images_file = self.selected_recording
    global_mean = np.load(images_file[:-4]+'_mean.npy' )
    print "... data_mean.shape: ", global_mean.shape 

    main_dir = os.path.dirname(os.path.dirname(self.selected_recording)[:-1])   #Strip file name and 'tif_files' directory 
    rec_name = self.selected_sort[:-5].replace(main_dir,'').replace('/tsf_files/','')   #Use the sort name - NOT the recording name 

    fs = np.loadtxt(main_dir+'/img_rate.txt')     #imaging rate
    #Check for filtered version of imaging data w. current filtering params
    lowcut = float(self.lowcut.text()); highcut = float(self.highcut.text())
    print "... frame rate: ", fs, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

    #*****************************************************************
    #************* LOAD FILTERED AND UNFILTERED IMAGES ***************
    #*****************************************************************
    #Load unfiltered and filtered data using mmap
    if self.selected_dff_filter!='nofilter': self.filtered_file = images_file[:-4]+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz.npy'
    else: self.filtered_file = images_file[:-4]+'.npy'   #Load unfiltered file and use it as "filtered_file"
    self.images_filtered = np.load(self.filtered_file,  mmap_mode='r+')
    print "...filtered_images_file.shape: ",  self.images_filtered.shape
    
    self.unfiltered_file = images_file[:-4]+'.npy'
    self.images_unfiltered = np.load(self.unfiltered_file,  mmap_mode='r+')

    
    #Check if DFF previously done
    print "... selected filter: ", self.selected_dff_filter
    if self.selected_dff_filter == 'nofilter':
        self.stm_filename = images_file[:-4]+"_"+ self.n_sec_window.text()+"sec_"+ self.dff_method+ '_' + self.selected_dff_filter + ".npy"
    else:
        self.stm_filename = images_file[:-4]+"_"+ self.n_sec_window.text()+"sec_"+ self.dff_method+ '_' + self.selected_dff_filter + '_' + str(lowcut)+"hz_"+str(highcut)+"hz.npy"
    
    print self.stm_filename
    
    if os.path.exists(self.stm_filename): 
        print "... DFF already computed ...skiping processing..."
        return

    #*****************************************************************
    #************** LOAD CAMERA ON/OFF AND IMG_RATE ******************
    #*****************************************************************
    #Load ON/OFF light and compute interpolated img rate based on # of frames aquired divided by ephys time stamps for start and end
    imaging_onoff_file = images_file[:-4].replace('tif_files','camera_files')+'_camera_onoff.npy'
    onoff_pulse = np.load(imaging_onoff_file)
    indexes = np.where(onoff_pulse==1)[0]
    img_start = indexes[0]/25000.         #DON"T HARDWIRE THE AQUISITION RATE; READ FROM RAW .TSF HEADER!!!
    img_end = indexes[-1]/25000.         #DON"T HARDWIRE THE AQUISITION RATE 
    
    self.reclength = img_end - img_start
    print "...imaging start offset, img_end, rec_length: ", img_start, img_end, self.reclength
    
    self.n_images=len(self.images_filtered)

    #Check img_rate; if incorrect, exit
    session_img_rate = self.n_images/self.reclength
    print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate

    if abs(session_img_rate-float(fs))<0.01:         #Compare computed session img_rate w. experimentally set img_rate
        print "Correct img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
        np.save(images_file.replace('_aligned.npy','')+'_img_rate', session_img_rate)
    else:
        print "***Incorrect img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
        np.save(images_file.replace('_aligned.npy','')+'_img_rate', session_img_rate)
        return

    #*****************************************************************
    #************************ LOAD EVENT FILE*************************
    #*****************************************************************
    if '.ptcs' in self.selected_sort:
        Sort = Ptcs(self.selected_sort)
        events = np.float32(Sort.units[int(self.selected_unit)])/Sort.samplerate         
        if 'compressed' in self.selected_sort: 
            events = events*compress_factor 
        events = events - img_start                                 #Align to imaging times by removing pre-imaging period
    else:
        events = np.loadtxt(self.selected_sort)                     #Manual time points are relative to ophys start; 

    print "... original events: ", events[:10], events[-10:]

    
    if self.selected_control =="yes":       #Compute control DFF; select random time points same # as 
        events = np.random.random(len(events))*events[-1]
    
        print "...control events: ", events[:10], events[-10:]
    
    
    
    #*****************************************************************
    #******************* MAKE IMAGING TRIGGER TIMES ******************
    #*****************************************************************
    #Find nearest frame indexes for each event
    trigger_times = events
    frame_times = np.linspace(0, self.reclength, self.n_images)             #Divide reclength by number of images; seconds
    img_frame_triggers = []
    self.window = float(self.n_sec_window.text()) * session_img_rate  #Window width in # frames
    for i in range(len(trigger_times)):
        #img_frame_triggers.append(self.find_previous(frame_times, trigger_times[i])) 
        nearest_frame = find_nearest(frame_times, trigger_times[i])  #Two different flags in the function; CHECK PERIODICALLY 
        if (nearest_frame < (2*self.window)) or (nearest_frame>(self.n_images-float(self.window))): continue  #Skip if too close to start/end; 
        img_frame_triggers.append(nearest_frame)     
        
    print "...# img_frame_triggers...", len(img_frame_triggers)
    

    #*****************************************************************
    #**************************** COMPUTE DFF ************************
    #*****************************************************************
    print "...computing DF/F..."
    n_pixels = len(self.images_filtered[0])
    data_stm = np.zeros((int(self.window*2),n_pixels,n_pixels), dtype=np.float32)
    for trigger in img_frame_triggers:
        print "...frame0: ", trigger

        data_chunk = np.float32(self.images_filtered[int(trigger-self.window):int(trigger+self.window)])[:int(self.window*2)]
        
        if self.dff_method == 'globalAverage':
            if self.selected_dff_filter!='nofilter':  data_stm+=data_chunk/global_mean    #Only need to divide by global mean as original data_chunk did not have mean img added in
            else: data_stm+=(data_chunk-global_mean)/global_mean
            
        elif self.dff_method == 'slidingWindow':            #Use baseline -2*window .. -window
            baseline = np.average(np.float32(self.images_unfiltered[int(trigger-2*self.window):int(trigger-self.window)]), axis=0)
            
            if self.selected_dff_filter!='nofilter': data_stm+=data_chunk/baseline      #ignore subtracting baseline because it was never added back in 
            else: data_stm+=(data_chunk-baseline)/baseline
            
    average_stm = data_stm/len(img_frame_triggers)

    #Save average of event triggered STMs
    print "Saving trial DFF...",

    if self.selected_dff_filter !='nofilter':
        stm_file_name = main_dir + '/stm_files/img_avg_' + rec_name+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window'
    else:
        stm_file_name = main_dir + '/stm_files/img_avg_' + rec_name+'_'+self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window'
        
    
    if self.selected_control =="yes":       stm_file_name = stm_file_name + '_control'

    np.save(stm_file_name, average_stm)

    print ''
    self.images_filtered = 0;     #Set aligned_images to empty 
    
    

def compute_dff_mouse_lever(self):
    print "\n\n... dff computation..."

    #Load average frame
    self.rec_filename = self.selected_session.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')

    images_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_aligned.npy'
    data_mean = np.load(images_file[:-4]+'_mean.npy' )
    print "... data_mean.shape: ", data_mean.shape #; plt.imshow(data_mean); plt.show()
    
    #Check for filtered version of imaging data w. current filtering params
    self.lowcut = float(self.parent.filter_low.text()); self.highcut = float(self.parent.filter_high.text())
    fs = self.parent.animal.img_rate
    print "... frame rate: ", fs, "  low_cutoff: ", self.lowcut, "  high_cutoff: ", self.highcut
    
    print "...self.locs_44threshold: "; print self.locs_44threshold
    print "...self.code_44threshold: "; print self.code_44threshold
    print "...self.selected_code: ", self.selected_code
    
    indexes = np.where(self.code_44threshold==self.selected_code)[0]
    print "...indexes: "; print indexes
    self.code_44threshold_selected = self.code_44threshold[indexes]
    self.locs_44threshold_selected = self.locs_44threshold[indexes]
    
    compute_DFF_function(self)


def compute_DFF_function(self):

    self.parent.n_sec = float(self.n_sec_window.text())
    #Check if already done
    print "... selected filter: ", self.selected_dff_filter
    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+str(self.lowcut)+"hz_"+str(self.highcut)+"hz_"+str(self.selected_code)+"code_traces.npy"
    
    if os.path.exists(self.traces_filename): 
        print "... DFF already computed ...skiping processing..."
        return

    #Load aligned/filtered data and find ON/OFF light; #*******************SKIP THIS IF ALREADY SAVED ON FILE***************
    images_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_aligned.npy'
    self.aligned_images = np.load(images_file)
    blue_light_threshold = 400  #Intensity threshold; when this value is reached - imaging light was turned on
    start_blue = 0; end_blue = len(self.aligned_images)
    
    if np.average(self.aligned_images[0])> blue_light_threshold:    #Case #1: imaging starts with light on; need to remove end chunk; though likely bad recording
        for k in range(len(self.aligned_images)):
            if np.average(self.aligned_images[k])< blue_light_threshold:
                #self.aligned_images = self.aligned_images[k:]
                end_blue = k
                break
    else:                                                           #Case #2: start with light off; remove starting and end chunks;
        #Find first light on
        for k in range(len(self.aligned_images)):
            if np.average(self.aligned_images[k])> blue_light_threshold:
                #self.aligned_images = self.aligned_images[k:]
                start_blue = k
                break

        #Find light off - count backwards from end of imaging data
        for k in range(len(self.aligned_images)-1,0,-1):
            if np.average(self.aligned_images[k])> blue_light_threshold:
                #self.aligned_images = self.aligned_images[:k]
                end_blue= k
                break
                
                
    self.lowcut = float(self.parent.filter_low.text())
    self.highcut = float(self.parent.filter_high.text())
        
    if self.selected_dff_filter == 'nofilter':
        pass; #already loaded nonfiltered self.aligned_images above
    else:
        filtered_filename = images_file[:-4]+'_'+self.selected_dff_filter+'_'+str(self.lowcut)+'hz_'+str(self.highcut)+'hz.npy'
        if os.path.exists(filtered_filename): self.aligned_images = np.load(filtered_filename)
        else: 
            print filtered_filename
            print "***filtered file does not exist*** "
            return
    
    self.aligned_images = self.aligned_images[start_blue:end_blue]
    print "...images_file.shape: ",  self.aligned_images.shape
    
    
    self.n_images=len(self.aligned_images)


    #Check img_rate; First, need to load session, so find index of tif_file name in tif_files.npy
    temp_tif_files = np.load(self.parent.animal.home_dir+self.parent.animal.name+"/tif_files.npy")
    temp_event_files = np.load(self.parent.animal.home_dir+self.parent.animal.name+"/event_files.npy")
    for k in range(len(temp_tif_files)):
        if self.rec_filename in temp_tif_files[k]:
            index = k; break

    self.reclength = self.parent.animal.load_reclength(temp_event_files[index])
    
    session_img_rate = self.n_images/self.reclength
    print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate

    if abs(session_img_rate-float(self.parent.animal.img_rate))<0.01:         #Compare computed session img_rate w. experimentally set img_rate
        print "Correct img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
        np.save(images_file.replace('_aligned.npy','')+'_img_rate', session_img_rate)
    else:
        print "***Incorrect img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
        np.save(images_file.replace('_aligned.npy','')+'_img_rate', session_img_rate)
        return


    #Find times of triggers from lever pull threshold times
    trigger_times = self.locs_44threshold_selected
    frame_times = np.linspace(0, self.reclength, self.n_images)             #Divide up reclength in number of images
    img_frame_triggers = []
    for i in range(len(trigger_times)):
        #img_frame_triggers.append(self.find_previous(frame_times, trigger_times[i])) 
        img_frame_triggers.append(find_nearest(frame_times, trigger_times[i]))     #Two different functions possible here; 
    print "...img_frame_triggers...", img_frame_triggers
    
    
    #BASELINE FOR GLOBAL BASLINE REMOVAL
    mean_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_aligned_mean.npy'
    global_mean = np.load(mean_file)

    self.abstimes = np.load(self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_abstimes.npy')
    self.abspositions = np.load(self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_abspositions.npy')

    print "...computing DF/F..."
    data_stm = []; traces = []; locs = []; codes = []
    counter=-1
    self.window = self.parent.n_sec * session_img_rate
    print "Trigger frame: ", 
    for trigger in img_frame_triggers:
        counter+=1
        print trigger,
        #NB: Ensure enough space for the sliding window; usually 2 x #frames in window
        if trigger < (2*self.window) or trigger>(self.n_images-self.window): 
            continue  #Skip if too close to start/end

        #add locs and codes
        locs.append(self.locs_44threshold_selected[counter])
        codes.append(self.code_44threshold_selected[counter])

        print "...self.dff_method: ", self.dff_method
        #dff_list = ['Global Average', 'Sliding Window: -6s..-3s']
    
        data_chunk = self.aligned_images[int(trigger-self.window):int(trigger+self.window)]

        if self.dff_method == 'globalAverage':
            data_stm.append((data_chunk-global_mean)/global_mean)    #Only need to divide by global mean as original data_chunk did not have mean img added in
            
        elif self.dff_method == 'slidingWindow':            #Use baseline -2*window .. -window
            baseline = np.average(self.aligned_images[int(trigger-2*self.window):int(trigger-self.window)], axis=0)
            data_stm.append((data_chunk-baseline)/baseline)
        
        #***PROCESS TRACES - WORKING IN DIFFERENT TIME SCALE
        lever_window = 120*self.parent.n_sec    #NB: Lever window is computing in real time steps @ ~120Hz; and discontinuous;
        t = np.linspace(-lever_window*0.0082,lever_window*0.0082, lever_window*2)
        lever_position_index = find_nearest(np.array(self.abstimes), self.locs_44threshold[counter])
        
        lever_trace = self.abspositions[int(lever_position_index-lever_window):int(lever_position_index+lever_window)]

        if len(lever_trace)!=len(t):    #Extraplote missing data
            print "...missing lever trace data ... extrapolating..."
            lever_trace = np.zeros(lever_window*2,dtype=np.float32)
            for k in range(-lever_window,lever_window,1):
                lever_trace[k+lever_window] = self.abspositions[k+lever_window]     #Double check this...

        traces.append(lever_trace)

    #Save traces, and 44 threshold locations and codes for trials within boundaries
    np.save(self.traces_filename, traces)
    #np.save(self.tif_file[:-4]+'_locs44threshold', locs)        #MAY WISH TO OVERRIDE ORIGINAL TIMES AND LOCS - FEW TRIALS OUTOF BOUNDS MAY AFFECT OTHER ANALYSIS
    #np.save(self.tif_file[:-4]+'_code44threshold', codes)       #i.e. some of the original values may be out-of-bounds

    #Save individual trial time dynamics
    data_stm = np.float16(data_stm)
    print "Saving trial DFF...",
    np.save(self.traces_filename.replace('_traces.npy','')+'_stm', data_stm)

    print ''
    #Set aligned_images to empty 
    self.aligned_images = []


def view_static_stm_events(self):
    
    block_save = int(self.block_save.text())
    
    control_text=''
    if self.selected_control=="yes": control_text='_control'
    
    print "...control flag: ", self.selected_control, control_text
    
    main_dir = os.path.dirname(os.path.dirname(self.selected_recording)[:-1])   #Strip file name and 'tif_files' directory 
    rec_name = self.selected_sort[:-5].replace(main_dir,'').replace('/tsf_files/','')

    if self.selected_dff_filter !='nofilter':
        stm_file_name = main_dir + '/stm_files/img_avg_' + rec_name+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window'+control_text+'.npy'
    else:
        stm_file_name = main_dir + '/stm_files/img_avg_' + rec_name+'_'+self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window'+control_text+'.npy'
        
    print stm_file_name
    
    data = np.float32(np.load(stm_file_name))
    print data.shape
    
    temp_array = data

    plt.close()
    ax = plt.subplot(1,1,1)
    img_rate = float(np.loadtxt(main_dir+'/img_rate.txt'))
    start_time = float(self.window_start.text()); end_time = float(self.window_end.text()); window_len = float(self.n_sec_window.text())
    
    img_out = []
    for i in range(int(img_rate*(window_len+start_time)),int(img_rate*(window_len+end_time)), block_save):
        img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
    
    #Mask data
    img_out = quick_mask_event(main_dir+'/genericmask.txt', img_out, int(self.midline_mask.text()))
    img_out = np.ma.hstack((img_out))

    #Compute max/min values for non-control runs
    #if self.selected_control=='no':
    if (self.vmin_default.text()=='0.0') and (self.vmax_default.text()=='0.0'):
        self.v_abs = max(np.nanmax(img_out),-np.nanmin(img_out))
        v_min = -self.v_abs; v_max = self.v_abs
    else:
        v_min = float(self.vmin_default.text()); v_max = float(self.vmax_default.text())

    plt.imshow(img_out, vmin = v_min, vmax=v_max)

    plt.ylabel(str(round(self.v_abs*100,2))+"%", fontsize=14)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([])

    #if ctr==(len(select_units)-1): 
    plt.xlabel("Time from spike (sec)", fontsize=25)
    old_xlabel = np.linspace(0,img_out.shape[1], 11)
    new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=18)

    plt.title(stm_file_name)
    #plt.suptitle(animal.ptcsName)
    plt.show()




def view_static_stm(self):
    
    self.parent.n_sec = float(self.n_sec_window.text())

    block_save = int(self.block_save.text())

    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

    filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'
    print "...stm_name: ", filename

    data = np.load(filename)
    print data.shape
    
    #Mask data
    print "...selected trial for stm: ", self.selected_trial
    temp_array = quick_mask(self, data[int(self.selected_trial)])

    plt.close()
    ax = plt.subplot(1,1,1)
    img_rate = self.parent.animal.img_rate
    start_time = -self.parent.n_sec; end_time = self.parent.n_sec
    
    img_out = []
    #for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
    for i in range(0,int(2*img_rate*self.parent.n_sec), block_save):
        print i
        img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
    img_out = np.ma.hstack((img_out))
    
    v_abs = max(np.nanmax(img_out),-np.nanmin(img_out))
    plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)

    plt.ylabel(str(round(v_abs*100,2))+"%", fontsize=14)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([])

    #if ctr==(len(select_units)-1): 
    plt.xlabel("Time from '44' threshold crossing (sec)", fontsize=25)
    old_xlabel = np.linspace(0,img_out.shape[1], 11)
    new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=18)

    plt.title(self.traces_filename)
    #plt.suptitle(animal.ptcsName)
    plt.show()

    
      

def view_video_stm(self):

    self.parent.n_sec = float(self.n_sec_window.text())

    block_save = int(self.block_save.text())

    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

    filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'
    print "...stm_name: ", filename

    data = np.load(filename)
    print data.shape
    
    #Mask data
    print "...selected trial for stm: ", self.selected_trial
    vid_array = quick_mask(self, data[int(self.selected_trial)])
    
    start_time = -self.parent.n_sec; end_time = self.parent.n_sec

    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    
    fig = plt.figure()
    #fig.tight_layout()

    ax = plt.subplot(1,1,1)
    v_max = np.ma.max(np.ma.abs(vid_array)); v_min = -v_max
    print v_max, v_min
    print vid_array.shape
    
    ax.get_xaxis().set_visible(False)
    ax.yaxis.set_ticks([])
    plt.ylabel(str(int(v_max*100))+ ".."+str(int(v_min*100)), labelpad=-1, fontsize=6)
    
    im = plt.imshow(vid_array[0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
        
    #function to update figure
    def updatefig(j):
        print j
        plt.suptitle("  Frame: "+str(j)+"  " +str(format(float(j)/self.parent.animal.img_rate+start_time,'.2f'))+"sec", fontsize = 15)

        # set the data in the axesimage object
        im.set_array(vid_array[j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array)), interval=100, blit=False, repeat=True)

    if True:
    #if save_animation:
        ani.save(filename.replace('tif_files', 'video_files').replace(self.selected_session,'')+'_.mp4', writer=writer)

    plt.show()

    #quit() 

def quick_mask_event(generic_mask_file, data, midline_mask_n_pixels):
        
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        generic_mask_indexes=np.zeros((128,128))
        for i in range(len(generic_coords)):
            generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True
    else:
        print "...generic mask not found..."
        return
    
    #Load midline mask
    for i in range(int(midline_mask_n_pixels)):
        generic_mask_indexes[:,64+int(int(midline_mask_n_pixels)/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    n_pixels = 128
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    for i in range(0, len(data),1):
        temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes, fill=np.nan)
    
    return temp_array
    


def quick_mask(self, data):
        
    generic_mask_file = self.parent.animal.home_dir+self.parent.animal.name + '/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        print "...generic mask not found..."
        return
    
    #Load generic mask
    generic_mask_indexes=np.zeros((128,128))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    #Load midline mask
    for i in range(int(self.midline_mask.text())):
        generic_mask_indexes[:,64+int(int(self.midline_mask.text())/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    n_pixels = 128
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    for i in range(0, len(data),1):
        temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes, fill=np.nan)
    
    return temp_array
    
      
def synchrony_index(data, SampleFrequency, si_limit):

    """Calculate an LFP synchrony index, using potentially overlapping windows of width
    and tres, in sec, from the LFP spectrogram, itself composed of bins of lfpwidth and
    lfptres. Note that for power ratio methods (kind: L/(L+H) or L/H), width and tres are
    not used, only lfpwidth and lfptres. Options for kind are:
    'L/(L+H)': fraction of power in low band vs total power (Saleem2010)
    'L/H': low to highband power ratio (Li, Poo, Dan 2009)
    'cv': coefficient of variation (std / mean) of all power
    'ncv': normalized CV: (std - mean) / (std + mean)
    'nstdmed': normalized stdmed: (std - med) / (std + med)
    'n2stdmean': normalized 2stdmean: (2*std - mean) / (2*std + mean)
    'n3stdmean': normalized 3stdmean: (3*std - mean) / (3*std + mean)
    relative2t0 controls whether to plot relative to t0, or relative to start of ADC
    clock. lim2stim limits the time range only to when a stimulus was presented, i.e. to
    the outermost times of non-NULL din.
    """

    ts = np.arange(0,len(data),1E3)   #set of timestamps, in sec
    t0i, t1i = int(ts[0]), int(ts[-1])
    t0 = t0i
    t1 = t1i
    
    x = data[t0i:t1i] / 1e3 # slice data, convert from uV to mV
    #x = filter.notch(x)[0] # remove 60 Hz mains noise

    lfpwidth=30
    lfptres=5     #Time resolution: bin width for analysis in seconds
    
    lowband = [0.1, 4]; highband = [15,100]     #Martin's suggestions; Saleem 2010?
    #lowband = [0.1, 5]; highband = [15,100]     #variations
    #lowband = [0.1, 4]; highband = [20,100]     #variations

    f0, f1 = lowband
    f2, f3 = highband

    assert lfptres <= lfpwidth
    NFFT = intround(lfpwidth * SampleFrequency)
    noverlap = intround(NFFT - lfptres * SampleFrequency)
    # t is midpoints of timebins in sec from start of data. P is in mV^2?:
    P, freqs, t = mpl.mlab.specgram(x, NFFT=NFFT, Fs=SampleFrequency, noverlap=noverlap)

    # keep only freqs between f0 and f1, and f2 and f3:
    f0i, f1i, f2i, f3i = freqs.searchsorted([f0, f1, f2, f3])
    lP = P[f0i:f1i] # nsubfreqs x nt
    hP = P[f2i:f3i] # nsubfreqs x nt
    lP = lP.sum(axis=0) # nt
    hP = hP.sum(axis=0) # nt
    
    # calculate some metric of each column, ie each width:
    kind = 'L/(L+H)'
    si = lP/(lP + hP) #Saleem 2010

    #Insert initial time point into si index
    si = np.insert(si,0,si[0])
    #si = np.insert(si,len(si),si[-1])
    t = np.insert(t,0,0)
    #t = np.insert(t,len(t),len(t))

    #**** Find periods of synchrony > si_limit = 0.7:
    #si_limit = 0.7
    print "...finding sync periods..."
    sync_periods = []
    in_out = 0
    for i in range(len(t)):
        if si[i]>si_limit:
            if in_out==0:
                ind1=t[i]
                in_out=1
        else:
            if in_out==1:
                ind2=t[i-1]
                sync_periods.append([ind1,ind2])
                in_out=0
    if in_out==1:
        sync_periods.append([ind1,t[-1]])

    return si, t, sync_periods # t are midpoints of bins, from start of ADC clock

def ncs_to_tsf(self, ncs_files):
    """ Neuralynx data conversion """
    from ncs import loadNcs

    print "... .ncs to .tsf ..."
    
    for k in range(len(ncs_files)):  print ncs_files[k]
    
    tsf = Object_empty()
    tsf.header = 'Test spike file '
    tsf.iformat = 1002
    tsf.n_electrodes = len(ncs_files)
    tsf.n_cell_spikes = 0

    tsf.Siteloc = np.zeros((tsf.n_electrodes*2), dtype=np.int16) #Read as 1D array
    for i in range (tsf.n_electrodes):
        tsf.Siteloc[i*2]=0
        tsf.Siteloc[i*2+1]=i*50 #GUESSING each tetrode is 50um apart vertically
            
    tsf.ec_traces = []
    min_samples = 1E14
    for ctr, file_name in enumerate(ncs_files):
        print "...loading: ", file_name
        data = loadNcs(file_name)
        tsf.vscale_HP = float(data[2][14].replace('-ADBitVolts ',''))*1E6
        
        tsf.SampleFrequency =float(data[2][12].replace('-SamplingFrequency ' ,''))
        print "SampleFrequency = ", tsf.SampleFrequency,
        
        tsf.n_vd_samples = len(data[0])
        print "... #samples: ", tsf.n_vd_samples
        if tsf.n_vd_samples<min_samples: min_samples = tsf.n_vd_samples

        tsf.ec_traces.append(data[0]) 
        plt.plot(data[0][:1000000])
        plt.show()
    
    #Trunkate extra voltage values (some chs in Neuralynx recs have more/less values than others)
    tsf.n_vd_samples = min_samples 
    for k in range(len(tsf.ec_traces)):
        tsf.ec_traces[k]=tsf.ec_traces[k][:min_samples]
        
    #tsf.ec_traces = np.array(tsf.ec_traces, dtype=np.int16)
    
    
    #******************SAVE HIGH PASS RECORD******************
    if self.parent.make_hp.text()=='True':
        print '\n...saving alltrack _hp_fromlfp.tsf...'
        file_name = ncs_files[0][:-4]+"_alltrack_hp_fromlfp.tsf"
        save_tsf_single(tsf, file_name)
    else:
        print "...skipping hp save..."
    
    #*************SAVE LOW PASS RECORD @ 1KHZ***************
    print '... processing low pass record...'
    #Wavelet filter record first
    #tsf.ec_traces = wavelet(tsf.ec_traces)
    
    temp_traces = []
    lowcut = 0.1; highcut=110; fs=1000
    for k in range(tsf.n_electrodes):
        temp = np.array(butter_bandpass_filter(tsf.ec_traces[k][::int(tsf.SampleFrequency/1000)], lowcut, highcut, fs, order = 2), dtype=np.float32)
        #temp = np.array(tsf.ec_traces[k][::int(tsf.SampleFrequency/1000.)])
        #Apply 60Hz Notch filter
        temp_traces.append(Notch_Filter(temp))

    #tsf.ec_traces = np.int16(temp_traces)
    tsf.ec_traces = np.array(temp_traces)*tsf.vscale_HP
    tsf.ec_traces = np.int16(tsf.ec_traces)
    
    #saving low pass record
    tsf.SampleFrequency = 1000
    tsf.vscale_HP = 1.0
    print ''; print "...saving alltrack _lp.tsf..."
    file_name = ncs_files[0][:-4]+"_alltrack_lp.tsf"
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    save_tsf_single(tsf, file_name)


def ntt_to_tsf(self, ntt_files):
    """ Neuralynx data conversion """
    from ncs import loadNtt

    print "... .ncs to .tsf ..."
    
    for k in range(len(ntt_files)):  print ntt_files[k]
    
    tsf = Object_empty()
    max_samples = 0
    spike_times = []
    spike_data = []
    for ctr, file_name in enumerate(ntt_files):
        data = loadNtt(file_name)
        spike_times.append(np.int32(data[0]/1.E6*data[2]))

        temp_data = data[1].swapaxes(0,2).swapaxes(1,2)
        spike_data.append(temp_data)    #Transpose data so each row is a channel
        if np.max(spike_times[ctr])>max_samples: max_samples = np.max(spike_times[ctr])
        
    tsf.SampleFrequency = data[2]

    tsf.ec_traces=np.zeros((len(ntt_files)*4, max_samples+32), dtype=np.int16)
 
    print tsf.ec_traces.shape
    for tetr in range(len(spike_data)):
        for ch in range(4):
            print "... processing tetrode: ", tetr, " ch: ", ch
            for s in range(len(spike_data[tetr][ch])):
                #print "tetr: ", tetr, " spk: ", s, " ch: ", ch, "spk time: ", spike_times[tetr][s]
                tsf.ec_traces[tetr*4+ch][spike_times[tetr][s]:spike_times[tetr][s]+32] = spike_data[tetr][ch][s]*10. 
    

    #******************SAVE HIGH PASS RECORD******************
    tsf.header = 'Test spike file '
    tsf.vscale_HP = 0.1     #USE ADperBit TO INCREASE RESOLUTIN
    tsf.iformat = 1002
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    tsf.n_electrodes = len(ntt_files)*4
    tsf.n_cell_spikes = 0

    tsf.Siteloc = np.zeros((tsf.n_electrodes*2), dtype=np.int16) #Read as 1D array
    for i in range (tsf.n_electrodes):
        tsf.Siteloc[i*2]=0
        tsf.Siteloc[i*2+1]=i*50 #GUESSING each tetrode is 50um apart
            
    #Wavelet filter record first
    #tsf.ec_traces = wavelet(tsf.ec_traces)

    print ''; print "...saving alltrack _hp.tsf..."
    file_name = ntt_files[0][:-4]+"_alltrack_hp.tsf"
    save_tsf_single(tsf, file_name)

    

def concatenate_tsf(self):
    """ Function doc """
    
    print "...concatenate multiple .tsf..."
    
    total_n_vd_samples = 0
    for k in range(len(self.tsf_files)):  
        tsf = Tsf_file(self.tsf_files[k])
        total_n_vd_samples += tsf.n_vd_samples
    
    print "...original layout: ",     
    print tsf.Siteloc
    print tsf.Readloc
    
    
    print "...total length of recs: ", total_n_vd_samples
    tsf.n_vd_samples = total_n_vd_samples #Set total tsf file # samples 
    tsf.n_cell_spikes = 0
    
    #Initialize ec_traces total
    tsf.ec_traces = np.zeros((tsf.n_electrodes, total_n_vd_samples), dtype=np.int16)
    print tsf.ec_traces.shape
    
    #Load each tsf data file 
    tsf_index = 0
    for ctr, file_name in enumerate(self.tsf_files):
        print "... loading: ", file_name
        temp_tsf = Tsf_file(file_name)
        temp_tsf.read_ec_traces()
        
        for ch in range(len(temp_tsf.ec_traces)):
            tsf.ec_traces[ch,tsf_index:tsf_index+len(temp_tsf.ec_traces[ch])] = temp_tsf.ec_traces[ch]
        
        tsf_index+=len(temp_tsf.ec_traces[ch])
        
    print ''; print "...saving alltrack .tsf..."

    file_name = self.tsf_files[0][:-4]+"_alltrack.tsf"
    save_tsf_single(tsf, file_name)
        

def concatenate_lfp_zip(self):
    """ Function doc """

    print "...concatenate lfp.zip files..."
    
    #This is only for Nick, Martin cat data; for intan data, can make lfp files from raw .tsf files directly
    for ctr, dir_name in enumerate(self.parent.animal.tsf_files):
        
        file_name = dir_name+dir_name[dir_name.rfind('/'):]+'.lfp.zip'
        print file_name
        if ctr==0:
            self.parent.animal.load_lfp_all(file_name)
            tsf = self.parent.animal.tsf
            tsf.iformat = 1002
            tsf.header = 'Test spike file '
        else:
            self.parent.animal.load_lfp_all(file_name)
            tsf_temp = self.parent.animal.tsf
            
            temp_ec_traces=[]
            for ch in range(tsf.n_electrodes):
                temp_ec_traces.append(np.append(tsf.ec_traces[ch],tsf_temp.ec_traces[ch]))
            
            tsf.ec_traces=np.int16(temp_ec_traces)
            tsf.n_vd_samples += tsf_temp.n_vd_samples
    
    tsf.n_cell_spikes = 0
    tsf.Siteloc=np.ravel(tsf.Siteloc[tsf.chans])        #Channel locations saved as flattened x,y coords
    
    print tsf.iformat
    print tsf.SampleFrequency
    print tsf.n_electrodes
    print tsf.n_vd_samples, tsf.n_vd_samples/float(tsf.SampleFrequency)
    print tsf.vscale_HP
    print tsf.chans
    print tsf.Siteloc
    
    print tsf.ec_traces

    tsf.ec_traces = np.int16(tsf.ec_traces*tsf.vscale_HP)
    tsf.vscale_HP = 1.0
    
    print ''; print "...saving alltrack .tsf..."
    file_name = self.parent.animal.tsf_files[0] + "_alltrack_lp.tsf"
    save_tsf_single(tsf, file_name)    

def compress_lfp(self):
    
    print "...making compressed lfp files ..."
    
    print self.parent.animal.tsf_file
    compression_factor = int(self.parent.compress_factor.text())
    print "...compressed factor: ", compression_factor
    
    tsf = Tsf_file(self.parent.animal.tsf_file)
    tsf.read_ec_traces()
    
    #Leave ADC convertion intact it possible
    #tsf.ec_traces= np.int16(tsf.ec_traces*tsf.vscale_HP)
    #tsf.vscale_HP = 1.0 
    
    if tsf.SampleFrequency != 1000:
        print "...lfp frequency not = 1000Hz... exiting ..."
        return
        
    #*********** SAVE COMPRESSED LOW PASS .TSF FILE *********
    #Save compression file name
    file_out = self.parent.animal.tsf_file[:-4]+'_'+str(compression_factor)+'compressed.tsf'
    print "Saving LFP : ", file_out
    
    #DON"T USE SUBSAMPLING - CAUSES PROBLEMS LATER
    #traces_out = []
    #for k in range(len(tsf.ec_traces)):
    #    #traces_out.append(tsf.ec_traces[k][::int(compression_factor/25)])      
    #    traces_out.append(tsf.ec_traces[k])
    #tsf.ec_traces = np.array(traces_out)

    tsf.SampleFrequency = compression_factor*1000     
    tsf.subsample = 1.0
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    tsf.save_tsf(file_out)
    
    

               
def load_lfpzip(file_name):     #Nick/Martin data has different LFP structure to their data.
    
    tsf = Object_empty()
    tsf.file_name = file_name
    
    #Set defaults
    tsf.header = 'Test spike file '
    tsf.iformat = 1002
    tsf.n_cell_spikes = 0
    
    #Load tsf.data
    data_in = np.load(file_name)
    
    tsf.SampleFrequency = data_in['tres']
    tsf.chans = data_in['chans']
    tsf.n_electrodes = len(tsf.chans)
    tsf.Siteloc = data_in['chanpos'][tsf.chans]
    tsf.vscale_HP = data_in['uVperAD']
    tsf.ec_traces = data_in['data']

    #Apply notch filter
    offset = np.zeros(int(data_in['t0']*1E-3), dtype=np.int16)      #Convert microsecond offset to miliseconds; Martin didn't save offset period as zeros... so add it back in
    temp_traces = []
    for k in range(tsf.n_electrodes):
        tsf.ec_traces[k] = Notch_Filter(tsf.ec_traces[k])
        temp_traces.append(np.append(offset, tsf.ec_traces[k]))
    
    tsf.ec_traces = np.int16(temp_traces)
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    
    return tsf
        
def rhd_to_tsf(filenames):
    '''Read .rhd files, convert to correct electrode mapping and save to .tsf file
    NB: There are 2 possible mapping depending on the insertion of the AD converter 
    TODO: implement a wavelet high pass filter directly to avoid SpikeSorter Butterworth filter artifacts
    '''
    
    print "...reading amp data..."

    probe = Probe()

    for file_name in filenames:
        print file_name
        #Delete previous large arrays; Initialize arrays; IS THIS REDUNDANT?
        ec_traces = 0.; ec_traces_hp = 0.; data=0.
        
        #file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_hp.tsf'
        file_out = file_name[:-4].replace('rhd_files','tsf_files')
        if os.path.exists(file_out)==True: continue

        print "Processing: \n", file_name

        data = read_data(file_name)
        ec_traces = data['amplifier_data'] #*10       #Multiply by 10 to increase resolution for int16 conversion
        ec_traces*=10.

        SampleFrequency = int(data['frequency_parameters']['board_adc_sample_rate']); print "SampleFrequency: ", SampleFrequency

        header = 'Test spike file '
        iformat = 1002
        n_vd_samples = len(ec_traces[0]); print "Number of samples: ", n_vd_samples
        vscale_HP = 0.1                             #voltage scale factor
        n_cell_spikes = 0

        print "Converting data to int16..."
        ec_traces = np.array(ec_traces, dtype=np.int16)

        #print "...plotting..."
        #plt.plot(ec_traces[30])
        #plt.show()


        #SAVE RAW DATA - ******NB:  SHOULD CLEAN THIS UP: the write function should be shared by all, just data is changing so no need to repeat;
        if True:
            print "Writing raw data ..."
            #print "CHANGE THIS TO WORK THROUGH FUNCTION WITHOUT REPEATING"
            fout = open(file_out+'_raw.tsf', 'wb')
            fout.write(header)
            fout.write(struct.pack('i', 1002))
            fout.write(struct.pack('i', SampleFrequency))
            fout.write(struct.pack('i', probe.n_electrodes))
            fout.write(struct.pack('i', n_vd_samples))
            fout.write(struct.pack('f', vscale_HP))
            
            for i in range (probe.n_electrodes):
                fout.write(struct.pack('h', probe.Siteloc[i][0]))
                fout.write(struct.pack('h', probe.Siteloc[i][1]))
                fout.write(struct.pack('i', i+1))

            for i in range(probe.n_electrodes):
                print "...writing ch: ", i
                ec_traces[probe.layout[i]].tofile(fout)  #Frontside

            fout.write(struct.pack('i', n_cell_spikes))
            fout.close()
            
        #SAVE HIGH PASS WAVELET FILTERED DATA
        if True:
            print "Writing hp data ..."
            fout = open(file_out+'_hp.tsf', 'wb')
            fout.write(header)
            fout.write(struct.pack('i', 1002))
            fout.write(struct.pack('i', SampleFrequency))
            fout.write(struct.pack('i', probe.n_electrodes))
            fout.write(struct.pack('i', n_vd_samples))
            fout.write(struct.pack('f', vscale_HP))
            
            for i in range (probe.n_electrodes):
                fout.write(struct.pack('h', probe.Siteloc[i][0]))
                fout.write(struct.pack('h', probe.Siteloc[i][1]))
                fout.write(struct.pack('i', i+1))

            print "Wavelet filtering..."
            ec_traces_hp = wavelet(ec_traces, wname="db4", maxlevel=6)
            print ec_traces_hp.shape

            for i in range(probe.n_electrodes):
                print "...writing ch: ", i
                ec_traces_hp[probe.layout[i]].tofile(fout)  #Frontside

            fout.write(struct.pack('i', n_cell_spikes))
            fout.close()


def rhd_digital_save(file_name):
    '''Read .rhd files, and save digital channels.
    NB: there can be 2, 4 or 6 digital channels inside Intan file
    chs 1 and 2 are laser meta data and laser pulse times (these are off for other experiments)
    chs 3 and 4 are camera pulse times and on/off times from clampx computer (these are chs 1 and 2 usually as laser chs are off)
    chs 5 and 6 are vis stim pulse times and meta data (these are usually 3 and 4 as not recorded w. laser on)
    '''
    
    print "...reading digital amp data..."

    
    #if os.path.exists(camera_onoff_filename+'.npy')==True: continue

    data = read_data(file_name)

    print "...# digital channels: ", len(data['board_dig_in_data'])

    SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
    print "SampleFrequency: ", SampleFrequency
    
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure

    plt.close()
    fig = Figure(figsize=(5,5), dpi=100)
    ax1 = fig.add_subplot(111)

    for ch in range(len(data['board_dig_in_data'])):
        ax1.plot(data['board_dig_in_data'][ch][:100000])
        plt.show()

        response = raw_input("Please enter filename extension: ")
        np.save(file_name[:-4].replace('rhd_files','camera_files')+'_'+response, data['board_dig_in_data'][ch])
        print file_name[:-4].replace('rhd_files','camera_files')+'_'+response
            


def parse_camera_pulses(self):
    
    upstates = np.where(self.camera_pulses==1)[0] #Find all index values where camera_pulse = 1
    
    #plt.plot(self.camera_pulses)
    #plt.show()
    
    triggers = []; triggers.append(upstates[0])
    triggers_length = []; marker = upstates[0]
    for k in range(len(upstates)-1):
        if (upstates[k+1] - upstates[k])>1:
            triggers.append(upstates[k+1])
            triggers_length.append(upstates[k]-marker); marker = upstates[k+1]
    
    self.triggers = np.array(triggers)
    self.triggers_length = np.array(triggers_length)


def compute_lfp_triggered_template(self):
    
    print "...excluded trials: ", self.excluded_trials.text()
    excluded_trials = self.excluded_trials.text()
    
    
    tsf = Tsf_file(self.selected_recording)
    tsf.read_ec_traces()
    
    print self.triggers
    print tsf.n_vd_samples

    font_size = 30
    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./1.) #int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())
    #compression = 50.   #Need this to convert from compressed sample points to realtime
    
    print self.selected_sort

    #Remove bottom power..
    if False:
    #if self.low_cutoff.text()!='0.0':
        for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
            print "...filtering ch: ", k
            self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
    #load single units
    #Sort = Ptcs(self.selected_sort)

    ax = plt.subplot(1,1,1)
    t = np.arange(-n_samples,n_samples+1,1)
    
    """ Spikes are saved in # of sample points so no need to scale them up from compressed .ptcs sort file to uncompressed lfp file.
    """
 
    traces = []
    for k in range(0, tsf.n_electrodes, electrode_rarifier):
        traces.append([])
    
    gs = gridspec.GridSpec(2,12)
    for ctr,trigger in enumerate(self.triggers):
        if str(ctr) in excluded_trials: continue
        print int(ctr/6), ctr%6
        ax = plt.subplot(gs[int(ctr/6),ctr%6])
        #ax = plt.subplot(2,6,ctr+1)
        print "...plotting event: ", ctr
        print "...trigger: ", trigger, " time: ", float(trigger)/tsf.SampleFrequency
        #for spike in Sort.units[int(self.selected_unit.text())]:
        for k in range(0, tsf.n_electrodes, electrode_rarifier):
            trace_out = tsf.ec_traces[k][int(trigger-n_samples):int(trigger+n_samples+1)]*tsf.vscale_HP
            traces[k].append(trace_out)
            
            plt.plot(t, trace_out-voltage_scaling*tsf.Siteloc[k*2+1], color='black', linewidth=1, alpha=.5)
            plt.yticks([])
            
            
    #Plot average
    ax = plt.subplot(gs[:,6:])
    for k in range(tsf.n_electrodes):
        trace_ave = np.average(np.array(traces[k]), axis=0)
        trace_std = np.std(np.array(traces[k]), axis=0)
       
        #offset = -voltage_scaling*tsf.Siteloc[k*2+1]
        plt.plot(t, trace_ave-voltage_scaling*tsf.Siteloc[k*2+1], color='black', linewidth=3)
        ax.fill_between(t, trace_ave+trace_std-voltage_scaling*tsf.Siteloc[k*2+1], trace_ave-trace_std-voltage_scaling*tsf.Siteloc[k*2+1], color='blue', alpha=0.25)

    #plt.plot([t[-1]+10,t[-1]+10], [-250, 0 ], color='black', linewidth=3)

    ##Set ylabel
    #old_ylabel = -voltage_scaling*np.linspace(0, np.max(tsf.Siteloc), 5)
    #new_ylabel = np.int16(np.linspace(0, np.max(tsf.Siteloc), 5))
    #plt.locator_params(axis='y',nbins=5)
    #plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    #plt.ylabel("Depth (um)", fontsize=font_size)

        ##Set xlabel
        #old_xlabel = np.linspace(t[0],t[-1],3)
        #new_xlabel = np.linspace(float(t[0])/tsf.SampleFrequency*1E3, float(t[-1])/tsf.SampleFrequency*1E3, 3)
        #plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)

        #plt.xlabel("Time (ms)", fontsize = font_size)
        #plt.tick_params(axis='both', which='both', labelsize=font_size)
        #plt.locator_params(axis='x',nbins=10)

        #plt.ylim(old_ylabel[-1],old_ylabel[0])
        #plt.xlim(old_xlabel[0], old_xlabel[-1])
    
    plt.show()

    

def tsf_to_lfp(filenames):
    '''Read .tsf files - subsample to 1Khz, save as *_lp.tsf
    '''
    
    print "...making low-pass tsf files (1Khz sample rates)..."

    for file_name in filenames:
        
        file_out = file_name[:-4]+'_lp.tsf'
        if os.path.exists(file_out)==True: continue

        print "Processing: \n", file_name

        tsf = Tsf_file(file_name)
        tsf.read_ec_traces()
        print tsf.Siteloc.shape
        
        n_vd_samples = len(tsf.ec_traces[0]); print "Number of samples: ", n_vd_samples
        
        print "...converting raw to .lfp (1Khz) sample rate tsf files ..."
        temp_traces = []
        lowcut = 0.1; highcut=110; fs=1000
        for k in range(tsf.n_electrodes):
            #Butter band pass and subsample to 1Khz simultaneously
            temp = np.array(butter_bandpass_filter(tsf.ec_traces[k][::int(tsf.SampleFrequency/1000)], lowcut, highcut, fs, order = 2), dtype=np.int16)

            #Apply 60Hz Notch filter
            temp_traces.append(Notch_Filter(temp))
        
        tsf.ec_traces = np.int16(temp_traces)
        tsf.n_vd_samples = len(tsf.ec_traces[0])
        tsf.SampleFrequency = fs
        
        #Save data to .tsf file
        tsf.save_tsf(file_name[:-4]+'_lp.tsf')


def lfp_to_lptsf(lfpzip_file):
    '''Read .tsf files - subsample to 1Khz, save as *_lp.tsf
    '''
    
    tsf = load_lfpzip(lfpzip_file)

    print tsf.Siteloc
    temp_chs = []
    for k in range(len(tsf.Siteloc)):
        temp_chs.extend(tsf.Siteloc[k])
    tsf.Siteloc = np.array(temp_chs)
    print tsf.Siteloc

    #Save data to .tsf file
    save_tsf_single(tsf, lfpzip_file.replace('.lfp.zip','_lp.tsf'))
    


def Notch_Filter(data, fs=1000, band=.5, freq=60., ripple=10, order=2, filter_type='ellip'):
    """Using iirfilter instead of Martin's LFP work
    """
    from scipy.signal import iirfilter, lfilter
    #fs   = 1/time
    nyq  = fs/2.0
    low  = freq - band/2.0
    high = freq + band/2.0
    low  = low/nyq
    high = high/nyq
    
    b, a = iirfilter(order, [low, high], rp=ripple, rs=50, btype='bandstop',
                     analog=False, ftype=filter_type)
    filtered_data = lfilter(b, a, data)
    return filtered_data

def view_templates(self):
    print "..."

    font_size = 30
    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())
    #compression = 50.   #Need this to convert from compressed sample points to realtime
    
    print self.selected_sort

    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
            print "...filtering ch: ", k
            self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
    #load single units
    Sort = Ptcs(self.selected_sort)

    ax = plt.subplot(1,1,1)
    t = np.arange(-n_samples,n_samples+1,1)
    
    """ Spikes are saved in # of sample points so no need to scale them up from compressed .ptcs sort file to uncompressed lfp file.
    """
    for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
        print "...plotting ch: ", k
        traces = []
        for spike in Sort.units[int(self.selected_unit.text())]:
            trace_out = self.tsf.ec_traces[k][int(spike-n_samples):int(spike+n_samples+1)]*self.tsf.vscale_HP
            traces.append(trace_out)
            #plt.plot(t, trace_out-voltage_scaling*self.tsf.Siteloc[k*2+1], color='black', linewidth=1, alpha=.1)
        
        trace_ave = np.average(traces, axis=0)
        trace_std = np.std(traces, axis=0)
       
        offset = -voltage_scaling*self.tsf.Siteloc[k*2+1]
        plt.plot(t, trace_ave+offset, color='black', linewidth=3)
        ax.fill_between(t, trace_ave+trace_std+offset, trace_ave-trace_std+offset, color=self.selected_colour, alpha=0.4)

    plt.plot([t[-1]+10,t[-1]+10], [-250, 0 ], color='black', linewidth=3)

    #Set ylabel
    old_ylabel = -voltage_scaling*np.linspace(0, np.max(self.tsf.Siteloc), 5)
    new_ylabel = np.int16(np.linspace(0, np.max(self.tsf.Siteloc), 5))
    plt.locator_params(axis='y',nbins=5)
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    plt.ylabel("Depth (um)", fontsize=font_size)

    #Set xlabel
    old_xlabel = np.linspace(t[0],t[-1],3)
    new_xlabel = np.linspace(float(t[0])/self.tsf.SampleFrequency*1E3, float(t[-1])/self.tsf.SampleFrequency*1E3, 3)
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)

    plt.xlabel("Time (ms)", fontsize = font_size)
    plt.tick_params(axis='both', which='both', labelsize=font_size)
    plt.locator_params(axis='x',nbins=10)

    plt.ylim(old_ylabel[-1],old_ylabel[0])
    plt.xlim(old_xlabel[0], old_xlabel[-1]*5)
    
    plt.show()



def view_traces(self):
    """ Display raw traces
    """
    
    font_size = 30
    voltage_scaling = float(self.voltage_scale.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    
    t0 = int(float(self.start_time.text())*self.tsf.SampleFrequency)
    t1 = int(float(self.end_time.text())*self.tsf.SampleFrequency)
    
    t = np.linspace(float(self.start_time.text()), float(self.end_time.text()), t1-t0)/60.
    
    ax = plt.subplot(1,1,1)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)

    print self.tsf.Siteloc
    print self.tsf.SampleFrequency
    print float(self.low_cutoff.text())

    for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
        trace_out = self.tsf.ec_traces[k][t0:t1]
        if float(self.low_cutoff.text())!=0.0: trace_out = butter_bandpass_filter(trace_out, float(self.low_cutoff.text()), 110., fs=1000, order = 2)
        
        plt.plot(t, trace_out-voltage_scaling*self.tsf.Siteloc[k*2+1], color='black', linewidth=1)

    plt.plot(t, -1500*self.camera_pulses[t0:t1], color='blue')
    #Set labels
    depth_offset = float(self.probe_penentration.text()) #Percent of probe inserted into brain
    old_ylabel = -voltage_scaling*np.linspace(np.max(self.tsf.Siteloc) - depth_offset*np.max(self.tsf.Siteloc),np.max(self.tsf.Siteloc), 5)
    new_ylabel = np.int16(np.linspace(0, depth_offset*np.max(self.tsf.Siteloc), 5))
    plt.locator_params(axis='y',nbins=5)
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    plt.ylabel("Depth (um)", fontsize=font_size)

    #Set xlabel
    plt.xlabel("Time (min)", fontsize = font_size)
    plt.tick_params(axis='both', which='both', labelsize=font_size)
    plt.locator_params(axis='x',nbins=10)

    plt.ylim(old_ylabel[-1],old_ylabel[0])
    plt.xlim(t[0], t[-1])

    plt.show()

def view_all_csd(self): 

    import imp

    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())

    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        if os.path.exists(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')==False:
            for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
                print "...filtering ch: ", k
                self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
            np.save(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut', self.tsf.ec_traces)
        else: 
            self.tsf.ec_traces = np.load(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')
            
    
    #Raw traces
    tsf = self.tsf

    #load single units
    Sort = Ptcs(self.selected_sort)

    CSD = []
    #Load spiketimes as event_trigges 
    for p in range(len(Sort.units)):
        event_times = Sort.units[p]

        lfp_ave = np.zeros((len(tsf.ec_traces),2*n_samples),dtype=np.float32)
        for ch in range(tsf.n_electrodes):
            ctr=0
            for event in event_times:
                trace_temp = tsf.ec_traces[ch][int(event-n_samples):int(event+n_samples)]
                if len(trace_temp)==(n_samples*2):
                    lfp_ave[ch]+= trace_temp
                    ctr+=1
            lfp_ave[ch]=lfp_ave[ch]/ctr

        lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())]
            
        #********Compute CSD
        if tsf.n_electrodes >10: 
            print "...loading every other channel, NN A64 probe ..."
            probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())][::2]
            lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())][::2]
        else:
            probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())]
            lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())]

        z = probe_layout*1E-3 #Convert to mm size

        csdmod = imp.load_source('csd_est_funds','csd_est_funcs.py')
        
        sigma = 0.3  # extracellular conductivity (mS/mm)
        b = 1.0   # assumed radius of the column (mm)
        SNR = float(self.snr_value.text())    # average ratio of signal to noise on each channel

        [A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
        [W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
        [W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check

        CSD.append(np.dot(W,lfp_ave))   # units: mS/mm^3*mV = uA/mm^3
        
    #v_max = np.max(np.abs(CSD)); v_min = -v_max
    
    for p in range(len(Sort.units)):
        ax = plt.subplot(4,4,p+1)
        t = np.linspace(-n_samples, n_samples, 6)

        #ax.imshow(CSD[p], vmin = v_min, vmax=v_max, extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
        ax.imshow(CSD[p], extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
        plt.plot([0,0], [0,np.max(tsf.Siteloc)*1E-3], 'r--', linewidth=3, color='black', alpha=0.6)

        plt.ylim(np.max(tsf.Siteloc)*1E-3,0)
        plt.xlim(t[0], t[-1])

        plt.tick_params(axis='both', which='major', labelsize=15)
        #plt.ylabel("Depth along probe (mm)", fontsize=20)
        #plt.xlabel("Time (msec)",fontsize=20)

        plt.suptitle("Cluster: " + str(p), fontsize=25)
    plt.show()

def compute_csd_event_triggered(self):
    
    print "...excluded trials: ", self.excluded_trials.text()
    excluded_trials = self.excluded_trials.text()
        
    print self.triggers

    font_size = 30
    n_samples = int(self.n_sample_pts.text())
    
    print self.selected_sort

    #Remove bottom power..
    #if self.low_cutoff.text()!='0.0':
        #if os.path.exists(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')==False:
            #for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
                #print "...filtering ch: ", k
                #self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
            #np.save(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut', self.tsf.ec_traces)
        #else: 
            #self.tsf.ec_traces = np.load(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')
    
    #Raw traces
    tsf = Tsf_file(self.selected_recording)
    tsf.read_ec_traces()
    
    #load single units
    #Sort = Ptcs(self.selected_sort)

    #Load spiketimes as event_trigges 
    event_times = self.triggers

    lfp_ave = np.zeros((tsf.n_electrodes,2*n_samples),dtype=np.float32)
    for ch in range(tsf.n_electrodes):
        ctr=0
        for idx, event in enumerate(event_times):
            if str(idx) in excluded_trials: continue
            trace_temp = tsf.ec_traces[ch][int(event-n_samples):int(event+n_samples)]
            if len(trace_temp)==(n_samples*2):
                lfp_ave[ch]+= trace_temp
                ctr+=1
        lfp_ave[ch]=lfp_ave[ch]/ctr


    #********Compute CSD
    if tsf.n_electrodes >10: 
        print "...loading every other channel, NeuroNexus A64 probe ..."
        probe_layout = tsf.Siteloc[1::2][::2]
        lfp_ave=lfp_ave[::2]
    else:
        probe_layout = tsf.Siteloc[1::2]
        lfp_ave=lfp_ave

    print probe_layout
    
    z = probe_layout*1E-3 #Convert to mm size


    import imp
    csdmod = imp.load_source('csd_est_funds','csd_est_funcs.py')
    
    sigma = 0.3  # extracellular conductivity (mS/mm)
    b = 1.0   # assumed radius of the column (mm)
    SNR = float(self.snr_value.text())    # average ratio of signal to noise on each channel

    [A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
    [W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
    [W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check


    CSD=np.dot(W,lfp_ave)   # units: mS/mm^3*mV = uA/mm^3
    
    t = np.linspace(-n_samples, n_samples, 6)

    #ax.imshow(CSD[p], vmin = v_min, vmax=v_max, extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
    plt.imshow(CSD, extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
    plt.plot([0,0], [0,np.max(tsf.Siteloc)*1E-3], 'r--', linewidth=3, color='black', alpha=0.6)

    plt.ylim(np.max(tsf.Siteloc)*1E-3,0)
    plt.xlim(t[0], t[-1])

    plt.tick_params(axis='both', which='major', labelsize=15)
    #plt.ylabel("Depth along probe (mm)", fontsize=20)
    old_xlabel = np.float32(np.linspace(t[0], t[-1], 6), decimals=2)
    new_xlabel = old_xlabel/tsf.SampleFrequency*1E3
    plt.xticks(old_xlabel, new_xlabel, fontsize=18)
    plt.xlabel("Time (msec)",fontsize=20)
            
    plt.show()



def view_csd(self):
    
    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())
    
    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        if os.path.exists(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')==False:
            for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
                print "...filtering ch: ", k
                self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
            np.save(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut', self.tsf.ec_traces)
        else: 
            self.tsf.ec_traces = np.load(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')
    
    #Raw traces
    tsf = self.tsf
    
    #load single units
    Sort = Ptcs(self.selected_sort)

    #Load spiketimes as event_trigges 
    event_times = Sort.units[int(self.selected_unit.text())]

    lfp_ave = np.zeros((len(tsf.ec_traces),2*n_samples),dtype=np.float32)
    for ch in range(tsf.n_electrodes):
        ctr=0
        for event in event_times:
            trace_temp = tsf.ec_traces[ch][int(event-n_samples):int(event+n_samples)]
            if len(trace_temp)==(n_samples*2):
                lfp_ave[ch]+= trace_temp
                ctr+=1
        lfp_ave[ch]=lfp_ave[ch]/ctr

    #********Compute CSD
    print '.......testing....'
    print tsf.Siteloc
    print tsf.Siteloc[1::2]
    print tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())]
    print tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())][::2]
    
    if tsf.n_electrodes >10: 
        print "...loading every other channel, NN A64 probe ..."
        probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())][::2]
        lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())][::2]
    else:
        probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())]
        lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())]

    print probe_layout
    
    z = probe_layout*1E-3 #Convert to mm size
    

    import imp
    csdmod = imp.load_source('csd_est_funds','csd_est_funcs.py')
    
    sigma = 0.3  # extracellular conductivity (mS/mm)
    b = 1.0   # assumed radius of the column (mm)
    SNR = float(self.snr_value.text())    # average ratio of signal to noise on each channel

    [A,A0] = csdmod.forward_operator(z,b,sigma) # compute the forward operator: A -dimensional (mm^3/mS) and A0 -dimensionless operators 
    [W,R] = csdmod.inverse_tikhonov(A,SNR) # compute the inverse operator, units (mS/mm^3)
    [W0,R0] = csdmod.inverse_tikhonov(A0,SNR)  # the zeros are dimensionless operators which we do not use but display for sanity check


    CSD=np.dot(W,lfp_ave)   # units: mS/mm^3*mV = uA/mm^3
    
    t = np.linspace(-n_samples, n_samples, 6)

    #ax.imshow(CSD[p], vmin = v_min, vmax=v_max, extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
    plt.imshow(CSD, extent=[t[0],t[-1],z[-1],z[0]], cmap='jet', aspect='auto', interpolation='none'); 
    plt.plot([0,0], [0,np.max(tsf.Siteloc)*1E-3], 'r--', linewidth=3, color='black', alpha=0.6)

    plt.ylim(np.max(tsf.Siteloc)*1E-3,0)
    plt.xlim(t[0], t[-1]*5)

    plt.tick_params(axis='both', which='major', labelsize=15)
    #plt.ylabel("Depth along probe (mm)", fontsize=20)
    #plt.xlabel("Time (msec)",fontsize=20)

    plt.show()




def Specgram_syncindex(self):
    
    channel=int(self.parent.specgram_ch.text())
    
    #if self.parent.exp_type=="mouse":
        #self.parent.animal.load_channel(self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+
                        #'_lp.tsf', channel) #Loads single channel as animal.tsf
                        
    #elif self.parent.exp_type=="cat":
        #temp_name = self.parent.animal.recName.replace('.ptcs','.lfp.zip').replace('.tsf','.lfp.zip')
        #self.parent.animal.load_lfp_channel(temp_name, channel) #Loads single channel as animal.tsf
    #elif self.parent.exp_type=='rat':
        #temp_name = self.parent.animal.recName.replace('_hp.tsf','_lp.tsf')
        #self.parent.animal.load_channel(temp_name, channel) #Loads single channel as animal.tsf
        
        ##print "exp_type unknown"
        ##return
    
    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    if '.tsf' in self.parent.recName:
        tsf = Tsf_file(self.parent.recName)
        tsf.read_ec_traces()
        for k in range(len(tsf.Siteloc)):
            print k, tsf.Siteloc[k]
            
    elif '.lfp.zip' in self.parent.recName:
        tsf = load_lfpzip(self.parent.recName)
        for k in range(len(tsf.Siteloc)):
            print k, tsf.Siteloc[k]

    
    samp_freq = tsf.SampleFrequency
    print "rec length: ", len(tsf.ec_traces[channel])/float(tsf.SampleFrequency), " sec."

    ax = plt.subplot(1,1,1)
    font_size = 30
    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    #Compute Specgram
    print "computing specgram..."
    data_in = tsf.ec_traces[channel]
    f0 = 0.1; f1 = 100
    p0 = float(self.parent.specgram_db_clip.text())
    P, extent = Compute_specgram_signal(data_in, samp_freq, f0, f1, p0)
    plt.imshow(P, extent=extent, aspect='auto')

    #Compute sync index
    si_limit=0.7
    si, t, sync_periods = synchrony_index(data_in, samp_freq, si_limit)
    
   # print "t: ", t
    
    sync_0 = -30
    sync_1 = -10
    sync_scale = 20
    
    plt.plot(t, si*sync_scale+sync_0, linewidth=6, color='black')
    plt.plot([0,max(t)],[-sync_scale*.3+sync_1,-sync_scale*.3+sync_1], 'r--', color='black', linewidth = 3, alpha=0.8)
    plt.plot([0,max(t)],[sync_1,sync_1], color='black', linewidth = 2, alpha=0.8)
    plt.plot([0,max(t)],[sync_0,sync_0], color='black', linewidth = 2, alpha=0.8)
    
    xx = np.linspace(0,max(t),5)
    x_label = np.round(np.linspace(0, max(t)/60.,5))
    plt.xticks(xx, x_label, fontsize=20)
       
    ax.set_xlim((0,P.shape[1]/2))    
    ax.set_ylim((sync_0-1,f1))    
    
    old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    
    ax.tick_params(axis='both', which='both', labelsize=font_size)

    plt.ylabel("Synchrony Index             Specgram Frequency (Hz)      ", fontsize=font_size-5)           
    plt.xlabel("Time (mins)", fontsize = font_size)
    plt.title(self.parent.recName, fontsize=font_size-10)
    plt.show()

def Compute_specgram_signal(data, SampleFrequency, f0=0.1, f1=110, p0=-40):

    t0=None
    t1=None
    #f0=0.1
    #f1=110
    #p0         #clipping bottom of specgram
    p1=None
    chanis=-1
    width=2
    tres=.5
    cm=None
    colorbar=False
    title=True
    figsize=(20, 6.5)
    
    P0, P1 = None, None
    chanis = -1

    sampfreq=SampleFrequency #1KHZ LFP SAMPLE RATE for Nick's data; Otherwise full sample rates;

    #NFFT = intround(width * sampfreq)
    #NOVERLAP = intround(NFFT - tres * SAMPFREQ)

    length = len(data)

    ts = np.arange(0,len(data),1.0)/sampfreq

    if t0 == None:
        t0, t1 = ts[0], ts[-1] # full duration
    #if t1 == None:
    #    t1 = t0 + 10 # 10 sec window
    if width == None:
        width = uns['LFPWIDTH'] # sec
    if tres == None:
        tres = uns['LFPTRES'] # sec
    assert tres <= width

    NFFT = intround(width * sampfreq)
    noverlap = intround(NFFT - tres * sampfreq)

    t0i, t1i = ts.searchsorted((t0, t1))

    #data = filter.notch(data)[0] # remove 60 Hz mains noise

    print "Computing regular fft specgram"
    P, freqs, t = mpl.mlab.specgram(data/1e3, NFFT=NFFT, Fs=sampfreq, noverlap=noverlap)
    
    # convert t to time from start of acquisition:
    t += t0
    # keep only freqs between f0 and f1:
    if f0 == None:
        f0 = freqs[0]
    if f1 == None:
        f1 = freqs[-1]
    lo, hi = freqs.searchsorted([f0, f1])
    P, freqs = P[lo:hi], freqs[lo:hi]
    #print P
    
    # check for and replace zero power values (ostensibly due to gaps in recording)
    # before attempting to convert to dB:
    zis = np.where(P == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        P[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float  #CAT: This can probably be unhacked using nanmax or masked arrays
        minnzval = P.min() # get minimum nonzero value
        P[zis] = minnzval # replace with min nonzero values
    P = 10. * np.log10(P) # convert power to dB wrt 1 mV^2?

    # for better visualization, clip power values to within (p0, p1) dB
    if p0 != None:
        P[P < p0] = p0
    if p1 != None:
        P[P > p1] = p1

    extent = ts[0], ts[-1], freqs[0], freqs[-1]

    return P[::-1], extent

def intround(n):
    """Round to the nearest integer, return an integer. Works on arrays.
    Saves on parentheses, nothing more"""
    if iterable(n): # it's a sequence, return as an int64 array
        return np.int64(np.round(n))
    else: # it's a scalar, return as normal Python int
        return int(round(n))

def iterable(x):
    """Check if the input is iterable, stolen from numpy.iterable()"""
    try:
        iter(x)
        return True
    except:
        return False


def compute_sta(self, ptcs_file):

    #self.animal.root_dir = self.animal_name.text()
    self.animal.name = self.animal_name.text()

    #PARAMETERS
    overwrite = False   #overwrite previous results
    n_procs=5
    window = 3          #number of secs frames pre- and post- spiking to compute STMTDs for
    min_spikes = 0  #Min number of spikes in unit 
    spike_mode = 'all'     #Trigger off all spkes
    compress_factor = 1.    #for normal .ptcs file sampled at 25Khz
    images_file_name = self.animal.ptcsName.replace('tsf_files','tif_files').replace('_hp.ptcs','').replace('_lp_compressed.ptcs','')+'.npy'
    rec_name = images_file_name.replace('.npy','').replace(self.animal.home_dir+self.animal.name+'/tif_files/','')
    file_extension = self.animal.ptcsName.replace('tsf_files','tif_files').replace(images_file_name[:-4],'').replace('.ptcs','')
    
    if 'compress' in self.animal.ptcsName:    
        if '07_11' in self.animal.ptcsName: compress_factor = 40.  #for july 11th timecompressed files
        else: compress_factor = 50.  #for all other files should be correct.

    print "compression: ", compress_factor
    #LOAD SORT
    print self.animal.ptcsName
    if (os.path.exists(self.animal.ptcsName)==False): return
    Sort = Ptcs(self.animal.ptcsName) #Auto load flag for Nick's data
    n_units = len(Sort.units); print "# units: ", n_units
    print np.float32(Sort.units[0])/(float(Sort.samplerate)/compress_factor)

    #LOAD CAMERA ON/OFF TIMES
    camera_file_name = self.animal.ptcsName.replace('tsf_files','camera_files').replace(file_extension+'.ptcs','')+'_camera_onoff.npy'
    print camera_file_name
    camera_onoff = np.load(camera_file_name)
    print len(camera_onoff)

    #FIND START/END OF IMAGING (from ephys data; find where camera val goes to '1')
    indexes = np.where(camera_onoff==1)[0]
    start_offset = float(indexes[0])/25000      #NB: THIS NEEDS TO BE LOADED FROM tsf.SampleFrequency or highpass Sort.samplerate
    end_offset = float(indexes[-1])/25000       #NB: compressed data samplerate is 50Khz or other so DO NOT FIX.
    print "start/end: ", start_offset, end_offset

    #LOAD RAW IMAGES; FILTERED OR NOT
    if False:
        print "...loading: ", images_file_name[:-4]+'_filtered.npy'
        images_raw = np.load(images_file_name[:-4]+'_filtered.npy')    #Zero out images before loading - might free up memory first;
    else:
        print "...loading: ", images_file_name
        images_raw = np.load(images_file_name)    #Zero out images before loading - might free up memory first;
        print images_raw.shape
        #print "...converting int16 to float32..."; images_raw = images_raw.astype(np.float32)#, copy=False)

    print images_raw.shape
    if len(images_raw.shape)<3:
        print "..reshaping..."
        images_raw = images_raw.reshape(-1, 128, 128)
        print images_raw.shape
       
    
    #COMPUTE AVE FRAME - FOR BASELINE
    print "... skipping baseline ..."
    #baseline = np.average(images_raw, axis=0)
    
    #make_vids([(images_raw[1000:2000]-baseline)/baseline], images_file_name[:-4]+'_filtered.npy')
    #quit()

    
    #COMPUTE IMG RATE
    if ('2016_07_20_vsd' in self.animal.ptcsName) and ('2nd' not in self.animal.ptcsName):        #SPECIAL CASE: Forgot to turn off imaging and it maxed out... was working on rig while recording...
        print "acquistion was corrupted...returning..."
        return

    else: #Compute img_rate and img_times normally:
        img_rate = float(len(images_raw))/float(end_offset-start_offset)
        img_times = np.linspace(0, float(end_offset-start_offset), len(images_raw))

    print "img_rate: ", img_rate
    print "# Interpolated img_times: ", len(img_times)
    
        
    for unit in range(n_units):
        print "\n\n****************"
        print "Processing Unit: ", unit
        #REMOVE TIME TO CAMERA ON FROM EPHYS TIMES; Even for triggered data, there is still ~1sec of buffered ephys data saved
        #spikes = np.array(Sort.units[unit])/Sort.samplerate - start_offset 
        spikes = np.float32(Sort.units[unit])/(float(Sort.samplerate)/compress_factor) - start_offset 
        #if unit != 1: continue

        channel = Sort.maxchan[unit]
    
        Compute_spike_triggered_average( unit, channel, spikes, images_raw, window, img_times, self.animal.home_dir+self.animal.name+'/', overwrite, img_rate, n_procs, spike_mode, self.animal.ptcsName, rec_name+file_extension)

    ##OPTIONAL: DON'T ERASE
    ##Load camera triggers to determine precision of interpolated frame times to saved pulse times.
    #camera_pulses = np.int32(np.load(animal.filenames[rec].replace('rhd_files','camera_files')+'_camera_pulses.npy'))
    #print len(camera_pulses) 
    #camera_pulses = camera_pulses[indexes[0]:indexes[-1]] #Look at pulses only between ON and OFF times...

    #if (os.path.exists(self.home_dir + camera_pulse_times)==False):
        #ctr = 0
        #ontimes = []
        #betweentimes = []
        #for k in range(len(camera_pulses)-1):
            #if ((camera_pulses[k] - camera_pulses[k+1])==-1): 
                #ontimes.append(k)
                #betweentimes.append(k-ctr)
                #ctr = k
        #np.save(self.home_dir + camera_pulse_times[:-4], ontimes)
        #np.save(self.home_dir + camera_pulse_times[:-4]+'_betweentimes',betweentimes)

    #ontimes = np.load(self.home_dir + camera_pulse_times)
    ##betweentimes = np.load(main_dir + camera_pulse_times[:-4]+'_betweentimes.npy')
    #print "Number of camera_pulses: ", len(ontimes)

    #max_len = len(ontimes)
    #if len(img_times)< max_len: max_len=len(img_times)
    #print "len(ontimes): ", len(ontimes), "  len(img_times): ", len(img_times)

def wavelet(data, wname="db4", maxlevel=6):
    """Perform wavelet multi-level decomposition and reconstruction (WMLDR) on data.
    See Wiltschko2008. Default to Daubechies(4) wavelet"""
    import pywt

    data = np.atleast_2d(data)  #TODO***************************PARALLELIZE THIS***********

    for i in range(len(data)):
        print "Wavelet filtering channel: ", i
        # decompose the signal:
        
        
        c = pywt.wavedec(data[i], wname, level=maxlevel)
        # destroy the appropriate approximation coefficients:
        c[0] = None
        # reconstruct the signal:
        data[i] = pywt.waverec(c, wname)
    
    return data
    
def Compute_spike_triggered_average( unit, channel, spikes, images_raw, window, img_times,  main_dir, overwrite, img_rate, n_procs, spike_mode, ptcs_file, rec_name):
    '''Computes average frame values from t=-window .. +window (usually 180 frames for 3sec window and 30Hz sampling) '''
    
    global n_pixels, images_temp, temp_window, temp_img_rate

    #Check to see if images already loaded and saved as .npy
    #stm_file_name = main_dir + 'stm_files/' + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
    stm_file_name = main_dir + 'stm_files/img_avg_' + rec_name+'_unit'+str(unit).zfill(3)+'_ch'+str(channel).zfill(3)+'_'+spike_mode+'_'+str(window)+'sec_window_'
    
    if (overwrite) or (len(glob.glob(stm_file_name+"*"))==0):
        
        n_pixels = len(images_raw[0])
        temp_window = window
        temp_img_rate = img_rate
        print "No. of processors: ", n_procs

        #Remove spikes outside of imaging frames
        print "Total no. spikes: ", len(spikes)
       
        temp0 = np.where(np.logical_and(spikes>=img_times[0]+2*window, spikes<=img_times[-1]-window))[0]    #Exclude spikes too close to beginning or end of recordings.
        spikes = spikes[temp0]

        temp3 = []
        if spike_mode =='all': #Use all spikes
            for spike in spikes:        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                temp3.append(np.where(np.logical_and(img_times>=spike-window, img_times<=spike+window))[0][0:2*int(window*img_rate)]) #Fixed this value or could be off +/-1 frame
        elif spike_mode =='lockout': 
            for s in range(1, len(spikes)-1,1):        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                if ((spikes[s]-window)>=spikes[s-1]) and ((spikes[s]+window)<=spikes[s+1]):
                    array_temp = np.where(np.logical_and(img_times>=spikes[s]-window, img_times<=spikes[s]+window))[0][0:2*int(window*img_rate)]
                    temp3.append(array_temp) #Fixed this value or could be off +/-1 frame
        elif spike_mode =='burst': 
            for s in range(1, len(spikes),1):        #temp3 contains all the frame indexes from img_times for each spike in raster; e.g. 180 frames for each spike automatically aligned
                if ((spikes[s]-window)>=spikes[s-1]):
                    array_temp = np.where(np.logical_and(img_times>=spikes[s]-window, img_times<=spikes[s]+window))[0][0:2*int(window*img_rate)]
                    temp3.append(array_temp) #Fixed this value or could be off +/-1 frame

        print "No. spks in img period: ", len(spikes), " rate: ", float(len(spikes))/(img_times[-1]-img_times[0]), " Hz, no. spks pass criteria: ", len(temp3)
            
        if len(temp3)==0: return 

        images_temp = images_raw

        #Compute all frames based on image index;
        print "... stm processing in parallel for window: ", window, " secs ..."
        images_triggered_temp=[]
        temp4 = []  #Contains only first spikes from bursts in cell

        #*********** USE SINGLE CORE ************
        if True: 
            temp3 = np.array(temp3)
            images_processed = np.zeros((len(temp3[0]),n_pixels,n_pixels), dtype=np.float32)
            
            for i in range(0,len(temp3),1):
                #if (i == 0) or (i==2500): 
                    print "...averaging spike: ", i, " / ", len(temp3), "  time: ", spikes[i]
                    
                    ##Konnerth method; Dongsheng paper method
                    temp_img = np.array(images_temp[temp3[i]]) #Remove avg of 1st half of images
                    baseline = np.average(images_temp[np.arange(temp3[i][0]-int(temp_img_rate)*3,temp3[i][0],1)], axis=0)  #Go back 3 sec prior to temp3[0] value...
                    temp_frame = (temp_img - baseline)/baseline
                    images_processed += temp_frame

                    #Allen, Jeff method
                    #temp_img = np.array(images_temp[temp3[i]]) #Remove avg of 1st half of images
                    #temp_frame = (temp_img - baseline)/baseline
                    #images_processed += temp_frame
                    
                    #Baseline already removed, just add images
                    #images_processed+= images_raw[temp3[i]]
                    
            images_processed = images_processed/max(len(temp3),1)

            #for k in range(0,300,3):
                #ax = plt.subplot(10,10,k/3+1)
                #plt.imshow(images_processed[300+k])
                ##plt.imshow(img1[400+k])
            #plt.title("Spike: "+str(i))
            #plt.show()

                
        #************ PARALLEL CODE  **************
        else:
            if len(temp3) < 30:
                pool = mp.Pool(1)
                temp4.append(temp3)            
            else:
                pool = mp.Pool(n_procs)
                chunks = int(len(temp3)/n_procs) #Break up the temp3 array into n_procs that are "chunk" long each
                for i in range(n_procs):
                    temp4.append(temp3[i*chunks:(i+1)*chunks])
                    
                #DISCARD RESIDUE: May wish to recapture this eventually.
                
            #print "Removing average of all pre spike frames - (time: -", window, "sec .. 0sec)"
            images_triggered_temp.extend(pool.map(Spike_averages_parallel_prespike_3sec, temp4))
            
            pool.close()
            print "... done "

            #Sum over all spikes
            print "Summing Number of chunks: ", len(images_triggered_temp)

            temp_images = np.zeros((int(window*img_rate)*2, n_pixels, n_pixels), dtype=np.float16)
            for i in range(len(images_triggered_temp)):
                temp_images += images_triggered_temp[i]
            
            #DIVIDE BY NUMBER OF CHUNKS; Note used to be divided by number of spikes; also residue is being thrown out...
            images_processed = temp_images/float(len(images_triggered_temp))

                
        #npy_file_name = file_dir + file_name + '/img_avg_' + file_name+ '_unit'+str(unit).zfill(2)+'_ch'+str(channel).zfill(2)+'_'+plot_string+'_'+str(window)+'sec_window_'+str(len(spikes)).zfill(5)+"_spikes"
        stm_file_name= stm_file_name+str(len(temp3)).zfill(5)+"_spikes"
        
        np.save(stm_file_name, images_processed)

    else: 
        print "Skipping processing of images ... loading from file"
        #images_processed = np.load(glob.glob(stm_file_name+"*")[0])

    #print "Shape mean: ", images_processed.shape



def sta_maps(self):

    animal = self.animal
    print animal.name
    print animal.recName
    Sort = Ptcs(animal.ptcsName)
    print "... n_units: ", len(Sort.units)
    print "... start time: ", self.start_time.text(), "  end time: ", self.end_time.text()
    print "... start cell: ", self.start_cell.text(), "  end cell: ", self.end_cell.text()
 
    block_save = int(self.block_save.text())

    select_units = np.arange(min(int(self.start_cell.text()),len(Sort.units)),min(int(self.end_cell.text()),len(Sort.units)),1)

    start_time = float(self.start_time.text())  #time before t=0 (secs)
    end_time = float(self.end_time.text())

    #**** LOAD GENERIC MASK
    generic_mask_file = animal.home_dir+animal.name + '/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        fname = glob.glob(animal.filenames[0].replace('rhd_files/','tif_files/')[:animal.filenames[0].find('rhd_files/')+10]+"*std*")[0]
        images_temp = np.load(fname)
        #Define_generic_mask(images_temp, animal.home_dir)
        Define_generic_mask_single_frame(images_temp, animal)
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        
    generic_mask_indexes=np.zeros((128,128))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True


    min_spikes = 0
    plot_img = []
    for ctr, k in enumerate(select_units):
        if len(Sort.units[k])<min_spikes: continue
        print "Loading saved .npy files: ", k
        temp_name = glob.glob(animal.ptcsName.replace('tsf_files/','stm_files/img_avg_').replace('.ptcs','')+'_unit'+str(k).zfill(3)+"*")[0]
        STM = np.load(temp_name)
        
        #Apply mask
        n_pixels=128
        temp_array = np.ma.array(np.zeros((len(STM),n_pixels,n_pixels),dtype=np.float32), mask=True)
        for i in range(len(STM)):
            temp_array[i] = np.ma.array(STM[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        
        ax = plt.subplot(len(select_units), 1, ctr+1)
        img_out = []
        #img_rate = float(len(STM))/6.
        img_rate = self.animal.img_rate
        print img_rate
        #for i in range(0,len(STM),block_save):
        for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
            img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
        img_out = np.ma.hstack((img_out))
        
        v_abs = max(np.ma.max(img_out),-np.ma.min(img_out))
        plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
        #plt.imshow(img_out)

        plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k])), rotation='horizontal',  labelpad=50, fontsize=10)
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])

        if ctr==(len(select_units)-1): 
            plt.xlabel("Time from spike (sec)", fontsize=25)
            old_xlabel = np.linspace(0,img_out.shape[1], 11)
            new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
            plt.xticks(old_xlabel, new_xlabel, fontsize=18)
        
        #plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k]))+"\nDepth: "+str(Sort.chanpos[Sort.maxchan[k]][1])+"um", fontsize=15)
    
    plt.suptitle(animal.ptcsName)
    plt.show()
    #quit()




#def sta_maps_lfp(animal, k):

    #start_time = -.2  #1 sec before t=0
    #end_time = +.2   #3 seconds after t=0

    ##**** LOAD GENERIC MASK
    #generic_mask_file = animal.home_dir + 'genericmask.txt'
    #if (os.path.exists(generic_mask_file)==True):
        #generic_coords = np.int32(np.loadtxt(generic_mask_file))
    #else:
        #images_temp = np.load(glob.glob(animal.filenames[s].replace('rhd_files/','stm_files/img_avg_') +'_unit'+str(0).zfill(3)+'_ch'+"*")[0])
        #Define_generic_mask(images_temp, animal.home_dir)
        #generic_coords = np.int32(np.loadtxt(generic_mask_file))
        
    #generic_mask_indexes=np.zeros((128,128))
    #for i in range(len(generic_coords)):
        #generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    ##**** MAKE STATIC FIGS
    ##Sort = Ptcs(animal.filenames[s].replace('rhd_files','tsf_files')+'_hp.ptcs')
    #Sort = Ptcs('/media/cat/12TB/in_vivo/tim/cat/2016_07_11_vsd/tsf_files/track2_150Hz_iso1.0_spontaneous_lp_compressed.ptcs')
    ##select_units = [0]
    #select_units = np.arange(0,len(Sort.units),1)

    #min_spikes = 50
    #plot_img = []
    #for ctr, k in enumerate(select_units):
        #if len(Sort.units[k])<min_spikes: continue
        #print "Loading saved .npy files: ", k
        #channel = Sort.maxchan[k]
        ##print glob.glob(animal.filenames[k].replace('rhd_files/','stm_files/img_avg_') +'_unit'+str(0).zfill(3)+"*")
        #temp_name = glob.glob(animal.filenames[s].replace('rhd_files/','stm_files/img_avg_') +'_unit'+str(k).zfill(3)+"*")[0]
        #STM = np.load(temp_name)
        
        ##Apply mask
        #n_pixels=128
        #temp_array = np.ma.array(np.zeros((len(STM),n_pixels,n_pixels),dtype=np.float32), mask=True)
        #for i in range(len(STM)):
            #temp_array[i] = np.ma.array(STM[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True)
        
        #ax = plt.subplot(len(select_units), 1, ctr+1)
        #img_out = []
        #block_save = 1
        #img_rate = float(len(STM))/6.
        ##for i in range(0,len(STM),block_save):
        #for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
            #img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
        #img_out = np.ma.hstack((img_out))
        
        #v_abs = max(np.ma.max(img_out),-np.ma.min(img_out))
        #plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
        ##plt.imshow(img_out)

        #plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k])), rotation='horizontal',  labelpad=50, fontsize=10)
        #ax.yaxis.set_ticks([])
        #ax.xaxis.set_ticks([])

        #if ctr==(len(select_units)-1): 
            #plt.xlabel("Time from spike (sec)", fontsize=25)
            #old_xlabel = np.linspace(0,img_out.shape[1], 11)
            #new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
            #plt.xticks(old_xlabel, new_xlabel, fontsize=18)
        
        ##plt.ylabel("U: " + str(k) + ", #spk: " + str(len(Sort.units[k]))+"\nDepth: "+str(Sort.chanpos[Sort.maxchan[k]][1])+"um", fontsize=15)
        
    #plt.show()
    ##quit()

def plot_rasters(self):
    
    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    #UPDATE PARAMS FROM CURRENT WIDGET TEXTBOXES;  IS THIS NECESSARY? YES, BECAUSE ACTIVELY CHANGING PARAMS ON SCREEN
    
    #print self.animal.name
    #print self.animal.recName
    #self.animal.name = self.animal_name.text()
    #self.animal.recName = self.root_dir+self.animal.name+'/rhd_files/'+self.rec_name.text()
        
    #Load SUA Sort
    sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    n_spikes = []
    for k in range(len(Sort_sua.units)):
        n_spikes.append(len(Sort_sua.units[k]))
    n_spikes = np.array(n_spikes)
    indexes = np.argsort(n_spikes)
    print indexes

    #Load LFP Sort
    lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(lfp_file) #Auto load flag for Nick's data

    ax = plt.subplot(1, 1, 1)
    y = []
    for i in indexes: #range(len(Sort_sua.units)):
        x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)/float(Sort_sua.samplerate) #float(Sort1.samplerate)*2.5

        ymin=np.zeros(len(Sort_sua.units[indexes[i]]))
        ymax=np.zeros(len(Sort_sua.units[indexes[i]]))
        ymin+=i+0.4
        ymax+=i-0.4

        plt.vlines(x, ymin, ymax, linewidth=1, color='black') #colors[mod(counter,7)])

        y.append(x)

    #Plot LFP spike
    y = []
    for i in range(len(Sort_lfp.units)):
        x = np.array(Sort_lfp.units[i],dtype=np.float32)/float(Sort_sua.samplerate)*50#float(Sort1.samplerate)*2.5
        
        ymin=np.zeros(len(Sort_lfp.units[0]))
        ymax=np.zeros(len(Sort_lfp.units[0]))
        
        ymin+=-i*2-5.
        ymax+=-i*2-7.
        
        plt.vlines(x, ymin, ymax, linewidth=3, color=colors[i%9]) #colors[mod(counter,7)])
    
    plt.xlabel('Time (seconds)',fontsize=35, weight='bold')
    plt.ylabel('Single Unit ID',multialignment='center', fontsize=35)

    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    plt.xlim(0,300)

    plt.ylabel('LFP Cluster Raster         Single Unit IDs',multialignment='center', fontsize=35, weight='bold')
    
    ax.tick_params(axis='both', which='major', labelsize=30)

    plt.show()
    
#Function convolves spike times with a 20ms gaussian; 
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def drift_movies(self):

    '''Make sequence order movies
    '''

    #****************************************************************************************
    #**************************** COMPUTE DEPTH, TEMPLATES, CELL TYPING *********************
    #****************************************************************************************

    plotting = True

    from scipy.interpolate import InterpolatedUnivariateSpline

    min_spikes = float(self.min_spikes.text())
    min_fire_rate = float(self.min_fire_rate.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
    
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    #for k in range(total_units): 
    peaks=[]
    troughs=[]
    for k in range(total_units): 
        print "... unit: ", k , "    #spikes: ", len(Sort_sua.units[k]), Sort_sua.xpos[k], Sort_sua.ypos[k]

        waveform = Sort_sua.wavedata[k][Sort_sua.maxchan[k]]

        #Interpolate max chan template to find fwhm values and do cell typing
        xi = np.arange(0, len(waveform), 1)
        yi = waveform

        x = np.linspace(0, len(waveform), 1000)
        order = 4
        s = InterpolatedUnivariateSpline(xi, yi, k=order)
        yinterp = s(x)

        #find peak width
        max_loc = np.argmax(yinterp); max_val = yinterp[max_loc]
        t2 = len(yinterp)-1
        for j in range(max_loc, len(yinterp),1):
            if yinterp[j]<(max_val/2.):
                t2 = j; break
        t1 = 0
        for j in range(max_loc, 0, -1):
            if yinterp[j]<(max_val/2.):
                t1 = j; break
        peak = x[t2]-x[t1]

        #find trough width
        min_loc = np.argmin(yinterp); min_val = yinterp[min_loc]
        t4 = len(yinterp)-1
        for j in range(min_loc, len(yinterp),1):
            if yinterp[j]>(min_val/2.):
                t4 = j; break
        t3 = 0
        for j in range(min_loc, 0, -1):
            if yinterp[j]>(min_val/2.):
                t3 = j; break
        trough = x[t4]-x[t3]

        #Compute cross product to find which side of line points lie on
        a = [0.4, 0]; b = [0, 0.3]
        c = [peak*0.04, trough*0.04]
        x_product = (b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0])

        if (x_product>0): color='blue'; m = 'o'
        else: color='red'; m = (3, 0) 

        if plotting: 
            #Plot scatter depth
            ax = plt.subplot(1,3,1)
            plt.scatter(Sort_sua.xpos[k], -Sort_sua.ypos[k], s=100, marker = m, color=color, alpha=.9)
            plt.plot([-100,100],[0,0], 'r--', color='black', linewidth=2, alpha=.7)
            plt.plot([-50,-50],[0,-2000], 'r--',color='black', linewidth=2, alpha=.7)
            plt.plot([50,50],[0,-2000], 'r--', color='black', linewidth=2, alpha=.7)
            plt.xlim(-100,100)
            plt.ylim(-2000,100)


            #Plot curve
            ax = plt.subplot(1,3,2)
            plt.plot(x*0.04-x[len(x)/2]*0.04, yinterp, linewidth=2, color=color, alpha=.9)
            plt.xlim(-0.5, 0.5)

            #Plot scatter plot distributions
            ax = plt.subplot(1,3,3)
            plt.scatter(peak*0.04, trough*0.04, color=color)
            plt.ylim(0,.5)
            plt.xlim(0,.6)

    if plotting: 
            ax = plt.subplot(1,3,1)
            plt.axhspan(-2000, 0, facecolor='black', alpha=0.1)
            plt.show()        


    #****************************************************************************************
    #********************************* COMPUTE DRIFT VALUES *********************************
    #****************************************************************************************

    selected_unit = int(self.starting_cell.text())

    #Find recording length; needed for plotting distributions
    if os.path.exists(self.parent.sua_file.replace('.ptcs','.tsf')):
        tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
        rec_length = tsf.n_vd_samples/float(tsf.SampleFrequency)
    else:
        rec_length = 0
        for k in range(len(Sort_sua.units)):
            if np.max(Sort_sua.units[k])>rec_length: rec_length = np.max(Sort_sua.units[k])
        rec_length = rec_length*1.E-6
    
    #Load LFP Sort
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data
    lfp_cluster = int(self.parent.lfp_cluster.text())

    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events              
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    
    n_units = int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**********************************************************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ********************************************************
    #**********************************************************************************************************************

    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    #temp_chunks.append(Sort_lfp.units[lfp_cluster][-1])

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]

    #Loop over chunks
    sig = float(self.sigma_width.text())

    #Make list to hold control peaks
    control_array = []
    for k in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1): control_array.append([])

    for ctr in range(len(time_chunks)):
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins."
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        
        for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):
            locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
            all_spikes = np.sort(np.hstack(locked_spikes))*1E-3

            if float(len(all_spikes))/(float(time_chunk[1])*1E-6-float(time_chunk[0])*1E-6) >= min_fire_rate:    #**********NB: only include chunk for cell if sufficient spikes

                fit_even = np.zeros(2000, dtype=np.float32)
                fit_odd = np.zeros(2000, dtype=np.float32)
                
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(all_spikes)):
                    mu = int(all_spikes[g])
                    if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                    else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

                control_array[unit-int(self.starting_cell.text())].append([x[950+np.argmax(fit_even[950:1050])],x[950+np.argmax(fit_odd[950:1050])]])  #******NB: LIMITING SEARCH TO 100ms window
            else:
                control_array[unit-int(self.starting_cell.text())].append([50, 50]) 
            
    control_array = np.float32(control_array)
    
    plotting = False
    if plotting: 
        fig, ax = plt.subplots()
        ax.ticklabel_format(useOffset=False, style='plain')

    #Compute distance of control to y=x line
    distances = []
    for unit in range(len(control_array)):
        distances.append([])
        for p in range(len(control_array[unit])):
            if (control_array[unit][p][0]!=50.) and (control_array[unit][p][1]!=50.):          #**********Exclude chunks/epochs with insufficient spikes
                distances[unit].append(abs(float(control_array[unit][p][0]-control_array[unit][p][1]))/np.sqrt(2))
        print unit, distances[unit]

    pts_array = []
    for unit in range(len(control_array)):
        fire_rate = len(Sort_sua.units[unit])/float(rec_length)
        control_ave = np.average(distances[unit], axis=0)
       
        if control_ave<1: color='blue'
        elif control_ave<2: color='green'
        else: 
            color = 'black'
            pts_array.append([])
            continue

        print "...unit: ", unit, "  ave MSL error: ", control_ave, "   fire rate: ", fire_rate

        pts = [];  ctr = 0
        for k in range(len(control_array[unit])):
            if (control_array[unit][k][0]==50) or (control_array[unit][k][1]==50):          #**********Exclude chunks/epochs with insufficient spikes
                pass
            else:
                pts.append([k, np.average(control_array[unit][k])])
                ctr+=1
        pts_array.append(pts)
        if plotting: 
            for k in range(len(pts)-1):
                plt.scatter(pts[k][0],pts[k][1], s=500, color=color)
                plt.scatter(pts[k+1][0],pts[k+1][1], s=500, color=color)
                plt.plot([pts[k][0],pts[k+1][0]], [pts[k][1], pts[k+1][1]], linewidth=5, color=color, alpha=.35)
    
    #plt.yscale('symlog', linthreshx=(-1E0,1E0))
    #plt.yscale('log')
    
    xtick_lbls = []
    for k in range(len(time_chunks)):
        xtick_lbls.append(int(time_chunks[k][1]*1E-6/60.))
    
    old_xlabel = np.arange(0, len(time_chunks), 1)
    plt.xticks(old_xlabel, xtick_lbls, fontsize=20) #,rotation='vertical')

    plt.xlabel("Time (mins)", fontsize=30)
    plt.ylabel("Phase Lock in Each Epoch (ms)", fontsize=30)

    plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(len(Sort_lfp.units[lfp_cluster]))+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
                '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)
    
    if plotting: 
        ax.tick_params(axis='both', which='major', labelsize=30)
        
        plt.ylim(-1E2,1E2)
        plt.xlim(0.1,len(time_chunks)+0.1)
        plt.show()


    #****************************************************************************************
    #************************************* MAKE DRIFTING MOVIES *****************************
    #****************************************************************************************

    ##REMOVE THIS AFTER COMPLETING MOVIES
    #for unit in range(len(control_array)):
        #fire_rate = len(Sort_sua.units[unit])/float(rec_length)
        #control_ave = np.average(distances[unit], axis=0)
       
        #if control_ave<1: color='blue'
        #elif control_ave<2: color='green'
        #else: color = 'black'; continue     #Skip cells with > 2ms control error values

        #print "...unit: ", unit, "  ave MSL error: ", control_ave, "   fire rate: ", fire_rate

        #pts = [];  ctr = 0
        #for k in range(len(control_array[unit])):
            #if (control_array[unit][k][0]==50) or (control_array[unit][k][1]==50):   #**********Exclude chunks/epochs with insufficient spikes
                #pass
            #else:
                #pts.append([k, np.average(control_array[unit][k])])
                #ctr+=1
        
        #for k in range(len(pts)-1):
            #plt.scatter(pts[k][0],pts[k][1], s=500, color=color)
            #plt.scatter(pts[k+1][0],pts[k+1][1], s=500, color=color)
            #plt.plot([pts[k][0],pts[k+1][0]], [pts[k][1], pts[k+1][1]], linewidth=5, color=color, alpha=.35)


    vid_array = np.zeros((len(control_array), time_chunks[len(time_chunks)-1][1]*1E-6/60., 2), dtype=np.float32)
    print vid_array.shape
    #pts_array contains [epoch, MSL time] pairs - or empty; len(control_array) = # epochs
    for p in range(len(control_array)):
        for k in range(len(pts_array[p])-1):
            t1 = time_chunks[pts_array[p][k][0]][0]*1E-6/60.
            t2 = time_chunks[pts_array[p][k][0]][1]*1E-6/60.
            print p, k, pts_array[p][k], pts_array[p][k+1], t1, t2
            
            #Interpolate positions over epoch beginnning and ends;
            interp_drift = np.linspace (pts_array[p][k][1], pts_array[p][k+1][1], int(t2)-int(t1))
            print interp_drift
            ctr = 0
            for q in range(int(t1), int(t2)):
                vid_array[p, q] = [interp_drift[ctr], 1]; ctr+=1
                
        print vid_array[p]
        plt.close()
        plt.plot(vid_array[p])
        plt.show()
        return
    return
        #for k in range(len(pts)-1):
            


    #
    xtick_lbls = []
    for k in range(len(time_chunks)):
        xtick_lbls.append(int(time_chunks[k][1]*1E-6/60.))
    
    old_xlabel = np.arange(0, len(time_chunks), 1)
    plt.xticks(old_xlabel, xtick_lbls, fontsize=20) #,rotation='vertical')
    
    
    
    




def drift_trends(self):
    
    '''Plot single cell drift location relative MSL
    '''
    
    
    min_spikes = float(self.min_spikes.text())
    min_fire_rate = float(self.min_fire_rate.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
    
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    for k in range(total_units): print "... unit: ", k , "    #spikes: ", len(Sort_sua.units[k])
    selected_unit = int(self.starting_cell.text())
    
    #Find recording length; needed for plotting distributions
    if os.path.exists(self.parent.sua_file.replace('.ptcs','.tsf')):
        tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
        rec_length = tsf.n_vd_samples/float(tsf.SampleFrequency)
    else:
        rec_length = 0
        for k in range(len(Sort_sua.units)):
            if np.max(Sort_sua.units[k])>rec_length: rec_length = np.max(Sort_sua.units[k])
        rec_length = rec_length*1.E-6
    
    #Load LFP Sort
    #lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    #end_lfp = min(int(self.parent.end_lfp.text()),len(Sort_lfp.units)-1)
    #start_lfp = min(int(self.parent.start_lfp.text()), len(Sort_lfp.units))
    #end_lfp = min(int(self.parent.end_lfp.text()), len(Sort_lfp.units))
    lfp_cluster = int(self.parent.lfp_cluster.text())

    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    
    n_units = int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    #temp_chunks.append(Sort_lfp.units[lfp_cluster][-1])

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]


    #Load chunks into data
    #chunk_index = []
    #for chk in self.chunks_to_plot.text().split(','):
        #chunk_index.append(int(chk))
    
    #print "...chunk_index: ", chunk_index
    
    #**************************************************************************
    #************************ LOOP OVER CHUNKS ********************************
    #**************************************************************************


    sig = float(self.sigma_width.text())
    #for chunk_ctr in range(int(self.chunks_to_plot.text())):

    #Make list to hold control peaks
    control_array = []
    for k in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1): control_array.append([])

    for ctr in range(len(time_chunks)):
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins."
        #ax = plt.subplot(1, len(time_chunks), ctr+1)
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        #print temp3
        
        for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):
            locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
            all_spikes = np.sort(np.hstack(locked_spikes))*1E-3

            if float(len(all_spikes))/(float(time_chunk[1])*1E-6-float(time_chunk[0])*1E-6) >= min_fire_rate:    #**********NB: only include chunk for cell if sufficient spikes

                fit_even = np.zeros(2000, dtype=np.float32)
                fit_odd = np.zeros(2000, dtype=np.float32)
                
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(all_spikes)):
                    mu = int(all_spikes[g])
                    if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                    else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

                #t = np.linspace(-1000, 1000, 2000)
                control_array[unit-int(self.starting_cell.text())].append([x[950+np.argmax(fit_even[950:1050])],x[950+np.argmax(fit_odd[950:1050])]])  #******NB: LIMITING SEARCH TO 100ms window
            else:
                control_array[unit-int(self.starting_cell.text())].append([50, 50]) 
            
    control_array = np.float32(control_array)
                
                
    fig, ax = plt.subplots()
    ax.ticklabel_format(useOffset=False, style='plain')

    #Compute distance of control to y=x line
    distances = []
    for unit in range(len(control_array)):
        distances.append([])
        for p in range(len(control_array[unit])):
            if (control_array[unit][p][0]!=50.) and (control_array[unit][p][1]!=50.):          #**********Exclude chunks/epochs with insufficient spikes
                distances[unit].append(abs(float(control_array[unit][p][0]-control_array[unit][p][1]))/np.sqrt(2))
        print unit, distances[unit]

    highest_rate = 0
    lowest_rate = 1000
    for unit in range(len(control_array)):
        fire_rate = len(Sort_sua.units[unit])/float(rec_length)
        control_ave = np.average(distances[unit], axis=0)
       
        if control_ave<1: color='blue'
        elif control_ave<2: color='green'
        else: color = 'black'; continue

        print "...unit: ", unit, "  ave MSL error: ", control_ave, "   fire rate: ", fire_rate

        pts = [];  ctr = 0
        for k in range(len(control_array[unit])):
            if (control_array[unit][k][0]==50) or (control_array[unit][k][1]==50):          #**********Exclude chunks/epochs with insufficient spikes
                pass
            else:
                pts.append([k, np.average(control_array[unit][k])])
                ctr+=1
        
        for k in range(len(pts)-1):
            plt.scatter(pts[k][0],pts[k][1], s=500, color=color)
            plt.scatter(pts[k+1][0],pts[k+1][1], s=500, color=color)
            plt.plot([pts[k][0],pts[k+1][0]], [pts[k][1], pts[k+1][1]], linewidth=5, color=color, alpha=.35)
    
    
    #plt.plot([lowest_rate,highest_rate], [1,1], 'r--', color='blue', linewidth=4, alpha=.7)
    #plt.plot([lowest_rate,highest_rate], [2,2], 'r--', color='green', linewidth=4, alpha=.7)

    #plt.yscale('symlog', linthreshx=(-1E0,1E0))
    #plt.yscale('log')
    ax.tick_params(axis='both', which='major', labelsize=30)
    
    plt.ylim(-1E2,1E2)
    plt.xlim(0.1,len(time_chunks)+0.1)
    
    xtick_lbls = []
    for k in range(len(time_chunks)):
        #print time_chunks[k]
        xtick_lbls.append(int(time_chunks[k][1]*1E-6/60.))
    
    old_xlabel = np.arange(0, len(time_chunks), 1)
    #new_xlabel = np.arange(-50, 51, 25)
    plt.xticks(old_xlabel, xtick_lbls, fontsize=20) #,rotation='vertical')
    #plt.xticks(np.arange(0,len(time_chunks),1), fontsize=30) #,rotation='vertical')    
    


    plt.xlabel("Time (mins)", fontsize=30)
    plt.ylabel("Phase Lock in Each Epoch (ms)", fontsize=30)
    

    plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(len(Sort_lfp.units[lfp_cluster]))+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
                '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)
    plt.show()

    

def all_cell_msl_stats(self):
    

    '''Plot single cell drift location relative MSL
    '''

    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    for k in range(total_units): print "... unit: ", k , "    #spikes: ", len(Sort_sua.units[k])
    selected_unit = int(self.starting_cell.text())


    #Find recording length; needed for plotting distributions
    if os.path.exists(self.parent.sua_file.replace('.ptcs','.tsf')):
        tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
        rec_length = tsf.n_vd_samples/float(tsf.SampleFrequency)
    else:
        rec_length = 0
        for k in range(len(Sort_sua.units)):
            if np.max(Sort_sua.units[k])>rec_length: rec_length = np.max(Sort_sua.units[k])
        rec_length = rec_length*1.E-6


    min_spikes = float(self.min_spikes.text())
    min_fire_rate = float(self.min_fire_rate.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
    

    
    
    #Load LFP Sort
    #lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    #end_lfp = min(int(self.parent.end_lfp.text()),len(Sort_lfp.units)-1)
    #start_lfp = min(int(self.parent.start_lfp.text()), len(Sort_lfp.units))
    #end_lfp = min(int(self.parent.end_lfp.text()), len(Sort_lfp.units))
    lfp_cluster = int(self.parent.lfp_cluster.text())

    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    
    n_units = int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    #temp_chunks.append(Sort_lfp.units[lfp_cluster][-1])

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]


    #Load chunks into data
    #chunk_index = []
    #for chk in self.chunks_to_plot.text().split(','):
        #chunk_index.append(int(chk))
    
    #print "...chunk_index: ", chunk_index
    
    #**************************************************************************
    #************************ LOOP OVER CHUNKS ********************************
    #**************************************************************************


    sig = float(self.sigma_width.text())
    #for chunk_ctr in range(int(self.chunks_to_plot.text())):

    #Make list to hold control peaks
    control_array = []
    for k in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1): control_array.append([])

    for ctr in range(len(time_chunks)):
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins."
        #ax = plt.subplot(1, len(time_chunks), ctr+1)
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        #print temp3
        
        for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):
            locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
            all_spikes = np.sort(np.hstack(locked_spikes))*1E-3

            if float(len(all_spikes))/(float(time_chunk[1])*1E-6-float(time_chunk[0])*1E-6) >= min_fire_rate:    #**********NB: only include chunk for cell if sufficient spikes

                fit_even = np.zeros(2000, dtype=np.float32)
                fit_odd = np.zeros(2000, dtype=np.float32)
                
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(all_spikes)):
                    mu = int(all_spikes[g])
                    if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                    else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

                t = np.linspace(-1000, 1000, 2000)
                control_array[unit].append([t[950+np.argmax(fit_even[950:1050])],t[950+np.argmax(fit_odd[950:1050])]])  #******NB: LIMITING SEARCH TO 100ms window
            else:
                control_array[unit].append([50, 50]) 
            
                
    #fig, ax = plt.subplots()
    #ax.ticklabel_format(useOffset=False, style='plain')

    #Compute distance of control to y=x line
    distances = []
    for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):
        distances.append([])
        if len(Sort_sua.units[unit])<min_spikes: continue
        for p in range(len(control_array[unit])):
            
            if control_array[unit][p]!=[50, 50]:          #**********Exclude chunks/epochs with insufficient spikes
                distances[unit].append(abs(float(control_array[unit][p][0]-control_array[unit][p][1]))/np.sqrt(2))
    

    highest_rate = 0
    lowest_rate = 1000
    for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):
        if len(Sort_sua.units[unit])<min_spikes: continue
        
        fire_rate = len(Sort_sua.units[unit])/float(rec_length)
        control_ave = np.average(distances[unit], axis=0)
        
        print "...unit: ", unit, "  ave MSL error: ", control_ave, "   fire rate: ", fire_rate
        
        color = 'black'
        if control_ave<1: color='blue'
        elif control_ave<2: color='green'
        
        
        plt.scatter(fire_rate, control_ave, s=100, color=color)
        
        if fire_rate> highest_rate: highest_rate = fire_rate
        if fire_rate< lowest_rate: lowest_rate = fire_rate
    
    plt.plot([lowest_rate,highest_rate], [1,1], 'r--', color='blue', linewidth=4, alpha=.7)
    plt.plot([lowest_rate,highest_rate], [2,2], 'r--', color='green', linewidth=4, alpha=.7)

    plt.yscale('symlog', linthreshx=1E0)
    #plt.yscale('log')
    plt.xscale('log')
    
    plt.ylim(0,1E2)
    plt.xlim(lowest_rate,highest_rate)
    
    plt.xlabel("Firing Rate (Hz)", fontsize=20)
    plt.ylabel("Control Error Averages (ms)", fontsize=20)
    

    plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(len(Sort_lfp.units[lfp_cluster]))+",  sigma: " + str(sig) +"(ms)" + ", min rate: "+ str(min_fire_rate) + " hz" +\
                '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=10)

    ax = plt.gca()
    #ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.tick_params(axis='both', which='major', labelsize=20)
    from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))


    plt.show()

    

def cell_msl_drift(self):
    
    '''Plot single cell drift location relative MSL
    '''
    
    min_spikes = float(self.min_spikes.text())
    min_fire_rate = float(self.min_fire_rate.text())
    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
    
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    for k in range(total_units): print "... unit: ", k , "    #spikes: ", len(Sort_sua.units[k])
    selected_unit = int(self.starting_cell.text())
    
    
    #Load LFP Sort
    #lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    #end_lfp = min(int(self.parent.end_lfp.text()),len(Sort_lfp.units)-1)
    #start_lfp = min(int(self.parent.start_lfp.text()), len(Sort_lfp.units))
    #end_lfp = min(int(self.parent.end_lfp.text()), len(Sort_lfp.units))
    lfp_cluster = int(self.parent.lfp_cluster.text())

    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    print type(pop_spikes[0])
    
    n_units = int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    #temp_chunks.append(Sort_lfp.units[lfp_cluster][-1])

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]



    #Load chunks into data
    #chunk_index = []
    #for chk in self.chunks_to_plot.text().split(','):
        #chunk_index.append(int(chk))
    
    #print "...chunk_index: ", chunk_index
    
    #**************************************************************************
    #************************ LOOP OVER CHUNKS ********************************
    #**************************************************************************


    sig = float(self.sigma_width.text())
    #for chunk_ctr in range(int(self.chunks_to_plot.text())):

    msl_array = []
    control_array = []

    for ctr in range(len(time_chunks)):
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins.",
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]

        unit = int(self.starting_cell.text())
        locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
        all_spikes = np.sort(np.hstack(locked_spikes))*1E-3
        
           
        print " ...unit: ", unit, " #spikes locked: ", len(all_spikes), " / ", len(Sort_sua.units[unit]), \
        "   duplicates: ", len(all_spikes) - len(np.unique(all_spikes))
        
        if float(len(all_spikes))/(float(time_chunk[1])*1E-6-float(time_chunk[0])*1E-6) > min_fire_rate:    #**********NB: only include chunk for cell if sufficient spikes

            fit_even = np.zeros(2000, dtype=np.float32)
            fit_odd = np.zeros(2000, dtype=np.float32)
            
            x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
            sig_gaussian = np.float32(gaussian(x, 0, sig))
            for g in range(len(all_spikes)):
                #print "...spike: ", all_spikes[g]
                mu = int(all_spikes[g])
                if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

            #plt.plot(t, fit_even/np.max(fit_even)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='blue', linewidth=6, alpha=.6)
            #plt.plot(t, fit_odd/np.max(fit_odd)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='red', linewidth=6, alpha=.6)
            
            #bin_width = sig*1E3   # histogram bin width in usec
            #y = np.histogram(all_spikes, bins = np.arange(-1000000,1000000,bin_width))
            #plt.bar(y[1][:-1]*1E-3, np.float32(y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
            #plt.bar(y[1][:-1], y[0], bin_width, color='blue')
        
            fit_sum= (fit_even+fit_odd)/2.
            
            msl_array.append(fit_sum[950:1050]/np.max(fit_sum[950:1050]))
            t = np.linspace(-1000, 1000, 2000)
            control_array.append([t[950+np.argmax(fit_even[950:1050])],t[950+np.argmax(fit_odd[950:1050])]])  #******NB: LIMITING SEARCH TO 100ms window
        else:
            msl_array.append(np.zeros(100, dtype=np.float32))
            control_array.append([50, 50]) 
        
        print unit, control_array[ctr]

    ax = plt.subplot(1,2,1)

    #Plot jet colored MSL image
    max_locs = []
    max_y = []
    for k in range(len(msl_array)):
        if control_array[k]!=[50,50]:
            max_locs.append(np.argmax(msl_array[k]))        #SKipping periods that do not have sufficient spikes
            max_y.append(k)

    msl_array = np.vstack(msl_array)
    ax.imshow(msl_array, aspect='auto', cmap=plt.get_cmap('jet'), interpolation='none', alpha=0.4)

    #Plot location of maxima and connect them with lines
    for k in range(len(max_y)-1):
        plt.plot([max_locs[k],max_locs[k+1]],[max_y[k],max_y[k+1]], linewidth=10, color='black')
        plt.scatter([max_locs[k],max_locs[k+1]], [max_y[k],max_y[k+1]], s=400, color='black')

    plt.xlim([0,101])
    plt.ylim([-0.5, len(msl_array)-0.5])
    
    ax.set_yticks([])
    ax.set_ylabel("Early <-------Time--------> Later",  fontsize = 40, labelpad=0)
    
    old_xlabel = np.arange(0, 101, 25)
    new_xlabel = np.arange(-50, 51, 25)
    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #,rotation='vertical')
    plt.xlabel("Time from event (ms)", fontsize=30)
    
    plt.plot([50,50],[-0.5,len(msl_array)-0.5], 'r--', color='black', linewidth=3, alpha=.9)
    
    
    #Plot controls
    ax = plt.subplot(1,2,2)
    for k in range(len(control_array)):
        if control_array[k]!=[50,50]:           #Skipping perdios with insufficient spikes
            plt.scatter(control_array[k][0], control_array[k][1], s=500, color='black')
            print control_array[k][0], control_array[k][1]
        
    plt.plot([-50, 50], [-50, 50], 'r--', linewidth = 5, color='black', alpha=0.5)
    #min_limit = min(np.min(control_array[:][0]), np.min(control_array[:][1]))
    #max_limit = max(np.max(control_array[:][0]), np.max(control_array[:][1]))
    
    plt.xticks(np.arange(-50,51,25), fontsize=30) #,rotation='vertical')
    plt.yticks(np.arange(-50,51,25), fontsize=30) #,rotation='vertical')
    
    
    plt.xlim([-51,51])
    plt.ylim([-51,51])
    plt.tick_params(axis='both', which='both', labelsize=30)
    plt.xlabel("Even Peaks (ms)", fontsize=30)
    plt.ylabel("Odd Peaks (ms)", fontsize=30, labelpad=-20)



    #t_max = t[np.argmax(fit_sum_even)]
    #plt.plot([t_max,t_max], [0,5000], 'r--', linewidth = 6, color='black')
    #plt.title("#spks: "+ str(len(all_spikes)) + "  Lock time: "+str(round(t_max,1)))

    ##plt.xticks(list(plt.xticks()[0]) + [round(t_max,1)])
    #ax.xaxis.set_ticks([50, 50, round(t_max,1), 0, 10])

    #plt.plot([0,0], [0, n_lfp_spikes*n_units], 'r--', linewidth=3, color='black', alpha=.8)

    #plt.xlim(lock_window_start-1,lock_window_end+1)
    #plt.ylim(0, n_lfp_spikes*n_units)
    
    ##old_ylabel = np.linspace(n_lfp_spikes/2.,n_lfp_spikes*n_units-n_lfp_spikes/2.,n_units)
    ##new_ylabel = np.arange(1, 11,1)
    ##plt.yticks(old_ylabel, new_ylabel, fontsize=30) #,rotation='vertical')
    
    ##old_xlabel = np.linspace(-100000, 100000, 101)
    ##new_xlabel = np.linspace(-100, 100,101)
    ##plt.xticks(old_xlabel, new_xlabel, fontsize=30) #,rotation='vertical')
    #plt.yticks([])
    ##ax.xaxis.set_major_locator(MaxNLocator(9)) 

    ##start, end = ax.get_xlim()
    ##ax.xaxis.set_ticks(np.linspace(start, end, 9))
    ##ax.locator_params(nbins=9, axis='x')
    ##ax.locator_params(nbins=9, axis='y')


    #plt.ylabel("Time Chunk: "+str(int(time_chunk[0]*1E-6/60.)) + ".."+str(int(time_chunk[1]*1E-6/60.))+" mins", fontsize=30,  labelpad=-2)
    #plt.xlabel("Time from event (msec)", fontsize = 30)
    #plt.tick_params(axis='both', which='both', labelsize=30)

    plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(len(Sort_lfp.units[lfp_cluster]))+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
                '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)
    plt.show()



    

def peth_scatter_plots(self):
    
    '''Plot peri-event time histograms using LFP events
    '''
    
    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
    
    ##Loading .tsf to compute length of record
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    #time_chunks=[]
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print temp_chunks
    
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    for k in range(total_units): print "... unit: ", k , "    #spikes: ", len(Sort_sua.units[k])

    #Load LFP Sort
    #lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    #end_lfp = min(int(self.parent.end_lfp.text()),len(Sort_lfp.units)-1)
    #start_lfp = min(int(self.parent.start_lfp.text()), len(Sort_lfp.units))
    #end_lfp = min(int(self.parent.end_lfp.text()), len(Sort_lfp.units))
    lfp_cluster = int(self.parent.lfp_cluster.text())

    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    #original_n_popspikes = len(pop_spikes)
    #pop_spikes=np.sort(np.unique(pop_spikes))       
    #print " ... # LFP events: ", len(pop_spikes), "   #duplicates: ", original_n_popspikes-len(pop_spikes)
    print type(pop_spikes[0])
    
    
    n_units = int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    temp_chunks.append(Sort_lfp.units[lfp_cluster][-1]*compress_factor)

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]
    

    chunk_index = []
    for chk in self.chunks_to_plot.text().split(','):
        chunk_index.append(int(chk))
    print "...chunk_index: ", chunk_index
    
    
    sig = float(self.sigma_width.text())
    ax = plt.subplot(1,1,1)
    for ctr, chunk_ctr in enumerate(chunk_index):
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[chunk_ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins."
        ax = plt.subplot(1, len(chunk_index), ctr+1)
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        #print temp3
        
        #for unit in range(total_units):
        for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):

            locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
            print type(locked_spikes[0])
            
               
            print "... chunk: ", time_chunk, " ...unit: ", unit, " #spikes locked: ", len(np.hstack(locked_spikes)), " / ", len(Sort_sua.units[unit]), \
            "   duplicates: ", len(np.hstack(locked_spikes)) - len(np.unique(np.hstack(locked_spikes)))

            #Plot rasters
            for event in range(len(locked_spikes)):
                ymin=np.zeros(len(locked_spikes[event]))
                ymax=np.zeros(len(locked_spikes[event]))
                ymin+=offset
                ymax+=offset+10
                offset+=1

                plt.vlines(np.float64(locked_spikes[event])*1E-3, ymin, ymax, linewidth=5, color='black',alpha=1) #colors[mod(counter,7)])
            plt.plot([-100,100], [offset,offset], linewidth=3, color='black', alpha=.8)
            
            #Don't convert xx1 into stack until after plotting rasters
            n_lfp_spikes= len(locked_spikes)

            #if len(xx1)>0: xx1 = np.hstack(xx1)
            #cell_rasters[chunk_ctr].append(xx1)
            
            all_spikes = np.sort(np.hstack(locked_spikes))
            if len(all_spikes)>0:
                #print all_spikes[:30]
                #Convolve w. gaussian
                fit_sum_even = np.zeros(2000000, dtype=np.float32)
                fit_sum_odd = np.zeros(2000000, dtype=np.float32)
                
                x = np.linspace(-1000,1000, 2000000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(all_spikes)):
                    print "...spike: ", all_spikes[g]
                    mu = int(all_spikes[g])
                    if g%2==0: fit_sum_even += np.roll(sig_gaussian, mu)
                    else: fit_sum_odd += np.roll(sig_gaussian, mu , axis=0)
                

                t = np.linspace(-1000, 1000, 2000000)
                
                if True: 
                    plt.plot(t, fit_sum_even/np.max(fit_sum_even)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='blue', linewidth=6, alpha=.6)
                    plt.plot(t, fit_sum_odd/np.max(fit_sum_odd)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='red', linewidth=6, alpha=.6)
                else:
                    fit_sum_even = (fit_sum_even+fit_sum_odd)/2.
                    plt.plot(t, fit_sum_even/np.max(fit_sum_even)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='blue', linewidth=7, alpha=.9)

                    
                #bin_width = sig*1E3   # histogram bin width in usec
                #y = np.histogram(all_spikes, bins = np.arange(-1000000,1000000,bin_width))
                #plt.bar(y[1][:-1]*1E-3, np.float32(y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
                #plt.bar(y[1][:-1], y[0], bin_width, color='blue')

        t_max = t[np.argmax(fit_sum_even)]
        plt.plot([t_max,t_max], [0,5000], 'r--', linewidth = 6, color='black')
        plt.title("#spks: "+ str(len(all_spikes)) + "  Lock time: "+str(round(t_max,1)))

        plt.xticks(list(plt.xticks()[0]) + [round(t_max,1)])
        ax.xaxis.set_ticks([50, 50, round(t_max,1), 0, 10])

        plt.plot([0,0], [0, n_lfp_spikes*n_units], 'r--', linewidth=3, color='black', alpha=.8)

        plt.xlim(lock_window_start-1,lock_window_end+1)
        plt.ylim(0, n_lfp_spikes*n_units)
        
        #old_ylabel = np.linspace(n_lfp_spikes/2.,n_lfp_spikes*n_units-n_lfp_spikes/2.,n_units)
        #new_ylabel = np.arange(1, 11,1)
        #plt.yticks(old_ylabel, new_ylabel, fontsize=30) #,rotation='vertical')
        
        #old_xlabel = np.linspace(-100000, 100000, 101)
        #new_xlabel = np.linspace(-100, 100,101)
        #plt.xticks(old_xlabel, new_xlabel, fontsize=30) #,rotation='vertical')
        plt.yticks([])
        #ax.xaxis.set_major_locator(MaxNLocator(9)) 

        #start, end = ax.get_xlim()
        #ax.xaxis.set_ticks(np.linspace(start, end, 9))
        #ax.locator_params(nbins=9, axis='x')
        #ax.locator_params(nbins=9, axis='y')


        plt.ylabel("Time Chunk: "+str(int(time_chunk[0]*1E-6/60.)) + ".."+str(int(time_chunk[1]*1E-6/60.))+" mins", fontsize=30,  labelpad=-2)
        plt.xlabel("Time from event (msec)", fontsize = 30)
        plt.tick_params(axis='both', which='both', labelsize=30)

    plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(n_lfp_spikes)+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
                ",  sigma: " + str(sig) +"(ms)" + '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)
    plt.show()



def msl_plots(self):
    '''Align msl latencies by lfp cluster 
    '''
    
    #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'

    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1000000  # window width in usec
    
    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())


    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())
        


    
    #UPDATE PARAMS FROM CURRENT WIDGET TEXTBOXES
    #self.parent.name = self.animal_name.text()
    #self.parent.recName = self.root_dir+self.animal.name+'/rhd_files/'+self.rec_name.text()
        
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    print type(Sort_sua.units[0][-1])

    #Load LFP Sort
    #lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    #end_lfp = min(int(self.parent.end_lfp.text()),len(Sort_lfp.units)-1)
    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardcoded to: ", compress_factor
    #try:
    #    self.subsample
    #except NameError:
    #    self.subsample = 1.0

    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())
      
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    #print Sort_lfp.units[lfp_cluster]
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    original_n_popspikes = len(pop_spikes)
    pop_spikes=np.sort(np.unique(pop_spikes))       
    print " ... # LFP events: ", len(pop_spikes)
    print type(pop_spikes[0])

    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    

    #OPTION 2: Divide into chunks of LFP events
    n_spikes = len(Sort_lfp.units[lfp_cluster])
    temp_chunks=[]
    chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    for t in range(0, n_spikes, chunk_width):
        temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
    temp_chunks.append(Sort_lfp.units[lfp_cluster][-1]*compress_factor)

    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    #int(self.starting_cell.text())      


    #OPTION 3: Divide into chunks of single unit spikes; DO LOCKED SPIKES VS ALL SPIKES
    #n_spikes = len(Sort_sua.units[selected_unit])
    #temp_chunks=[]
    #chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
    #for t in range(0, n_spikes, chunk_width):
        #temp_chunks.append(Sort_sua.units[selected_unit][t])
    #temp_chunks.append(Sort_sua.units[selected_unit][-1])

    #time_chunks = []
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print time_chunks[:10]

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    if os.path.exists(cell_rasters_filename+'.npy')==False:
        lfp_ctr=0
        
        print "LFP Cluster # : ", lfp_cluster, " / ", len(Sort_lfp.units)
          

        ##Compute periods of synchrony from si index                    #***********************************REIMPLEMENT ASAP
        #data_in = lfp.data[rec_index][9]
        #si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)
        #temp_list = []
        #for p in range(len(sync_periods)):
        #    indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
        #    temp_list.extend(pop_spikes[indexes])
        #pop_spikes=np.array(temp_list)


        #CELL RASTER LISTS TO HOLD SPIKES FOR EACH LFP EPOCH
        cell_rasters = []
        for unit in range(total_units):         #PREPROCESS ALL UNIT MSLs

            #Load unique track-wide unit id 
            #unique_unit = Sort_sua.uid[unit]  #NOT USED FOR TRACK WIDE CONCATENATED SORTS
            
            #Load sua spikes during synch periods (secs); use default unit
            spike_array = Sort_sua.units[unit]  #This loads into usec
            original_nspikes = len(spike_array)
            spike_array = np.sort(np.unique(spike_array))       
            print " ... Unit # : ", unit, " total # spikes: ", len(spike_array), "   #duplicates: ", original_nspikes-len(spike_array)
            spike_array = np.int64(spike_array)
            
            if len(spike_array)<min_spikes:
                cell_rasters.append([]) 
                continue     #If too few spikes in cell exclude - BUT SHOULD ADD BLANK SPACE---!!!!!!!!!!!!!

            #temp3 = np.where(np.logical_and(spike_array>=time_chunk[0], spike_array<=time_chunk[1]))[0]
            #spike_array= spike_array[temp3]
            
            #print "...unit: ", unit, " uid: ", unique_unit,  " #spikes: ", len(spike_array), " / ", len(Sort_sua.units[unit])
            

            #SKIP SYNC PERIOD FOR NOW                        #***********************************REIMPLEMENT ASAP
            #temp_list = []
            #for p in range(len(sync_periods)):
            #    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
            #    temp_list.extend(spike_array[indexes])
            #spike_array=np.array(temp_list)


            xx_even=[]          #collect even spikes                    
            xx_odd=[]           #collect odd spikes                     
            locked_spikes=[]              #collect all spikes for KS stat test
            for j in range(len(pop_spikes)):
                #Skip pop spikes that occur w/in 100ms of each other
                if j<(len(pop_spikes)-1):
                    if (pop_spikes[j+1]-pop_spikes[j])<100000: 
                        locked_spikes.append([])                    #*************************NB: ADDING NO LOCKING FOR CLOSE LFP CLUSTERS; OVERDOING IT
                        print "... lfp events too close, skipping..."
                        continue

                #find spikes that fall w/in +/- window of pop event
                temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                x=spike_array[temp2]-pop_spikes[j]  #Offset time of spike to t_0; stay in usec units
                locked_spikes.append(x)

            #Save cell rasters for each epoch chunk
            cell_rasters.append(locked_spikes)             #************** IS THIS MAINTAINING THE SAME CELL ORDER ACROSS CHUNKS?! CHECK!!!!!!!!
       
        np.save( cell_rasters_filename, cell_rasters)
        print cell_rasters[0][0]
        print type(cell_rasters[0][0])
    
    else:
        cell_rasters = np.load(cell_rasters_filename+".npy")


    

    chunk_index = []
    for chk in self.chunks_to_plot.text().split(','):
        chunk_index.append(int(chk))
    
    print "...chunk_index: ", chunk_index
    
    ##Loading .tsf to compute length of record
    #self.parent.tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    #temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    #time_chunks=[]
    #for t in range(len(temp_chunks)-1):
        #time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    #print temp_chunks
    



    sig = float(self.sigma_width.text())
    for ctr, chunk_ctr in enumerate(chunk_index):
        print len(time_chunks), chunk_ctr
        time_chunk = time_chunks[chunk_ctr]
        
        print time_chunk
        fit_sum = np.zeros((total_units,2000), dtype=np.float32)

        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        
        for unit in range(total_units):
        #for unit in range(1):
            if len(Sort_sua.units[unit])<min_spikes: continue
            locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])
            
            n_spikes_pre = len(locked_spikes)
            locked_spikes = np.unique(np.sort(locked_spikes))
               
            print "... chunk: ", time_chunk, " ...unit: ", unit, " #spikes locked: ", len(locked_spikes), " / ", len(Sort_sua.units[unit]), \
            "   duplicates: ", n_spikes_pre - len(locked_spikes)

            #NEW METHOD: JUST COMPUTE GAUSSIAN ONCE THEN TRANSLATE THE ARRAY
            if len(locked_spikes)>0 : 
                #for g in range(len(locked_spikes)):
                #    mu = np.array(locked_spikes[g])*1E-3        #Convert back to milliseconds; No need for MSL to be in usec
                #    fit_sum[unit] += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)
            
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(locked_spikes)):
                    mu = int(locked_spikes[g]*1E-3)
                    fit_sum[unit] += np.roll(sig_gaussian, mu)
                                
        img=[]
        lock_time=[]
        for unit in range(total_units):
            if np.max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end])>0:
                lock_time.append(np.argmax(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
                img.append(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]/max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
            else:
                lock_time.append(lock_window_end*4)
                img.append(np.zeros(lock_window_end-lock_window_start, dtype=np.float32))
                
        
        #ORDER MSL IMG BY LOCK TIME OF FIRST EPOCH
        if (ctr ==0): inds = np.array(lock_time).argsort()
        print "Order: ", inds

        
        img=np.array(img)[inds]
        temp_img = []
        for i in range(len(img)):
            #if np.max(img[i])!=0:
                temp_img.append(img[i])
        img=np.array(temp_img)
        
        
        #********** PLOTING ROUTINES **************
        #ax=plt.subplot(1,int(self.chunks_to_plot.text()),chunk_ctr+1)
        ax = plt.subplot(1, len(chunk_index), ctr+1)


        im = ax.imshow(img, origin='upper', extent=[0,(lock_window_end-lock_window_start), len(img),0], aspect='auto', interpolation='none')


        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if ctr ==0:
            xx = np.arange(0,(lock_window_end-lock_window_start)+1,(lock_window_end-lock_window_start)/2.)
            x_label = np.arange(lock_window_start,lock_window_end+1,(lock_window_end-lock_window_start)/2.)
            plt.xticks(xx,x_label, fontsize=25)
            plt.xlabel("Time from LFP event (ms)", fontsize=30)

            yy = np.arange(0,len(img),20)
            y_label = np.arange(0,len(img),20)
            plt.yticks(yy, y_label, fontsize=30)

            plt.ylabel("Cell #", fontsize=30)

        plt.title(str(int(time_chunk[0]/(60.*1E6)))+".."+str(int(time_chunk[1]/(60.*1E6)))+" mins", fontsize=15)

        plt.ylim(len(inds),0)

        
        #plt.ylabel("Neuron (lock order)", fontsize=30)
        plt.plot([-lock_window_start,-lock_window_start],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.xaxis.labelpad = 0

    plt.suptitle(self.parent.sua_file.replace('.ptcs','')+",  sigma: " + str(sig) +"(ms)", fontsize=20)
    plt.show()

def compute_msl_pvals(self):
    
    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']


    #COMPUTE KS P-VALUES FOR CURRENT LFP CLUSTER
    from scipy import stats

    tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))
    Sort_lfp = Ptcs(self.parent.lfp_event_file) 
    Sort_sua = Ptcs(self.parent.sua_file) 
    
    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    #Load pop events during synch periods (secs)
    compress_factor = 50.
    print "... compress_factor is hardcoded to: ", compress_factor


    n_spikes = int(self.min_spikes.text()) #Min number of spikes in clusters considered

    print "... loading data from file..."
    
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+'.npy')        #array dimensions (no. cells, no. lfp events/triggers)
    
    #temp_rasters = []
    #for k in range(len(cell_rasters)):
        #temp = cell_rasters[k]
        #print temp
        #temp = np.unique(temp)
        #temp = np.sort(temp)
        #temp_rasters.append(temp)
        ##print "...duplicate spikes in epoch/unit: ", k, p, "   = ", original_len - len(cell_rasters[k][p]), "  /  ", original_len
    #cell_rasters = temp_rasters

    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = np.float32(Sort_lfp.units[lfp_cluster])/Sort_lfp.samplerate*compress_factor#*1E-3  
    original_n_popspikes = len(pop_spikes)
    pop_spikes=np.sort(np.unique(pop_spikes))       
    print " ... # LFP events: ", len(pop_spikes), "   #duplicates: ", original_n_popspikes-len(pop_spikes)
     
    #time_chunks = np.load(self.parent.sua_file.replace('.ptcs','')+"_time_chunks.npy")
    rec_length = tsf.n_vd_samples/float(tsf.SampleFrequency)
    time_chunks = [[0.0,3600.], [rec_length-3600, rec_length]]  #Look at 1st and last hour only.
    #time_chunks = [[0.0,3600.], [3600, 7200]]  #Look at 1st and last hour only.

    locked_spikes = []
    sig = float(self.sigma_width.text())
    for chunk_ctr in range(len(time_chunks)): 
    #for chunk_ctr in range(len(self.chunks_to_plot.text())):
        time_chunk = time_chunks[chunk_ctr]
        #print time_chunk
        locked_spikes.append([])
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        
        for unit in range(len(cell_rasters)):
            temp = np.hstack(cell_rasters[unit][temp3])
            temp_len = len(temp)
            temp = np.unique(temp)
            print "...duplicate spikes in epoch: ", time_chunk, "  unit: ", unit, "   = ", temp_len - len(temp), "  /  ", temp_len
            locked_spikes[chunk_ctr].append(temp)
    

    ax = plt.subplot(1,1,1)
    p_val_array = []
    for p in range(len(Sort_sua.units)):
        p_val_array.append([])
        if len(Sort_sua.units[p])<n_spikes: continue
        for c in range(1,len(time_chunks),1):
            print p, c
            if (len(locked_spikes[0][p])>0) and (len(locked_spikes[c][p])>0):
                KS, p_val = stats.ks_2samp(locked_spikes[0][p], locked_spikes[c][p])
                p_val_array[p].append(p_val)
            else:
                p_val_array[p].append(1.1)
    
    
    for p in range(len(Sort_sua.units)):
        if len(Sort_sua.units[p])<n_spikes: continue
        plt.scatter(p_val_array[p], -np.array([p]*len(p_val_array[p])), s=75, color=colors[p%9])
        #plt.scatter([p]*len(p_val_array[p]), p_val_array[p], color='black')
    
    plt.plot([5E-2, 5E-2], [0 , -len(Sort_sua.units)], 'r--', color='black', linewidth=3, alpha=0.7)
    plt.plot([1E-2, 1E-2], [0 , -len(Sort_sua.units)], 'r--', color='black', linewidth=4, alpha=0.7)
    plt.plot([1E-3, 1E-3], [0 , -len(Sort_sua.units)], 'r--', color='black', linewidth=5, alpha=0.7)
    
    plt.xlim(-1E-7,1)
    plt.tick_params(axis='both', which='both', labelsize=25)

    plt.xlabel("Kolmogorov-Smirnov 2-Sample P-Value", fontsize=30)
    plt.ylabel("Cell #", fontsize=30)
    
    plt.xscale('symlog', linthreshx=1E-7)
    
    plt.ylim(-len(Sort_sua.units),1)
    plt.show()

    #quit()




def sua_lock_percentage(self):

    self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'

    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    

    lock_window = int(self.parent.lock_window.text())
    
    #Loading low pass LFP; read just header and then required channel for sync index computation below
    tsf = Tsf_file(self.parent.sua_file.replace('_hp.ptcs','_lp.tsf')) 
    tsf.read_trace(int(self.specgram_ch.text()))
    print tsf.header
    print tsf.iformat
    print tsf.SampleFrequency
    print tsf.n_vd_samples
    print tsf.vscale_HP
    print len(tsf.ec_traces)
        
        
    print tsf.n_electrodes
    
    
    data_in = tsf.ec_traces
    if len(tsf.ec_traces)==0: print "... channel incorrect..."; return
        
    #Load SUA Sort
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    #Load LFP Sort
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data

    start_lfp = min(int(self.parent.start_lfp.text()), len(Sort_lfp.units))
    end_lfp = min(int(self.parent.end_lfp.text()), len(Sort_lfp.units))
    
   
    
    #******************* OLD CODE STARTS ****************

    window=1
    small_window = 100 # ms of window for luczak plots
    large_window = window*1E3
    #start_lfp = 0; end_lfp = len(Sorts_lfp[0])# 

    plotting=True
    plotting_sync = False
    start_window = -0.04                            #*********************************** CODE THESE INTO THE GUI AS TEXT BOXES
    end_window = 0.01
    min_lfp_isi = 0.100   #Lockout lfp events more 
    
    si_limit = 0.7                                  #*********************************** MAYBE ALSO THIS!?
    
    
    ##Convert rec name from string to relative index in concatenated data; 
    ##also compute # recordings w. minimum required lfp events
    #rec_indexes=[]
    #n_lfprecs = 0
    #for rec in recs:
        #rec_indexes.append(lfp.rec.index(rec))  #Function .index provides the index of (rec) in lfp.rec
        #if len(Sorts_lfp[lfp.rec.index(rec)][start_lfp])>5: n_lfprecs+=1
    #if n_lfprecs==0:
        #print "No lfp_event recordings!"
        #return


    #compute track recording length:
    track_length = tsf.n_vd_samples/float(tsf.SampleFrequency)
    print "...rec length: ", track_length
    

    #Search all SUA sorts to find max # units 
    total_units = len(Sort_sua.units)
    #for rec_index in rec_indexes:
    #    if max(Sorts_sua[rec_index].uid)>total_units: total_units = max(Sorts_sua[rec_index].uid)
    #total_units +=1 #Adjust for zero based indexes
    
    #Collect all locked spikes for each LFP cluster
    total_locked_spikes_allrecs=np.zeros((end_lfp-start_lfp, total_units),dtype=np.float32)
    
    #Make array to collect all single unit spikes:
    sua_allspikes = np.zeros(total_units, dtype=np.float32)


    #********************* COMPUTE SYNC PERIODS AND ONLY SELECT POP SPIKES DURING THOSE PERIODS *************************
    samp_freq = 1000 #Sample frequency
    print "... computing sync index..."
    filename_sync_values = self.parent.sua_file.replace('_hp.ptcs','_lp_syncindex_ch'+self.specgram_ch.text()+'.npz')
    if os.path.exists(filename_sync_values):
        data = np.load(filename_sync_values)
        si = data['arr_0']; t = data['arr_1']; sync_periods = data['arr_2'] #Load default array names; if time, save proper variable names in .npz file
    else:
        si, t, sync_periods = synchrony_index(data_in, samp_freq, si_limit)
        np.savez(filename_sync_values[:-4], si, t, sync_periods)

    if plotting_sync: 
        plt.plot(t, si*10-10, linewidth=3, color='black')
        plt.plot([t[0],t[-1]],[0,0], linewidth=2, color='black')
        plt.plot([t[0],t[-1]],[-10,-10], linewidth=2, color='black')
        plt.plot([t[0],t[-1]],[-10*si_limit,-10*si_limit], 'r--', linewidth=2, color='red')
        
        print sync_periods
        
        print "computing specgram..."
        f0 = 0.1; f1 = 110
        p0 = -60. #float(self.parent.specgram_db_clip.text())
        P, extent = Compute_specgram_signal(data_in, samp_freq, f0, f1, p0)
        plt.imshow(P, extent=extent, aspect='auto')
        plt.show()

    #Save total length of sync periods; DO IT ONCE ONLY!
    sync_period_total_time = 0 #total period of synchronization across all recs (secs)
    for k in range(len(sync_periods)):
        sync_period_total_time += sync_periods[k][1]-sync_periods[k][0]


    #*********** LOOP OVER POP SPIKES
    all_pop_spikes = 0 #add up all pop spikes over all recordings to use as control for final plot
    
    for ss in range(start_lfp, end_lfp, 1):
        print ""
        print ""
        print "LFP Cluster # : ", ss+1, " / ", len(Sort_lfp.units)

        cluster_pop_spikes = 0  #Keep track of all pop spikes during sync periods for each LFP cluster
       
        #* LIMIT ANALYSIS TO SYNC PERIODS - HARRIS SPECIAL!!!
        #Make lists to hold # unique spikes that lock to lfp events
        n_lock_spikes = [[] for x in xrange(total_units)]

        #Load pop events during synch periods (secs)
        pop_spikes = np.array(Sort_lfp.units[ss])*1E-3      #Convert from 1Khz sampling rate to seconds
        print "... total lfp events: ", len(pop_spikes)
        temp_list = []
        for p in range(len(sync_periods)):
            indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
            temp_list.extend(pop_spikes[indexes])
        pop_spikes=np.array(temp_list)
        
        print "...no. of sync period lfp events: ", len(pop_spikes)
        
        
        #************************TRACK POP SPIKES *****************************
        #Track all pop_spikes for individual clusters
        cluster_pop_spikes+= len(pop_spikes)
        
        #Track cumulative total of pop_spikes
        all_pop_spikes+= len(pop_spikes)
        
        Sorts_sua_sync_spikes = np.zeros(total_units, dtype=np.float32)
        
        #Loop over all single units for each recording
        for unit in range(len(Sort_sua.units)):
            #Load unique track-wide unit id 
            #unique_unit = Sort_sua.uid[unit]   #NOT NEEDED FOR CONCATENATED SORTS

            #Load sua spikes during synch periods (secs); use default unit
            spike_array = np.array(Sort_sua.units[unit],dtype=np.float32)/float(Sort_sua.samplerate)  #This converts to seconds
            temp_list = []
            for p in range(len(sync_periods)):
                indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                temp_list.extend(spike_array[indexes])
            
            spike_array=np.array(temp_list)
            
            #Track total # spikes for each unit during synch periods   #NB: COUNT ONLY ONCE DURING MULTIPLE LFP LOOPS!!!
            if (ss-start_lfp)==0: sua_allspikes[unit]+=len(spike_array)
            
            #Save # of spikes during sync period; use sequential unit id - only used w/in a recording
            Sorts_sua_sync_spikes[unit]=len(spike_array)
            
            #NB: May wish to remove LARGE SPIKE TIMES due to bug 1E+8 for secs should do it.
            xx1=[]              #collect all spikes for KS stat test
            
            for j in range(len(pop_spikes)):

                #Skip pop spikes that occur w/in 100ms of each other
                if j<(len(pop_spikes)-1):
                    if (pop_spikes[j+1]-pop_spikes[j])<min_lfp_isi: continue
                
                #find spikes that fall w/in +/- window of pop event
                temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                xx1.append(x)

                #Add # spikes occuring w/in start_window..end_window (~50ms) of lfp event
                n_lock_spikes[unit].extend(spike_array[np.where(np.logical_and(spike_array>=pop_spikes[j]+start_window, spike_array<=pop_spikes[j]+end_window))[0]])


        #**********************************************************************************************
        #******************************** PERCENTAGE LOCK PLOTS ***************************************
        #**********************************************************************************************
        
        percent_lock = []
        for unit in range(len(Sort_sua.units)):
            n_spikes_unique = len(np.unique(n_lock_spikes[unit]))*100
            percent_lock.append(float(n_spikes_unique)/Sorts_sua_sync_spikes[unit])
            
            #Save number of unique spikes locked for each lfp cluster, each unit, and each recording
            total_locked_spikes_allrecs[ss-start_lfp][Sort_sua.uid[unit]] += n_spikes_unique*1E-2          #Add # locked spikes to total for each unique unit id
        
        percent_lock = np.array(percent_lock)
        

        #First, normalize the number of spikes for each unit:
        for u in range(total_units):
            if sua_allspikes[u]>0: #make sure recording looped over contain spikes from unit - not all recs do
                total_locked_spikes_allrecs[ss-start_lfp][u] = total_locked_spikes_allrecs[ss-start_lfp][u] / sua_allspikes[u]
            else:
                total_locked_spikes_allrecs[ss-start_lfp][u] = 0
        
        #Plot bar graphs by depth of cell
        ax1 = plt.subplot(2,1,1)
        indexes = np.arange(0,len(sua_allspikes),1)
        if (ss-start_lfp)==0: #Plot first bar graphs
            x_count = 0
            for u in indexes: 
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss])
                x_count+=1
                
        #Plot cumulative bar graphs for additional lfp clusters
        else:
            cumulative=np.zeros(len(indexes),dtype=np.float32)
            for u in indexes:
                for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                    cumulative[u] += total_locked_spikes_allrecs[p][u]
            x_count=0
            for u in indexes:
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
                x_count+=1

        #Plot bar graphs by firing rate order 
        ax1 = plt.subplot(2,1,2)
        
        #set indexes in order of firing rate; use master list of sua_allspikes, otherwise incorrect;
        indexes = sua_allspikes.argsort()
        
        #Plot lfp cluster 0
        if (ss-start_lfp)==0: #Plot first bar graphs
            x_count = 0
            for u in indexes: 
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss])
                x_count+=1
        
        #Plot cumulative bar graphs for additional lfp clusters
        else:
            cumulative=np.zeros(len(indexes),dtype=np.float32)
            for u in indexes:
                for p in range(0,ss-start_lfp):                         #Sum all previous % up to current LFP cluster ss
                    cumulative[u] += total_locked_spikes_allrecs[p][u]
            x_count=0
            for u in indexes:
                plt.bar(x_count, total_locked_spikes_allrecs[ss-start_lfp][u]*100., 1, color=colors[ss%10], bottom=cumulative[u]*100.)
                x_count+=1
                


    #Compute expected lock for allc ells also
    exp_lock = all_pop_spikes*(end_window-start_window)/sync_period_total_time*1E2
    print "total # pop spikes: ", all_pop_spikes
    print "track_length: ", track_length
    print "track_length - sync periods only: ", sync_period_total_time
    print "expected % lock: ", exp_lock


    #********** LABELS FOR DEPTH HISTOGRAMS
    ax1 = plt.subplot(2,1,1)
    plt.plot([0,len(sua_allspikes)],[exp_lock,exp_lock], 'r--', color='cyan', linewidth=2) #Plot cyan expected lock line
    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.tick_params(axis='x', which='both', labelsize=8)
    plt.ylabel("%locked")
    plt.xlabel("Depth of cell (relative Units)", fontsize=15)
    

    #********* LABELS FOR FIRE RATE HISTOGRAMS
    ax1 = plt.subplot(2,1,2)
    plt.plot([0,len(sua_allspikes)],[exp_lock,exp_lock], 'r--', color='cyan', linewidth=2) #Plot cyan expected lock line

    #Compute ratio of expected lock to actual lock and plot on top of bar graphs
    cumulative=np.zeros(len(indexes),dtype=np.float32)
    for u in indexes:
        for p in range(0,ss-start_lfp+1):                         #Sum all locking percentages over all LFP spikes
            cumulative[u] += total_locked_spikes_allrecs[p][u]

    x_count=0
    for u in indexes:
        #text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)
        text_temp = str(int(round(cumulative[u]/exp_lock*100.)))
        plt.text(x_count, cumulative[u]*100.+2 , text_temp, fontsize=7)
        x_count+=1

    #Use firing rates for each cell to label axes; NB: use sua_allspikes which contains only spikes during sync periods
    x_label = []
    sua_allspikes = sua_allspikes[indexes] #Reorder all units by firing rate
    for k in range(len(sua_allspikes)):
        x_label.append(str(round(sua_allspikes[k]/track_length,3)))
    xx = np.linspace(0,len(sua_allspikes)-1,len(sua_allspikes))+0.5

    plt.xlim(0,len(sua_allspikes))
    plt.ylim(0,100.0)
    plt.ylabel("%locked \n(%exp: "+str(round(exp_lock,2))+")")

    plt.xticks(xx, x_label, fontsize=8)
    plt.xticks(rotation=90)
    plt.xlabel("Firing rates - during sync periods (Hz)", fontsize=15)
    plt.suptitle(self.parent.sua_file+"\n # Clusters: "+ str(len(Sort_lfp.units))+", si_limit: "+str(si_limit)+ ", "+str(int(start_window*1E3))+
    "ms.."+str(int(end_window*1E3))+"ms.", fontsize=25)
    
        
    #SHOW PLOT
    plt.show()
    


def cell_count_matrix(self):
    print "count matrix"
    
    #n_spikes = int(self.n_spikes.text()) #Min number of spikes in clusters considered
    #cell_rasters = self.cell_rasters
    #time_chunks = self.t_chunks
       
    #Load SUA Sort
    if self.exp_type=='mouse': sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    elif self.exp_type =='cat': sua_file = self.animal.recName.replace('.tsf','.ptcs')

    print sua_file
    
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    for u in range(len(Sort_sua.units)):
        print u, Sort_sua.uid[u], len(Sort_sua.units[u])

    cell1 = int(self.start_cell.text())
    cell2 = int(self.end_cell.text())
    cell3 = int(self.other_cell.text())
    
    #for s in range(len(Sort_sua.units[start_cell])):
    spike_array1 = np.array(Sort_sua.units[cell1])/Sort_sua.samplerate #Work in seconds
    spike_array2 = np.array(Sort_sua.units[cell2])/Sort_sua.samplerate
    spike_array3 = np.array(Sort_sua.units[cell3])/Sort_sua.samplerate

    #w = 0.020       #size of bin in sec
    window =1   #size of window to search in sec
    count_matrix = np.zeros((window*4000,window*4000), dtype=np.int32)
    
    print len(spike_array1), len(spike_array2), len(spike_array3)
    for s1 in spike_array1:
        temp2 = np.where(np.logical_and((spike_array2-s1)>=-window, (spike_array2-s1)<=window))[0]
        #print "s1: ", s1
        for s2 in spike_array2[temp2]:
            #print "s2: ", s2
            temp3 = np.where(np.logical_and((spike_array3-s2)>=-window, (spike_array3-s2)<=window))[0]
            #print "s3: ", spike_array3[temp3]

            if len(temp3)>0: 
                t2_t1 = int((s2-s1)*1000)
                t3_t1 = np.int32((spike_array3[temp3]-s1)*1000)
            
                #print t2_t1, t3_t1
                #print t2_t1+2000, t3_t1+2000
                count_matrix[t2_t1+2000, t3_t1+2000]+=1
            
            #quit()
            
    #Plot matrix
    binning = int(self.isi_binning.text())     #ms per bin   
    print "... binning results..." 
    if binning == 1: 
        ave_count_matrix = np.array(count_matrix)
    else:
        ave_count_matrix = []
        for k in range(0, len(count_matrix), binning):
            temp_array = []
            for p in range(0, len(count_matrix[k]), binning):
                temp_array.extend([np.sum(count_matrix[k:k+binning,p:p+binning])])
            ave_count_matrix.append(temp_array)
        
        ave_count_matrix=np.vstack((ave_count_matrix))

    print ave_count_matrix.shape
    
    size_window = int(self.zoom.text()) #ms of view screen
    count_window = size_window  / binning
    print count_window
    len_matrix = len(ave_count_matrix)

    
    plt.imshow(ave_count_matrix[len_matrix/2 - count_window:len_matrix/2 + count_window
                                ,len_matrix/2 - count_window:len_matrix/2 + count_window], 
                                
                                extent=[0, count_window*2, count_window*2, 0],
                                interpolation=self.interpolation.text())
    
    old_xlabel = np.linspace(0,count_window*2,int(self.zoom.text())/10+1)
    new_xlabel = np.around(np.linspace(-size_window,size_window, int(self.zoom.text())/10+1), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=12,rotation='vertical')
    plt.yticks(old_xlabel, new_xlabel, fontsize=12)
    
    plt.ylim(-old_xlabel[0],old_xlabel[-1])
    plt.xlim(-old_xlabel[0],old_xlabel[-1])
    
    plt.show()
    
def cell_count_matrix_chunks(self):
    print "count matrix"
    
    #n_spikes = int(self.n_spikes.text()) #Min number of spikes in clusters considered
    #cell_rasters = self.cell_rasters
    #time_chunks = self.t_chunks
       
    #Load SUA Sort
    sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    for u in range(len(Sort_sua.units)):
        print u, len(Sort_sua.units[u])

    cell1 = int(self.start_cell.text())
    cell2 = int(self.end_cell.text())
    cell3 = int(self.other_cell.text())
    
    #Divide rasters into n chunks
    time_chunks = int(self.time_chunks.text())
    rec_len = self.animal.rec_length
    print time_chunks, rec_len

    spikes1 = np.array(Sort_sua.units[cell1])/Sort_sua.samplerate #Work in seconds
    spikes2 = np.array(Sort_sua.units[cell2])/Sort_sua.samplerate
    spikes3 = np.array(Sort_sua.units[cell3])/Sort_sua.samplerate
    spike_array1 = []; spike_array2 = []; spike_array3 = []
    
    
    
    t_chunks = np.linspace(0, rec_len, time_chunks+1)
    print t_chunks
    
    for t in range(len(t_chunks)-1):
        print t
        spike_array1.append(spikes1[np.where(np.logical_and(spikes1>=t_chunks[t], spikes1<=t_chunks[t+1]))[0]])
        spike_array2.append(spikes2[np.where(np.logical_and(spikes2>=t_chunks[t], spikes2<=t_chunks[t+1]))[0]])
        spike_array3.append(spikes3[np.where(np.logical_and(spikes3>=t_chunks[t], spikes3<=t_chunks[t+1]))[0]])
        
    #print spike_array1
    #print spike_array2
    #print spike_array3

    for c in range(time_chunks):
        window =1   #size of window to search in sec
        count_matrix = np.zeros((window*4000,window*4000), dtype=np.int32)
        print "TIME CHUNK: ", c
        for s1 in spike_array1[c]:
            #print s1
            temp2 = np.where(np.logical_and((spike_array2[c]-s1)>=-window, (spike_array2[c]-s1)<=window))[0]
            #print "s1: ", s1
            if len(temp2)==0: continue
            #print spike_array2[c][temp2]
            for s2 in spike_array2[c][temp2]:

                #print "s2: ", s2
                temp3 = np.where(np.logical_and((spike_array3[c]-s2)>=-window, (spike_array3[c]-s2)<=window))[0]
                #print "s3: ", spike_array3[temp3]

                if len(temp3)>0: 
                    t2_t1 = int((s2-s1)*1000)
                    t3_t1 = np.int32((spike_array3[c][temp3]-s1)*1000)
                
                    #print t2_t1, t3_t1
                    #print t2_t1+2000, t3_t1+2000
                    count_matrix[t2_t1+2000, t3_t1+2000]+=1
                
                #quit()
                
        #Plot matrix
        binning = int(self.isi_binning.text())     #ms per bin   
        print "... binning results..." 
        if binning == 1: 
            ave_count_matrix = np.array(count_matrix)
        else:
            ave_count_matrix = []
            for k in range(0, len(count_matrix), binning):
                temp_array = []
                for p in range(0, len(count_matrix[k]), binning):
                    temp_array.extend([np.sum(count_matrix[k:k+binning,p:p+binning])])
                ave_count_matrix.append(temp_array)
            
            ave_count_matrix=np.vstack((ave_count_matrix))

        size_window = int(self.zoom.text()) #ms of view screen
        count_window = size_window  / binning
        len_matrix = len(ave_count_matrix)

        ax = plt.subplot(2,2,c+1)
        plt.imshow(ave_count_matrix[len_matrix/2 - count_window:len_matrix/2 + count_window,
                                    len_matrix/2 - count_window:len_matrix/2 + count_window], 
                                    extent=[0, count_window*2, count_window*2, 0],
                                    interpolation=self.interpolation.text())
        
        old_xlabel = np.linspace(0, count_window*2, int(self.zoom.text())/10+1)
        new_xlabel = np.around(np.linspace(-size_window,size_window, int(self.zoom.text())/10+1), decimals=2)
        plt.xticks(old_xlabel, new_xlabel, fontsize=12, rotation='vertical')
        plt.yticks(old_xlabel, new_xlabel, fontsize=12)
        
        plt.title(str(int(t_chunks[c]/60.))+" .. "+str(int(t_chunks[c+1]/60.))+" (mins) "+
                        "(#spks:   "+str(len(spike_array1[c]))+', '+str(len(spike_array2[c])) +', '+str(len(spike_array3[c]))+')', fontsize=14)
        plt.ylim(-old_xlabel[0],old_xlabel[-1])
        plt.xlim(-old_xlabel[0],old_xlabel[-1])
    
    plt.suptitle(self.animal.name + '/'+self.animal.recName+'\n'+str(cell1)+' '+str(cell2)+' '+str(cell3))
    plt.show()
    


def all_cell_count_matrix(self):
    print "compute count matrix for all cells"
    
    #Load SUA Sort
    sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    for u in range(len(Sort_sua.units)):
        print u, len(Sort_sua.units[u])
    
    print self.animal.home_dir 
    
    from scipy import sparse, io
    import cPickle as pickle


    for p in range(len(Sort_sua.units)):
        spike_array1 = np.array(Sort_sua.units[p])/Sort_sua.samplerate #Work in seconds
        
        for q in range(p+1, len(Sort_sua.units), 1):
            spike_array2 = np.array(Sort_sua.units[q])/Sort_sua.samplerate
            
            for r in range(q+1, len(Sort_sua.units), 1):
                msl_file = self.animal.recName.replace('rhd_files','msl_files/data').replace('.rhd','')+str(p).zfill(3)+'_'+str(q).zfill(3)+'_'+str(r).zfill(3)
                if (os.path.exists(msl_file)==True): continue

                spike_array3 = np.array(Sort_sua.units[r])/Sort_sua.samplerate

                #w = 0.020       #size of bin in sec
                window =1   #size of window to search in sec
                count_matrix = np.zeros((window*4000,window*4000), dtype=np.int8)
                
                #print len(spike_array1), len(spike_array2), len(spike_array3)
                for s1 in spike_array1:
                    temp2 = np.where(np.logical_and((spike_array2-s1)>=-window, (spike_array2-s1)<=window))[0]
                    #print "s1: ", s1
                    for s2 in spike_array2[temp2]:
                        #print "s2: ", s2
                        temp3 = np.where(np.logical_and((spike_array3-s2)>=-window, (spike_array3-s2)<=window))[0]
                        #print "s3: ", spike_array3[temp3]

                        if len(temp3)>0: 
                            t2_t1 = int((s2-s1)*1000)
                            t3_t1 = np.int32((spike_array3[temp3]-s1)*1000)
                        
                            #print t2_t1, t3_t1
                            #print t2_t1+2000, t3_t1+2000
                            count_matrix[t2_t1+2000, t3_t1+2000]+=1
           
               
                count_matrix = sparse.csr_matrix(count_matrix)

                with open(msl_file, 'wb') as outfile:
                    pickle.dump(count_matrix, outfile, pickle.HIGHEST_PROTOCOL)


def view_all_cell_count_matrix(self):
    
    from scipy import sparse, io
    import cPickle as pickle

    print self.animal.home_dir 

    #Load SUA Sort
    sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    nspikes = int(self.cmatrix_nspikes.text())     #min # spikes in cell to be considered   
    
    for k in range(len(Sort_sua.units)):
        print k, len(Sort_sua.units[k])
    
    ctr=0
    for p in range(len(Sort_sua.units)):
        if len(Sort_sua.units[p])<nspikes:continue
        for q in range(p+1, len(Sort_sua.units), 1):
            if len(Sort_sua.units[q])<nspikes:continue
            for r in range(q+1, len(Sort_sua.units), 1):
                if len(Sort_sua.units[r])<nspikes:continue
                print ctr
                                
                ax = plt.subplot(5,5,ctr+1)
                msl_file = self.animal.recName.replace('rhd_files','msl_files/data').replace('.rhd','')+'_'+str(p).zfill(3)+'_'+str(q).zfill(3)+'_'+str(r).zfill(3)
                if (os.path.exists(msl_file)==False): #USE OLD FORMAT WITH MISSING UNDERSCORE
                    msl_file = self.animal.recName.replace('rhd_files','msl_files/data').replace('.rhd','')+str(p).zfill(3)+'_'+str(q).zfill(3)+'_'+str(r).zfill(3)
                    
                
                with open(msl_file, 'rb') as infile:    data = pickle.load(infile)

                data = sparse.csr_matrix(data)
                count_matrix = data.toarray()

                ax.set_xticklabels([])
                ax.set_yticklabels([])
                #plt.imshow(data[1700:2300, 1700:2300], interpolation='none')

                #BINNING RESULTS
                binning = int(self.isi_binning.text())     #ms per bin   
                if binning == 1: 
                    ave_count_matrix = np.array(count_matrix)
                else:
                    ave_count_matrix = []
                    for k in range(0, len(count_matrix), binning):
                        temp_array = []
                        for m in range(0, len(count_matrix[k]), binning):
                            temp_array.extend([np.sum(count_matrix[k:k+binning,m:m+binning])])
                        ave_count_matrix.append(temp_array)
                    
                    ave_count_matrix=np.vstack((ave_count_matrix))

                size_window = int(self.zoom.text()) #ms of view screen
                count_window = size_window  / binning
                len_matrix = len(ave_count_matrix)
                
                plt.ylabel(str(p).zfill(3)+'_'+str(q).zfill(3)+'_'+str(r).zfill(3), fontsize=8)
                plt.imshow(ave_count_matrix[len_matrix/2 - count_window:len_matrix/2 + count_window
                                            ,len_matrix/2 - count_window:len_matrix/2 + count_window], 
                                            extent=[0, count_window*2, count_window*2, 0],
                                            interpolation=self.interpolation.text())

                ctr+=1
                if ctr>24:
                    #plt.show()
                    plt.savefig(msl_file.replace('data','')+".png")
                    ctr=0
                    #return p, q, r
        
    #plt.show()
    plt.savefig(msl_file.replace('data','')+".png")
    ctr=0
    plt.close()
    #return p, q, r           
        
def count_matrix(self):
    print "count matrix"         

def sta_movies(self):
    import matplotlib.animation as animation
    from math import sqrt

    animal = self.animal
    Sort = Ptcs(animal.ptcsName)
    print "... n_units: ", len(Sort.units)
    print "... start time: ", self.start_time.text(), "  end time: ", self.end_time.text()
    print "... start cell: ", self.start_cell.text(), "  end cell: ", self.end_cell.text()

    select_units = np.arange(min(int(self.start_cell.text()),len(Sort.units)),min(int(self.end_cell.text()),len(Sort.units)),1)

    #Set parameters
    start_time = float(self.start_time.text())  #time before t=0 (secs)
    end_time = float(self.end_time.text())
    img_rate = self.animal.img_rate
    block_save = int(self.block_save.text())

    #**** LOAD GENERIC MASK
    generic_mask_file = animal.home_dir+animal.name + '/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        fname = glob.glob(animal.filenames[0].replace('rhd_files/','tif_files/')[:animal.filenames[0].find('rhd_files/')+10]+"*std*")[0]
        images_temp = np.load(fname)
        #Define_generic_mask(images_temp, animal.home_dir)
        Define_generic_mask_single_frame(images_temp, animal)
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        
    generic_mask_indexes=np.zeros((128,128))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    print "... # units loading: ", select_units

    min_spikes = 0

    vid_array = []
    spikes_array = []
    for k in select_units:
        if len(Sort.units[k])<min_spikes: continue
        #if k !=4: continue
        print "Loading saved .npy files: ", k
        channel = Sort.maxchan[k]
        temp_name = glob.glob(animal.ptcsName.replace('tsf_files/','stm_files/img_avg_').replace('.ptcs','')+'_unit'+str(k).zfill(3)+"*")[0]
        print temp_name

        criteria_spikes = int(temp_name[temp_name.find('spikes')-6:temp_name.find('spikes')-1])

        spikes_array.append(criteria_spikes)
        STM = np.load(temp_name)

        for v in range(len(STM)):
            STM[v]= ndimage.gaussian_filter(STM[v], sigma=2)
            
        #Apply mask
        n_pixels=128
        #temp_array = np.ma.array(np.zeros((len(STM),n_pixels,n_pixels),dtype=np.float32), mask=True)
        temp_array=[]
        for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
        #for i in range(len(STM)):
            temp_array.append(np.ma.array(STM[i], mask=generic_mask_indexes, fill_value = 0., hard_mask = True))
        
        vid_array.append(temp_array)

    #if os.path.exists(animal.ptcsName.replace('rhd_files','movie_files')+'_units_'+str(len(vid_array))+'.mp4'): return

    #**** INITIALIZE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure() # make figure

    print "... generating animation..." 

    #Fix size of concatenated animations
    extra_row = 0
    x = int(sqrt(len(vid_array)))
    if x*int(x*1.4)< len(vid_array): 
        x+=1
        if (x-1)*int(x*1.4)>=len(vid_array): extra_row=-1

    im=[]
    for k in range(len(vid_array)):
        ax = plt.subplot(2,len(vid_array),k+1)
        
        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        ax.yaxis.labelpad = 0
        #ax.set_ylabel("",fontsize=6)
        ax.set_title(str(spikes_array[k])+", "+str(round(np.ma.min(vid_array[k])*100,2))+"%.."+str(round(np.ma.max(vid_array[k])*100,2))+"%", fontsize=5)

        Ev_max = np.max(np.abs(vid_array[k]))#*.9
        #v_min = -v_max 
        print "Max/min DFF (%): ", Ev_max*100
        
        im.append([])
        im[k] = plt.imshow(vid_array[k][0], cmap=plt.get_cmap('jet'), vmin=-Ev_max, vmax=Ev_max, interpolation='none')#, vmin=0, vmax=v_max)
        #im[k] = plt.imshow(vid_array[k][0], cmap=blue_red1, vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)

    #function to update figure
    def updatefig(j):
        print j,
        plt.suptitle(animal.ptcsName.replace(animal.home_dir,'') +"\nFrame: "+str(j)+"  " +str(round((float(j)/img_rate+start_time),2))+"sec")

        # set the data in the axesimage object
        for k in range(len(vid_array)):
            im[k].set_array(vid_array[k][j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array[0])), interval=100, blit=False, repeat=True)

    if True:
    #if save_animation:
        ani.save(animal.ptcsName.replace('tsf_files','movie_files')+'_units_'+str(len(vid_array))+'.mp4', writer=writer)
        print " DONE! "

    #plt.show()



def dim_reduction(mouse, method):
    
    matrix_in = []
    for k in range(len(mouse.traces)):
        matrix_in.extend(mouse.traces[k])
    matrix_in = np.array(matrix_in)
    print matrix_in.shape
    
    methods = ['MDS - SMACOF', 't-SNE', 'PCA', 'Sammon']
    
    print "Computing dim reduction, size of array: ", np.array(matrix_in).shape
    
    if method==0:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        if os.path.exists(mouse.home_dir+mouse.name+'/MDS.npy')==False:
            print "... MDS-SMACOF..."
            print "... pairwise dist ..."
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing MDS ..."
            mds_clf = manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
            results = mds_clf.fit(adist)
            Y = results.embedding_ 

            np.save(mouse.home_dir+mouse.name+'/MDS', Y)
        else:
            Y = np.load(mouse.home_dir+mouse.name+'/MDS.npy')
                
    elif method==1:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        if os.path.exists(mouse.home_dir+mouse.name+'/tSNE.npy')==False:
            print "... tSNE ..."
            print "... pairwise dist ..."
            
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
            
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing tSNE ..."
            model = manifold.TSNE(n_components=3, init='pca', random_state=0)
            Y = model.fit_transform(adist)
            #Y = model.fit(adist)
        
            np.save(mouse.home_dir+mouse.name+'/tSNE', Y)
        
        else:
            Y = np.load(mouse.home_dir+mouse.name+'/tSNE.npy')

    elif method==2:

        Y, X = PCA(matrix_in, 3)
        np.save(mouse.home_dir+mouse.name+'/PCA', Y)

        if False: 
            if os.path.exists(mouse.home_dir+mouse.name+'/PCA.npy')==False:
                print "...computing PCA..."
                Y, X = PCA(matrix_in, 3)

                np.save(mouse.home_dir+mouse.name+'/PCA', Y)
            else:
                Y = np.load(mouse.home_dir+mouse.name+'/PCA.npy')
            
                
    elif method==4:

        if os.path.exists(mouse.home_dir+mouse.name+'/tSNE_barnes_hut.npy')==False:
            print "... computing Barnes-Hut tSNE..."
            Y = bh_sne(np.array(matrix_in))
        
            np.save(mouse.home_dir+mouse.name+'/tSNE_barnes_hut', Y)
        else:
            Y = np.load(mouse.home_dir+mouse.name+'/tSNE_barnes_hut.npy')
    
    
    elif method==5: 
        print "NOT IMPLEMENTED"
        
        sammon(matrix_in)    

    return Y



def dim_reduction_stack(self):
    """ Input is 2D stack: (samples, dimension)
    """
    import sklearn

    
    matrix_in = self.stack
    method = self.selected_dim_red
    
    file_out = self.filtered_file[:-4]+'_'+method+'_'+self.starting_frame.text()+"start_"+self.number_frame.text()+"frames"
    

    methods = ['PCA', 'MDS', 'tSNE', 'tSNE_Barnes_Hut']
    
    print "Computing dim reduction, size of array: ", matrix_in.shape
    
    if method==methods[1]:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        if os.path.exists(file_out+'.npy')==False:
            print "... MDS-SMACOF..."
            print "... pairwise dist ..."
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing MDS ..."
            mds_clf = sklearn.manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
            results = mds_clf.fit(adist)
            Y = results.embedding_ 

            np.save(file_out, Y)
        else:
            Y = np.load(file_out+'.npy')
                
    elif method==methods[2]:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        if os.path.exists(file_out+'.npy')==False:
            print "... tSNE ..."
            print "... pairwise dist ..."
            
            dists = sklearn.metrics.pairwise.pairwise_distances(matrix_in)
            
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing tSNE ..."
            model = sklearn.manifold.TSNE(n_components=3, init='pca', random_state=0)
            Y = model.fit_transform(adist)
            #Y = model.fit(adist)
        
            np.save(file_out, Y)
        
        else:
            Y = np.load(file_out+'.npy')

    elif method==methods[0]:

        Y, X = PCA(matrix_in, 3)
        np.save(file_out+'.npy', Y)

        if False: 
            if os.path.exists(mouse.home_dir+mouse.name+'/PCA.npy')==False:
                print "...computing PCA..."
                Y, X = PCA(matrix_in, 3)

                np.save(file_out, Y)
            else:
                Y = np.load(file_out+'.npy')
            
                
    elif method==methods[3]:

        if os.path.exists(file_out+'.npy')==False:
            print "... computing Barnes-Hut tSNE..."
            Y = bh_sne(np.array(matrix_in))
        
            np.save(file_out, Y)
        else:
            Y = np.load(file_out+'.npy')

    print "...DONE..."
    
    self.dim_reduction_out = Y
    
    
def plot_3D_distribution_new(self):

    #main_widget = self.parent
    
    #Load saved dim_red data
    temp_file = self.parent.root_dir+self.selected_animal+"/tif_files/"+self.selected_recording+'.npy'
    filename = temp_file[:-4]+'_'+self.selected_filter+'_'+self.parent.filter_low.text()+'hz_'+self.parent.filter_high.text()+'hz_'+ self.selected_dim_red+'_'+self.starting_frame.text()+"start_"+self.number_frame.text()+"frames.npy"

    print filename

    data = np.load(filename)
    data=data*float(self.scaling_factor.text())
    print data.shape
    
    #Start window
    self.parent.glwindow = ClusterWindow(pos=(0, 0), size=(500, 500))
    self.parent.glwindow.show()
    self.parent.glwindow.glWidget.points_selected=[]    #Empty list to gather all selected points.
    
    #Load data as points
    X = data
    X = np.float32(X)
    sids = np.arange(len(data))
    nids = sids % len(CLUSTERCOLOURSRGB)

    #Load data as triangles
    points, colours = points_to_triangles(data)
    self.parent.glwindow.glWidget.points_pyramids = points
    self.parent.glwindow.glWidget.colours_pyramids = colours   

    #Load data as lines
    points, colours = points_to_lines(data)
    self.parent.glwindow.glWidget.points_lines = points
    self.parent.glwindow.glWidget.colours_lines = colours   
    self.parent.glwindow.glWidget.points_selected=[]
    
    self.parent.glwindow.plot(X, sids, nids)

def points_to_lines(data):
    """ Make lines connecting each point in the data """

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.Normalize(vmin=0, vmax=len(data)) #Set maximum for color scheme to length of data array;
    cmap = cm.jet
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    
    points=[]
    colors=[]
    for ctr in range(len(data[:-1])):
        colors.extend([[x*255 for x in m.to_rgba(ctr)[:3]]]*2) # uint8   
        
        points.append(data[ctr])
        points.append(data[ctr+1])

    points = np.vstack([ points ])
    colours = np.uint8(np.vstack([ colors ]))
    
    return points, colours
    
def points_to_triangles(data):
    
    from math import sqrt

    import matplotlib as mpl
    import matplotlib.cm as cm

    norm = mpl.colors.Normalize(vmin=0, vmax=len(data))
    cmap = cm.jet
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    

    points=[]
    colors=[]
    for ctr, point in enumerate(data):
        soma_xyz=point

        size=3.  #Size of cell soma; CAN IMPLEMENT THIS TO BE CELL SPECIFIC EVENTUALLY
        
        colors.extend([[x*255 for x in m.to_rgba(ctr)[:3]]]*12) # uint8   
        #print [int(x*255) for x in m.to_rgba(ctr)[:3]]
        #return
        
        #Start making tetraheadrons;
        #NB: Need to offset the start point to centre of tetrahaedron which is sqrt(1/6) x size higher
        centre_offset = sqrt(1./6.)*size
        
        floats=[]
        floats.append(soma_xyz[0])
        floats.append(soma_xyz[1])
        floats.append(soma_xyz[2])

        #************************************************************************************
        # Coordinate #1; 1st coordinate = top of tetraheadron; 1st side triangle
        points.append(floats[:3])
        floats[0]=floats[0]+size/2          
        floats[1]=floats[1]-sqrt(size**2-sqrt(size**2-(size/2)**2)) 
        floats[2]=floats[2]+sqrt(size**2-(size/2)**2)/2
        points.append(floats[:3])
        floats[0]=floats[0]-size
        floats[1]=floats[1]
        floats[2]=floats[2]
        points.append(floats[:3])

        # Coordinate #2; Reset location first; 1st crd = top of tetraheadron; 2nd side 
        floats[0]=floats[0]+size/2
        floats[1]=floats[1]+sqrt(size**2-sqrt(size**2-(size/2)**2))
        floats[2]=floats[2]-sqrt(size**2-(size/2)**2)/2
        points.append(floats[:3])
        floats[0]=floats[0]+size/2
        floats[1]=floats[1]-sqrt(size**2-sqrt(size**2-(size/2)**2))
        floats[2]=floats[2]+sqrt(size**2-(size/2)**2)/2
        points.append(floats[:3])
        floats[0]=floats[0]-size/2
        floats[1]=floats[1]
        floats[2]=floats[2]-sqrt(size**2-(size/2)**2)
        points.append(floats[:3])

        # Coordinate #3; Reset location first; 1st coord = top of tetraheadron; 3rd side
        floats[0]=floats[0]
        floats[1]=floats[1]+sqrt(size**2-sqrt(size**2-(size/2)**2))
        floats[2]=floats[2]+sqrt(size**2-(size/2)**2)/2
        points.append(floats[:3])
        floats[0]=floats[0]
        floats[1]=floats[1]-sqrt(size**2-sqrt(size**2-(size/2)**2))
        floats[2]=floats[2]-sqrt(size**2-(size/2)**2)/2
        points.append(floats[:3])
        floats[0]=floats[0]-size/2
        floats[1]=floats[1]
        floats[2]=floats[2]+sqrt(size**2-(size/2)**2)
        points.append(floats[:3])

        # Coordinate #4; Use last location; bottom triangle
        points.append(floats[:3])
        floats[0]=floats[0]+size/2
        floats[1]=floats[1]
        floats[2]=floats[2]-sqrt(size**2-(size/2)**2)
        points.append(floats[:3])
        floats[0]=floats[0]+size/2
        floats[1]=floats[1]
        floats[2]=floats[2]+sqrt(size**2-(size/2)**2)
        points.append(floats[:3])

    points = np.vstack([ points ])
    colours = np.uint8(np.vstack([ colors ]))
    
    return points, colours
    
def Define_generic_mask_single_frame(images_processed, animal):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed
    
    fig, ax = plt.subplots()

    if (os.path.exists(animal.home_dir+animal.name + '/genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed)#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click_single_frame)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed)):
            for j in range(len(images_processed)):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed)):
            for j in range(len(images_processed)):
                if mask[counter] == False:
                    images_processed[i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed)
        plt.show()
       
        genericmask_file = animal.home_dir+animal.name + '/genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    #else:
    #    print "Loading saved general mask"
        
        
    #if (os.path.exists(main_dir + 'bregmamask.txt')==False):
        #bregma_coords = []
        #print "Making Bregma mask"
        #ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        #ax.set_title("Compute bregma mask")
        ##figManager = plt.get_current_fig_manager()
        ##figManager.window.showMaximized()
        #cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        #plt.show()

       
        #bregmamask_file = main_dir + 'bregmamask.txt'
        #np.savetxt(bregmamask_file, bregma_coords)

        #print "Finished Bregma Mask"

    #else:
    #    print "Loading saved bregma mask"
        
    #return generic_coords
#*******************************




def Define_generic_mask(images_processed, main_dir):

    global coords, images_temp, ax, fig, cid
    
    images_temp = images_processed
    
    fig, ax = plt.subplots()

    if (os.path.exists(main_dir + 'genericmask.txt')==False):
        coords=[]

        ax.imshow(images_processed[len(images_processed)/2])#, vmin=0.0, vmax=0.02)
        ax.set_title("Compute generic (outside the brain) mask")
        #figManager = plt.get_current_fig_manager()
        #figManager.window.showMaximized()
        cid = fig.canvas.mpl_connect('button_press_event', on_click)
        plt.show()

        #******* MASK AND DISPLAY AREAS OUTSIDE GENERAL MASK 
        #Search points outside and black them out:
        all_points = []
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                all_points.append([i,j])

        all_points = np.array(all_points)
        vertixes = np.array(coords) 
        vertixes_path = Path(vertixes)
        
        mask = vertixes_path.contains_points(all_points)
        counter=0
        coords_save=[]
        for i in range(len(images_processed[0][0])):
            for j in range(len(images_processed[0][0])):
                if mask[counter] == False:
                    images_processed[len(images_processed)/2][i][j]=0
                    coords_save.append([i,j])
                counter+=1

        fig, ax = plt.subplots()
        ax.imshow(images_processed[len(images_processed)/2])
        plt.show()
       
        genericmask_file = main_dir + '/genericmask.txt'
        np.savetxt(genericmask_file, coords_save)

        print "Finished Making General Mask"

    #else:
    #    print "Loading saved general mask"
        
        
    #if (os.path.exists(main_dir + 'bregmamask.txt')==False):
        #bregma_coords = []
        #print "Making Bregma mask"
        #ax.imshow(images_processed[100])#, vmin=0.0, vmax=0.02)
        #ax.set_title("Compute bregma mask")
        ##figManager = plt.get_current_fig_manager()
        ##figManager.window.showMaximized()
        #cid = fig.canvas.mpl_connect('button_press_event', remove_bregma)
        #plt.show()

       
        #bregmamask_file = main_dir + 'bregmamask.txt'
        #np.savetxt(bregmamask_file, bregma_coords)

        #print "Finished Bregma Mask"

    #else:
    #    print "Loading saved bregma mask"
        
    #return generic_coords
#*******************************
def make_vids(data, file_):
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=50, metadata=dict(artist='Me'), bitrate=6400)

    fig = plt.figure() # make figure

    im = []
    for k in range(len(data)):
        ax = plt.subplot(1,1,k+1)

        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        ax.yaxis.labelpad = 0
        Ev_max = np.max(np.abs(data[k]))*.9
        #v_min = -v_max 
        
        im.append([])
        im[k] = plt.imshow(data[k][0], cmap=plt.get_cmap('gray'), vmin = -Ev_max, vmax = Ev_max, interpolation='none')#, vmin=0, vmax=v_max)

    #function to update figure
    def updatefig(j):
        print "... frame: ", j
        plt.suptitle("Frame: "+str(j)+"\n"+str(round(float(j)/150,2))+"sec")

        # set the data in the axesimage object
        for k in range(len(data)):
            im[k].set_array(data[k][j])

        # return the artists set
        return im
        
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(data[0])), interval=10, blit=False, repeat=True)
    ani.save(file_+'.mp4', writer=writer)
    plt.show()


#def on_click_single_frame_light(event):
    
    
def on_click_single_frame(event):
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp)
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        ax.imshow(images_temp)
        fig.canvas.draw()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)
        


def on_click(event):
    
    global coords, images_temp, ax, fig, cid
    
    n_pix = len(images_temp[0])
    
    if event.inaxes is not None:
        coords.append((event.ydata, event.xdata))
        for j in range(len(coords)):
            for k in range(3):
                for l in range(3):
                    images_temp[len(images_temp)/2][min(n_pix,int(coords[j][0])-1+k)][min(n_pix,int(coords[j][1])-1+l)]=0

        ax.imshow(images_temp[len(images_temp)/2])
        #plt.show()
        fig.canvas.draw()
                    #figManager = plt.get_current_fig_manager()
                    #figManager.window.showMaximized()
    else:
        print 'Exiting'
        plt.close()
        fig.canvas.mpl_disconnect(cid)



      
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y





def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

def find_previous(array,value):
    temp = (np.abs(array-value)).argmin()
    if array[temp]>value: return temp-1
    else: return temp
    


def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T #THIS IS REDUNDANT... REDUCE IT
