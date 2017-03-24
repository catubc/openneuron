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
from scipy import stats
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from sklearn import decomposition
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

import multiprocessing

#Smooth and convolve original data to look for flow:
import cv2
import scipy
import scipy.ndimage

from load_intan_rhd_format import *

sys.path.append('/home/cat/code/')
import TSF.TSF as TSF
import PTCS.PTCS as PTCS


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


#class Tsf_file(object):

    #def __init__(self, file_name):
        
        #self.read_header(file_name)
        
    #def read_header(self, file_name):
        
        #self.fin = open(file_name, "rb")
        
        #self.header = self.fin.read(16)
        #self.iformat = struct.unpack('i',self.fin.read(4))[0] 
        #self.SampleFrequency = struct.unpack('i',self.fin.read(4))[0] 
        #self.n_electrodes = struct.unpack('i',self.fin.read(4))[0] 
        #self.n_vd_samples = struct.unpack('i',self.fin.read(4))[0] 
        #self.vscale_HP = struct.unpack('f',self.fin.read(4))[0] 

        #if self.iformat==1001:
            #self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            #self.Siteloc = struct.unpack(str(2*self.n_electrodes)+'h', self.fin.read(2*self.n_electrodes*2))
        #if self.iformat==1002:
            #self.Siteloc = np.zeros((2*self.n_electrodes), dtype=np.int16)
            #self.Readloc = np.zeros((self.n_electrodes), dtype=np.int32)
            #for i in range(self.n_electrodes):
                #self.Siteloc[i*2] = struct.unpack('h', self.fin.read(2))[0]
                #self.Siteloc[i*2+1] = struct.unpack('h', self.fin.read(2))[0]
                #self.Readloc[i] = struct.unpack('i', self.fin.read(4))[0]

    #def read_ec_traces(self):
        #print " ... reading data, #chs: ", self.n_electrodes, " nsamples: ", self.n_vd_samples, " len: ", float(self.n_vd_samples)/float(self.SampleFrequency), " sec."
        #self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        #self.ec_traces.shape = self.n_electrodes, self.n_vd_samples

        #self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 
        
        ##print "No. ground truth cell spikes: ", self.n_cell_spikes
        #if (self.n_cell_spikes>0):
            #if (self.iformat==1001):
                #self.vertical_site_spacing = struct.unpack('i',self.fin.read(4))[0] 
                #self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 

            #self.fake_spike_times =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            #self.fake_spike_assignment =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            #self.fake_spike_channels =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
        
        #self.fin.close()

    #def read_trace(self, channel):
        ##Load single channel 

        #indent = 16+20+self.n_electrodes*8

        #self.fin.seek(indent+channel*2*self.n_vd_samples, os.SEEK_SET)         #Not 100% sure this indent is correct.
        #self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_vd_samples)
        #self.fin.close()
    
    #def save_tsf(self, file_name):
        
        #fout = open(file_name, 'wb')
        #print "...saving: ",  file_name
        #fout.write(self.header)
        #fout.write(struct.pack('i', self.iformat))
        #fout.write(struct.pack('i', self.SampleFrequency))
        #fout.write(struct.pack('i', self.n_electrodes))
        #fout.write(struct.pack('i', self.n_vd_samples))
        #fout.write(struct.pack('f', self.vscale_HP))
        #for i in range (self.n_electrodes):
            #fout.write(struct.pack('h', self.Siteloc[i*2]))
            #fout.write(struct.pack('h', self.Siteloc[i*2+1]))
            #fout.write(struct.pack('i', i+1))                 #CAREFUL, SOME FILES MAY USE ReadLoc values..

        #self.ec_traces.tofile(fout)

        #fout.write(struct.pack('i', self.n_cell_spikes))

        ##try:
            ##self.subsample
        ##except NameError:
            ##self.subsample = 1.0

        ##fout.write(struct.pack('i', self.subsample))

        #fout.close()


    #def read_footer(self):
        
        ##Header indent
        #indent = 16+20+self.n_electrodes*8

        ##Voltage indent
        #self.fin.seek(indent+self.n_electrodes*2*self.n_vd_samples, os.SEEK_SET)         #Not 100% sure this indent is correct.
        
        ##Read last record
        #self.n_cell_spikes =  np.fromfile(self.fin, dtype=np.int32, count=1)
        #print "... n spikes: ", self.n_cell_spikes
        

        ##print "No. ground truth cell spikes: ", self.n_cell_spikes
        #if (self.n_cell_spikes>0):
            #self.fake_spike_times =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            #self.fake_spike_assignment =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
            #self.fake_spike_channels =  np.fromfile(self.fin, dtype=np.int32, count=self.n_cell_spikes)
        
        ##Save meta data into files.
        
        ##CHECK TO SEE IF 
        #self.n_files =  np.fromfile(self.fin, dtype=np.int32, count=1)
        #if len(self.n_files)==0:
            #print "... reached end of file ... older version of .tsf does not contain original file names or original # of samples.."
            #return

        #self.file_names = []
        #self.n_samples = []
        #self.n_digital_chs = []
        #self.digital_channels = []
        
        #print "... n files: ", self.n_files                                 #THIS IS FOR GENERAL CASE where > 1 .tsf file saved.
        #for k in range(len(self.n_files)):
            #self.file_names.append(self.fin.read(256))
            #print "... original file name: ", self.file_names[k]

            #self.n_samples.append(np.fromfile(self.fin, dtype=np.int32, count=1))
            #print "... original n_samples: ", self.n_samples[k]
            
            ##Load digital channels
            #self.n_digital_chs.append(np.fromfile(self.fin, dtype=np.int32, count=1))

            #if len(self.n_digital_chs[0])==0:
                #print "... reached end of file ... older version of .tsf does not contain digital channel information..."
                #return

            #print "... # of digital chs: ", self.n_digital_chs
            
            #temp_chs = []
            #for ch in range(self.n_digital_chs[k]):
                #temp_chs.append(np.fromfile(self.fin, dtype=bool, count=self.n_samples[k]))         #Load the original #samples saved to disk, NOT n_vd_samples

            #self.digital_channels.append(temp_chs)


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
        self.Siteloc = np.zeros((self.n_electrodes*2), dtype=np.int16) #Read as 1D array
        for i in range (self.n_electrodes):
            self.Siteloc[i*2]=30*(i%2)
            self.Siteloc[i*2+1]=i*23


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

    fout.write(tsf.header)
    fout.write(struct.pack('i', tsf.iformat))
    fout.write(struct.pack('i', tsf.SampleFrequency))
    fout.write(struct.pack('i', tsf.n_electrodes))
    fout.write(struct.pack('i', tsf.n_vd_samples))
    fout.write(struct.pack('f', tsf.vscale_HP))
    
    for i in range(tsf.n_electrodes):
        fout.write(struct.pack('h', tsf.Siteloc[i*2]))
        fout.write(struct.pack('h', tsf.Siteloc[i*2+1]))
        fout.write(struct.pack('i', i+1))

    for i in range(tsf.n_electrodes):
        np.int16(tsf.ec_traces[tsf.layout[i]]).tofile(fout)  #Save ec_traces for each channel while converting to int16; the channel order comes from probe
        #tsf.ec_traces[i].tofile(fout)  #Save ec_traces for each channel while converting to int16; the channel order comes from probe

    fout.write(struct.pack('i', tsf.n_cell_spikes))

    #Save additional spike information if value non-zero
    if tsf.n_cell_spikes!=0:
        np.int32(tsf.fake_spike_times).tofile(fout)
        np.int32(tsf.fake_spike_assignment).tofile(fout)
        np.int32(tsf.fake_spike_channels).tofile(fout)

    #Footer information
    #number of files saved in tsf; 
    n_files = len(tsf.file_names)       #Save # of files first
    np.int32(n_files).tofile(fout)
    
    #Save: name of each file; number of digital channels; digital channels in boolean format
    for fname, n_samples, n_dig_chs, dig_chs in zip(tsf.file_names, tsf.n_samples, tsf.n_digital_chs, tsf.digital_chs):
        fname = fname+" "*(256-len(fname))
        print "...saving file_name: ", fname
        fout.write(fname)
        fout.write(struct.pack('i', n_samples))

        #Save digital channels
        print "... # of digital chs: ", n_dig_chs
        np.int32(n_dig_chs).tofile(fout)        #Save # of digital channels
        print "... number of actual saved data streams: ", len(dig_chs)
        for ch in range(len(dig_chs)):
            print "...xaving ch: ", ch
            dig_chs[ch].tofile(fout)
    
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
    #length = vid.get(cv2.cv.CV_CAP_PROP_FRAME_COUNT)
    #width  = vid.get(cv2.cv.CV_CAP_PROP_FRAME_WIDTH)
    #height = vid.get(cv2.cv.CV_CAP_PROP_FRAME_HEIGHT)
    #fps    = vid.get(cv2.cv.CV_CAP_PROP_FPS)

    #NB: CAN ALSO INDEX DIRECTLY INTO .WMV FILES:
    #time_length = 30.0
    #fps=25
    #frame_seq = 749
    #frame_no = (frame_seq /(time_length*fps))

    ##The first argument of cap.set(), number 2 defines that parameter for setting the frame selection.
    ##Number 2 defines flag CV_CAP_PROP_POS_FRAMES which is a 0-based index of the frame to be decoded/captured next.
    ##The second argument defines the frame number in range 0.0-1.0
    #cap.set(2,frame_no);

   #print length, width, height, fps
    #if length==0: print "... no movie file... returning..."; return
        
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

def behavioural_stack(self):
    
    #Load behaviour camera and annotations data - if avialable
    if os.path.exists(self.vid_rate_filename):

        #********** Load Annotations ********
        print "... making stacks of annotation arrays ...",
        areas = ['_lick', '_lever'] #, '_pawlever', '_lick', '_snout', '_rightpaw', '_leftpaw', '_grooming'] 
        annotation_arrays = []
        
        for ctr, area in enumerate(areas):
            annotation_arrays.append([])
            for k in range(len(self.movie_data)): annotation_arrays[ctr].append([])
            
            #data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+area+'_clusters.npz')
            data = np.load(self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+area+'_clusters.npz')
            cluster_indexes=data['cluster_indexes'] 
            cluster_names=data['cluster_names']
            
            for k in range(len(cluster_indexes)):
                for p in range(len(cluster_indexes[k])):
                    annotation_arrays[ctr][cluster_indexes[k][p]] = cluster_names[k]

            annotation_arrays[ctr] = np.array(annotation_arrays[ctr])[self.movie_indexes]

        self.annotation_arrays=[]
        for k in range(len(annotation_arrays)):                  #***************************** SAME DUPLICATION AS ABOVE; Video is 15Hz, imaging is 30Hz
            self.annotation_arrays.append([])
            for p in range(len(annotation_arrays[k])):
                self.annotation_arrays[k].append(annotation_arrays[k][p])
                self.annotation_arrays[k].append(annotation_arrays[k][p])

        self.annotation_stacks = []
        print "...making annotation_stacks..."
        areas = ['_Tongue', '_Lever'] #, '_pawlever', '_lick', '_snout', '_rightpaw', '_leftpaw', '_grooming'] 
        for k in range(len(self.annotation_arrays[0])):
            #print k
            fig = Figure()
            canvas = FigureCanvas(fig)
            ax = fig.gca()

            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            #ax.tick_params(axis='both', which='both', labelsize=30)
            #plt.title(self.annotation_arrays[0][k]+'_'+self.annotation_arrays[1][k], fontsize=40)
            for p in range(len(self.annotation_arrays)):
                ax.text(-10, p*20+30, areas[p][1:]+':  '+self.annotation_arrays[p][k], fontsize=40, fontweight='bold')
                #print areas[p][1:]+':  '+self.annotation_arrays[p][k]
            
            ax.set_ylim(0,len(self.annotation_arrays)*20+20)
            ax.set_xlim(0,120)
            ax.axis('off')

            #ax.plot(x_val[:k],y_val[:k], linewidth=3)
            #ax.set_ylim(0,120)
            #ax.set_xlim(0,len(self.movie_stack))
            canvas.draw()       # draw the canvas, cache the renderer
            
            data = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
            data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            self.annotation_stacks.append(data)
        
        print "...behavioural stack done..."
        
        plt.close()
    
    else:
        print "... behavioural video data doesn't exist ... "
        
        self.annotation_stacks = []

        self.movie_stack = np.zeros((len(self.ca_stack[0]), 30, 40), dtype=np.int8)

    filesave = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+ \
                self.selected_session+'_'+self.selected_code+'_'+self.selected_trial+'_annotation_stack'
               
    np.save(filesave, self.annotation_stacks)
    
    #conn.send(self.annotation_stacks)
    #conn.close()

def calcium_stack(self):
    
    
    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_stm.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_stm.npy"
    print "...stm_name: ", self.traces_filename
    
    data = np.load(self.traces_filename, mmap_mode='c')[int(self.selected_trial)]  #Load only selected trial
    print data.shape


    self.data_norm = []
    self.ca_stack = []
    n_stacks = 2
    for k in range(n_stacks): 
        self.ca_stack.append([])
    
    self.ca_stack[0] = data
    
    #Normalize data; opencv functions largely work on unit8 data types
    
    #img_rate = 30.0
    #start_time = float(self.stm_start_time.text()); end_time = float(self.stm_end_time.text())
    #for k in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), 1):
    for k in range(len(data)):
        data_norm = ((data[k]-np.min(data[k]))/(np.max(data[k])-np.min(data[k]))*255).astype(np.uint8)      #Normalize data to gray scale 0..255
        #data_norm = (data[k]-np.min(data[k]))/(np.max(data[k])-np.min(data[k]))*2. - 1.   #Normalize to: -1.. 0 .. +1 scale
        self.data_norm.append(data_norm)
    
    self.data_norm = np.array(self.data_norm)
    self.data_norm_ave = np.average(self.data_norm, axis = 0)
    
    
    filter_power = float(self.mask_power.text())
    
    for k in range(len(self.data_norm)): 
        stack_ctr = 1
        
        sigma=2
        img_arctan = np.ma.arctanh((self.data_norm[k]-128.)/128.)
        
        data_neg = -1. * np.ma.clip(img_arctan, -100., 0)
        data_neg = -np.ma.power(data_neg, filter_power)
        data_pos = np.ma.clip(img_arctan, 0, 100.)
        data_pos = np.ma.power(data_pos, filter_power)
        data_total = np.float32(data_neg + data_pos)
        
        #img_uint8 = ((data_total-np.min(data_total))/(np.max(data_total)-np.min(data_total))*255.).astype(np.uint8)
        img_gaussian = ndimage.gaussian_filter(data_total, sigma=sigma)    #Recentre image around zero
        #negatives = np.clip(img_gaussian, -1E-6, 1E-6)/1E-6
        #data_out = np.power(np.ma.abs(img_gaussian), filter_power)*negatives

        self.ca_stack[stack_ctr].append(img_gaussian);  stack_ctr+=1


        #kernel_5pix = np.ones((5,5),np.float32)/5.
        #data_out = cv2.filter2D(self.data_norm[k],-1,kernel_5pix)
        #self.ca_stack[stack_ctr].append(data_out);  stack_ctr+=1

        #sigma = 10
        #data_out = ndimage.gaussian_filter(self.data_norm[k], sigma=sigma)
        #self.ca_stack[stack_ctr].append(data_out);  stack_ctr+=1
        
        
        #sx = scipy.ndimage.sobel(data_out.astype(np.uint8), axis=0, mode='nearest')
        #sy = scipy.ndimage.sobel(data_out.astype(np.uint8), axis=1, mode='nearest')
        #data_out = np.hypot(sx, sy)
        #self.ca_stack[stack_ctr].append(data_out);  stack_ctr+=1
    
    
    #Remove averages from the smoothed stacks;
    #self.ca_stack[2] = self.ca_stack[2] - np.average(self.ca_stack[2], axis=0)
    #self.ca_stack[2] = self.ca_stack[2] - np.average(self.ca_stack[2], axis=0)
    
   
    #np.save(self.traces_filename[:-4] + "_Ca_stacks" , self.ca_stack)   #SAVE ARRAYS BEFORE MASKING
    temp_stack = self.ca_stack   #Save for loading below
    
    #Apply Generic Mask
    print "...selected trial for stm: ", self.selected_trial
    for k in range(len(self.ca_stack)):
        self.ca_stack[k] = quick_mask(self, self.ca_stack[k])
    

    #Last and apply Motion Mask
    filename_motion_mask = self.traces_filename.replace('_traces.npy','')[:-4]+'_motion_mask.npy'
    print filename_motion_mask
    if os.path.exists(filename_motion_mask)==False:
        print "...motion mask missing..."
        return
    else: 
        motion_mask = np.load(filename_motion_mask)
        for k in range(len(self.ca_stack)):
            self.ca_stack[k] = self.ca_stack[k] * motion_mask
    
    self.ca_stack = np.ma.array(self.ca_stack)
    print "...self.ca_stack.shape: ", self.ca_stack.shape


    #************************************************************************************************************
    #**************************************** PCA SPACE CALCIUM IMAGING PANEL ***********************************
    #************************************************************************************************************
    
    data = temp_stack[0]
    subsampled_array = []
    for k in range(len(data)):
        subsampled_array.append(scipy.misc.imresize(data[k], .9, interp='bilinear', mode=None))

    methods = ['MDS', 'tSNE', 'PCA', 'Sammon']
    method = methods[2]
    print "... computing original dim reduction ..."

    X = []
    for k in range(len(subsampled_array)):
        X.append(np.ravel(subsampled_array[k]))

    X = PCA_reduction(X, n_components=3)

    cm = plt.cm.get_cmap('jet')
    colors = range(len(X))
    
    self.pca_stack = []
    print "... making pca_stack..."
    for k in range(len(data)):
        #print k
        fig = Figure()
        canvas = FigureCanvas(fig)
        #ax = fig.gca()
        fig.set_size_inches(10, 10)

        #fig = plt.figure(1, figsize=(4, 3))
        plt.clf()
        #ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
        ax = Axes3D(fig, rect=[0, 0, 1, 1], elev=48, azim=134)

        #ax.scatter(X[:, 0], X[:, 1], X[:, 2], c = colors, cmap=plt.cm.spectral)
        ax.scatter(X[:k+1, 0], X[:k+1, 1], X[:k+1, 2], s =200, c = colors[:k+1], cmap=cm)
        if k>1: ax.plot3D (X[:k, 0], X[:k, 1], X[:k,2])

        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        
        ax.set_xlim(np.min(X[:, 0])-1, np.max(X[:, 0])+1)
        ax.set_ylim(np.min(X[:, 1])-1, np.max(X[:, 1])+1)
        ax.set_zlim(np.min(X[:, 2])-1, np.max(X[:, 2])+1)


        canvas.draw()       # draw the canvas, cache the renderer
        
        data = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        self.pca_stack.append(data)
    
    print "...done pca_stack..."
    
    plt.close()

    filesave = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+ \
               self.selected_session+'_'+self.selected_code+'_'+self.selected_trial
               
    np.save(filesave+'_ca_stack', np.ma.filled(self.ca_stack, np.nan))
    
    #self.ca_stack.dump(filesave+'_ca_stack')
    
    np.save(filesave+'_pca_stack', self.pca_stack)
    #self.ca_stack.dump(filesave+'_ca_stack')

    
def lever_position_stack(self):

    lever_file = self.temp_file+'_abspositions.npy'
    times_file = self.temp_file+'_abstimes.npy'
    
    lever_data = np.load(lever_file)
    times_data = np.load(times_file)
    
    start_index = find_nearest(times_data, self.selected_locs_44threshold)
    #print start_index, times_data[start_index]
    
    #Initialize lever_stack and find indexes in lever_data @~120Hz that match 30Hz sampling rate; both img_rates should be saved to disk so can use exact vals
    #self.lever_stack = np.zeros((self.movie_stack.shape[0], 120, self.movie_stack.shape[0]), dtype=np.float32)+255
    self.lever_stack = []
    
    #Make plots and convert to img stack
    x_val = []
    y_val = []
    for k in range(len(self.movie_stack)):
        x_val.append(k)
        #y_val.append(lever_data[int(start_index+k*4-self.movie_stack.shape[0]/2*4)])
        y_val.append(abs(lever_data[int(start_index+k*(120./self.img_rate)-self.movie_stack.shape[0]/2*(120./self.img_rate))]))

    print y_val

    print "... making stacks of lever pull panels..."
    for k in range(len(self.movie_stack)):
        #print k
        fig = Figure()
        canvas = FigureCanvas(fig)
        ax = fig.gca()
        
        ax.tick_params(axis='both', which='both', labelsize=45)
        old_xlabel = np.linspace(0,len(self.movie_stack), 2*int(self.n_sec_window.text()))
        new_xlabel = np.around(np.linspace(-int(self.n_sec_window.text()), int(self.n_sec_window.text()), 2*int(self.n_sec_window.text())), decimals=2)

        ax.set_xticks(old_xlabel)
        ax.set_xticklabels(new_xlabel)

        ax.plot([0,len(self.movie_stack)], [10, 10], color = 'black', linewidth=2, alpha = 0.8)
        ax.plot([0,len(self.movie_stack)], [34, 34], 'r--', color = 'blue', linewidth=3, alpha = 0.8)
        ax.plot([0,len(self.movie_stack)], [60, 60], color = 'blue', linewidth=2, alpha = 0.8)
        
        ax.plot(x_val[:k],y_val[:k], linewidth=6)
        ax.set_ylim(0,120)
        ax.set_xlim(0,len(self.movie_stack))
        
        canvas.draw()       # draw the canvas, cache the renderer
        
        data = np.fromstring(canvas.tostring_rgb(), dtype='uint8')
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        self.lever_stack.append(data)
    
    print "...done lever stack..."
    
    plt.close()
    
    filesave = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+ \
               self.selected_session+'_'+self.selected_code+'_'+self.selected_trial+'_lever_stack'
               
    np.save(filesave, self.lever_stack)
    
    
def event_triggered_movies_single_Ca(self):
    """ Load [Ca] imaging and behavioural camera data and align to selected trial"""

    self.parent.n_sec = float(self.n_sec_window.text())
    self.start_time = -self.parent.n_sec; self.end_time = self.parent.n_sec
    self.temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session    
    self.abstimes = np.load(self.temp_file+'_abstimes.npy')

    self.img_rate = np.load(self.temp_file+'_img_rate.npy') #imaging rate

    #Process reward triggered data
    if (self.selected_code =='02') or (self.selected_code =='04') or (self.selected_code =='07'):
        self.locs_44threshold = np.load(self.temp_file+'_locs44threshold.npy')
        self.code_44threshold = np.load(self.temp_file+'_code44threshold.npy')
        
        indexes = np.where(self.code_44threshold==self.selected_code)[0]
        print "...indexes: "; print indexes

        self.code_44threshold = self.code_44threshold[indexes]  #Select only indexes that match the code selected
        self.locs_44threshold = self.locs_44threshold[indexes]

    #Process behaviour triggered data;
    else:
        load_behavioural_annotation_data(self)
        
        print len(self.code_44threshold)
        print len(self.locs_44threshold)
    
    
    print "...self.selected_code: ", self.selected_code
    print "...self.selected_trial: ", self.selected_trial
    

    self.selected_locs_44threshold = self.locs_44threshold[int(self.selected_trial)]        #selected_locs should already have been selected above
    self.selected_code_44threshold = self.code_44threshold[int(self.selected_trial)]
    
    print self.selected_locs_44threshold
    print self.selected_code_44threshold

    #****************************** Load behaviour video ******************************
    #Load behaviour camera and annotations data - if avialable
    self.vid_rate_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_vid_rate.npy'
    print "...loading behavioural camera data..."
    if os.path.exists(self.vid_rate_filename):
        
        self.vid_rate = np.loadtxt(self.vid_rate_filename)

        #Load original movie data and index only during blue_light_frames
        movie_data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy', mmap_mode='c')
        self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+self.selected_session+'_blue_light_frames.npy'
        self.movie_data = movie_data[np.load(self.blue_light_filename)]
        
        #Find movie frame corresponding to lever pull trigger
        movie_times = np.linspace(0, self.abstimes[-1], self.movie_data.shape[0])
        self.movie_04frame_locations = find_nearest(movie_times, self.selected_locs_44threshold)
        print "... frame event triggers: ", self.movie_04frame_locations

        #Make movie stack
        self.movie_indexes = np.arange(self.movie_04frame_locations+int(-self.parent.n_sec*self.vid_rate-1), self.movie_04frame_locations+int(self.parent.n_sec*self.vid_rate+1), 1)
        self.movie_stack = self.movie_data[self.movie_indexes]
        
        #self.movie_stack = self.movie_data[self.movie_04frame_locations+int(-self.parent.n_sec*self.vid_rate-1): self.movie_04frame_locations+int(self.parent.n_sec*self.vid_rate+1)]
        #print len(self.movie_stack)
        #quit()

        #Duplicate movie stack to match [Ca] imaging rate     #******************NB INTERPOLATION IS KIND OF HARDWIRED TO 30HZ & 15HZ.... PERHAPS THIS CAN SKIP FRAME SOMETIMES !?!
        new_stack = []
        for frame in range(len(self.movie_stack)):
            new_stack.append(self.movie_stack[frame])
            #new_stack.append((np.int16(self.movie_stack[frame])+np.int16(self.movie_stack[frame+1]))/2.)        #Interpolation, maybe not use it.
            new_stack.append(self.movie_stack[frame])
            
        #new_stack.append(self.movie_stack[-1]);  new_stack.append(self.movie_stack[-1])
        self.movie_stack = np.uint8(new_stack)
        
        print self.movie_stack.shape
        
    else:
        print "... behavioural video data doesn't exist ... "
        
        self.movie_stack = np.zeros((len(self.ca_stack[0]), 30, 40), dtype=np.int8)


    #*************************************************************************************************************
    #Process stacks in parallel - save data to disk
    procs=[]
    procs.append(multiprocessing.Process(target=behavioural_stack, args=(self,)))
    procs.append(multiprocessing.Process(target=calcium_stack, args=(self,)))
    procs.append(multiprocessing.Process(target=lever_position_stack, args=(self,)))

    map(lambda x: x.start(), procs)
    map(lambda x: x.join(), procs)
    
    
    #*******************************************************************************************************
    #Load processed data from disk
    filesave = self.parent.root_dir+self.parent.animal.name+"/tif_files/"+self.selected_session+'/'+ \
               self.selected_session+'_'+self.selected_code+'_'+self.selected_trial
               
    self.annotation_stacks = np.load(filesave+'_annotation_stack.npy')
    self.ca_stack = np.load(filesave+'_ca_stack.npy', allow_pickle=True)
    self.pca_stack = np.load(filesave+'_pca_stack.npy')
    self.lever_stack = np.load(filesave+'_lever_stack.npy')    
    
    print "... len pca_stack: ", len(self.pca_stack)
    print "... len ca_stack: ", len(self.ca_stack)
    
    #Make movies
    make_movies_ca(self)


def make_movies_ca(self):
    
    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=15000)
  
    fig = plt.figure()
    im = []
    
    #gs = gridspec.GridSpec(2,len(self.ca_stack)*2)
    gs = gridspec.GridSpec(4,6)
    
    #[Ca] stacks
    titles = ["GCamp6s Activity", "GCamp6s (z-Tranformed)"]
    for k in range(len(self.ca_stack)):    
        ax = plt.subplot(gs[0:2,k*2:k*2+2])
        plt.title(titles[k], fontsize = 12)

        v_max = np.nanmax(np.ma.abs(self.ca_stack[k])); v_min = -v_max
        ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
        im.append(plt.imshow(self.ca_stack[k][0], vmin=v_min, vmax = v_max, cmap=plt.get_cmap('jet'), interpolation='none'))

    #PCA stack
    ax = plt.subplot(gs[0:2,4:])
    plt.title("PCA State Space", fontsize = 12)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.pca_stack[0], cmap=plt.get_cmap('gray'), interpolation='none'))

    #Camera stack
    ax = plt.subplot(gs[2:4,0:4])
    plt.title("Behavioural Camera", fontsize = 12)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.movie_stack[0], cmap=plt.get_cmap('gray'), interpolation='none'))

    #Lever position trace
    ax = plt.subplot(gs[2:3,4:])
    plt.title("Lever Position", y = .95 , fontsize = 10)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.lever_stack[0], cmap=plt.get_cmap('gray'), interpolation='none'))
    
    #Annotation Stck
    ax = plt.subplot(gs[3:4,4:])
    plt.title("Annotations", y = .95 , fontsize = 10)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    im.append(plt.imshow(self.annotation_stacks[0], cmap=plt.get_cmap('gray'), interpolation='none'))
    
    #Loop to combine all video insets into 1
    print "...making final video..."
    def updatefig(j):
        print "...frame: ", j
        #plt.suptitle(self.selected_dff_filter+'  ' +self.dff_method + "\nFrame: "+str(j)+"  " +str(format(float(j)/self.img_rate-self.parent.n_sec,'.2f'))+"sec  ", fontsize = 15)
        plt.suptitle("Time: " +str(format(float(j)/self.img_rate-self.parent.n_sec,'.2f'))+"sec  Frame: "+str(j), fontsize = 15)

        # set the data in the axesimage object
        ctr=0
        for k in range(len(self.ca_stack)): 
            im[ctr].set_array(self.ca_stack[k][j]); ctr+=1
        
        im[ctr].set_array(self.pca_stack[j]); ctr+=1
        im[ctr].set_array(self.movie_stack[j]); ctr+=1
        im[ctr].set_array(self.lever_stack[j]); ctr+=1
        im[ctr].set_array(self.annotation_stacks[j]); ctr+=1

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.movie_stack)), interval=100, blit=False, repeat=True)
    #ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.ca_stack[1])), interval=100, blit=False, repeat=True)

    if True:
        #ani.save(self.parent.root_dir+self.parent.animal.name+"/movie_files/"+self.selected_session+'_'+str(len(self.movie_stack))+'_'+str(self.selected_trial)+'trial.mp4', writer=writer, dpi=300)
        ani.save(self.parent.root_dir+self.parent.animal.name+"/movie_files/"+self.selected_session+'_'+str(self.selected_code)+"_"+str(self.selected_trial)+'trial.mp4', writer=writer, dpi=600)
    plt.show()

def annotate_movies(self):
    ''' Make annotated movies
    '''
    #Constants for processing
    n_frames = 21000        #Number of frames of video to process
    video_rate = 15.        #Video rate in Hz.
    
    self.parent.n_sec = float(self.n_sec_window.text())
    self.start_time = -self.parent.n_sec; self.end_time = self.parent.n_sec
    self.temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session   
    self.abstimes = np.load(self.temp_file+'_abstimes.npy')
    self.img_rate = np.load(self.temp_file+'_img_rate.npy') #imaging rate


    #****************** LOAD RAW VIDEO ***********************
    #movie_raw = np.load(self.temp_file+'.m4v')
    #movie_raw = np.load(self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/movie.npy', mmap_mode='r')
    #print movie_raw.shape
    movie_stack = []
    n_movie_frames = 20000
    if True: 
        camera = cv2.VideoCapture(self.temp_file+'.m4v')

        #Find 200th frame in video: #Save cropped raw image into .npy array
        ctr = 0
        while ctr<n_movie_frames: 
            (grabbed, frame) = camera.read()
            ctr+=1
            movie_stack.append(cv2.cvtColor(frame, cv2.COLOR_RGB2GRAY))
            print ctr
        #image_original = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        #image_original_gray = cv2.cvtColor(frame, cv2.COLOR_RGB2GRAY)

    
    #************************ PROCESS ANNOTATIONS ********************* 

    #Process reward triggered data
    self.locs_44threshold = np.load(self.temp_file+'_locs44threshold.npy')
    self.code_44threshold = np.load(self.temp_file+'_code44threshold.npy')
       
    indexes = np.where(self.code_44threshold=="02")[0]
    self.locs_02 = self.locs_44threshold[indexes]
    indexes = np.where(self.code_44threshold=="04")[0]
    self.locs_04 = self.locs_44threshold[indexes]
    indexes = np.where(self.code_44threshold=="07")[0]
    self.locs_07 = self.locs_44threshold[indexes]

    print self.locs_02[:10],"...", self.locs_02[-10:]
    #print self.locs_04
    #print self.locs_07

    
    self.locs_annotated = []
    for annotated_cluster in self.annotated_clusters:
        self.locs_annotated.append(load_behavioural_annotation_data_all(self, annotated_cluster))   
    
    
    print self.locs_annotated[0][:10],"...", self.locs_annotated[0][-10:]
    print self.locs_annotated[1][:10],"...", self.locs_annotated[1][-10:]
    print self.locs_annotated[2][:10],"...", self.locs_annotated[2][-10:]


    #Make matrix for plotting
    annotated_matrix = np.zeros((6,n_frames), dtype=np.float32)
    
    #Shift time of behaviours based on blue_light data
    blue_light_index = np.load(self.temp_file+'_blue_light_frames.npy')[0]  #load number of video frames at which excitation light starts

    #Extend reward codes over 1second; mostly onwards from time of code; save for about 8 behavioural imaging frames ~=500ms
    for k in range(0,8,1):
        annotated_matrix[0][np.int32(self.locs_02*video_rate)+k+blue_light_index] = 1
        annotated_matrix[1][np.int32(self.locs_04*video_rate)+k+blue_light_index] = 2
        annotated_matrix[2][np.int32(self.locs_07*video_rate)+k+blue_light_index] = 3
    
    for k in range(len(self.locs_annotated)):
        annotated_matrix[k+3][self.locs_annotated[k]+blue_light_index]=4+k

    annotated_stack = []
    for k in range(n_movie_frames):
        annotated_stack.append(annotated_matrix[:, k:k+80])
    
    annotated_stack = np.array(annotated_stack)
    print annotated_stack.shape

    #ax = plt.subplot(111)
    #Make discrete colour map:
    cmap = mpl.colors.ListedColormap(['w','r', 'g', 'c','m','y','k'])
    
    #Plot data
    

    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=15000)
  
    fig = plt.figure()
    gs = gridspec.GridSpec(6,3)
    
        
    #PLOT MOVIES
    im = [] 
    ax = plt.subplot(gs[0:2,0:3])
    #ax = plt.subplot(2,1,1)
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    #v_max1 = np.nanmax(np.ma.abs(self.stack)); v_min1 = -v_max1; print v_min1, v_max1
    im.append(plt.imshow(movie_stack[0], cmap=plt.get_cmap('gray'), interpolation='none'))


    #PLOT ANNOTATIONS
    #ax = plt.subplot(2,1,2)
    ax = plt.subplot(gs[2:6,0:3])
    ax.get_xaxis().set_visible(False); ax.yaxis.set_ticks([]); ax.yaxis.labelpad = 0
    #Save label data
    old_ylabel = np.linspace(0,len(annotated_matrix),len(annotated_matrix)+1) #+.5
    new_ylabel = ['02','04','07']
    for k in self.annotated_clusters:
        if k == "paw_to_mouth":
            new_ylabel.append("pawing")
        else:
            new_ylabel.append(k)
    #new_ylabel = reversed(new_ylabel)
    plt.yticks(old_ylabel, new_ylabel, fontsize=10)    
    
    im.append(plt.imshow(annotated_stack[0], cmap=cmap,  aspect='auto', interpolation='none'))


    def updatefig(j):
        print "...making frame: ", j
        plt.suptitle("Frame: "+str(j)+"  " +str(format(float(j)/15.,'.2f'))+"sec", fontsize = 15)

        # set the data in the axesimage object
        im[0].set_array(movie_stack[j])
        im[1].set_array(annotated_stack[j])
        #im[1].set_array(self.stack_mean[j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(movie_stack)), interval=100, blit=False, repeat=True)

    if True:
        ani.save(self.temp_file+'.mp4', writer=writer)
    plt.show()    
    
    

    
    return






    ax = plt.subplot(111)
    #Make discrete colour map:
    cmap = mpl.colors.ListedColormap(['w','r', 'g', 'c','m','y','k'])
    
    #Plot data
    plt.imshow(annotated_matrix,  extent=[0,n_frames/video_rate,0,len(annotated_matrix)], aspect='auto', cmap=cmap)
    
    #Save label data
    old_ylabel = np.linspace(0,len(annotated_matrix),len(annotated_matrix)+1)+.5
    new_ylabel = ['02','04','07']
    for k in self.annotated_clusters:
        new_ylabel.append(k)
    new_ylabel = reversed(new_ylabel)
    plt.yticks(old_ylabel, new_ylabel, fontsize=18)    


    plt.tick_params(axis='both', which='both', labelsize=25)
    plt.xlabel("Time (sec)", fontsize=25)

    
    
    
    
    plt.show()

    
def find_isolated(self):
    
    
    
    #************* FIND ISOLATED BEHAIVOURS - 3 SEC GAPS *******************
    #LOOP OVER ALL ARRAYS
    for k in range(len(self.locs_annotated)):
        for frame in self.locs_annotated[k]:
            ctr=0
            print_array = ["behaviour: "+str(k)+" time: "+str(frame/video_rate)]
            for p in range(len(self.locs_annotated)):
                if p == k: break    #Don't check array against itself
                if abs(frame - self.locs_annotated[p][find_nearest(self.locs_annotated[p], frame)]) < 30:       #If nearest other behaviour is less than 30 frames away exit loop
                    break
                else:
                    #print "...behaviour: ", k, " frame: ", frame, "   nearest frame: ", self.locs_annotated[p][find_nearest(self.locs_annotated[p], frame)], \
                    #        "   compared behaviour: ", p
                    print_array.append("nearste behaviour: "+str(p)+ " time: "+ str(self.locs_annotated[p][find_nearest(self.locs_annotated[p], frame)]/video_rate))
                    ctr+=1
            if ctr==(len(self.locs_annotated)-1):
                print "... isolated behaviour..."
                print print_array
                print ""
   


def show_trial_locations(self):
    
    #Constants for processing
    n_frames = 21000        #Number of frames of video to process
    video_rate = 15.        #Video rate in Hz.
    
    self.temp_file = self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session    
    self.img_rate = np.load(self.temp_file+'_img_rate.npy') #imaging rate

    #Process reward annotated data
    self.locs_44threshold = np.load(self.temp_file+'_locs44threshold.npy')
    self.code_44threshold = np.load(self.temp_file+'_code44threshold.npy')
       
    indexes = np.where(self.code_44threshold=="02")[0]
    self.locs_02 = self.locs_44threshold[indexes]
    indexes = np.where(self.code_44threshold=="04")[0]
    self.locs_04 = self.locs_44threshold[indexes]
    indexes = np.where(self.code_44threshold=="07")[0]
    self.locs_07 = self.locs_44threshold[indexes]

    #Process behaviourally annotated data
    self.locs_annotated = []
    for annotated_cluster in self.annotated_clusters:
        self.locs_annotated.append(load_behavioural_annotation_data_all(self, annotated_cluster))   

    #Make matrix for plotting
    annotated_matrix = np.zeros((6,n_frames), dtype=np.float32)
    
    #Shift time of behaviours based on blue_light data
    blue_light_index = np.load(self.temp_file+'_blue_light_frames.npy')[0]  #load number of video frames at which excitation light starts

    #Extend reward codes over 1second; mostly onwards from time of code; save for about 8 behavioural imaging frames ~=500ms
    for k in range(0,8,1):
        annotated_matrix[0][np.int32(self.locs_02*video_rate)+k+blue_light_index] = 1
        annotated_matrix[1][np.int32(self.locs_04*video_rate)+k+blue_light_index] = 2
        annotated_matrix[2][np.int32(self.locs_07*video_rate)+k+blue_light_index] = 3
    
    for k in range(len(self.locs_annotated)):
        annotated_matrix[k+3][self.locs_annotated[k]+blue_light_index]=4+k
        
    ax = plt.subplot(111)
    #Make discrete colour map:
    cmap = mpl.colors.ListedColormap(['w','r', 'g', 'c','m','y','k'])
    
    #Plot data
    plt.imshow(annotated_matrix,  extent=[0,n_frames/video_rate,0,len(annotated_matrix)], aspect='auto', cmap=cmap)
    
    #Save label data
    old_ylabel = np.linspace(0,len(annotated_matrix),len(annotated_matrix)+1)+.5
    new_ylabel = ['02','04','07']
    for k in self.annotated_clusters:
        new_ylabel.append(k)
    new_ylabel = reversed(new_ylabel)
    
    print new_ylabel
    plt.yticks(old_ylabel, new_ylabel, fontsize=18)    
    plt.tick_params(axis='both', which='both', labelsize=25)
    plt.xlabel("Time (sec)", fontsize=25)
    
    
    plt.show()

    
def filter_data(self):
    """ Filter _aligned.npy files for lever_pull analysis.  
        NB: mean value of stack is added back into filtered data - so it isn't a purely filtered 
    """
    
    plotting = False
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
        print "...data already filtered...\n\n\n"
        return

    #Load aligned images
    print "... loading aligned imgs..."
    if os.path.exists(images_file):
        images_aligned = np.load(images_file)
    else:
        print "...images file does not exist (likely session is incomplete)...\n\n\n"
        return
    
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
        
    data= np.load(filtered_file,  mmap_mode='c')
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
        
    data= np.load(filtered_file,  mmap_mode='c')
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
        
    data= np.load(filtered_file,  mmap_mode='c')
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
        
    data= np.load(filtered_file,  mmap_mode='c')
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
        
    data= np.load(self.filtered_file,  mmap_mode='c')
    print data.shape

    #Load stack and mean of filtered data
    self.stack = data[start_frame:start_frame+n_frames]
    print self.stack.shape
    self.stack = self.stack.reshape(self.stack.shape[0],-1)
    print self.stack.shape
    
    dim_reduction_stack(self)



def PCA_reduction(X, n_components):


    plt.cla()
    #pca = decomposition.SparsePCA(n_components=3, n_jobs=1)
    pca = decomposition.PCA(n_components=n_components)

    print "... fitting PCA ..."
    pca.fit(X)
    
    print "... pca transform..."
    return pca.transform(X)
        

def compute_dff_events(self):

    compress_factor = 50.   #Needed to uncompress the LFP compressed sorts

    print "\n\n... dff computation event triggers..."
    print "... control computation: ", self.selected_control

    images_file = self.selected_recording
    print "... images file: ", images_file

    #*****************************************************************
    #************* LOAD FILTERED AND UNFILTERED IMAGES ***************
    #*****************************************************************
    #Load unfiltered and filtered data using mmap
    if self.selected_dff_filter!='nofilter': 
        self.filtered_file = images_file[:-4]+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz.npy'
    else: 
        self.filtered_file = images_file[:-4]+'.npy'   #Load unfiltered file and use it as "filtered_file"
    self.images_filtered = np.load(self.filtered_file,  mmap_mode='c')
    print "...images_file.shape: ",  self.images_filtered.shape
    
    self.unfiltered_file = images_file[:-4]+'.npy'
    self.images_unfiltered = np.load(self.unfiltered_file,  mmap_mode='c')
    
    if os.path.exists(images_file[:-4]+'_mean.npy'):
        global_mean = np.load(images_file[:-4]+'_mean.npy')
    else:
        global_mean = np.average(self.images_unfiltered, axis=0)
        np.save(images_file[:-4]+'_mean.npy', global_mean)
    print "... data_mean.shape: ", global_mean.shape 



    main_dir = os.path.dirname(os.path.dirname(self.selected_recording)[:-1])   #Strip file name and 'tif_files' directory 
    rec_name = self.selected_sort[:-5].replace(main_dir,'').replace('/tsf_files/','')   #Use the sort name - NOT the recording name 

    fs = np.loadtxt(main_dir+'/img_rate.txt')     #imaging rate
    #Check for filtered version of imaging data w. current filtering params
    lowcut = float(self.lowcut.text()); highcut = float(self.highcut.text())
    print "... frame rate: ", fs, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut



    
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
    

    
def compute_dff_events_mcd_all(self):

    print "\n\n... dff computation event triggers..."
    print "... control computation: ", self.selected_control

    images_file = self.selected_recording
    print "... images file: ", images_file


    #*****************************************************************
    #************* LOAD FILTERED AND UNFILTERED IMAGES ***************
    #*****************************************************************
    #Load unfiltered and filtered data using mmap
    if self.selected_dff_filter!='nofilter': 
        self.filtered_file = images_file[:-4]+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz.npy'
    else: 
        self.filtered_file = images_file[:-4]+'.npy'   #Load unfiltered file and use it as "filtered_file"
    self.images_filtered = np.load(self.filtered_file,  mmap_mode='c')
    print "...images_file.shape: ",  self.images_filtered.shape
    n_pixels = len(self.images_filtered[0])

    self.unfiltered_file = images_file[:-4]+'.npy'
    self.images_unfiltered = np.load(self.unfiltered_file,  mmap_mode='c')
    
    if os.path.exists(images_file[:-4]+'_mean.npy'):
        global_mean = np.load(images_file[:-4]+'_mean.npy')
    else:
        global_mean = np.average(self.images_unfiltered, axis=0)
        np.save(images_file[:-4]+'_mean.npy', global_mean)
    print "... data_mean.shape: ", global_mean.shape 
    

    #*****************************************************************
    #************** LOAD CAMERA ON/OFF AND IMG_RATE ******************
    #*****************************************************************
    #Load ON/OFF light and compute interpolated img rate based on # of frames aquired divided by ephys time stamps for start and end
    path_dir = os.path.dirname(images_file)
    epochs_file = path_dir+"/epochs.txt"

    #Check to see if multi-epoch file
    if os.path.exists(epochs_file[0]):  
        print "... multi-epoch recording, loading epochs and rec index..."
        epochs = np.loadtxt(epochs_file)
        rec_index = int(np.loadtxt(path_dir+'/rec_index.txt'))
        img_start, img_end = epochs[rec_index]
    else: 
        mcd_file = glob.glob(path_dir+"/*.mcd")
        print mcd_file

        if len(mcd_file)>1: print "... TOO MANY .MCD FILES..."; return
        imaging_onoff_file = mcd_file[0][:-4]+'_imagingonoff.txt'
        
        if os.path.exists(imaging_onoff_file)==False:
            #Find .mcd file, make sure only a single file:
            MCD_read_imagingtimes(mcd_file[0])

        onoff_pulse = np.loadtxt(imaging_onoff_file)
        img_start = onoff_pulse[0]  
        img_end = onoff_pulse[1] 
    
    self.reclength = img_end - img_start
    print "...imaging start offset, img_end: ", img_start, img_end
    
    self.n_images=len(self.images_filtered)
    session_img_rate = self.n_images/self.reclength
    print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate
    
    np.savetxt(path_dir+'/img_rate.txt', [session_img_rate])

    #Check for filtered version of imaging data w. current filtering params
    lowcut = float(self.lowcut.text()); highcut = float(self.highcut.text())
    print "... session img rate: ", session_img_rate, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

    #Check if DFF previously done
    print "... selected filter: ", self.selected_dff_filter
    if self.selected_dff_filter == 'nofilter':
        self.stm_filename = images_file[:-4]+"_"+ self.n_sec_window.text()+"sec_"+ self.dff_method+ '_' + self.selected_dff_filter + ".npy"
    else:
        self.stm_filename = images_file[:-4]+"_"+ self.n_sec_window.text()+"sec_"+ self.dff_method+ '_' + self.selected_dff_filter + '_' + str(lowcut)+"hz_"+str(highcut)+"hz.npy"
    
    print self.stm_filename
    
    

    #Load units:
    Sort = Ptcs(self.selected_sort)

    for unit in range(len(Sort.units)):

        #*****************************************************************
        #************************ LOAD EVENT FILE*************************
        #*****************************************************************
        if '.ptcs' in self.selected_sort:
            Sort = Ptcs(self.selected_sort)
            spikes = np.float32(Sort.units[unit])/1.E6    #Convert int64 usec to float32 seconds
            
            #compute the window of the recording window; 
            spike_indexes = np.where(np.logical_and(spikes>=img_start, spikes<=img_end))[0]    #Exclude spikes too close to beginning or end of recordings.
            events = spikes[spike_indexes]-img_start           #Align to imaging times by removing pre-imaging period
          
        else:
            events = np.loadtxt(self.selected_sort)                     #Manual time points are relative to ophys start; 

        print "... first 10 events: ", events[:10]
        
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
        self.img_frame_triggers = img_frame_triggers


        #import parmap
        ##motion_mask_array = parmap.map(do_filter, pixels, b, a, processes=30)
        #data_stm = parmap.map(parallel_dff, img_frame_triggers, self.images_filtered, self.window, self.dff_method, self.selected_dff_filter, processes=10)


        data_stm = np.zeros((1000, int(self.window*2),n_pixels,n_pixels), dtype=np.float32) #Make STMs for only 1000 spikes; 
        #data_stm = []
        for ctr, trigger in enumerate(img_frame_triggers):
            
            if ctr==1000: break #Only consider 1000 spikes in train
            print "...frame0: ", trigger, "  event #: ", ctr, " / ", len(img_frame_triggers)

            data_chunk = np.float32(self.images_filtered[int(trigger-self.window):int(trigger+self.window)])[:int(self.window*2)] #Fix # of frames
            
            if self.dff_method == 'globalAverage':
                if self.selected_dff_filter != 'nofilter':
                    temp_data = data_chunk/global_mean
                    data_stm[ctr] = temp_data    #Only need to divide by global mean as original data_chunk did not have mean img added in
                else: 
                    temp_data =(data_chunk-global_mean)/global_mean
                    data_stm[ctr] = temp_data
                
            elif self.dff_method == 'slidingWindow':            #Use baseline -2*window .. -window
                baseline = np.average(np.float32(self.images_unfiltered[int(trigger-2*self.window):int(trigger-self.window)]), axis=0)
                if self.selected_dff_filter != 'nofilter': 
                    temp_data= data_chunk/baseline
                    data_stm[ctr] = temp_data      #ignore subtracting baseline because it was never added back in 
                else:
                    temp_data= (data_chunk-baseline)/baseline
                    data_stm[ctr] = temp_data
        
        #plt.imshow(data_stm[0][90])
        #plt.show()
        #average_stm = data_stm/len(img_frame_triggers)

        #Save average of event triggered STMs
        print "Saving trial DFF...",

        if self.selected_dff_filter !='nofilter':
            stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz_'+\
            self.dff_method+'_unit'+str(unit).zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window'
        else:
            stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
            self.dff_method+'_unit'+str(unit).zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window'
        
        if self.selected_control =="yes": stm_file_name = stm_file_name + '_control'


        #Save average and variance data
        data_var = np.var(data_stm, axis=0)
        np.save(stm_file_name+'_var', data_var)
        
        data_ave = np.average(data_stm, axis=0)
        np.save(stm_file_name+'_mean', data_ave)

        #STD = sqrt(variance)
        #data_std = np.std(data_stm, axis=0)
        #np.save(stm_file_name+'_std', data_std)


        #Save all original data stack
        #np.save(images_file[:-4], data_stm)

        print ''
        self.images_filtered = np.load(self.filtered_file,  mmap_mode='c')
        
        
        

        
def compute_dff_events_mcd(self):

    print "\n\n... dff computation event triggers..."
    print "... control computation: ", self.selected_control

    images_file = self.selected_recording
    print "... images file: ", images_file


    #*****************************************************************
    #************* LOAD FILTERED AND UNFILTERED IMAGES ***************
    #*****************************************************************
    #Load unfiltered and filtered data using mmap
    if self.selected_dff_filter!='nofilter': 
        self.filtered_file = images_file[:-4]+'_'+self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz.npy'
    else: 
        self.filtered_file = images_file[:-4]+'.npy'   #Load unfiltered file and use it as "filtered_file"
    self.images_filtered = np.load(self.filtered_file,  mmap_mode='c')
    print "...images_file.shape: ",  self.images_filtered.shape
    
    self.unfiltered_file = images_file[:-4]+'.npy'
    self.images_unfiltered = np.load(self.unfiltered_file,  mmap_mode='c')
    
    if os.path.exists(images_file[:-4]+'_mean.npy'):
        global_mean = np.load(images_file[:-4]+'_mean.npy')
    else:
        global_mean = np.average(self.images_unfiltered, axis=0)
        np.save(images_file[:-4]+'_mean.npy', global_mean)
    print "... data_mean.shape: ", global_mean.shape 
    

    #*****************************************************************
    #************** LOAD CAMERA ON/OFF AND IMG_RATE ******************
    #*****************************************************************
    #Load ON/OFF light and compute interpolated img rate based on # of frames aquired divided by ephys time stamps for start and end
    path_dir = os.path.dirname(images_file)
    epochs_file = path_dir+"/epochs.txt"
    
    #Check to see if multi-epoch file
    if os.path.exists(epochs_file):  
        print "... multi-epoch recording, loading epochs and rec index..."
        epochs = np.loadtxt(epochs_file)
        rec_index = int(np.loadtxt(path_dir+'/rec_index.txt'))
        img_start, img_end = epochs[rec_index]
    else: 
        mcd_file = glob.glob(path_dir+"/*.mcd")
        print mcd_file

        if len(mcd_file)>1: print "... TOO MANY .MCD FILES..."; return
        imaging_onoff_file = mcd_file[0][:-4]+'_imagingonoff.txt'
        
        MCD_read_imagingtimes(mcd_file[0])

        onoff_pulse = np.loadtxt(imaging_onoff_file)
        
        img_start = onoff_pulse[0]  
        img_end = onoff_pulse[1] 
    
    self.reclength = img_end - img_start
    print "...imaging start offset, img_end: ", img_start, img_end
    
    self.n_images=len(self.images_filtered)
    session_img_rate = self.n_images/self.reclength
    print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate
    
    np.savetxt(path_dir+'/img_rate.txt', [session_img_rate])

    #Check for filtered version of imaging data w. current filtering params
    lowcut = float(self.lowcut.text()); highcut = float(self.highcut.text())
    print "... session img rate: ", session_img_rate, "  low_cutoff: ", lowcut, "  high_cutoff: ", highcut

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
    #************************ LOAD EVENT FILE*************************
    #*****************************************************************
    if '.ptcs' in self.selected_sort:
        Sort = Ptcs(self.selected_sort)
        spikes = np.float32(Sort.units[int(self.selected_unit)])/1.E6    #Convert int64 usec to float32 seconds
        
        #compute the window of the recording window; 
        spike_indexes = np.where(np.logical_and(spikes>=img_start, spikes<=img_end))[0]    #Exclude spikes too close to beginning or end of recordings.
        events = spikes[spike_indexes]-img_start           #Align to imaging times by removing pre-imaging period
        
    else:
        events = np.loadtxt(self.selected_sort)                     #Manual time points are relative to ophys start; 

    print "... # events: ", events
    
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
    self.img_frame_triggers = img_frame_triggers


    #import parmap
    ##motion_mask_array = parmap.map(do_filter, pixels, b, a, processes=30)
    #data_stm = parmap.map(parallel_dff, img_frame_triggers, self.images_filtered, self.window, self.dff_method, self.selected_dff_filter, processes=10)


    data_stm = np.zeros((1000, int(self.window*2),n_pixels,n_pixels), dtype=np.float32) #Make 
    #data_stm = []
    for ctr, trigger in enumerate(img_frame_triggers):
        
        if ctr==1000: break #Only consider 1000 spikes in train
        print "...frame0: ", trigger, "  event #: ", ctr, " / ", len(img_frame_triggers)

        data_chunk = np.float32(self.images_filtered[int(trigger-self.window):int(trigger+self.window)])[:int(self.window*2)] #Fix # of frames
        
        if self.dff_method == 'globalAverage':
            if self.selected_dff_filter != 'nofilter':
                temp_data = data_chunk/global_mean
                data_stm[ctr] = temp_data    #Only need to divide by global mean as original data_chunk did not have mean img added in
            else: 
                temp_data =(data_chunk-global_mean)/global_mean
                data_stm[ctr] = temp_data
            
        elif self.dff_method == 'slidingWindow':            #Use baseline -2*window .. -window
            baseline = np.average(np.float32(self.images_unfiltered[int(trigger-2*self.window):int(trigger-self.window)]), axis=0)
            if self.selected_dff_filter != 'nofilter': 
                temp_data= data_chunk/baseline
                data_stm[ctr] = temp_data      #ignore subtracting baseline because it was never added back in 
            else:
                temp_data= (data_chunk-baseline)/baseline
                data_stm[ctr] = temp_data
    
    
    #plt.imshow(data_stm[0][90])
    #plt.show()
    #average_stm = data_stm/len(img_frame_triggers)

    #Save average of event triggered STMs
    print "Saving trial DFF...",

    if self.selected_dff_filter !='nofilter':
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+self.lowcut.text()+'hz_'+self.highcut.text()+'hz_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window'
    else:
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window'
    
    if self.selected_control =="yes": stm_file_name = stm_file_name + '_control'


    #Save average and variance data
    print "...computing variance..."
    data_var = np.var(data_stm, axis=0)
    np.save(stm_file_name+'_var', data_var)
    
    print "...computing mean..."
    data_ave = np.average(data_stm, axis=0)
    np.save(stm_file_name+'_mean', data_ave)

    #STD = sqrt(variance)
    #data_std = np.std(data_stm, axis=0)
    #np.save(stm_file_name+'_std', data_std)


    #Save all original data stack
    #np.save(images_file[:-4], data_stm)

    print ''
    self.images_filtered = 0;     #Set aligned_images to empty 
    

def parallel_dff(trigger, images_filtered, window, dff_method, selected_dff_filter):
    
    #print "...frame0: ", trigger, "  event #: ", ctr, " / ", len(img_frame_triggers)

    data_chunk = np.float32(images_filtered[int(trigger-window):int(trigger+window)])[:int(window*2)] #Fix # of frames
    
    if dff_method == 'globalAverage':
        if selected_dff_filter!='nofilter':
            temp_data = data_chunk/global_mean
            data_stm= temp_data    #Only need to divide by global mean as original data_chunk did not have mean img added in
        else: 
            temp_data =(data_chunk-global_mean)/global_mean
            data_stm = temp_data
        
    elif dff_method == 'slidingWindow':            #Use baseline -2*window .. -window
        baseline = np.average(np.float32(self.images_unfiltered[int(trigger-2*self.window):int(trigger-self.window)]), axis=0)
        if selected_dff_filter!='nofilter': 
            temp_data= data_chunk/baseline
            data_stm = temp_data      #ignore subtracting baseline because it was never added back in 
        else:
            temp_data= (data_chunk-baseline)/baseline
            data_stm = temp_data
    
    return data_stm


def load_behavioural_annotation_data(self):
    ''' Find all annotation fields called "*_clusters.npz" and add their annotations '''
    filenames = glob.glob(self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session+"*_clusters.npz")
    
    print "...img_rate: ", self.parent.animal.img_rate

    for k in range(len(filenames)):
        print filenames[k]
        data = np.load(filenames[k])
        
        for cluster, cluster_data in zip(data['cluster_names'], data['cluster_indexes']):
            
            if cluster==self.selected_code: 
                print "... loading: ", cluster
                #Save locations and ids of events
                cluster_indexes = [cluster]*len(cluster_data)
                print "...cluster_indexes: ", cluster_indexes[:10]
                print "...cluster_data: ", cluster_data[:10], cluster_data[-10:]
                print len(cluster_data)

                #SET BOTH OVERALL THRESHOLDS AND SELECTED THRESHOLDS TO THE SAME VALUES....
                self.code_44threshold= cluster_indexes 
                self.locs_44threshold= np.int32(cluster_data)/float(15.)       #Unclear if this is matching up to data correctly

                self.code_44threshold_selected= cluster_indexes 
                self.locs_44threshold_selected= np.int32(cluster_data)/float(15.)       #Unclear if this is matching up to data correctly
                                                #************************CHECK THIS**************************
                return    
    
        
    
    #cluster_1 = data['cluster_indexes'][1]; cluster_2 = data['cluster_indexes'][0]  #Licking info is in 2nd cluster
    #t = np.arange(0,max(np.max(cluster_1), np.max(cluster_2)),1)/15.

    #clrs = []
    #locs = []
    #last_loc = 0
    #for k in range(len(t)):
        #if k in cluster_1: 
            #locs.append(1.)
            #last_loc=k
        #else:
            #if (k-last_loc)<2:
                #locs.append(1)
            #else:
                #locs.append(0)

    #plt.plot(t, locs, color='red', alpha=0.8)#, color=clrs)
    #ax.fill_between(t, np.zeros(len(locs)), locs, color='red', alpha=0.2)

        
    #indexes = np.where(self.code_44threshold==self.selected_code)[0]
    #print "...indexes: "; print indexes

    #self.code_44threshold_selected = self.code_44threshold[indexes]
    #self.locs_44threshold_selected = self.locs_44threshold[indexes]
    
    
    
def load_behavioural_annotation_data_all(self, annotated_cluster):
    
    ''' Find all annotation fields called "*_clusters.npz" and add their annotations '''
    filenames = glob.glob(self.parent.root_dir + self.parent.animal.name + '/tif_files/'+self.selected_session+'/'+self.selected_session+"*_clusters.npz")
    
    for k in range(len(filenames)):
        data = np.load(filenames[k])
        
        for cluster, cluster_data in zip(data['cluster_names'], data['cluster_indexes']):
            
            if cluster==annotated_cluster: 
                #Save locations and ids of events
                cluster_indexes = [cluster]*len(cluster_data)

                #SET BOTH OVERALL THRESHOLDS AND SELECTED THRESHOLDS TO THE SAME VALUES....
                self.code_44threshold= cluster_indexes 
                self.locs_44threshold= np.int32(cluster_data)       #LEAVE DATA HERE IN FRAME TIME INTEGERS UNTIL READY TO USE

                self.code_44threshold_selected= cluster_indexes 
                self.locs_44threshold_selected= np.int32(cluster_data)
                                                                                        #************************CHECK THIS**************************
                return self.locs_44threshold_selected
    
    
    

def compute_dff_mouse_lever(self):
    print "\n\n... dff computation..."

    #Load average frame
    self.rec_filename = self.selected_session.replace(self.parent.root_dir+self.parent.animal.name+"/tif_files/",'')

    images_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_aligned.npy'
    data_mean = np.load(images_file[:-4]+'_mean.npy')
    print "... data_mean.shape: ", data_mean.shape #; plt.imshow(data_mean); plt.show()
    
    #Check for filtered version of imaging data w. current filtering params
    self.lowcut = float(self.parent.filter_low.text()); self.highcut = float(self.parent.filter_high.text())
    fs = self.parent.animal.img_rate
    print "... frame rate: ", fs, "  low_cutoff: ", self.lowcut, "  high_cutoff: ", self.highcut
    
    print "...self.selected_code: ", self.selected_code

    #Process reward triggered data
    if (self.selected_code =='02') or (self.selected_code =='04') or (self.selected_code =='07'):
        print "...self.locs_44threshold: "; print self.locs_44threshold
        print "...self.code_44threshold: "; print self.code_44threshold
        
        indexes = np.where(self.code_44threshold==self.selected_code)[0]
        print "...indexes: "; print indexes

        self.code_44threshold_selected = self.code_44threshold[indexes]
        self.locs_44threshold_selected = self.locs_44threshold[indexes]

    #Process behaviour triggered data;
    else:
        print self.code_44threshold_selected[:10] 
        print self.locs_44threshold_selected[:10]
    
    
    print len(self.code_44threshold_selected)
    print len(self.locs_44threshold_selected)
    
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
    temp_filearray = []
    for p in range(len(temp_tif_files)): 
        temp_filearray.append(temp_tif_files[p].replace("12TB", self.parent.replacement_dir))
    temp_tif_files = temp_filearray

    temp_event_files = np.load(self.parent.animal.home_dir+self.parent.animal.name+"/event_files.npy")
    temp_filearray = []
    for p in range(len(temp_event_files)): 
        temp_filearray.append(temp_event_files[p].replace("12TB", self.parent.replacement_dir))
    temp_event_files = temp_filearray
        
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
    
    
    #BASELINE FOR GLOBAL BASELINE REMOVAL
    mean_file = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_aligned_mean.npy'
    global_mean = np.load(mean_file)

    self.abstimes = np.load(self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_abstimes.npy')
    self.abspositions = np.load(self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.rec_filename+'/'+self.rec_filename+'_abspositions.npy')

    print "...computing DF/F..."
    data_stm = []; traces = []; locs = []; codes = []
    counter=-1
    self.window = self.parent.n_sec * session_img_rate      #THIS MAY NOT BE GOOD ENOUGH; SHOULD ALWAYS GO BACK AT LEAST X SECONDS EVEN IF WINDOW IS ONLY 1SEC or 0.5sec...
                                                            #Alternatively: always compute using at least 3sec window, and then just zoom in
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
        #lever_position_index = find_nearest(np.array(self.abstimes), self.locs_44threshold[counter])
        lever_position_index = find_nearest(np.array(self.abstimes), self.locs_44threshold_selected[counter])
        
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


def view_mean_stm_events(self):
    
    block_save = int(self.block_save.text())
    
    stm_type = self.selected_stm_type

    images_file = self.selected_recording
    self.images_filtered = np.load(self.selected_recording, mmap_mode='c')

    #*****************************************************************
    #************** LOAD CAMERA ON/OFF AND IMG_RATE ******************
    #*****************************************************************
    #Load ON/OFF light and compute interpolated img rate based on # of frames aquired divided by ephys time stamps for start and end
    path_dir = os.path.dirname(images_file)
    epochs_file = path_dir+"/epochs.txt"
    
    #Check to see if multi-epoch file
    if os.path.exists(epochs_file):  
        print "... multi-epoch recording, loading epochs and rec index..."
        epochs = np.loadtxt(epochs_file)
        rec_index = int(np.loadtxt(path_dir+'/rec_index.txt'))
        img_start, img_end = epochs[rec_index]
    else: 
        mcd_file = glob.glob(path_dir+"/*.mcd")
        print mcd_file

        if len(mcd_file)>1: print "... TOO MANY .MCD FILES..."; return
        imaging_onoff_file = mcd_file[0][:-4]+'_imagingonoff.txt'
        
        MCD_read_imagingtimes(mcd_file[0])

        onoff_pulse = np.loadtxt(imaging_onoff_file)
        
        img_start = onoff_pulse[0]  
        img_end = onoff_pulse[1] 
    
    self.reclength = img_end - img_start
    print "...imaging start offset, img_end: ", img_start, img_end
    
    self.n_images=len(self.images_filtered)
    session_img_rate = np.loadtxt(path_dir+'/img_rate.txt')

    print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate
    
    
    #******************************************************************
    #**************** LOAD DATA ***************************************
    #******************************************************************

    path_dir = os.path.dirname(images_file)
    rate_file = path_dir+"/img_rate.txt"
    img_rate = np.loadtxt(rate_file)   
    

    if self.selected_dff_filter !='nofilter':
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window_mean'

    else:
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(self.parent.n_sec)+'sec_window_mean'
    
    
    path_dir = os.path.dirname(images_file)
    temp_file = os.path.split(path_dir)[1]

    search_string = path_dir+'/img_avg_'+temp_file+'_unit'+str(self.selected_unit.zfill(2))+"*"+stm_type+"*"
    stm_file_name = glob.glob(search_string)
    if len(stm_file_name)==1:
        stm_file_name = stm_file_name[0]
    else:
        print "... file not found: ", search_string
        return
    
    
    data = np.float32(np.load(stm_file_name))
    print data.shape
    
    temp_array = data
    main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/'

    plt.close()
    ax = plt.subplot(1,1,1)
    #img_rate = float(np.loadtxt(main_dir+'/img_rate.txt'))


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
        v_min = float(self.vmin_default.text()); v_max = float(self.vmax_default.text()); self.v_abs = max(v_min, v_max)

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
    plt.suptitle("STM - Mean   spiking type: " + stm_type)
    plt.show()



def view_var_stm_events(self):
    
    block_save = int(self.block_save.text())
    
    images_file = self.selected_recording
    path_dir = os.path.dirname(images_file)
    rate_file = path_dir+"img_rate.txt"
    img_rate = np.loadtxt(rate_file)   
    
    if self.selected_dff_filter !='nofilter':
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window_var'
    else:
        stm_file_name = images_file[:-4] + '_img_avg_' + self.selected_dff_filter+'_'+\
        self.dff_method+'_unit'+self.selected_unit.zfill(3)+'_'+str(int(self.parent.n_sec))+'sec_window_var'
    
    data = np.float32(np.load(stm_file_name+'.npy'))
    print data.shape
    
    temp_array = data
    main_dir = '/media/cat/8TB/in_vivo/tim/dongsheng/'

    plt.close()
    ax = plt.subplot(1,1,1)
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
        v_min = float(self.vmin_default.text()); v_max = float(self.vmax_default.text()); self.v_abs = max(v_min, v_max)

    
    print v_min, v_max

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
    plt.suptitle("STM - Variance")
    plt.show()




def sigmoid_function(x, a, b):

    return np.clip(a*(np.ma.log(x) - np.ma.log(1 - x))+b, 0, 1)       #Compute sigmoid and cut off values below 0 and above 1
    
    
def mangle(width, x, img_temp, maxval, power, val999):
    
    mu = 0 #Select approximate midline as centre of gaussian
    sig = width
    
    a = .005       #The steepness of the sigmoid function
    b = val999        #% of maxval to cutoff


    #Normalize img_temp for sigmoid to work properly
    #img_temp_norm = (img_temp-np.min(img_temp))/(np.max(img_temp) - np.min(img_temp))
    img_temp_norm = (img_temp-np.min(img_temp))/(np.max(img_temp) - np.min(img_temp))


    #Original root function
    #return -np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) * (abs(pix_val/maxval)**(1./power))
    return -np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) * sigmoid_function(img_temp_norm, a, b)



def regularization(maxval, val999, input_val):
    
    return -input_val*((input_val-val999)/(maxval-val999))      #Return a negative value up to -1.0



def motion_mask_parallel(img_temp, maxval, val999, width, power):
    '''Parallel computation of mask
    '''
   
    y_array = []
    for x in range(len(img_temp)):
        y_array.append(np.arange(0,len(img_temp), 1))
        
    y_array = np.vstack(y_array)
    motion_mask = img_temp*mangle(width, np.abs(64-y_array), img_temp, maxval, power, val999)
    

    return motion_mask


def view_static_stm(self):
    
    self.parent.n_sec = float(self.n_sec_window.text())
    width = int(self.mask_width.text()) #40
    power = float(self.mask_power.text()) #4.
    mask_percentile = float(self.mask_percentile.text()) #4.
    mask_power = float(self.mask_power.text()) #4.
    
    block_save = int(self.block_save.text())

    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

    filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'
    print "...stm_name: ", filename
    filename_motion_mask = filename[:-4]+"_motion_mask.npy"

    data = np.load(filename)
    print data.shape
    
    #Mask data
    print "...selected trial for stm: ", self.selected_trial
    temp_array = quick_mask(self, data[int(self.selected_trial)])


    #Load motion mask
    #motion_mask_array = np.load('/media/cat/12TB/in_vivo/tim/yuki/IA1/IA1_motion_mask.npy')
    #motion_mask = (motion_mask_array[10]-np.min(motion_mask_array[10]))/(np.max(motion_mask_array[10]) - np.min(motion_mask_array[10]))
    plt.close()


    #Make single trial motion mask:
    img_temp = temp_array
    maxval = np.max(img_temp)
    data_1d = img_temp.ravel()
    val999 = np.percentile(data_1d, mask_percentile)               #Mark stroke as the 97.5 percentile and higher values; 
    print "...masked percentile value: ", val999

    #Loop over every frame and scale data down
    motion_mask = np.zeros((len(img_temp), len(img_temp[0]), len(img_temp[0])), dtype=np.float32)

    motion_mask_array = []
    mask_start_frame = int(self.mask_start_frame.text());  mask_end_frame = int(self.mask_end_frame.text())
    if mask_start_frame != mask_end_frame:
        for k in range(int(len(img_temp)/2.)+mask_start_frame, int(len(img_temp)/2.) + mask_end_frame, 1):
            motion_mask_array.append(img_temp[k])
        #for k in range(0, len(img_temp), 1):
        #    motion_mask_array.append(img_temp[k])
            
        #Use parmap;  DO MASKING BEFORE BLOCK AVERAGING
        import parmap
        #motion_mask_array = parmap.map(do_filter, pixels, b, a, processes=30)
        motion_mask_array = parmap.map(motion_mask_parallel, motion_mask_array, maxval, val999, width, power, processes=10)
        

        #motion_file = self.parent.animal.home_dir + self.parent.animal.name + "/" + self.parent.animal.name+ '_motion_mask'
        #np.save(motion_file, motion_mask)

        #EITHER AVERAGE THE MASK, OR GRAB LARGEST VALUES IN MASK
        motion_mask_ave = np.ma.average(motion_mask_array, axis=0)
        motion_mask = np.power((motion_mask_ave-np.min(motion_mask_ave))/(np.max(motion_mask_ave) - np.min(motion_mask_ave)), np.zeros(motion_mask_ave.shape, dtype=np.float32)+mask_power)

    else:
        print "... no mask applied..."
        motion_mask = np.zeros(img_temp[0].shape)+1

    np.save(filename_motion_mask, np.ma.filled(motion_mask, 0))


    ax = plt.subplot(3,1,1)
    vabs = np.max(np.abs(motion_mask))
    plt.imshow(-motion_mask, vmin = -vabs, vmax=vabs)
    #plt.show()
    

    img_rate = self.parent.animal.img_rate
    #start_time = -self.parent.n_sec; end_time = self.parent.n_sec
    start_time = float(self.stm_start_time.text()); end_time = float(self.stm_end_time.text())

    img_out = []
    img_out_original = []
    for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
    #for i in range(0,int(2*img_rate*self.parent.n_sec), block_save):
        print i
        img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0)*motion_mask)
        img_out_original.append(np.ma.average(temp_array[i:i+block_save], axis=0))

    img_out = np.ma.hstack((img_out))

    
    #v_abs = max(np.nanmax(img_out),-np.nanmin(img_out))
    #plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
    v_abs = np.max(img_out)


    ax = plt.subplot(3,1,2)
    plt.imshow(img_out)

    plt.ylabel(str(round(v_abs*100,2))+"%", fontsize=14)
    ax.yaxis.set_ticks([])
    ax.xaxis.set_ticks([])

    #if ctr==(len(select_units)-1): 
    plt.xlabel("Time from '44' threshold crossing (sec)", fontsize=25)
    old_xlabel = np.linspace(0,img_out.shape[1], 11)
    new_xlabel = np.around(np.linspace(start_time,end_time, 11), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=18)

    plt.title(self.traces_filename)


    #SEE ORIGINAL DATA ALSO 
    img_out = np.ma.hstack((img_out_original))
    
    #v_abs = max(np.nanmax(img_out),-np.nanmin(img_out))
    #plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
    v_abs = np.max(img_out)
    
    ax = plt.subplot(3,1,3)

    plt.imshow(img_out)

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

    
    
def make_stm_motion_mask(self):
    ''' Average all/some STMs together to see if midline and other artifacts are substantial
    '''
    
    self.parent.n_sec = float(self.n_sec_window.text())

    block_save = int(self.block_save.text())

    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

    filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'
    print "...stm_name: ", filename                                         #Filename containing STMs for all trials for that particular code 

    #print "...resetting stm_name to: ", filename

    if os.path.exists(filename)==True: 
        data = np.load(filename,  mmap_mode='c')
    else:
        print "...data not yet processed ******************"
        return
    
    data_array = np.zeros(data.shape[1:], dtype=np.float32)
    for trial in range(len(data)): 
        data_array+= data[trial]
   
    data = data_array/len(data)
    
    #Mask data
    temp_array = quick_mask(self, data)
    
    block_save = 1

    plt.close()
    ax = plt.subplot(2,1,1)
    img_rate = self.parent.animal.img_rate
    #start_time = -self.parent.n_sec; end_time = self.parent.n_sec
    start_time = float(self.stm_start_time.text()); end_time = float(self.stm_end_time.text())

    img_temp = []
    for i in range(int(img_rate*(3+start_time)),int(img_rate*(3+end_time)), block_save):
    #for i in range(0,int(2*img_rate*self.parent.n_sec), block_save):
        print i
        #img_out.append(np.ma.average(temp_array[i:i+block_save], axis=0))
        img_temp.append(np.ma.average(temp_array[i:i+block_save], axis=0))

    img_out = np.ma.hstack((img_temp))
    
    #v_abs = max(np.nanmax(img_out),-np.nanmin(img_out))
    #plt.imshow(img_out, vmin = -v_abs, vmax=v_abs)
    v_abs = np.max(img_out)
    plt.imshow(img_out)

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

    #***************************************************************************************
    #PLOT PIXELS OVER 95 PERCENTILE
    ax = plt.subplot(2,1,2)
    
    x_shape = img_out.shape[0]
    y_shape = img_out.shape[1]
    
    data_1d = img_out.reshape(x_shape*y_shape)
    val999 = np.percentile(data_1d, 98)               #Mark stroke as the 97.5 percentile and higher values; 
    print val999
    maxval = np.max(data_1d)
    
    print maxval
    
    data_1d[data_1d > val999] = val999
    
    width = 60
    
    #Loop over every frame and scale data down
    motion_mask = np.zeros((len(img_temp), len(img_temp[0]), len(img_temp[0])), dtype=np.float32)
    for k in range(len(img_temp)/2.-3,len(img_temp)/2.+3, 1):
        print "... frame: ", k
        for x in range(len(img_temp[k])):
            for y in range(len(img_temp[k][x])):
                #if img_temp[k][x][y]>val999:     #OR JUST ALL PIXELS UNDER THE GAUSSIAN
                #if abs(64-y)<15:
                #if val999>img_temp[k][x][y]:
                if img_temp[k][x][y]> 0: 
                    motion_mask[k][x][y] = img_temp[k][x][y]*inverted_gaussian(width, abs(64-y), img_temp[k][x][y], maxval)
                    img_temp[k][x][y] = img_temp[k][x][y] + img_temp[k][x][y]*inverted_gaussian(width, abs(64-y), img_temp[k][x][y], maxval)
                
                    #img_temp[k][x][y] += regularization(maxval, val999, img_temp[k][x][y])

    motion_file = self.parent.animal.home_dir + self.parent.animal.name + "/" + self.parent.animal.name+ '_motion_mask'
    np.save(motion_file, motion_mask)

    motion_file_single = motion_file+"_single"
    motion_mask_single = np.average(motion_mask, axis=0)
    np.save(motion_file_single, motion_mask_single)
    
    img_out = np.ma.hstack((img_temp))

    #data_1d = img_out.reshape(img_out.shape[0]*img_out.shape[1])
    #data_1d[data_1d > val999] = 0
    
    #img_out = data_1d.reshape(x_shape, y_shape)

    v_abs = np.max(img_out)
    plt.imshow(img_out)

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
    start_time = -self.parent.n_sec; end_time = self.parent.n_sec

    block_save = int(self.block_save.text())

    if self.selected_dff_filter == 'nofilter':
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_"+ self.selected_dff_filter+'_' +self.dff_method+'_'+str(self.selected_code)+"code_traces.npy"
    else:
        self.traces_filename = self.parent.animal.home_dir+self.parent.animal.name+'/tif_files/'+self.selected_session+'/'+self.selected_session+"_"+ \
            str(self.parent.n_sec)+"sec_" + self.selected_dff_filter + "_"+self.dff_method+'_'+self.parent.filter_low.text()+"hz_"+self.parent.filter_high.text()+"hz_"+str(self.selected_code)+"code_traces.npy"

    #Load single trial data
    filename = self.traces_filename.replace('_traces.npy','')+'_stm.npy'
    print "...stm_name: ", filename

    data = np.load(filename)
    print data.shape
    
    #Apply Generic Mask 
    print "...selected trial for stm: ", self.selected_trial
    vid_array = quick_mask(self, data[int(self.selected_trial)])

   
    #Apply motion Mask
    filename_motion_mask = filename[:-4]+ self.selected_session+"_motion_mask.npy"
    if os.path.exists(filename_motion_mask)==False:
        print "...motion mask missing..."
        return
    
    motion_mask = np.load(filename_motion_mask)
    
    print motion_mask
    vid_array = vid_array * motion_mask
    

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
        #im.set_array(vid_array[j])
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
    
    n_pixels = len(data[0])
    
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
        generic_mask_indexes=np.zeros((n_pixels,n_pixels))
        for i in range(len(generic_coords)):
            generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True
    else:
        print "...generic mask not found..."
        return
    
    #Load midline mask
    for i in range(int(midline_mask_n_pixels)):
        generic_mask_indexes[:,n_pixels/2+int(int(midline_mask_n_pixels)/2)-i]=True

    #Apply full mask; probably FASTER METHOD
    n_pixels = n_pixels
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
    for i in range(0, len(data), 1):
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

    lfpwidth=10
    lfptres=5     #Time resolution: bin width for analysis in seconds
    
    lowband = [0.1, 4]; highband = [15,100]     #Martin's suggestions; Saleem 2010?
    #lowband = [0.1, 4]; highband = [4,100]     # Cat's implementation
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
    tsf.vscale_HP = 0.1     #USE ADperBit TO INCREASE RESOLUTIN
    tsf.iformat = 1002

    tsf.n_electrodes = len(ncs_files)
    tsf.n_cell_spikes = 0

    tsf.n_cell_spikes = 0
    tsf.layout = np.arange(tsf.n_electrodes)
    tsf.file_names = []
    tsf.n_samples = []
    tsf.n_digital_chs = []
    tsf.digital_chs = []
    
    tsf.Siteloc = np.zeros((tsf.n_electrodes*2), dtype=np.int16) #Read as 1D array
    for i in range (tsf.n_electrodes):
        tsf.Siteloc[i*2]=0
        tsf.Siteloc[i*2+1]=i*50 #GUESSING each tetrode is 50um apart
        
            
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
        #plt.plot(data[0][:1000000])
        #plt.show()
    
    #Trunkate extra voltage values (some chs in Neuralynx recs have more/less values than others)
    tsf.n_vd_samples = min_samples 
    for k in range(len(tsf.ec_traces)):
        tsf.ec_traces[k]=tsf.ec_traces[k][:min_samples]
        
    #tsf.ec_traces = np.array(tsf.ec_traces, dtype=np.int16)
    
    
    #******************SAVE HIGH PASS RECORD******************
    #if self.parent.make_hp.text()=='True':
    print '\n...saving alltrack _hp_fromlfp.tsf...'
    file_name = ncs_files[0][:-4]+"_alltrack_raw.tsf"
    save_tsf_single(tsf, file_name)
    #else:
    #    print "...skipping hp save..."
    
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
    tsf.layout = np.arange(tsf.n_electrodes)
    tsf.Siteloc = np.zeros((tsf.n_electrodes*2), dtype=np.int16) #Read as 1D array
    tsf.file_names = []
    tsf.n_samples = []
    tsf.n_digital_chs = []
    tsf.digital_chs = []
    
    for i in range (tsf.n_electrodes):
        tsf.Siteloc[i*2]=0
        tsf.Siteloc[i*2+1]=i*50 #GUESSING each tetrode is 50um apart
            
    #Wavelet filter record first
    #tsf.ec_traces = wavelet(tsf.ec_traces)

    print ''; print "...saving alltrack _hp.tsf..."
    file_name = ntt_files[0][:-4]+"_alltrack_hp.tsf"
    save_tsf_single(tsf, file_name)

def filter_ephys(self):
    
    tsf = TSF.TSF(self.tsf_filename)
    tsf.read_ec_traces()
    print tsf.SampleFrequency
    
    for ch in range(0, tsf.n_electrodes):

        trace_temp = tsf.ec_traces[ch] #[int(float(self.time_start.text())*tsf.SampleFrequency): int(float(self.time_end.text())*tsf.SampleFrequency)] 
        
        #trace_temp = butter_highpass_filter(trace_temp, float(self.low_cutoff.text()), tsf.SampleFrequency, 5)
        
        tsf.ec_traces[ch] = butter_lowpass_filter(trace_temp, float(self.high_cutoff.text()), tsf.SampleFrequency, 5)
    

    tsf.file_names=[]; tsf.n_samples=[]; tsf.n_digital_chs=[]; tsf.digital_chs=[]
    tsf.layout = np.arange(tsf.n_electrodes)

    save_tsf_single(tsf, self.tsf_filename[:-4]+"_low"+self.low_cutoff.text()+"_high"+self.high_cutoff.text()+".tsf")
    

def notch_ephys(self):
    
    tsf = TSF.TSF(self.tsf_filename)
    tsf.read_ec_traces()
    print tsf.SampleFrequency
    
    for ch in range(0, tsf.n_electrodes):

        trace_temp = tsf.ec_traces[ch] #[int(float(self.time_start.text())*tsf.SampleFrequency): int(float(self.time_end.text())*tsf.SampleFrequency)] 
        
        #trace_temp = butter_highpass_filter(trace_temp, float(self.low_cutoff.text()), tsf.SampleFrequency, 5)
        
        tsf.ec_traces[ch] = Notch_Filter(tsf.ec_traces[ch])

    tsf.file_names=[]; tsf.n_samples=[]; tsf.n_digital_chs=[]; tsf.digital_chs=[]
    tsf.layout = np.arange(tsf.n_electrodes)

  
    save_tsf_single(tsf, self.tsf_filename[:-4]+"_notch.tsf")
    
    

def Plot_rasters(self):

    #********************************************************************************************
    #*************************************** PLOT SPECGRAM FIRST ********************************
    #********************************************************************************************
    
    colors=['blue','red', 'green', 'violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    font_size = 30
    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    
    #********************************************************************************************
    #*************************************** PLOT RASTERS ***************************************
    #********************************************************************************************

    
    #Load LFP Sort
    #Sort_lfp = PTCS.PTCS(self.selected_sort_lfp) #Auto load flag for Nick's data

    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'

    Sort_sua = PTCS.PTCS(self.selected_sort_sua) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    if True: 
        n_spikes = []
        for k in range(len(Sort_sua.units)):
            n_spikes.append(len(Sort_sua.units[k]))
        n_spikes = np.array(n_spikes)
        indexes = np.argsort(n_spikes)
        print indexes
    else:
        indexes = np.arange(Sort_sua.n_units)

    offset = 10
    
    ax = plt.subplot(1, 1, 1)
    y = []
    for i in indexes: #range(len(Sort_sua.units)):
    #for i in range(10):
        print "... unit: ", i
    #for i in indexes[0:5]: #range(len(Sort_sua.units)):
        #x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)/float(Sort_sua.samplerate) #float(Sort1.samplerate)*2.5
        x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)*1E-6

        #x = spikes[np.where(np.logical_and(spikes>=int(self.time_start.text()), spikes<=int(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=offset+0.4
        ymax+=offset-0.4

        #plt.vlines(x-int(self.time_start.text()), ymin, ymax, linewidth=1, color='black') #colors[mod(counter,7)])
        plt.vlines(x, ymin, ymax, linewidth=1, color=colors[i%7], alpha=1) #colors[mod(counter,7)])

        y.append(x)
        
        offset=offset-1.0

    ##Plot LFP spike
    #offset = offset -10.
    #y = []
    #for i in range(len(Sort_lfp.units)):
        ##x = np.array(Sort_lfp.units[i],dtype=np.float32)/float(Sort_sua.samplerate)*50 #***************************** UNCOMPRESSSING LFP RASTERS
        #spikes = np.array(Sort_lfp.units[i],dtype=np.float32)*1E-6*50

        #x = spikes[np.where(np.logical_and(spikes>=int(self.time_start.text()), spikes<=int(self.time_end.text())))[0]]

        #ymin=np.zeros(len(x))
        #ymax=np.zeros(len(x))
        
        #ymin+=offset-5
        #ymax+=offset-7
        
        #plt.vlines(x-int(self.time_start.text()), ymin, ymax, linewidth=3, color=colors[i%9]) #colors[mod(counter,7)])
    
        #offset=offset-2



    ##********************************************************************************
    ##******************** LABELING ********************
    ##********************************************************************************
    
    #old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    #new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    #plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    

    #old_xlabel = np.linspace(0, float(self.time_end.text()) - float(self.time_start.text()), 5)
    #new_xlabel = np.round(np.linspace(float(self.time_start.text()), float(self.time_end.text()), 5), 1)
    #plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)


    #ax.tick_params(axis='both', which='both', labelsize=font_size)

    #plt.xlabel("Time (sec)", fontsize = font_size , weight = 'bold')        

        
    ##plt.xlabel('Time (seconds)',fontsize=35, weight='bold')
    ##plt.ylabel('Single Unit ID',multialignment='center', fontsize=35)

    ##ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    ##plt.ylabel('Single Unit Rasters     Synchrony Index   Specgram Frequency (Hz)', fontsize=font_size, weight = 'bold')           
    
    #plt.ylabel('LFP   Thresholded Events      Clustered Events', fontsize=font_size, weight = 'bold')           
    ##plt.ylabel('LFP Cluster Raster         Single Unit IDs',multialignment='center', fontsize=35, weight='bold')
    
    #plt.xlim(0, int(self.time_end.text())-int(self.time_start.text()))
    
    plt.show()
    
    

def concatenate_tsf(self):
    """ Function doc """
    
    print "...concatenate multiple .tsf..."
    
    tsf = TSF.TSF('')       #Initialize empty object; 
    
    total_n_vd_samples = 0
    temp_names = []
    temp_n_samples = []
    temp_n_digital_chs = []
    temp_digital_chs = []
    for k in range(len(self.tsf_files)):
        
        if '.tsf' in self.tsf_files[k]:                    #Concatenate standard .tsf files
            temp_tsf = TSF.TSF(self.tsf_files[k])

            print temp_tsf.Readloc
            
            print ">>>>>>>>> SITELOC ORDER CHANGED TO DEAL WITH NICK'S FILES: Make sure it still works for INTAN files <<<<<<<<<<<<<"
            ordered_indexes = np.argsort(temp_tsf.Siteloc[1::2])        #Is this ok for non cat data? ************************
            print ordered_indexes
            temp_tsf.Readloc = np.arange(1, temp_tsf.n_electrodes+1, 1)
            print temp_tsf.Readloc

            print temp_tsf.Siteloc
            
            temp_siteloc = []
            for p in ordered_indexes:
                print temp_tsf.Siteloc[p*2], temp_tsf.Siteloc[p*2+1]
                temp_siteloc.append(temp_tsf.Siteloc[p*2])
                temp_siteloc.append(temp_tsf.Siteloc[p*2+1])
            
            temp_tsf.Siteloc = temp_siteloc
            #print temp_tsf.Siteloc
            
            total_n_vd_samples += temp_tsf.n_vd_samples
            print "...len loaded data: ", temp_tsf.n_vd_samples
            
            temp_tsf.read_footer()                      #Reads the footer and original file_name if available;
            temp_names.extend(temp_tsf.file_names)
            temp_n_samples.extend(temp_tsf.n_samples)
            
            temp_n_digital_chs.extend(temp_tsf.n_digital_chs)   #I think this is a list also
            temp_digital_chs.extend(temp_tsf.digital_chs)


        elif ".lfp.zip" in self.tsf_files[k]:                #Concatenate martin's .lfp.zip files
            
            print ">>>>>>>>>>>>> ERROR  - USE LFP.ZIP FUNCTION INSTEAD <<<<<<<<<<<<<<<"
            temp_tsf = TSF.TSF('')
            temp_tsf.header = ''
            temp_tsf.iformat = 1002
            
            #Load .lfp data from Martin's files
            data = np.load(self.tsf_files[k])
            keys = ['t0', 't1', 'chanpos', 'uVperAD', 'chans', 'tres', 'data']
            
            print "... uVperAD: ", data['uVperAD']

            total_n_vd_samples += data['t1']*1E-3

            temp_names.append(self.tsf_files[k])
            temp_n_samples.append(int(data['t1']*1E-3))
            
            temp_n_digital_chs.extend([0])
            temp_digital_chs.extend([])

            print len(data['data'][0])
            print data['data'].shape

            temp_tsf.SampleFrequency = data['tres']
            temp_tsf.vscale_HP = data['uVperAD']
            temp_tsf.n_electrodes = len(data['chans'])
            temp_tsf.Siteloc = data['chanpos'][data['chans']].ravel()
            print temp_tsf.Siteloc
            ordered_indexes = np.argsort(temp_tsf.Siteloc[1::2])        #Get indexes for vertical ordering
            print ordered_indexes
            
            temp_tsf.Readloc = np.arange(1, temp_tsf.n_electrodes+1, 1)
            
    
    #Header
    tsf.header = temp_tsf.header
    tsf.iformat = temp_tsf.iformat
    tsf.SampleFrequency = temp_tsf.SampleFrequency
    tsf.vscale_HP = temp_tsf.vscale_HP
    tsf.n_vd_samples = int(total_n_vd_samples)
    tsf.Siteloc = temp_tsf.Siteloc
    tsf.n_electrodes = temp_tsf.n_electrodes
    tsf.n_cell_spikes = 0
    tsf.layout = temp_tsf.Readloc - 1   #The order of the channels is the same as Readloc; change to zero-based indices

    #Footer; tsf object doesn't exist above, otherwise could have assigend directly.
    tsf.file_names = temp_names     #These are the original filenames to be preserved;
    tsf.n_samples = temp_n_samples
    tsf.n_digital_chs = temp_n_digital_chs
    tsf.digital_chs = temp_digital_chs
    
    
    print "... channel order layout: ", tsf.layout
    
    #Initialize ec_traces total
    tsf.ec_traces = np.zeros((tsf.n_electrodes,  tsf.n_vd_samples), dtype=np.int16)
    print "... total rec length: ", tsf.ec_traces.shape
    
    #Load each tsf data file 
    tsf_index = 0
    for ctr, file_name in enumerate(self.tsf_files):
        
        print "... loading: ", file_name
        
        if ".tsf" in file_name: 
            temp_tsf = TSF.TSF(file_name)
            temp_tsf.read_ec_traces()
        
        elif ".lfp.zip" in file_name:                #Concatenate martin's .lfp.zip files
            
            data = np.load(self.tsf_files[k])
            keys = ['t0', 't1', 'chanpos', 'uVperAD', 'chans', 'tres', 'data']
                   
            temp_tsf = TSF.TSF('')  #Make empty file
            #temp_tsf.ec_traces = np.zeros((tsf.n_electrodes, int(data['t1']*1E-3)), dtype=np.int16)     #Make empty file; convert from usec to ms which is sample rate timeing.

            temp_tsf.ec_traces=[]
            for p in range(len(data['data'])):
                temp_tsf.ec_traces.append([data['data'][p][0]])     #This data is packed oddly, so need to pick the [0] element of the array            
        
        #for ch in range(len(temp_tsf.ec_traces)):
        for ch, index in enumerate(ordered_indexes):        #Use original order and reorder top down.
            tsf.ec_traces[ch,tsf_index:tsf_index+len(temp_tsf.ec_traces[ch])] = temp_tsf.ec_traces[index]
        
        tsf_index+=len(temp_tsf.ec_traces[ch])
        
    print "...saving alltrack .tsf..."
    
    
    #REFILTER DATA BEFORE SAVING:
    if False: 
        for ch in range(len(tsf.ec_traces)):
            print ".... bandpass filtering, lowcut: ", self.low_cutoff.text(), "   highcut: ", self.high_cutoff.text(), "   channel: ", ch
            temp_trace = tsf.ec_traces[ch]
            
            print tsf.SampleFrequency
            #plt.plot(temp_trace[0:10000], color='blue')
            #temp_trace = butter_bandpass_filter(temp_trace, float(self.low_cutoff.text()), float(self.high_cutoff.text()), tsf.SampleFrequency, order=5)
            temp_trace = butter_lowpass_filter(temp_trace, float(self.high_cutoff.text()), tsf.SampleFrequency, order=5)
            #print temp_trace[0:100]
            #plt.plot(temp_trace[0:10000], color='red')
            #plt.show()
            
            tsf.ec_traces[ch] = temp_trace

        file_name = self.tsf_files[0][:-4]+"_alltrack_lowcut"+self.low_cutoff.text()+"_highcut"+self.high_cutoff.text()+".tsf"
    
    else: 
        file_name = self.tsf_files[0][:-4]+"_alltrack.tsf"
    
  
    tsf.save_tsf(file_name)
        
    print "... done saving..."



    
    
def concatenate_lfp_zip(self):
    """ Function doc """

    print "...concatenate lfp.zip files..."
    
    #This is only for Nick, Martin cat data; for intan data, can make lfp files from raw .tsf files directly
    files_in = np.loadtxt(self.tsf_files[0], dtype=str)

    tsf = TSF.TSF('')       #Initialize empty tsf

    tsf.file_names = []
    tsf.n_samples = []
    tsf.n_digital_chs = []
    tsf.digital_chs = []
    tsf.n_vd_samples = 0

    tsf.ec_traces = [[]]*10
    for ctr, file_name in enumerate(files_in):
        tsf.file_names.append(file_name)
        
        tsf_temp = load_lfp_all(file_name)

        print tsf_temp.chans

        temp_ec_traces=[]
        for ch in range(tsf_temp.n_electrodes):
            temp_ec_traces.append(np.append(tsf.ec_traces[ch],tsf_temp.ec_traces[ch]))
        
        tsf.ec_traces=np.int16(temp_ec_traces)
        tsf.n_vd_samples += tsf_temp.n_vd_samples

        tsf.n_samples.append(tsf_temp.n_vd_samples)
        tsf.n_digital_chs.append(0)
        tsf.digital_chs.append([])

    tsf.iformat = tsf_temp.iformat
    tsf.header = tsf_temp.header    
    tsf.n_cell_spikes = 0
    tsf.Siteloc=np.ravel(tsf_temp.Siteloc[tsf_temp.chans])        #Channel locations saved as flattened x,y coords
    tsf.n_electrodes = tsf_temp.n_electrodes
    tsf.SampleFrequency = tsf_temp.SampleFrequency
    tsf.vscale_HP = tsf_temp.vscale_HP
    
    print tsf.iformat
    print tsf.SampleFrequency
    print tsf.n_electrodes
    print tsf.n_vd_samples, tsf.n_vd_samples/float(tsf.SampleFrequency)
    print tsf.vscale_HP
    print tsf.Siteloc
    tsf.Readloc = np.arange(1,tsf.n_electrodes+1,1)
    tsf.layout = tsf.Readloc-1
    
    tsf.ec_traces = np.int16(tsf.ec_traces*tsf.vscale_HP)
    tsf.vscale_HP = 1.0
    
    
    print ''; print "...saving alltrack .tsf..."
    file_name = files_in[0].replace('.lfp.zip', "_alltrack_lfp.tsf")
    save_tsf_single(tsf, file_name)    
    tsf.save_tsf
    print "... done saving..."

def compress_lfp(self):
    
    print "...making compressed lfp files ..."
    
    print self.tsf_file
    compression_factor = int(self.compress_factor.text())
    print "...compressed factor: ", compression_factor
    
    tsf = TSF.TSF(self.tsf_file)
    tsf.read_footer()

    tsf.read_ec_traces()
    
    tsf.layout = tsf.Readloc-1
    #Leave ADC convertion intact it possible
    #tsf.ec_traces= np.int16(tsf.ec_traces*tsf.vscale_HP)
    #tsf.vscale_HP = 1.0 
    
    if tsf.SampleFrequency != 1000:
        print "...lfp frequency not = 1000Hz... exiting ..."
        return
        
    #*********** SAVE COMPRESSED LOW PASS .TSF FILE *********
    #Save compression file name
    file_out = self.tsf_file[:-4]+'_'+str(compression_factor)+'compressed.tsf'
    print "Saving LFP : ", file_out
    
    #DON"T USE SUBSAMPLING - CAUSES PROBLEMS LATER
    #traces_out = []
    #for k in range(len(tsf.ec_traces)):
    #    #traces_out.append(tsf.ec_traces[k][::int(compression_factor/25)])      
    #    traces_out.append(tsf.ec_traces[k])
    #tsf.ec_traces = np.array(traces_out)

    tsf.SampleFrequency = compression_factor*tsf.SampleFrequency     
    tsf.subsample = 1.0
    #tsf.n_vd_samples = len(tsf.ec_traces[0])
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


def do_butter_highpass_filter(data, a, b):
    return filtfilt(b, a, data)
    
    
def ephys_to_tsf(filenames):
    '''Read .rhd files, convert to correct electrode mapping and save to .tsf file
    NB: There are 2 possible mapping depending on the insertion of the AD converter 
    TODO: implement a wavelet high pass filter directly to avoid SpikeSorter Butterworth filter artifacts
    '''
    
    print "...reading amp data..."
    
    probe = Probe()     #Load 64ch probe layout

    for file_name in filenames:
        #Make empty tsf object and then append to it. 
        tsf = TSF.TSF('')

        print file_name
        #Delete previous large arrays; Initialize arrays; IS THIS REDUNDANT?
        ec_traces = 0.; ec_traces_hp = 0.; data=0.
        
        file_out = file_name[:-4]

        #********** READ ALL DATA FROM INTAN HARDWARE ***********
        print "Processing: \n", file_name
        data = read_data(file_name)

        tsf.SampleFrequency = int(data['frequency_parameters']['board_adc_sample_rate'])

        #********** PROCESS DIGITAL CHANNELS FIRST ***************
        if 'board_dig_in_data' in data.keys():
            
            tsf.n_digital_chs = [len(data['board_dig_in_data'])]
            print "...# digital channels: ", len(data['board_dig_in_data'])
            temp_chs = []
            for ch in range(len(data['board_dig_in_data'])):
                temp_chs.append(np.array(data['board_dig_in_data'][ch],dtype=bool))     #Load digital channel data as boolean;
        
        else:
            print "... no digital channels recorded..."
            tsf.n_digital_chs=[0]
            temp_chs=[]
            
        tsf.digital_chs = []
        tsf.digital_chs.append(temp_chs)
        
        #****** SCALE EPHYS DATA *************
        tsf.ec_traces = data['amplifier_data']  #Data comes out as int32?  NOT SURE, NEED TO CHECK...
        tsf.ec_traces*=10.                      #*10       #Multiply by 10 to increase resolution for int16 conversion

        #***** META DATA FOR HEADER AND FOOTER
        tsf.n_electrodes = len(tsf.ec_traces)
        tsf.header = 'Test spike file '
        tsf.iformat = 1002
        if tsf.n_electrodes >32:
            tsf.layout = probe.layout
        else:
            tsf.layout = np.arange(tsf.n_electrodes)
        tsf.n_vd_samples = len(tsf.ec_traces[0])
        tsf.vscale_HP = 0.1                             #voltage scale factor set manually by scaling data...
        tsf.n_cell_spikes = 0       
        tsf.Siteloc = probe.Siteloc

        #Footer
        tsf.file_names = [file_name]
        tsf.n_samples = [tsf.n_vd_samples]              #Need this to save to footer
        
        #************************* THIS SHOULD BE DONE VIA FUNCTION ************************
        #SAVE RAW DATA
        if True:
            #tsf.ec_traces = tsf.ec_traces_raw
            tsf.save_tsf(file_out+'_raw.tsf')
        
        print tsf.n_electrodes
        print len(tsf.ec_traces)
        
        #SAVE LFP to 250Hz.
        if True: 
            tsf_lfp = TSF.TSF('')
            tsf_lfp.SampleFrequency = 1000   #LFP to be downsampled to 1khz from raw data.
            
            print "...converting raw to .lfp (1Khz) sample rate tsf files ..."
            lowpass_freq = 250      #250Hz for lowpass cutoff
            #lowpass_freq = 100      #250Hz for lowpass cutoff

            temp_traces = []
            temp_siteloc = []
            for k in range(tsf.n_electrodes):
               
                #Butter band pass and subsample to 1Khz simultaneously
                trace_out = tsf.ec_traces[k][::int(tsf.SampleFrequency/1000)]
                
                #print "...lowpass filter ch: ", k
                trace_out = butter_lowpass_filter(trace_out, lowpass_freq, fs=tsf_lfp.SampleFrequency, order = 2)

                #Apply 60Hz Notch filter - NOT ALWAYS REQUIRED; MAKE IT AN OPTION EVENTUALLY
                trace_out = Notch_Filter(trace_out)

                temp_traces.append(trace_out)
            
            print len(temp_traces)
            
            tsf_lfp.header = 'Test spike file '
            tsf_lfp.iformat = 1002
            tsf_lfp.vscale_HP = 0.1                             #voltage scale factor
            tsf_lfp.n_cell_spikes = 0 
            tsf_lfp.Siteloc = tsf.Siteloc

            tsf_lfp.ec_traces = temp_traces
            tsf_lfp.Siteloc = tsf.Siteloc
            tsf_lfp.n_electrodes = len(tsf_lfp.ec_traces)
            tsf_lfp.layout = tsf.layout
            
            tsf_lfp.n_vd_samples = len(tsf_lfp.ec_traces[0])

            #Footer
            tsf_lfp.file_names = tsf.file_names     #Same as original
            tsf_lfp.n_samples = [tsf.n_vd_samples]              #Need this to save to footer
            tsf_lfp.digital_chs = tsf.digital_chs 
            tsf_lfp.n_digital_chs = tsf.n_digital_chs
            
            print tsf_lfp.layout
            print len(tsf_lfp.ec_traces)
            
            tsf_lfp.save_tsf(file_name[:-4]+'_lfp_'+str(lowpass_freq)+'hz.tsf')

            
        ##SAVE HIGH PASS WAVELET FILTERED DATA
        #if False:
            
            #print "Wavelet filtering..."
            #wname='db4'
            #maxlevel = 4
            #tsf.ec_traces = wavelet(tsf.ec_traces_raw, wname, maxlevel)
            #print tsf.ec_traces.shape
            
            #print "Writing hp data wavelet filtered ..."
            #save_tsf_single(tsf, file_out+'_hp_wavelet.tsf')

        
        #SAVE BUTTER PASS FILTERED DATA
        if True:

            print "... highpass filtering data (parallel version)..."
            cutoff = 1000
            order = 2
            temp_traces = []

            nyq = 0.5 * tsf.SampleFrequency
            normal_cutoff = cutoff/nyq
            b, a = butter(order, normal_cutoff, btype='high', analog=False)
            
            #Use parmap
            import parmap
            #filtered_pixels = parmap.map(do_filter, pixels, b, a, processes=30)
            tsf.ec_traces = parmap.map(do_filter, tsf.ec_traces, b, a, processes=16)
            
            print "Writing butter pass filtered data ..."
            tsf.save_tsf(file_out+'_hp_butter.tsf')
                    
    print "... Done converting data to .tsf format ..."

    

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
        #file_out = file_name[:-4].replace('rhd_files','tsf_files')
        file_out = file_name[:-4]
        #if os.path.exists(file_out)==True: continue

        #********** READ ALL DATA FROM INTAN HARDWARE ***********
        print "Processing: \n", file_name
        data = read_data(file_name)

        SampleFrequency = int(data['frequency_parameters']['board_adc_sample_rate']); print "SampleFrequency: ", SampleFrequency

        #****** SCALE EPHYS DATA *************
        ec_traces = data['amplifier_data'] #*10       #Multiply by 10 to increase resolution for int16 conversion
        ec_traces*=10.
        print "...length original traces: ", len(ec_traces)
        
        #print "Converting data to int16..."
        #for k in range(len(ec_traces)):
        #    np.save(file_name+"_ch_"+str(k), np.int16(ec_traces[k]))
        
        n_electrodes = len(ec_traces)
        
        #ec_traces = []
        #for k in range(n_electrodes):
        #    ec_traces.append(np.load(file_name+"_ch_"+str(k)+'.npy'))

        #ec_traces = np.array(ec_traces)
        
        #print "...length reloaded traces: ", len(ec_traces)

        header = 'Test spike file '
        iformat = 1002
        n_vd_samples = len(ec_traces[0]); print "Number of samples: ", n_vd_samples
        vscale_HP = 0.1                             #voltage scale factor
        n_cell_spikes = 0

        #****** PROCESS DIGITAL CHANNELS
        n_digital_chs = len(data['board_dig_in_data'])
        print "...# digital channels: ", n_digital_chs
        
        #np.save(file_name[:-4].replace('rhd_files','camera_files')+'_'+response, data['board_dig_in_data'][ch])

        #SAVE RAW DATA
        if True:
            print "Writing raw data ..."
            fout = open(file_out+'.tsf', 'wb')
            fout.write(header)
            fout.write(struct.pack('i', 1002))
            fout.write(struct.pack('i', SampleFrequency))
            fout.write(struct.pack('i', probe.n_electrodes+n_digital_chs))
            fout.write(struct.pack('i', n_vd_samples))
            fout.write(struct.pack('f', vscale_HP))

            #Save ephys data location
            for i in range(probe.n_electrodes):
                fout.write(struct.pack('h', probe.Siteloc[i][0]))
                fout.write(struct.pack('h', probe.Siteloc[i][1]))
                fout.write(struct.pack('i', i+1))

            #Save additional channel locations
            for i in range(n_digital_chs):
                fout.write(struct.pack('h', 0))
                fout.write(struct.pack('h', 2000+100*i))
                fout.write(struct.pack('i', probe.n_electrodes+i+1))

            #Save Ephys data
            for i in range(probe.n_electrodes):
                print "...writing ch: ", i
                ec_traces[probe.layout[i]].tofile(fout)  #Frontside

            #Save digital channels
            for ch in range(n_digital_chs):
                np.save(file_name[:-4]+'_digitalchannel_'+str(ch), data['board_dig_in_data'][ch])
                temp_data = np.load(file_name[:-4]+'_digitalchannel_'+str(ch)+ '.npy')
                temp_data = np.int16(temp_data*10000.)
                #plt.plot(temp_data)
                #plt.show()
                
                temp_data.tofile(fout)  #Pack data into .tsf file
        
                temp_data.tofile(file_name[:-4]+'_digitalchannel_'+str(ch)+'.bin')  #save digital channel separately.

            fout.write(struct.pack('i', n_cell_spikes))
            fout.close()
            
        ##SAVE HIGH PASS WAVELET FILTERED DATA
        #if True:
            #print "Writing hp data ..."
            #fout = open(file_out+'_hp.tsf', 'wb')
            #fout.write(header)
            #fout.write(struct.pack('i', 1002))
            #fout.write(struct.pack('i', SampleFrequency))
            #fout.write(struct.pack('i', probe.n_electrodes))
            #fout.write(struct.pack('i', n_vd_samples))
            #fout.write(struct.pack('f', vscale_HP))
            
            #for i in range (probe.n_electrodes):
                #fout.write(struct.pack('h', probe.Siteloc[i][0]))
                #fout.write(struct.pack('h', probe.Siteloc[i][1]))
                #fout.write(struct.pack('i', i+1))

            #print "Wavelet filtering..."
            #ec_traces_hp = wavelet(ec_traces, wname="db4", maxlevel=6)
            #print ec_traces_hp.shape

            #for i in range(probe.n_electrodes):
                #print "...writing ch: ", i
                #ec_traces_hp[probe.layout[i]].tofile(fout)  #Frontside

            #fout.write(struct.pack('i', n_cell_spikes))
            #fout.close()

    print "... Done conversion ..."

def rhd_digital_save(file_names):
    '''Read .rhd files, and save digital channels.
    NB: there can be 2, 4 or 6 digital channels inside Intan file
    chs 1 and 2 are laser meta data and laser pulse times (these are off for other experiments)
    chs 3 and 4 are camera pulse times and on/off times from clampx computer (these are chs 1 and 2 usually as laser chs are off)
    chs 5 and 6 are vis stim pulse times and meta data (these are usually 3 and 4 as not recorded w. laser on)
    '''
    
    print "...reading digital amp data..."

    
    #if os.path.exists(camera_onoff_filename+'.npy')==True: continue
    for file_name in file_names:

        data = read_data(file_name)

        print "...# digital channels: ", len(data['board_dig_in_data'])

        SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
        print "SampleFrequency: ", SampleFrequency
        
        #from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
        #from matplotlib.figure import Figure

        #plt.close()
        #fig = Figure(figsize=(5,5), dpi=100)
        #ax1 = fig.add_subplot(111)
        

        
        for ch in range(len(data['board_dig_in_data'])):
            #plt.plot(np.array(data['board_dig_in_data'][ch][:1000000],dtype=bool))
            #plt.show()
            
            #ax1.plot(data['board_dig_in_data'][ch][:100000])
            #plt.show()

            #response = raw_input("Please enter filename extension: ")
            np.save(file_name[:-4]+'_channel_'+str(ch), np.array(data['board_dig_in_data'][ch],dtype=bool))
            


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
    
    
    tsf = TSF.TSF(self.selected_recording)
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

    

def tsf_to_lfp(self):
    '''Read .tsf files - subsample to 1Khz, save as *_lp.tsf
    '''

    lfp_samplerate = 1000   #Select LFP sample rate for export
    
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    
    print "...making low-pass tsf files (1Khz sample rates)..."

    for file_name in self.tsf_files:
        
        file_out = file_name[:-4]+'_lp.tsf'
        if os.path.exists(file_out)==True: continue

        print "Processing: \n", file_name

        tsf = TSF.TSF(file_name)
        tsf.read_ec_traces()
        print tsf.Siteloc.shape
        
        n_vd_samples = len(tsf.ec_traces[0]); print "Number of samples: ", n_vd_samples
        
        print "...converting raw to .lfp (1Khz) sample rate tsf files ..."
        temp_traces = []
        temp_siteloc = []
        for k in range(0, tsf.n_electrodes, electrode_rarifier):

            #Check depth of electrode; skip depth past 900
            if tsf.Siteloc[k*2+1]>1000: continue    #Skip electrodes deeper than 900um
           
            temp_siteloc.append(tsf.Siteloc[k*2]); temp_siteloc.append(tsf.Siteloc[k*2+1])
            
           
            #Butter band pass and subsample to 1Khz simultaneously
            trace_out = tsf.ec_traces[k][::int(tsf.SampleFrequency/1000)]
            
            if (float(self.low_cutoff.text())!=0) and (float(self.high_cutoff.text())!=0):
                print "...bandpass filter..."
                trace_out = np.int16(butter_bandpass_filter(trace_out, float(self.low_cutoff.text()), float(self.high_cutoff.text()), fs=lfp_samplerate, order = 2))

            if (float(self.low_cutoff.text())!=0) and (float(self.high_cutoff.text())==0):
                print "...highpass filter..."
                trace_out = np.int16(butter_highpass_filter(trace_out, float(self.low_cutoff.text()), fs=lfp_samplerate, order = 2))

            #Apply 60Hz Notch filter - NOT ALWAYS REQUIRED; MAKE IT AN OPTION EVENTUALLY
            #temp_traces.append(Notch_Filter(temp))
            temp_traces.append(trace_out)
            
        tsf.ec_traces = np.int16(temp_traces)
        tsf.Siteloc = np.int16(temp_siteloc)
        tsf.n_electrodes = len(tsf.ec_traces)

        tsf.n_vd_samples = len(tsf.ec_traces[0])
        tsf.SampleFrequency = lfp_samplerate
        
        #Save data to .tsf file
        tsf.save_tsf(file_name[:-4]+'_lfp_'+self.low_cutoff.text()+"hz_"+self.high_cutoff.text()+'hz.tsf')


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
    


def Notch_Filter(data, fs=1000, band=1., freq=60., ripple=10, order=4, filter_type='ellip'):
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
    
    #filtered_data = filtfilt(b,a,data)     #Try thsi!?
    
    return filtered_data


def view_templates(self):
    print "..."

    top_channel = int(np.loadtxt(os.path.split(os.path.split(self.selected_sort)[0])[0]+"/top_channel.txt") - 1)      #Load top channel for track; convert to 0-based ichannel values.
    
    bad_channels = np.loadtxt(os.path.split(os.path.split(self.selected_sort)[0])[0]+"/electrode_mask.txt", dtype=np.int32)
    
    if 0 in bad_channels:
        bad_channels = []
       

    font_size = 20
    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())

    
    #compression = 50.   #Need this to convert from compressed sample points to realtime
    
    print self.selected_sort

    self.tsf = TSF.TSF(self.selected_recording)
    self.tsf.read_ec_traces()

    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
            print "...filtering ch: ", k
            self.tsf.ec_traces[k] = butter_bandpass_filter(self.tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
    #load single units
    Sort = Ptcs(self.selected_sort)
    maxchan = Sort.maxchan[int(self.selected_unit.text())]
    print "... maxchan: ", maxchan
    
    ax = plt.subplot(1,2,1)
    t = np.arange(-n_samples,n_samples+1,1)
    
    """ Spikes are saved in # of sample points so no need to scale them up from compressed .ptcs sort file to uncompressed lfp file.
    """
    spikes = Sort.units[int(self.selected_unit.text())]*1E-6 * Sort.samplerate
    
    print self.tsf.Siteloc
    

    max_chan_traces = []
    for k in range(top_channel, self.tsf.n_electrodes, electrode_rarifier):
        print "...plotting ch: ", k, "  depth: ", -self.tsf.Siteloc[k*2+1]
        
        ch_offset = 0
        if k in bad_channels:   #SKIP BAD CHANNELS AND SELECT NEXT CHANNEL
            ch_offset=1
        
        if k !=maxchan:continue
        
        self.tsf.ec_traces[k+ch_offset] = butter_highpass_filter (self.tsf.ec_traces[k+ch_offset], 4.0, 1000., 5)

        traces = []
        for spike in spikes:
            #print spike
            trace_out = self.tsf.ec_traces[k+ch_offset][int(spike-n_samples):int(spike+n_samples+1)]*self.tsf.vscale_HP*voltage_scaling
            
            if len(trace_out)!=(n_samples*2+1):     #If end/begining of potentials;
                continue
            
            traces.append(trace_out)

            if (k+ch_offset) == maxchan:     #LFP STABILITY METRICS
                max_chan_traces.append(trace_out)
                
                #maxchan_peak.append(np.argmin(trace_out[n_samples-40:n_samples+40]))
                #maxchan_trough.append(np.argmin(trace_out[n_samples-40:n_samples+40]))
                #maxchan_trough.append(np.argmin(trace_out[n_samples:n_samples+60]))
                #maxchan_trough.append(np.argmin(trace_out[n_samples-100:n_samples]))
        
        print len(traces)
        
        trace_ave = np.average(traces, axis=0)
        trace_std = np.std(traces, axis=0)
       
        #offset = -voltage_scaling*self.tsf.Siteloc[k*2+1]
        offset = -self.tsf.Siteloc[k*2+1]
        
        if k == maxchan:
            plt.plot(t, trace_ave+offset, color='black', linewidth=5, alpha=1)
        else: 
            plt.plot(t, trace_ave+offset, color='black', linewidth=3, alpha=0.25)
        
        ax.fill_between(t, trace_ave+trace_std+offset, trace_ave-trace_std+offset, color=self.selected_colour, alpha=0.2)

        plt.plot([-n_samples,n_samples], [offset, offset], color='black', linewidth=2, alpha=.35)

    #maxchan_troughs = np.array(maxchan_trough)+20 #Centre data on exported times
    #print maxchan_troughs
    
    plt.plot([t[-1]+50,t[-1]+50], [-500, -400], color='black', linewidth=4)
    plt.plot([0,0], [150, -5000], 'r--', color='black', linewidth=2, alpha=.5)

    #Set ylabel
    #old_ylabel = -voltage_scaling*np.arange(0, np.max(self.tsf.Siteloc)+1, 100)
    #new_ylabel = np.int16(np.arange(0, np.max(self.tsf.Siteloc)+1, 100))
    #plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    #plt.locator_params(axis='y',nbins=5)
    
    plt.ylabel("Depth (um)", fontsize=font_size)

    #Set xlabel
    old_xlabel = np.linspace(t[0],t[-1],3)
    new_xlabel = np.linspace(float(t[0])/self.tsf.SampleFrequency*1E3, float(t[-1])/self.tsf.SampleFrequency*1E3, 3)
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)

    plt.xlabel("Time (ms)", fontsize = font_size)
    plt.tick_params(axis='both', which='both', labelsize=font_size)
    plt.locator_params(axis='x',nbins=10)

    plt.ylim(-2400, 150)
    plt.xlim(old_xlabel[0], old_xlabel[-1]*5)
    
    plt.title("#: "+ str(len(spikes)), fontsize = font_size)
    
    
    #***************************************************************************************************
    #************************************* PLOT LFP SHAPE STABILITY HISTOGRAMS ****************************
    #***************************************************************************************************


    #Plot trough histogram
    ax = plt.subplot(1,2,2)
    #plt.title("STD: "+ str(np.round(np.std(maxchan_troughs),2))+"ms.", fontsize = font_size)

    if True: 
        
        trace_ave = np.mean(max_chan_traces, axis=0)
        trough_loc = np.argmin(trace_ave)
        peak_loc = np.argmax(trace_ave)
        trough_time = []
        peak_time = []
        ptp_amplitude = []
       
        for k in range(len(max_chan_traces)): 
            trough_temp = np.argmin(max_chan_traces[k][trough_loc-25:trough_loc+25])+trough_loc-25
            trough_time.append(n_samples - trough_temp)
            
            peak_temp = np.argmax(max_chan_traces[k][peak_loc-25:peak_loc+25])+peak_loc-25
            peak_time.append(n_samples - peak_temp)
            
            ptp_amplitude.append(max_chan_traces[k][peak_temp]-max_chan_traces[k][trough_temp])
                
        bin_width = 5   # histogram bin width in usec
        print "...std time trough: ", np.std(trough_time)
        std_trough = np.std(trough_time)
        std_peak = np.std(peak_time)
        if std_trough<std_peak: 
            y = np.histogram(trough_time, bins = np.arange(-n_samples,n_samples,bin_width))
            plt.bar(y[1][:-1], y[0], bin_width, color=self.selected_colour, alpha=1)

        else: 
            print "...std time peak: ", np.std(peak_time)
            y = np.histogram(peak_time, bins = np.arange(-n_samples,n_samples,bin_width))
            plt.bar(y[1][:-1], y[0], bin_width, color=self.selected_colour, alpha=1)
        
        std_overall = np.round(min(np.std(trough_time),np.std(peak_time)),1)
        

        #if False: 
            #trough_time = []
            #peak_time = []
            #fwhm_trough = []
            #fwhm_peak = []
            
            #for k in range(len(max_chan_traces)): 
                #trough_temp = np.argmin(max_chan_traces[k][trough_loc-25:trough_loc+25])
                #trough_time.append(trough_temp)
                
                #peak_temp = np.argmax(max_chan_traces[k][peak_loc-25:peak_loc+25])
                #peak_time.append(peak_temp)
                
                ##FIND TROUGH FWHM
                #fwhm = []
                #for p in range(trough_temp, 0, -1):
                    #if (max_chan_traces[k][p]< max_chan_traces[k][trough_temp]/2.) and (max_chan_traces[k][p-1]> max_chan_traces[k][trough_temp]/2.):
                        #fwhm.append(p)
                        #break
                
                #for p in range(trough_temp, n_samples*2, 1):
                    #if (max_chan_traces[k][p]< max_chan_traces[k][trough_temp]/2.) and (max_chan_traces[k][p-1]> max_chan_traces[k][trough_temp]/2.):
                        #fwhm.append(p)
                        #break
                #print fwhm
                #if len(fwhm)>1:
                    #fwhm_trough.append(fwhm[1]-fwhm[0])
                
                ##FIND PEAK FWHM
                #fwhm = []
                #for p in range(peak_temp, 0, -1):
                    #if (max_chan_traces[k][p]> max_chan_traces[k][peak_temp]/2.) and (max_chan_traces[k][p-1]< max_chan_traces[k][peak_temp]/2.):
                        #fwhm.append(p)
                        #break
                
                #for p in range(peak_temp, n_samples*2, 1):
                    #if (max_chan_traces[k][p]> max_chan_traces[k][peak_temp]/2.) and (max_chan_traces[k][p+1]< max_chan_traces[k][peak_temp]/2.):
                        #fwhm.append(p)
                        #break
                #print fwhm
                
                #if len(fwhm)>1: 
                    #fwhm_peak.append(fwhm[1]-fwhm[0])
            
            #fwhm_trough = np.float32(fwhm_trough)
            #fwhm_peak = np.float32(fwhm_peak)
            
            #bin_width = 2   # histogram bin width in usec
            #print "...fwhm_peak std: ", np.std(fwhm_peak)
            #y = np.histogram(fwhm_peak, bins = np.arange(0,n_samples*2,bin_width))
            ##plt.bar(y[1][:-1]*1E-3, np.float32(y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
            #plt.bar(y[1][:-1], y[0], bin_width, color='green', alpha=.5)


            #bin_width = 2   # histogram bin width in usec
            #print "...fwhm_trough std: ", np.std(fwhm_trough)
            #y = np.histogram(fwhm_trough, bins = np.arange(0,n_samples*2,bin_width))
            ##plt.bar(y[1][:-1]*1E-3, np.float32(y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
            #plt.bar(y[1][:-1], y[0], bin_width, color='magenta', alpha=.5)




        plt.xlabel("Time (ms)", fontsize = font_size*1.5)
        plt.tick_params(axis='both', which='both', labelsize=font_size*1.5)
        plt.title("Unit: "+self.selected_unit.text()+ "  STD: " + str(std_overall)+ 'ms \n'+os.path.split(self.selected_recording)[1], fontsize=font_size)
        plt.yticks([])
        
        plt.xlim(-0,100)

    plt.show()
    
    return
    #PLOT 
    
    plt.plot(peak_time, color='red')
    plt.plot(trough_time, color='blue')
    plt.plot(ptp_amplitude, color='green', linewidth=3)
    plt.show()
    

def view_traces(self):
    """ Display raw traces
    """
    
    colors=['gray']
    colors=['blue', 'red', 'green', 'brown', 'pink', 'magenta', 'orange', 'black', 'cyan']
    
    compression_factor = 50
    font_size = 30
    voltage_scaling = float(self.voltage_scale.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    
    t0 = int(float(self.start_time.text())*1000)
    t1 = int(float(self.end_time.text())*1000)
    t = np.arange(t0, t1, 1)*1E-3

    #If loading compressed LFP file, uncompress the .tsf file for display in milisecond time steps:
    if "compressed" in self.selected_recording:
        self.tsf.SampleFrequency = 1000 
        #electrode_rarifier = 1              #Plot all electrodes
        
    #Sumbsample data if raw sampling rate being used
    elif self.tsf.SampleFrequency != 1000:
        
        subsample = int(self.tsf.SampleFrequency/1000)
        temp = []
        for k in range(self.tsf.n_electrodes):
            temp.append(self.tsf.ec_traces[k][::subsample])
    
        self.tsf.SampleFrequency = 1000
        self.tsf.ec_traces = np.array(temp)            
                
    print "... .tsf samplerate: ", self.tsf.SampleFrequency

    
    #Plot traces; subsample traces as per electrode-rarifier

    ax = plt.subplot(1,1,1)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    for k in range(0, self.tsf.n_electrodes, electrode_rarifier):
        if 'tim' in self.selected_recording:            #In 64 channel mouse recordings skip electrodes that are too deep
            if self.tsf.Siteloc[k*2+1]>1000: continue    

        trace_out = self.tsf.ec_traces[k][t0:t1]
        if (float(self.low_cutoff.text())!=0) and (float(self.high_cutoff.text())!=0):
            print "...bandpass filter..."
            trace_out = butter_bandpass_filter(trace_out, float(self.low_cutoff.text()), float(self.high_cutoff.text()), fs=1000, order = 2)

        if (float(self.low_cutoff.text())!=0) and (float(self.high_cutoff.text())==0):
            print "...highpass filter..."
            trace_out = butter_highpass_filter(trace_out, float(self.low_cutoff.text()), fs=1000, order = 2)
        
        print "...len(trace_out): ", len(trace_out), "t0: ", t0, "   t1: ", t1
        plt.plot(t, trace_out-voltage_scaling*self.tsf.Siteloc[k*2+1], color='black', linewidth=2)
        #plt.plot(trace_out-voltage_scaling*self.tsf.Siteloc[k*2+1], color='black', linewidth=2)


    #Plot 
    if True: 
        
        #colors=['gray']*100
        for k in range(len(self.Sort.units)):
            spikes = self.Sort.units[k]*1E-6 *compression_factor #X-axis is in seconds

            ymin = np.zeros(len(spikes),dtype=np.float32)+self.Sort.chanpos[self.Sort.maxchan[k]][1]*(-voltage_scaling)+600
            ymax = np.zeros(len(spikes),dtype=np.float32)+self.Sort.chanpos[self.Sort.maxchan[k]][1]*(-voltage_scaling)-600

            plt.vlines(spikes, ymin, ymax, color=colors[k], linewidth=10, alpha=0.35)
        
        
    #Set labels
    depth_offset = float(self.probe_penentration.text()) #Percent of probe inserted into brain
    
    #Setup new y-axis labels
    
    if 'tim' in self.selected_recording:            #In 64 channel mouse recordings skip electrodes that are too deep
        depth = 900
    else:
        depth = np.max(self.tsf.Siteloc)

    old_ylabel = -voltage_scaling*np.arange(depth - depth_offset*depth,depth+1, 100)
    #new_ylabel = np.int16(np.arange(0, depth_offset*np.max(self.tsf.Siteloc)+1, 100))
    new_ylabel = np.int16(np.arange(0, depth_offset*depth+1, 100))

    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)

    #plt.locator_params(axis='y',nbins=5)
    plt.ylabel("Depth (um)", fontsize=font_size)
    
        
    #Setup new x-axis labels
    old_xlabel = np.arange(t[0], t[-1], 0.5)
    new_xlabel = old_xlabel/50.
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)



    #Set xlabel
    plt.xlabel("Time (sec)", fontsize = font_size)
    plt.tick_params(axis='both', which='both', labelsize=font_size)
    plt.locator_params(axis='x',nbins=10)

    #Plot scale bar
    plt.plot([t[0], t[0]], [-50,-250], color='black', linewidth=10, alpha=1)

    #plot limits
    plt.ylim(-voltage_scaling*depth*1.1, -voltage_scaling*(-25))
    plt.xlim(t[0], t[-1])

    plt.title("Lowcut: " + self.low_cutoff.text() + "Hz     Highcut: " + self.high_cutoff.text()+"Hz.", fontsize=font_size)

    plt.show()



    #Plot location of imaging frames
    #plt.plot(t, -1500*self.camera_pulses[t0:t1], color='blue')


def view_all_csd_old(self): 

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
    tsf = TSF.TSF(self.selected_recording)
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


    tsf = TSF.TSF(self.selected_recording) 
    tsf.read_ec_traces()
    
    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        if os.path.exists(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')==False:
            for k in range(0, tsf.n_electrodes, electrode_rarifier):
                print "...filtering ch: ", k
                tsf.ec_traces[k] = butter_bandpass_filter(tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
            np.save(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut', tsf.ec_traces)
        else: 
            tsf.ec_traces = np.load(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')

    #load single units
    Sort = PTCS.PTCS(self.selected_sort)

    #Load spiketimes as event_triggers 
    event_times = Sort.units[int(self.selected_unit.text())]*1E-6 * Sort.samplerate #Convert from usec to sec back to sample-rate time
    print event_times
    
    lfp_ave = np.zeros((len(tsf.ec_traces),2*n_samples),dtype=np.float32)
    for ch in range(tsf.n_electrodes):
        ctr=0
        for event in event_times:
            trace_temp = tsf.ec_traces[ch][int(event-n_samples):int(event+n_samples)]*tsf.vscale_HP
            if len(trace_temp)==(n_samples*2):
                lfp_ave[ch]+= trace_temp
                ctr+=1
        lfp_ave[ch]=lfp_ave[ch]/ctr

    lfp_ave = lfp_ave[::2]      #Limit analysis to just one column of Neuronexus probe

    lfp_ave = lfp_ave[int(self.start_ch.text()):]      #Exclude top 10 channels

    np.save(self.selected_recording[:-4]+"_ave_lfp", lfp_ave)
    
    #*****************PLOT LFP 
    fig = plt.figure()
    ax = plt.subplot(1,2,1)
    vmax = np.max(np.abs(lfp_ave)); vmin=-vmax
    cax = ax.imshow(lfp_ave, vmin=vmin, vmax=vmax, aspect='auto', cmap='viridis_r', interpolation='none')

    old_ylabel = np.arange(0, len(lfp_ave), 5)
    new_ylabel = np.arange(0, len(lfp_ave)*46, 5*46)
    plt.yticks(old_ylabel, new_ylabel, fontsize=15)
    plt.ylabel("Depth (um)", fontsize = 20)
    
    old_xlabel = np.linspace(0, len(lfp_ave[0]), 5)
    new_xlabel = np.int32(np.linspace(-n_samples, n_samples, 5))
    plt.xticks(old_xlabel, new_xlabel, fontsize=15)
    plt.xlabel("Time (ms)", fontsize = 20)

    print "... lfp max: ", vmax

    cbar = fig.colorbar(cax)#, ticks=[x_ticks])
    #cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=15) 
    cbar.ax.set_ylabel("LFP (uV)", fontsize=15, labelpad=-6)

    #******************COMPUTE CSD
    lfp_ave_dif1 = np.gradient(lfp_ave)
    lfp_ave_dif1 = lfp_ave_dif1[0]
    
    lfp_ave_dif2 = np.gradient(lfp_ave_dif1)
    lfp_ave_dif2 = lfp_ave_dif2[0]*1E2      #***************CONVERSION TO CSD UNITS; matched from Gatue's CSD function
    
    lfp_ave_dif2 = -1. * lfp_ave_dif2     #Need to invert the currenst; matched from Gaute's CSD function
   
    #*****************PLOT CSD:
    ax = plt.subplot(1,2,2)
    
    vmax = np.max(np.abs(lfp_ave_dif2)); vmin=-vmax
    cax = ax.imshow(lfp_ave_dif2, vmin=vmin, vmax=vmax, aspect='auto',cmap='viridis_r', interpolation='none')
    
    old_xlabel = np.linspace(0, len(lfp_ave[0]), 5)
    new_xlabel = np.int32(np.linspace(-n_samples, n_samples, 5))
    plt.xticks(old_xlabel, new_xlabel, fontsize=15)
    plt.xlabel("Time (ms)", fontsize = 20)

    old_ylabel = np.arange(0, len(lfp_ave), 5)
    new_ylabel = np.arange(0, len(lfp_ave)*46, 5*46)
    plt.yticks(old_ylabel, new_ylabel, fontsize=15)
    

    print "... CSD max: ", vmax

    cbar = fig.colorbar(cax)#, ticks=[x_ticks])
    #cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=15) 
    cbar.ax.set_ylabel("CSD (A/m^3)", fontsize=15, labelpad=-6)

    plt.suptitle(os.path.split(self.selected_recording)[1][:-4] + ",  Cluster #: "+ self.selected_unit.text()+",  #events: "+str(len(Sort.units[int(self.selected_unit.text())])), fontsize=20)
    plt.show()
        
    

    print "... done..."
    return

    #****************** SELECTING PROBE CHANNELS ***********************
    if tsf.n_electrodes >10: 
        print "...loading every other channel, NN A64 probe ..."
        probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())][::2]
        lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())][::2]
    else:
        probe_layout = tsf.Siteloc[1::2][int(self.start_ch.text()):int(self.end_ch.text())]
        lfp_ave=lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())]

    print probe_layout
    
    plt.show()



    #************************ SERGEY'S CSD FUNCTIONS *****************************************
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


def view_all_csd(self):
    
    n_samples = int(self.n_sample_pts.text())
    electrode_rarifier = int(1./float(self.n_electrodes.text()))
    voltage_scaling = float(self.voltage_scale.text())


    tsf = TSF.TSF(self.selected_recording) 
    tsf.read_ec_traces()
    
    #Remove bottom power..
    if self.low_cutoff.text()!='0.0':
        if os.path.exists(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')==False:
            for k in range(0, tsf.n_electrodes, electrode_rarifier):
                print "...filtering ch: ", k
                tsf.ec_traces[k] = butter_bandpass_filter(tsf.ec_traces[k], float(self.low_cutoff.text()), 240., fs=1000, order = 2)
        
            np.save(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut', tsf.ec_traces)
        else: 
            tsf.ec_traces = np.load(self.selected_recording+'_'+self.low_cutoff.text()+'lowcut.npy')


    #Load top channels and channel masks
    path_dir, filename = os.path.split(self.selected_recording)
    root_dir, filename = os.path.split(path_dir)

    electrode_mask = np.loadtxt(root_dir+"/electrode_mask.txt", dtype='str')

    #= []
    #if temp_data != 0.0:
    #    for k in range(len(temp_data)):
    #        electrode_mask.append(temp_data[k])

    print electrode_mask

    top_channel = np.int16(np.loadtxt(root_dir+"/top_channel.txt"))-1       #Use 1-based indexes
    print top_channel
    
    #load single units
    Sort = PTCS.PTCS(self.selected_sort)
    fig = plt.figure()
    
    n_units_morethan_100spikes = 0
    for k in range(Sort.n_units):
        print len(Sort.units[k])
        if len(Sort.units[k])>100:
            n_units_morethan_100spikes= n_units_morethan_100spikes +1
          
    print n_units_morethan_100spikes
    #Load spiketimes as event_triggers 
    unit_ctr=-1
    for unit in range(Sort.n_units): 
        if len(Sort.units[unit])<100: 
            print "skipping unit"
            continue
            
        unit_ctr=unit_ctr+1
        print unit_ctr

        event_times = Sort.units[unit]*1E-6 * Sort.samplerate #Convert from usec to sec back to sample-rate time
    
        lfp_ave = np.zeros((len(tsf.ec_traces),2*n_samples),dtype=np.float32)
        for ch in range(top_channel, tsf.n_electrodes, 1):
            ctr=0
            
            print str(ch+1), electrode_mask
            if str(ch+1) in electrode_mask:
                correct_ec_traces = (tsf.ec_traces[ch-1]+tsf.ec_traces[ch+1])/2.
                print "masked channel"
            else:
                correct_ec_traces = tsf.ec_traces[ch]

            for event in event_times:
                trace_temp = correct_ec_traces[int(event-n_samples):int(event+n_samples)]*tsf.vscale_HP
                if len(trace_temp)==(n_samples*2):
                    lfp_ave[ch]+= trace_temp
                    ctr+=1
            lfp_ave[ch]=lfp_ave[ch]/ctr

        #Limit analysis to just one column of Neuronexus probe
        #if ('16_07_26' in self.selected_recording) or ('16_08_31' in self.selected_recording): #This probe had a lot of bad channels at the top - and they seem to be in the odd channel data;
        #    lfp_ave = lfp_ave[1::2]      
        #else:
        lfp_ave = lfp_ave[top_channel:]     #Only look at the bottom channels - if probe not inserted entire way;
        if "nick" in self.selected_recording:
            lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())]
        else:           #Select only a single column
            lfp_ave = lfp_ave[int(self.start_ch.text()):int(self.end_ch.text())][::2]

        np.save(self.selected_recording[:-4]+"_unit"+str(unit)+"_ave_lfp", lfp_ave)
        
        #*****************PLOT LFP 
        ax = plt.subplot(n_units_morethan_100spikes, 2, unit_ctr*2+1)
        vmax = np.max(np.abs(lfp_ave)); vmin=-vmax
        cax = ax.imshow(lfp_ave, vmin=vmin, vmax=vmax, aspect='auto', cmap='plasma', interpolation='none')

        old_ylabel = np.arange(0, len(lfp_ave), 5)
        new_ylabel = np.arange(0, len(lfp_ave)*46, 5*46)
        plt.yticks(old_ylabel, new_ylabel, fontsize=15)
        plt.ylabel(str(unit)+" #: "+ str(len(Sort.units[unit]))+"\nDepth (um)", fontsize = 15)
        
        if unit==(n_units_morethan_100spikes-1):
            old_xlabel = np.linspace(0, len(lfp_ave[0]), 5)
            new_xlabel = np.int32(np.linspace(-n_samples, n_samples, 5))
            plt.xticks(old_xlabel, new_xlabel, fontsize=15)
            plt.xlabel("Time (ms)", fontsize = 20)
        else:
            plt.xticks([])
        
        if unit==0: 
            plt.title("LFP", fontsize = 20)

        print "... lfp max: ", vmax

        cbar = fig.colorbar(cax)#, ticks=[x_ticks])
        #cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=15) 
        cbar.ax.set_ylabel("LFP (uV)", fontsize=15)

        #******************COMPUTE CSD
        lfp_ave_dif1 = np.gradient(lfp_ave)
        lfp_ave_dif1 = lfp_ave_dif1[0]
        
        lfp_ave_dif2 = np.gradient(lfp_ave_dif1)
        lfp_ave_dif2 = lfp_ave_dif2[0]*1E2      #***************CONVERSION TO CSD UNITS; matched from Gaute's CSD function
        
        lfp_ave_dif2 = -1. * lfp_ave_dif2     #***************************Need to invert the currents; matched from Gaute's CSD function
       
        #*****************PLOT CSD:
        ax = plt.subplot(n_units_morethan_100spikes,2,unit_ctr*2+2)
        
        vmax = np.max(np.abs(lfp_ave_dif2)); vmin=-vmax
        cax = ax.imshow(lfp_ave_dif2, vmin=vmin, vmax=vmax, aspect='auto',cmap='plasma', interpolation='none')
        
        if unit==(n_units_morethan_100spikes-1):
            old_xlabel = np.linspace(0, len(lfp_ave[0]), 5)
            new_xlabel = np.int32(np.linspace(-n_samples, n_samples, 5))
            plt.xticks(old_xlabel, new_xlabel, fontsize=15)
            plt.xlabel("Time (ms)", fontsize = 20)
        else:
            plt.xticks([])

        if unit==0: 
            plt.title("CSD", fontsize = 20)

        old_ylabel = np.arange(0, len(lfp_ave), 5)
        new_ylabel = np.arange(0, len(lfp_ave)*46, 5*46)
        plt.yticks(old_ylabel, new_ylabel, fontsize=15)
        

        print "... CSD max: ", vmax

        cbar = fig.colorbar(cax)#, ticks=[x_ticks])
        #cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
        cbar.ax.tick_params(labelsize=15) 
        cbar.ax.set_ylabel("CSD (A/m^3)", fontsize=15)


    rois = ['visual', 'barrel', 'auditory']
    for roi in rois:
        if roi in self.selected_recording:
            plt.suptitle(os.path.split(self.selected_recording)[1][:-4] + "    *"+roi+"*", fontsize=20)
            break

    plt.show()
        
    print "... done..."
    return


def Multitaper_specgram(time_series, sampfreq):
    
    from libtfr import * #hermf, dpss, trf_spec <- these are some of the fucntions; for som reason last one can't be explicitly imported

    
    plotting=False
    if plotting:
        ax2 = plt.subplot(1,1,1)
        ax2.autoscale(enable=True, tight=True)
    
    #**************************************
    
    s = time_series
    print "Length of recording: ", len(s)

    #******************************************************
    #Plot multi taper specgram

    #General parametrs
    f0=0.1
    f1=100
    sampfreq=1000

    #time parameters
    t0 = None
    t1 = None
    ts = np.arange(0,len(s),1.0)/sampfreq

    if t0 == None:
        t0, t1 = ts[0], ts[-1] # full duration
    if t1 == None:
        t1 = t0 + 10 # 10 sec window
            
    #Multi taper function parameters
    N = 512     #Number of tapers 
    NW = 40     #
    step = 10   #shift
    k = 6       #
    tm = 6.0    #time support of the 
    Np = 201    #

    w = np.hamming(N)

    h,Dh,Th = hermf(N, k, tm)
    E,V   = dpss(N, NW, k)
    
    print "Computing trf_specgram..."
    spec = tfr_spec(s, N, step, Np, k, tm)
    print "..."
    #OTHER TYPES OF FUNCTIONS... Did not check in detail
    #mpsd  = mtm_psd(s, NW)
    #J     = mtfft(s, NW)
    #spec  = stft(s, w, step)
    #mspec = mtm_spec(s, N, step, NW)
    #tspec_zoom = tfr_spec(s, N, step, Np, k, tm, fgrid=np.linspace(0.1,0.475,512))
    #tspec_log = tfr_spec(s, N, step, Np, k, tm, fgrid=log_fgrid(0.1, 0.45, 256))
    

    extent = t0, t1, f0, f1

    lo = 0
    hi = int(float(N)/1000. * f1)
    print lo, hi

    spec = spec[lo:hi]

    zis = np.where(spec == 0.0) # row and column indices where P has zero power
    if len(zis[0]) > 0: # at least one hit
        spec[zis] = np.finfo(np.float64).max # temporarily replace zeros with max float
        minnzval = spec.min() # get minimum nonzero value
        spec[zis] = minnzval # replace with min nonzero values
    spec = 10. * np.log10(spec) # convert power to dB wrt 1 mV^2?

    p0=-40
    p1=None
    if p0 != None:
        spec[spec < p0] = p0
    if p1 != None:
        spec[spec > p1] = p1

    if plotting:
        im = ax2.imshow(spec[::-1], extent=extent, aspect='auto', cmap=None) #, vmin=0, vmax=1e2)
        plt.show()

    return spec[::-1], extent


def plot_all_rasters(self):
    """ Function doc """
    

    colors=['blue','red', 'green', 'violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    channel=int(self.specgram_ch.text())
    top_channel = np.loadtxt(os.path.split(os.path.split(self.selected_sort_sua)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.

    tsf = TSF.TSF(self.selected_recording)
    tsf.read_ec_traces()
    print "... len of rec: ", tsf.n_vd_samples/tsf.SampleFrequency
        
    #Find deepest channel by parsing sorted units:
    ax = plt.subplot(111)
       
    sync_0 = 0
    sync_1 = -10
    sync_scale = 20
    font_size = 30
    
    #******************* PLOT SUA RASTERS ***************************
    #Load SUA Sort

    Sort_sua = PTCS.PTCS(self.selected_sort_sua) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    if False: 
        n_spikes = []
        for k in range(len(Sort_sua.units)):
            n_spikes.append(len(Sort_sua.units[k]))
        n_spikes = np.array(n_spikes)
        indexes = np.argsort(n_spikes)
        print indexes
    else:
        indexes = np.arange(Sort_sua.n_units)

    offset = sync_0-10
    
    y = []
    mua = []
    for i in indexes: #range(len(Sort_sua.units)):
        print "... unit: ", i
    #for i in indexes[0:5]: #range(len(Sort_sua.units)):
        #x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)/float(Sort_sua.samplerate) #float(Sort1.samplerate)*2.5
        spikes = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)*1E-6

        x = spikes[np.where(np.logical_and(spikes>=float(self.time_start.text()), spikes<=float(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=offset+0.2
        ymax+=offset-0.2

        plt.vlines(x-float(self.time_start.text()), ymin, ymax, linewidth=.5, color='black') #colors[mod(counter,7)])

        y.append(x)
        
        offset=offset-0.5

        mua.extend(x-float(self.time_start.text()))
        
    #******************** PLOT MUA HISTOGRAM *********************************
    if False: 
        #mua=[]
        #for i in range(len(Sort_sua.units)):
        #    mua.extend(np.array(Sort_sua.units[i])*1E-6)
        offset = offset - 15
        
        mua_bin_width = 0.020
        mua = np.array(mua)

        mua_plot=np.histogram(mua, bins = np.arange(0, float(self.time_end.text()) - float(self.time_start.text()),mua_bin_width))
        plt.plot(mua_plot[1][0:-1],mua_plot[0]*.5+offset-40, linewidth=3, color='blue', alpha=1)
    

    #**************************** PLOT LFP RASTERS *************************
    #Load LFP Sort
    print "...loading lfp sort..."
    Sort_lfp = PTCS.PTCS(self.selected_sort_lfp) #Auto load flag for Nick's data

    #Save rasters to csv files for test 
    for k in range(Sort_lfp.n_units):
        np.savetxt(self.selected_sort_lfp[:-4]+"_"+str(k), np.float32(Sort_lfp.units[k])*50.*1E-6)

    offset = offset - 2.
    y = []

        
    for i in range(len(Sort_lfp.units)):
        spikes = np.array(Sort_lfp.units[i],dtype=np.float32)*1E-6*50
        print "... unit: ", i, " colour: ", colors[i%9], "  # events: ", len(spikes)
        if len(spikes)<100: continue

        x = spikes[np.where(np.logical_and(spikes>=float(self.time_start.text()), spikes<=float(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))

        ymin+=offset
        ymax+=offset-5

        plt.vlines(x-float(self.time_start.text()), ymin, ymax, linewidth=2, color=colors[i%9]) #colors[mod(counter,7)])
    
        offset=offset-5
        
    
    #******************** LABELING ********************

   
    #old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    #new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    #plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    

    old_xlabel = np.linspace(0, float(self.time_end.text()) - float(self.time_start.text()), 9)
    #new_xlabel = np.round(np.linspace(float(self.time_start.text()), float(self.time_end.text()), 9), 1)
    new_xlabel = np.round(np.linspace(float(self.time_start.text())/60., float(self.time_end.text())/60., 9), 1)
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)

    plt.xlim(0, float(self.time_end.text()) - float(self.time_start.text()))
    ax.tick_params(axis='both', which='both', labelsize=font_size)

    #plt.xlabel("Time (sec)", fontsize = font_size , weight = 'bold')        
    plt.xlabel("Day 1  Day 2 (min)", fontsize = font_size , weight = 'bold')        


    plt.ylabel(' Multi-Laminar LFP   ', fontsize=font_size, weight = 'bold')           

    plt.ylim(-70, -10)

    plt.show()
    


    

def Specgram_syncindex_tfr(self):
    '''PLOT TFR Specgram + MUA histograms + LFP traces + SUA Rasters
    '''
    
    colors=['blue','red', 'green', 'violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    channel=int(self.specgram_ch.text())
    top_channel = np.loadtxt(os.path.split(os.path.split(self.selected_sort_sua)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.

    tsf = TSF.TSF(self.selected_recording)
    tsf.read_ec_traces()
        
    #Find deepest channel by parsing sorted units:
    ax = plt.subplot(111)
   
    #******************* PLOT TRF SPECGRAM ***************************
    if False: 
        data_in = tsf.ec_traces[channel][int(float(self.time_start.text())*tsf.SampleFrequency):int(float(self.time_end.text())*tsf.SampleFrequency)]

        #if (os.path.exists(fname+".npy")==False):
        img, extent = Multitaper_specgram(data_in, 1000)    #
        #    np.save(fname, img)
        #    np.save(fname+'_extent', extent)
        #else: 
        #    img = np.load(fname+'.npy')
        #    extent = np.load(fname+'_extent.npy')

        #extent = extent[0]+start_traces, extent[1]+start_traces, extent[2], extent[3]
        im = ax.imshow(img, extent=extent, aspect='auto') #, extent=tsf.extent, cmap=tsf.cm)
    
    sync_0 = 0
    sync_1 = -10
    sync_scale = 20
    font_size = 30
    
    #******************* PLOT SUA RASTERS ***************************
    #Load SUA Sort

    Sort_sua = PTCS.PTCS(self.selected_sort_sua) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    if False: 
        n_spikes = []
        for k in range(len(Sort_sua.units)):
            n_spikes.append(len(Sort_sua.units[k]))
        n_spikes = np.array(n_spikes)
        indexes = np.argsort(n_spikes)
        print indexes
    else:
        indexes = np.arange(Sort_sua.n_units)

    offset = sync_0-10
    
    y = []
    mua = []
    for i in indexes: #range(len(Sort_sua.units)):
        print "... unit: ", i
    #for i in indexes[0:5]: #range(len(Sort_sua.units)):
        #x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)/float(Sort_sua.samplerate) #float(Sort1.samplerate)*2.5
        spikes = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)*1E-6

        x = spikes[np.where(np.logical_and(spikes>=float(self.time_start.text()), spikes<=float(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=offset+0.8
        ymax+=offset-0.8

        plt.vlines(x-float(self.time_start.text()), ymin, ymax, linewidth=3, color='black') #colors[mod(counter,7)])

        y.append(x)
        
        offset=offset-1.0

        mua.extend(x-float(self.time_start.text()))
    #******************** PLOT MUA HISTOGRAM *********************************
    
    #mua=[]
    #for i in range(len(Sort_sua.units)):
    #    mua.extend(np.array(Sort_sua.units[i])*1E-6)
    offset = offset - 30
    
    mua_bin_width = 0.020
    mua = np.array(mua)

    mua_plot=np.histogram(mua, bins = np.arange(0, float(self.time_end.text()) - float(self.time_start.text()),mua_bin_width))
    plt.plot(mua_plot[1][0:-1],mua_plot[0]*.5+offset-40, linewidth=3, color='blue', alpha=1)
    

    #******************** PLOT LFP TRACES ******************************
    offset = offset - 100
    t = np.arange(0, float(self.time_end.text())-float(self.time_start.text()), 1./tsf.SampleFrequency)
    for ch in range(0, tsf.n_electrodes, 4):
        if ch<top_channel: continue
        trace_temp = tsf.ec_traces[ch][int(float(self.time_start.text())*tsf.SampleFrequency): int(float(self.time_end.text())*tsf.SampleFrequency)] 
        
        trace_temp = butter_highpass_filter(trace_temp, 1.0, 1000, 5)
        trace_temp = butter_lowpass_filter(trace_temp, 50.0, 1000, 5)
        
        plt.plot(t,trace_temp/150.-1+offset, color='black')
        
        offset=offset-50


    #**************************** PLOT LFP RASTERS *************************
    #Load LFP Sort
    print "...loading lfp sort..."
    Sort_lfp = PTCS.PTCS(self.selected_sort_lfp) #Auto load flag for Nick's data

    offset = offset - 50.
    y = []
    for i in range(len(Sort_lfp.units)):
        spikes = np.array(Sort_lfp.units[i],dtype=np.float32)*1E-6*50

        x = spikes[np.where(np.logical_and(spikes>=float(self.time_start.text()), spikes<=float(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))

        ymin+=offset
        ymax+=offset-50

        plt.vlines(x-float(self.time_start.text()), ymin, ymax, linewidth=3, color=colors[i%9]) #colors[mod(counter,7)])
    
        #offset=offset-50
        
    
    #******************** LABELING ********************

   
    #old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    #new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    #plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    

    old_xlabel = np.linspace(0, float(self.time_end.text()) - float(self.time_start.text()), 9)
    new_xlabel = np.round(np.linspace(float(self.time_start.text()), float(self.time_end.text()), 9), 1)
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)

    plt.xlim(0, float(self.time_end.text()) - float(self.time_start.text()))
    ax.tick_params(axis='both', which='both', labelsize=font_size)

    plt.xlabel("Time (sec)", fontsize = font_size , weight = 'bold')        


    plt.ylabel(' Multi-Laminar LFP   ', fontsize=font_size, weight = 'bold')           



    plt.show()
    




def Specgram_syncindex(self):
    
    channel=int(self.specgram_ch.text())
    
    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    if '.tsf' in self.selected_recording:
        tsf = TSF.TSF(self.selected_recording)
        tsf.read_ec_traces()
        #for k in range(0, len(tsf.Siteloc),2):
        #    print k, tsf.Siteloc[k], tsf.Siteloc[k+1]
            
    elif '.lfp.zip' in self.selected_recording:
        tsf = load_lfpzip(self.selected_recording)
        #for k in range(len(tsf.Siteloc)):
        #    print k, tsf.Siteloc[k]

    
    samp_freq = tsf.SampleFrequency
    print "rec length: ", len(tsf.ec_traces[channel])/float(tsf.SampleFrequency), " sec."

    ax = plt.subplot(1,1,1)
    font_size = 30
    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    #Compute Specgram
    print "computing specgram..."
    data_in = tsf.ec_traces[channel][int(self.time_start.text())*tsf.SampleFrequency:int(self.time_end.text())*tsf.SampleFrequency][::int(tsf.SampleFrequency/1000)]
    samp_freq = 1000
    data_in = Notch_Filter(data_in)  
    
    f0 = 0.1; f1 = 100
    p0 = float(self.specgram_db_clip.text())
    temp_file = self.selected_recording[:-4]+"_ch"+str(channel)+"specgram"
    #if os.path.exists(temp_file+'.npz'):
    #    data = np.load(temp_file+'.npz')
    #    P, extent = data['P'], data['extent']
    #else:
    P, extent = Compute_specgram_signal(data_in, samp_freq, f0, f1, p0)
    #    np.savez(temp_file, P=P, extent=extent)
        
    plt.imshow(P, extent=extent, aspect='auto')

    #Compute sync index
    si_limit=0.7
    si, t, sync_periods = synchrony_index(tsf.ec_traces[channel][int(self.time_start.text())*tsf.SampleFrequency:(int(self.time_end.text())+10)*tsf.SampleFrequency], samp_freq, si_limit)
    
   # print "t: ", t
    
    sync_0 = -30
    sync_1 = -10
    sync_scale = 20
    
    plt.plot(t, si*sync_scale+sync_0, linewidth=5, color='green')
    plt.plot([0,max(t)],[-sync_scale*.3+sync_1,-sync_scale*.3+sync_1], 'r--', color='black', linewidth = 3, alpha=0.8)
    plt.plot([0,max(t)],[sync_1,sync_1], color='black', linewidth = 2, alpha=0.8)
    plt.plot([0,max(t)],[sync_0,sync_0], color='black', linewidth = 2, alpha=0.8)
    
    #xx = np.linspace(0,max(t)+10,5)
    #x_label = np.round(np.linspace(0, max(t)+10,5))
    #plt.xticks(xx, x_label, fontsize=20)
       
    #ax.set_xlim((0,P.shape[1]/2))    
    ax.set_ylim((sync_0-1,f1))    
    
    old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    
    ax.tick_params(axis='both', which='both', labelsize=font_size)

    plt.ylabel("Synchrony Index             Specgram Frequency (Hz)      ", fontsize=font_size-5)           
    plt.xlabel("Time (sec)", fontsize = font_size)
    #plt.title(self.recName, fontsize=font_size-10)
    plt.show()

def Compute_specgram_signal(data, SampleFrequency, f0=0.1, f1=110, p0=-40):

    print "...computing spegram fft..."
    t0=None
    t1=None
    #f0=0.1
    #f1=110
    #p0         #clipping bottom of specgram
    p1=None
    chanis=-1
    width=1
    tres=.25
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

#def wavelet(data, wname="db4", maxlevel=4):    #MAXLEVEL 6 has lower distortion. 
def wavelet(data, wname, maxlevel):
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
    
    #********************************************************************************************
    #*************************************** PLOT SPECGRAM FIRST ********************************
    #********************************************************************************************

    channel=int(self.specgram_ch.text())
    
    top_channel = np.loadtxt(os.path.split(os.path.split(self.selected_sort_sua)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.
    
    colors=['blue','red', 'green', 'violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    if '.tsf' in self.selected_recording:
        tsf = TSF.TSF(self.selected_recording)
        tsf.read_ec_traces()
        #for k in range(0, len(tsf.Siteloc),2):
        #    print k, tsf.Siteloc[k], tsf.Siteloc[k+1]
            
    elif '.lfp.zip' in self.selected_recording:
        tsf = load_lfpzip(self.selected_recording)
        #for k in range(len(tsf.Siteloc)):
        #    print k, tsf.Siteloc[k]

    
    samp_freq = tsf.SampleFrequency
    print "rec length: ", len(tsf.ec_traces[channel])/float(tsf.SampleFrequency), " sec."

    ax = plt.subplot(1,1,1)
    font_size = 30
    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    #Compute Specgram
    print "computing specgram..."
    data_in = tsf.ec_traces[channel][int(self.time_start.text())*tsf.SampleFrequency:int(self.time_end.text())*tsf.SampleFrequency]
    #data_in = Notch_Filter(data_in)  
    
    f0 = 0.1; f1 = 100
    p0 = float(self.specgram_db_clip.text())
    temp_file = self.selected_recording[:-4]+"_ch"+str(channel)+"specgram"
    #if os.path.exists(temp_file+'.npz'):
    #    data = np.load(temp_file+'.npz')
    #    P, extent = data['P'], data['extent']
    #else:
    P, extent = Compute_specgram_signal(data_in, samp_freq, f0, f1, p0)
    #    np.savez(temp_file, P=P, extent=extent)
        
    plt.imshow(P, extent=extent, aspect='auto')

    #Compute sync index
    si_limit=0.7
    si, t, sync_periods = synchrony_index(tsf.ec_traces[channel][int(self.time_start.text())*tsf.SampleFrequency:(int(self.time_end.text())+10)*tsf.SampleFrequency], samp_freq, si_limit)
    
   # print "t: ", t
    
    sync_0 = -30
    sync_1 = -10
    sync_scale = 20
    
    plt.plot(t, si*sync_scale+sync_0, linewidth=3, color='red')
    plt.plot([0,max(t)],[-sync_scale*.3+sync_1,-sync_scale*.3+sync_1], 'r--', color='black', linewidth = 3, alpha=0.8)
    plt.plot([0,max(t)],[sync_1,sync_1], color='black', linewidth = 2, alpha=0.8)
    plt.plot([0,max(t)],[sync_0,sync_0], color='black', linewidth = 2, alpha=0.8)


    #ax.set_ylim((sync_0-1,f1))    
    
    #********************************************************************************************
    #*************************************** PLOT RASTERS ***************************************
    #********************************************************************************************

    
    #Load LFP Sort
    Sort_lfp = PTCS.PTCS(self.selected_sort_lfp) #Auto load flag for Nick's data

    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'

    Sort_sua = PTCS.PTCS(self.selected_sort_sua) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    
    if False: 
        n_spikes = []
        for k in range(len(Sort_sua.units)):
            n_spikes.append(len(Sort_sua.units[k]))
        n_spikes = np.array(n_spikes)
        indexes = np.argsort(n_spikes)
        print indexes
    else:
        indexes = np.arange(Sort_sua.n_units)

    offset = sync_0-10
    
    ax = plt.subplot(1, 1, 1)
    y = []
    for i in indexes: #range(len(Sort_sua.units)):
        print "... unit: ", i
    #for i in indexes[0:5]: #range(len(Sort_sua.units)):
        #x = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)/float(Sort_sua.samplerate) #float(Sort1.samplerate)*2.5
        spikes = np.array(Sort_sua.units[indexes[i]],dtype=np.float32)*1E-6

        x = spikes[np.where(np.logical_and(spikes>=int(self.time_start.text()), spikes<=int(self.time_end.text())))[0]]

        ymin=np.zeros(len(x))
        ymax=np.zeros(len(x))
        ymin+=offset+0.4
        ymax+=offset-0.4

        plt.vlines(x-int(self.time_start.text()), ymin, ymax, linewidth=1, color='black') #colors[mod(counter,7)])

        y.append(x)
        
        offset=offset-1.0

    ##Plot LFP spike
    #offset = offset -10.
    #y = []
    #for i in range(len(Sort_lfp.units)):
        ##x = np.array(Sort_lfp.units[i],dtype=np.float32)/float(Sort_sua.samplerate)*50 #***************************** UNCOMPRESSSING LFP RASTERS
        #spikes = np.array(Sort_lfp.units[i],dtype=np.float32)*1E-6*50

        #x = spikes[np.where(np.logical_and(spikes>=int(self.time_start.text()), spikes<=int(self.time_end.text())))[0]]

        #ymin=np.zeros(len(x))
        #ymax=np.zeros(len(x))
        
        #ymin+=offset-5
        #ymax+=offset-7
        
        #plt.vlines(x-int(self.time_start.text()), ymin, ymax, linewidth=3, color=colors[i%9]) #colors[mod(counter,7)])
    
        #offset=offset-2


    #********************************************************************************
    #******************** PLOT LFP TRACES ******************************
    #********************************************************************************
    offset = offset - 100
    t = np.arange(0, float(self.time_end.text())-float(self.time_start.text()), 1./tsf.SampleFrequency)
    ch_skip = 1
    if tsf.n_electrodes>10: ch_skip = 4
    for ch in range(0, tsf.n_electrodes, ch_skip):
        if ch<top_channel: continue
        trace_temp = tsf.ec_traces[ch][int(float(self.time_start.text())*tsf.SampleFrequency): int(float(self.time_end.text())*tsf.SampleFrequency)] 
        
        trace_temp = butter_highpass_filter(trace_temp, 4.0, 1000, 5)
        trace_temp = butter_lowpass_filter(trace_temp, 100.0, 1000, 5)
        
        plt.plot(t,trace_temp/5.-1+offset, color='black')
        
        offset=offset-50



    #********************************************************************************
    #******************** LABELING ********************
    #********************************************************************************
    
    old_ylabel = [sync_0, sync_0/2, sync_1, 0, f1/2, f1]
    new_ylabel = [0, 0.7, 1, 0, f1/2, f1]
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    

    old_xlabel = np.linspace(0, float(self.time_end.text()) - float(self.time_start.text()), 5)
    new_xlabel = np.round(np.linspace(float(self.time_start.text()), float(self.time_end.text()), 5), 1)
    plt.xticks(old_xlabel, new_xlabel, fontsize=font_size)


    ax.tick_params(axis='both', which='both', labelsize=font_size)

    plt.xlabel("Time (sec)", fontsize = font_size , weight = 'bold')        

        
    #plt.xlabel('Time (seconds)',fontsize=35, weight='bold')
    #plt.ylabel('Single Unit ID',multialignment='center', fontsize=35)

    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))

    #plt.ylabel('Single Unit Rasters     Synchrony Index   Specgram Frequency (Hz)', fontsize=font_size, weight = 'bold')           
    
    plt.ylabel('LFP   Thresholded Events      Clustered Events', fontsize=font_size, weight = 'bold')           
    #plt.ylabel('LFP Cluster Raster         Single Unit IDs',multialignment='center', fontsize=35, weight='bold')
    
    plt.xlim(0, int(self.time_end.text())-int(self.time_start.text()))
    
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

    plotting = False

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
    cell_type=[]
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
        else: color='red'; m = '^'

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

        peaks.append(peak)
        troughs.append(troughs)
        cell_type.append(m)
    if plotting: 
            ax = plt.subplot(1,3,1)
            plt.title("Unit location", fontsize=25)
            plt.axhspan(-2000, 0, facecolor='black', alpha=0.1)
            plt.xlabel("Location (um)", fontsize=25)

            plt.ylabel("Depth (um)", fontsize=25)
            old_ylabel = np.arange(-2000, 0, 500)
            new_ylabel = np.arange(2000, 0, -500)
            plt.yticks(old_ylabel, new_ylabel, fontsize=20) #,rotation='vertical')
            ax.tick_params(axis='both', which='major', labelsize=25)

            ax = plt.subplot(1,3,2)
            plt.title("Cell Types - Templates", fontsize=25)
            ax.tick_params(axis='both', which='major', labelsize=25)
            plt.ylabel("Template Amplitude (uV)", fontsize = 25)
            plt.xlabel("Time (ms)", fontsize = 25)

                


            plt.show()        


    #np.save( , peaks)
    #np.save( , troughs)

    #****************************************************************************************
    #********************************* COMPUTE DRIFT VALUES *********************************
    #****************************************************************************************

    selected_unit = int(self.starting_cell.text())

    #Find recording length; needed for plotting distributions
    if os.path.exists(self.parent.sua_file.replace('.ptcs','.tsf')):
        tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))
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
    
    
    if plotting: 
        xtick_lbls = []
        for k in range(len(time_chunks)):
            xtick_lbls.append(int(time_chunks[k][1]*1E-6/60.))
        
        old_xlabel = np.arange(0, len(time_chunks), 1)
        plt.xticks(old_xlabel, xtick_lbls, fontsize=20) #,rotation='vertical')

        plt.xlabel("Time (mins)", fontsize=30)
        plt.ylabel("Phase Lock in Each Epoch (ms)", fontsize=30)

        plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(len(Sort_lfp.units[lfp_cluster]))+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
                    '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)

        ax.tick_params(axis='both', which='major', labelsize=30)
        
        plt.ylim(-1E2,1E2)
        plt.xlim(0.1,len(time_chunks)+0.1)
        plt.show()
    
    #****************************************************************************************
    #************************************* MAKE DRIFTING MOVIES *****************************
    #****************************************************************************************

    plt.close()

    #Interpolate positions over epoch beginning and ends;

    vid_array = np.zeros((len(control_array), time_chunks[len(time_chunks)-2][1]*1E-6/60., 2), dtype=np.float32)
    print vid_array.shape
    #pts_array contains [epoch, MSL time] pairs - or empty; len(control_array) = # epochs
    for p in range(len(pts_array)):
        print "... cell: ", p
        
        temp_x = [] 
        temp_y = []
        for k in range(len(pts_array[p])):
            t1 = time_chunks[pts_array[p][k][0]][0]*1E-6/60.
            t2 = time_chunks[pts_array[p][k][0]][1]*1E-6/60.
            print p, k, pts_array[p][k], t1, t2
            
            temp_y.append(pts_array[p][k][1])
            temp_x.append(int(t1))
        
        
        #Interpolate positions over epoch beginnning and ends; do it every minute
        #Do linear interpolation if only 2 points
        if len(temp_x)==2: 
            #Interpolate positions over epoch beginnning and ends; do it every minute
            interp_drift = np.linspace (temp_y[0], temp_y[1], temp_x[1]-temp_x[0])
            #print interp_drift
            ctr = 0
            for q in range(int(temp_x[0]), int(temp_x[1]),1):
                vid_array[p, q] = [interp_drift[ctr], -p]
                ctr+=1

        #Use spline interpolation of higher order 3 or more points
        if len(temp_x)>2: 
            xi = np.hstack(temp_x)
            yi = np.hstack(temp_y)

            order = 4
            s = InterpolatedUnivariateSpline(xi, yi, k=min(order,len(temp_x)-1))
            x = np.arange(temp_x[0], temp_x[-1], 1)
            yinterp = s(x)


            ctr = 0
            for q in list(x):
                vid_array[p, q] = [yinterp[ctr], -p]
                ctr+=1
    print vid_array.shape
    vid_array = np.swapaxes(vid_array,0,1)
    print vid_array.shape
    #print vid_array[:,0]


    #**********************************************************************************
    #********************************** RUN MOVIES ************************************
    #**********************************************************************************

    #Initialize plot arrays
    plot_array = []

    out_x = []
    out_y = []
    out_color = []
    out_marker = []
    out_alpha = []
    out_edgecolor = []

    for k in range(len(vid_array)):
        temp_x = []
        temp_y = []
        temp_colour = []
        temp_alpha = []
        temp_edgecolour = []
        temp_marker = []
        for p in range(len(vid_array[k])):
            #print vid_array[p]
            if (vid_array[k,p,0]==0) and (vid_array[k,p,1]==0):             #If cell not active we save this data; 
                temp_colour.append([0,0,0,0])
                temp_edgecolour.append([0,0,0,0])
            else:
                if cell_type[p] == 'o': 
                    temp_colour.append([0,0,1,1])
                    temp_edgecolour.append([0,0,1,1])
                else: 
                    temp_colour.append([1,0,0,1])
                    temp_edgecolour.append([1,0,0,1])
            temp_x.append(vid_array[k,p][0])
            #temp_y.append(vid_array[k,p][1])
            temp_y.append(-Sort_sua.ypos[p])
            
            temp_marker.append(cell_type[p])
        
        out_x.append(np.vstack(temp_x))
        out_y.append(np.vstack(temp_y))
        out_color.append(np.vstack(temp_colour))
        out_edgecolor.append(np.vstack(temp_edgecolour))


    #**** INITIALIZE ANIMATIONS
    plt.close()


    from matplotlib import animation
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    plt.close('all') # close all previous plots

    fig = plt.figure(1)
    ax = plt.axes(xlim=(-50, 10), ylim=(-2000, 100))
    
    plt.plot([-50,10,],[0,0], 'r--', linewidth=3, color='black', alpha=.8)

    scat = ax.scatter([], [], s=60)
    plt.ylabel("Depth (um)", fontsize=20, labelpad=0)
    old_ylabel = np.arange(-2000, 0, 500)
    new_ylabel = np.arange(2000, 0, -500)
    plt.yticks(old_ylabel, new_ylabel, fontsize=20) #,rotation='vertical')
    ax.tick_params(axis='both', which='major', labelsize=15)

    plt.xlabel("Time from LFP event(ms)", fontsize=20)
    plt.plot([0,0],[-2000,1], 'r--', linewidth=3, color='black', alpha=.8)

    def init():
        scat.set_offsets([])
        return scat,

    def animate(i):
        print i, " / ", len(out_x)
        plt.suptitle("Time: "+str(i)+" min ", fontsize=20)

        scat.set_facecolor(out_color[i])
        scat.set_edgecolor(out_edgecolor[i])
        scat.set_offsets(np.hstack((out_x[i], out_y[i])))
        #scat.set_markers(temp_marker)
        return scat,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(out_x), interval=200, blit=False, repeat=False)

    print self.parent.sua_file.replace('tsf_files','movie_files').replace('.ptcs','')+'_LPF_'+str(lfp_cluster)
    anim.save(self.parent.sua_file.replace('tsf_files','movie_files').replace('.ptcs','')+'_LPF_'+str(lfp_cluster)+'.mp4', writer=writer)
    print " DONE! "
    plt.show()


    #Save lfp event times and order to a .npz file
    #out_x = the location of all cells at each interpolated time point (usually minute)
    #out_colour = the colour indicates whether cell passed control measures
    #

    out_file = self.parent.sua_file.replace('.ptcs','')+'_sequences_LPF_'+str(lfp_cluster)
    np.savez(out_file, sequences=out_x, controls=out_color, lfp_events = pop_spikes)


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
        tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))
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
        tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))
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



def Compute_LFP_rasters(self):
    ''' Saves spiketimes for each LFP event for each cell; 
        Also saves the gaussian convolved data, based on sigma provided
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
        
    #Load SUA Sort
    #sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = PTCS.PTCS(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)
    print type(Sort_sua.units[0][-1])

    #Load LFP Sort
    Sort_lfp = PTCS.PTCS(self.parent.lfp_event_file) #Auto load flag for Nick's data

    #start_lfp = min(int(self.parent.start_lfp.text()),len(Sort_lfp.units)-1)
    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardcoded to: ", compress_factor

    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())
      
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
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
    #if os.path.exists(cell_rasters_filename+'.npy')==False:
    if True:
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
                    if (pop_spikes[j+1]-pop_spikes[j])<100000:  #microsecond timing.
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
       
        np.save(cell_rasters_filename, cell_rasters)
   
def Compute_LFP_histograms(self):
    ''' Computes the histogram for each cell for each event; uses gaussian convolution
    '''

    sig = float(self.sigma_width.text())

    lfp_cluster = int(self.parent.lfp_cluster.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    print cell_rasters_filename
    cell_rasters = np.load(cell_rasters_filename+'.npy')
    print cell_rasters.shape
    
    
    x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
    sig_gaussian = np.float32(gaussian(x, 0, sig))

    
    cell_histograms = np.zeros((len(cell_rasters),len(cell_rasters[0]), 2000), dtype=np.float32)
    for unit in range(len(cell_rasters)):         #PREPROCESS ALL UNIT MSLs
        print "...processing cell: ", unit
        for lfp_event in range(len(cell_rasters[unit])):
            
            fit = np.zeros(2000, dtype=np.float32)

            if len(cell_rasters[unit][lfp_event])>0: 
                locked_spikes = np.int16(cell_rasters[unit][lfp_event]*1E-3)    #Convert to millseconds
                #Compute convolved gaussian
                for g in range(len(locked_spikes)):
                    mu = int(locked_spikes[g])
                    fit += np.roll(sig_gaussian, mu)

            #plt.plot(fit)
            #plt.show()

            cell_histograms[unit][lfp_event] = fit
        
    np.save(cell_rasters_filename+"_"+self.sigma_width.text()+'ms_histograms', cell_histograms)

    
    
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
    print "... selected lfp_cluster: ", lfp_cluster
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardwired to: ", compress_factor
    
    
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = Sort_lfp.units[lfp_cluster]*compress_factor#*1E-3  
    #original_n_popspikes = len(pop_spikes)
    #pop_spikes=np.sort(np.unique(pop_spikes))       
    #print " ... # LFP events: ", len(pop_spikes), "   #duplicates: ", original_n_popspikes-len(pop_spikes)
    
    
    n_units = 1 #int(self.ending_cell.text()) - int(self.starting_cell.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")


    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    #OPTION 1: Divide into chunks of recording length
    if False: 
        self.parent.tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))
        temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text())+1)
    
    #OPTION 2: Divide into chunks of LFP events
    else:
        n_spikes = len(Sort_lfp.units[lfp_cluster])
        temp_chunks=[]
        chunk_width = int(n_spikes/float(self.parent.time_chunks.text()))
        for t in range(0, n_spikes, chunk_width):
            temp_chunks.append(Sort_lfp.units[lfp_cluster][t]*compress_factor)
        temp_chunks.append(Sort_lfp.units[lfp_cluster][-1]*compress_factor)

    
    
    #LOAD temp chunks from above into time chunks 
    time_chunks = []
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    print time_chunks[:10]
    
    
    chunk_index = []
    for chk in range(int(self.parent.time_chunks.text())):
    #for chk in self.chunks_to_plot.text().split(','):
        chunk_index.append(int(chk))
    print "...chunk_index: ", chunk_index
    
    
    locked_spike_array = []
    f1 = plt.figure()

    sig = float(self.sigma_width.text())
    ax = plt.subplot(1,1,1)
    for ctr, chunk_ctr in enumerate(chunk_index):                   #******************* LOOP OVER EPOCHS **********************
        offset=0     #Used for plotting rasters from multiple cells
        time_chunk = time_chunks[chunk_ctr]
        print "...time chunk: ", time_chunk[0]*1E-6/60., time_chunk[1]*1E-6/60., "  mins."
        ax = plt.subplot(1, len(chunk_index), ctr+1)
        
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]

        locked_spike_array.append([])
        #for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):  #********************** LOOP OVER UNITS ********************
        unit = int(self.starting_cell.text())
        locked_spikes = cell_rasters[unit][temp3]       #Vertical stack of rasters during chunk
        
        if len(locked_spikes)>0:
            
            print "... chunk: ", time_chunk, " ...unit: ", unit, " #spikes locked: ", len(np.hstack(locked_spikes)), " / ", len(Sort_sua.units[unit]), \
            "   duplicates: ", len(np.hstack(locked_spikes)) - len(np.unique(np.hstack(locked_spikes)))

            locked_spike_array[ctr].append(np.hstack(locked_spikes))

            #Plot rasters
            for event in range(len(locked_spikes)):
                ymin=np.zeros(len(locked_spikes[event]))
                ymax=np.zeros(len(locked_spikes[event]))
                ymin+=offset
                ymax+=offset+1
                offset+=1

                plt.vlines(np.float64(locked_spikes[event])*1E-3, ymin, ymax, linewidth=5, color='black',alpha=1) #colors[mod(counter,7)])
                
            plt.plot([-100,100], [offset,offset], linewidth=3, color='black', alpha=.8)
            
            #Don't convert xx1 into stack until after plotting rasters
            n_lfp_spikes= len(locked_spikes)

            #if len(xx1)>0: xx1 = np.hstack(xx1)
            #cell_rasters[chunk_ctr].append(xx1)
            
            all_spikes = np.sort(np.hstack(locked_spikes))
        
            #Convolve w. gaussian
            fit_sum_even = np.zeros(2000000, dtype=np.float32)
            fit_sum_odd = np.zeros(2000000, dtype=np.float32)
            
            x = np.linspace(-1000,1000, 2000000)    #Make an array from -1000ms .. +1000ms with microsecond precision
            sig_gaussian = np.float32(gaussian(x, 0, sig))
            for g in range(len(all_spikes)):
                #print "...spike: ", all_spikes[g]
                mu = int(all_spikes[g])
                if g%2==0: fit_sum_even += np.roll(sig_gaussian, mu)
                else: fit_sum_odd += np.roll(sig_gaussian, mu , axis=0)
            

            t = np.linspace(-1000, 1000, 2000000)
            
            if True: 
                #plt.plot(t, fit_sum_even/np.max(fit_sum_even)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='blue', linewidth=6, alpha=.6)
                #plt.plot(t, fit_sum_odd/np.max(fit_sum_odd)*len(locked_spikes) +len(locked_spikes)*(unit-int(self.starting_cell.text())), color='red', linewidth=6, alpha=.6)
                plt.plot(t, fit_sum_even/np.max(fit_sum_even)*len(locked_spikes) , color='blue', linewidth=6, alpha=.6)
                plt.plot(t, fit_sum_odd/np.max(fit_sum_odd)*len(locked_spikes), color='red', linewidth=6, alpha=.6)
            else:
                fit_sum_even = (fit_sum_even+fit_sum_odd)/2.
                plt.plot(t, fit_sum_even/np.max(fit_sum_even)*len(locked_spikes) +len(locked_spikes), color='blue', linewidth=7, alpha=.9)


        else:
            print "... no spikes in epoch..."
            all_spikes = []
            locked_spike_array[ctr].append([])    
            
            fit_sum_even = 0

            #bin_width = sig*1E3   # histogram bin width in usec
            #y = np.histogram(all_spikes, bins = np.arange(-1000000,1000000,bin_width))
            #plt.bar(y[1][:-1]*1E-3, np.float32(y[0])/np.max(y[0])*len(locked_spikes), bin_width*1E-3, color='blue', alpha=0.2)
            #plt.bar(y[1][:-1], y[0], bin_width, color='blue')

        t_max = t[np.argmax(fit_sum_even)]
        #plt.plot([t_max,t_max], [0,5000], 'r--', linewidth = 6, color='black')
        plt.title("#spks: "+ str(len(all_spikes)) + "  Lock time: "+str(round(t_max,1)))
        print "...peak max: ", t_max
        
        #plt.xticks(list(plt.xticks()[0]) + [round(t_max,1)])
        #ax.xaxis.set_ticks([50, 50, round(t_max,1), 0, 10])
        #ax.xaxis.set_ticks

        plt.plot([0,0], [0, n_lfp_spikes*n_units], 'r--', linewidth=3, color='black', alpha=.8)

        #plt.xlim(lock_window_start-1,lock_window_end+1)
        #plt.xlim(-100,100)
        plt.xlim(int(self.parent.lock_window_start.text()), int(self.parent.lock_window_end.text()))
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


        #plt.ylabel("Time Chunk: "+str(int(time_chunk[0]*1E-6/60.)) + ".."+str(int(time_chunk[1]*1E-6/60.))+" mins", fontsize=30,  labelpad=-2)
        #plt.ylabel("Epoch: "+str(ctr+1), fontsize=30,  labelpad=-2)
        plt.title("E: "+str(ctr+1), fontsize=30)
        plt.xlabel("Time (ms)", fontsize = 30)
        plt.tick_params(axis='both', which='both', labelsize=30)

    #plt.suptitle("LFP Cluster: "+str(lfp_cluster)+ " # events: " + str(n_lfp_spikes)+",  Unit: "+self.starting_cell.text()+ " #spikes: " +str(len(Sort_sua.units[unit]))+ \
    #            ",  sigma: " + str(sig) +"(ms)" + '\n'+self.parent.sua_file.replace(self.parent.root_dir,''), fontsize=20)
    plt.suptitle(os.path.split(self.parent.sua_file.replace(self.parent.root_dir,''))[1]+"  "+ str(sig) +"(ms)", fontsize=20)


    #******************* PVALUE MATRIX ***************

    #COMPUTE KS P-VALUES FOR CURRENT LFP CLUSTER

    f2 = plt.figure()
    empty_row = np.zeros(len(locked_spike_array),dtype=np.float32)+1.0


    #for unit in range(int(self.starting_cell.text()), int(self.ending_cell.text()),1):  #********************** LOOP OVER UNITS ********************
    kstest_array = []
    for k in range(len(locked_spike_array)):
        print "... unit: ", int(self.starting_cell.text()), " epoch: ", k
        kstest_array.append([])    
        

        spk1 = np.hstack(locked_spike_array[k])*1E-3
        indexes = np.where(np.logical_and(spk1>=int(self.parent.lock_window_start.text()), spk1<=int(self.parent.lock_window_end.text())))[0]   #Work in milliseconds            #********************* ONLY USING SPIKES WITHIN 100ms of t=0 for KS TEST
        if len(indexes)==0:
            kstest_array[k] = empty_row             #DESYNCH STATES; APPEND DUMMY VALUE
            continue

        spk1 = spk1[indexes]
          
        for p in range(len(locked_spike_array)):
            spk2 = np.hstack(locked_spike_array[p])*1E-3
            #indexes = np.where(np.logical_and(spk2>=-100, spk2<=100))[0]   #Work in milliseconds
            indexes = np.where(np.logical_and(spk2>=int(self.parent.lock_window_start.text()), spk2<=int(self.parent.lock_window_end.text())))[0]   #Work in milliseconds            #********************* ONLY USING SPIKES WITHIN 100ms of t=0 for KS TEST

            if len(indexes)==0:
                kstest_array[k].append(1.0) #DESYNCH STATES; APPEND DUMMY VALUE
                continue 

            spk2 = spk2[indexes]
            #print spk1
            #print spk2
            KS, p_val = stats.ks_2samp(spk1, spk2)
            
            #print "... p_val: ", p_val
            kstest_array[k].append(p_val)
   
        kstest_array[k]=np.hstack(kstest_array[k])
   
    kstest_array = np.vstack(kstest_array)[::-1]
    kstest_array = np.log10(kstest_array)
    vmin_value = np.min(kstest_array)
    if vmin_value < -6: vmin_value = -6
    vmax_value = max(0,np.max(kstest_array))

    cax = plt.imshow(kstest_array, vmin=vmin_value, vmax=vmax_value, cmap=self.cmap.text(), interpolation='none')


    #cax = ax.imshow(ks_img, vmax=float(self.vmin_value.text()), vmin=float(self.vmax_value.text()), cmap=cm.jet_r, interpolation='none')
    #plt.ylim(len(ks_img),0)

    #ax = plt.subplot(2,3,4)
    
    ##Plot SUA Events
    #spikes = Sort_sua.units[unit]*1E-6/60.
    #ymin=np.zeros(len(spikes))
    #ymax=np.zeros(len(spikes))
    #ymin+=0
    #ymax+=10
    #plt.vlines(spikes, ymin, ymax, linewidth=1, color='black',alpha=.5) #colors[mod(counter,7)])
    
    plt.tick_params(axis='both', which='both', labelsize=30)
    plt.xlabel("Epochs", fontsize=30)
    plt.ylabel("Epochs", fontsize=30)
    
    tstep= 1
    old_xlabel = np.arange(0, len(kstest_array), tstep)
    new_xlabel = np.arange(1,len(kstest_array)+1,tstep)
    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    
    
    old_ylabel = np.arange(len(kstest_array)-1, -1, -tstep)
    new_ylabel = np.arange(1,len(kstest_array)+1,tstep)
    plt.yticks(old_ylabel, new_ylabel, fontsize=30) #, rotation='vertical')    
    
    #vmin_value = np.min(kstest_array)
    #vmax_value = 1.0
    print vmin_value, vmax_value
    #x_ticks = np.arange(float(self.vmin_value.text()), float(self.vmax_value.text())+1, 1)
    x_ticks = np.arange(vmax_value, vmin_value, -1)
    x_ticks = np.append(x_ticks, vmin_value)
    print x_ticks
    #x_ticks = [float(vmax_value), float(vmin_value)]
    cbar = f2.colorbar(cax, ticks=x_ticks)
    
    x_tick_labels = []
    for k in range(0,len(x_ticks)-1,1):
        x_tick_labels.append('$10^{'+str(int(x_ticks[k]))+"}$")
    
    x_tick_labels.append('$10^{'+str(round(x_ticks[-1],1))+"}$")
    

    print x_tick_labels
    cbar.ax.set_yticklabels(x_tick_labels)  # vertically oriented colorbar
    #cbar.ax.set_xticklabels([str(vmax_value), str(vmin_value)])  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=30) 
    
    #plt.suptitle('$log_{10}(P-Value)$', fontsize=30)
    plt.suptitle('P-Value', fontsize=30, fontweight='bold')

    
    plt.show()
   
    
    ##ks_img = np.ma.log10(np.clip(1-ks_img,0,1))
    #ks_img = np.ma.log10(ks_img)
    
    ##print ks_img
    
    #cax = ax.imshow(ks_img, vmax=float(self.vmin_value.text()), vmin=float(self.vmax_value.text()), cmap=cm.jet_r, interpolation='none')
    ##plt.ylim(len(ks_img),0)

    ##ax = plt.subplot(2,3,4)
    
    ###Plot SUA Events
    ##spikes = Sort_sua.units[unit]*1E-6/60.
    ##ymin=np.zeros(len(spikes))
    ##ymax=np.zeros(len(spikes))
    ##ymin+=0
    ##ymax+=10
    ##plt.vlines(spikes, ymin, ymax, linewidth=1, color='black',alpha=.5) #colors[mod(counter,7)])
    
    #plt.tick_params(axis='both', which='both', labelsize=30)
    #plt.xlabel("Time (min)", fontsize=30)
    #plt.ylabel("Time (min)", fontsize=30)
    
    #tstep= 50
    ##old_xlabel = np.arange(len(ks_img), 0, -tstep)
    #new_xlabel = np.arange(0,len(ks_img),tstep)
    #plt.xticks(new_xlabel, fontsize=30) #, rotation='vertical')    
    
    #old_ylabel = np.arange(len(ks_img), 0, -tstep)
    #new_ylabel = np.arange(0,len(ks_img),tstep)
    #plt.yticks(old_ylabel, new_ylabel, fontsize=30) #, rotation='vertical')    
        
    #x_ticks = np.arange(float(self.vmin_value.text()), float(self.vmax_value.text())+1, 1)
    #cbar = fig.colorbar(cax, ticks=[x_ticks])
    ##cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
    #cbar.ax.tick_params(labelsize=30) 

    #plt.suptitle("2-Sample KS-Test: P-values "+ " Unit: " + str(unit), fontsize=30)
    #plt.show()
    


    #plt.show()

def compute_nspikes_histograms(self):
    '''  Compute number of spikes from each cell locked to each lfp event
    '''

    Sort_sua = Ptcs(self.parent.sua_file)

    lfp_cluster = int(self.parent.lfp_cluster.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    data = np.load(cell_rasters_filename+".npy")

    for u in range(len(Sort_sua.units)):
        ax = plt.subplot(11,11,u+1)
        ax.tick_params(axis='both', which='both', labelsize=8)
        length_array = []
        for k in range(len(data[u])):
            length_array.append(len(data[u][k]))
        length_array = np.array(length_array)

        bin_width = 1   #100ms bins
        y = np.histogram(length_array, bins = np.arange(0,np.max(length_array),bin_width))
        if len(y[0])>1:
            #print y[0]
            plt.bar(y[1][:-1], y[0], bin_width, color='blue')
            plt.xlim(1,20)
            plt.ylim(0,np.max(y[0][1:]))
        if u!= 110: ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.ylabel(len(Sort_sua.units[u]), fontsize=10, labelpad=-2)
    plt.show()


def compute_isi_histograms(self):

    Sort_sua = Ptcs(self.parent.sua_file)

    lfp_cluster = int(self.parent.lfp_cluster.text())

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    data = np.load(cell_rasters_filename+".npy")

    for u in range(len(Sort_sua.units)):
    #for u in range(10):
        print "... unit: ", u
        ax = plt.subplot(11,11,u+1)
        ax.tick_params(axis='both', which='both', labelsize=8)

        isi_array = []
        for k in range(len(data[u])):
            #print "...lfp event: ", k
            for j in range(len(data[u][k])-1):
                #print data[u][k][j]
                isi_array.append((data[u][k][j+1]-data[u][k][j])*1E-3)
        isi_array = np.array(isi_array)
        #return
        #print isi_array
        #print len(isi_array)

        bin_width = 1   #100ms bins
        y = np.histogram(isi_array, bins = np.arange(0,50,bin_width))
        if len(y[0])>1:
            #print y[0]
            plt.bar(y[1][:-1], y[0], bin_width, color='blue')
            #plt.xlim(1,20)
            #plt.ylim(0,np.max(y[0][1:]))
        if u!= 110: ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.ylabel(len(Sort_sua.units[u]), fontsize=10, labelpad=-2)
    plt.show()

def plt_msl_discrete_single(self):
    

    sig = float(self.sigma_width.text())            #SHOULD EVENTUALLY SAVE THIS IN THE META DATA AS WELL; or save information as .npz file <- even better

    start_window = int(self.parent.lock_window_start.text())
    end_window = int(self.parent.lock_window_end.text())
    
    lfp_cluster = int(self.parent.lfp_cluster.text())
    compress_factor = 50
    

    Sort_lfp = PTCS.PTCS(self.parent.lfp_event_file)
    

    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    print cell_rasters_filename
    cell_rasters = np.load(cell_rasters_filename+'.npy')
    print cell_rasters.shape
    

    
    cell_histograms = np.load(cell_rasters_filename+"_"+ self.sigma_width.text() + "ms" + "_histograms.npy")
    print cell_histograms.shape
    
    img_out = cell_histograms[int(self.starting_cell.text())]
    #Normalize data
    temp_array = np.zeros(img_out.shape)
    maxloc_array = []
    for k in range(len(img_out)):
        if np.max(img_out[k])!=0:
            temp_array[k] = img_out[k]/np.max(img_out[k])

        if (1000-np.argmax(img_out[k]) < int(self.parent.lock_window_end.text())):
            maxloc_array.append([Sort_lfp.units[lfp_cluster][k]*1E-6*50, np.argmax(img_out[k])])
        #else:
        #    maxloc_array.append(123)
        
        
    print maxloc_array
    maxloc_array = np.int32(maxloc_array)
    
    np.savetxt(cell_rasters_filename+"_maxlocs", maxloc_array, delimiter = ',')
    quit()
        
    img_out = np.array(temp_array).T

    #DO same test on shuffled data
    img_out_shuffled = cell_histograms[int(self.starting_cell.text())]
    #print img_out_shuffled.shape
    shuffle_indexes = np.random.randint(len(img_out_shuffled), size=len(img_out_shuffled))
    #print shuffle_indexes
    img_out_shuffled = img_out_shuffled[shuffle_indexes]
    
    temp_array = np.zeros(img_out.shape)
    maxloc_array_shuffled = []
    temp_array = np.zeros(img_out_shuffled.shape)
    for k in range(len(img_out_shuffled)):
        if np.max(img_out_shuffled[k])!=0:
            temp_array[k] = img_out_shuffled[k]/np.max(img_out_shuffled[k])

        #if (1000-np.argmax(img_out_shuffled[k]) < int(self.parent.lock_window_end.text())):
        maxloc_array_shuffled.append(np.argmax(img_out_shuffled[k]))
        #else:
        #    #maxloc_array_shuffled.append(123)
        #    maxloc_array_shuffled.append(np.argmax(img_out_shuffled[k]))
            
    img_out_shuffled = np.array(temp_array).T

    #************************
    f1 = plt.figure()
    ax = plt.subplot(2, 1, 1)
    
    #print img_out.shape
    #window = 100

    img_out = img_out[1000+start_window:1000+end_window]
    plt.imshow(img_out, aspect='auto', interpolation='none')
    
    plt.plot([0,len(img_out[0])] , [end_window ,end_window ], 'r--', linewidth=2, color='white', alpha=0.8)
    

    #LABELS
    
    old_ylabel = np.linspace(0, end_window*2, 5)
    new_ylabel = np.linspace(start_window, end_window, 5)
    plt.yticks(old_ylabel, new_ylabel, fontsize=30) #, rotation='vertical')    
    
    #print Sort_lfp.units[lfp_cluster]
    #print len(Sort_lfp.units[lfp_cluster])
    
    old_xlabel = np.arange(0, len(Sort_lfp.units[lfp_cluster]), 10)
    new_xlabel = np.round(Sort_lfp.units[lfp_cluster][::10]*1E-6*compress_factor,1)
    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    

    

    plt.ylabel("Time (ms)", fontsize=30, fontweight='bold')

    plt.suptitle(os.path.split(self.parent.sua_file)[1]+" unit: "+self.starting_cell.text(), fontsize=25)
    
    
    ax = plt.subplot(2, 1, 2)
    
    img_out = img_out_shuffled[1000+start_window:1000+end_window]
    plt.imshow(img_out, aspect='auto', interpolation='none')
    
    
    #******************************** DISTRIBUTION PLOTS *******************************
    if False:
        f2 = plt.figure()
        bin_width = 1   #100ms bins

        #print maxloc_array
        #print len(maxloc_array)
        
        loc_distribution = []
        for k in range(len(maxloc_array)-1):
            if (maxloc_array[k]!=123) and (maxloc_array[k+1]!=123): 
                loc_distribution.append(maxloc_array[k+1]-maxloc_array[k])

        #print loc_distribution
        ax = plt.subplot(1,2,1)
        
        y = np.histogram(loc_distribution, bins = np.arange(0,100,bin_width))
        plt.bar(y[1][:-1], y[0], bin_width, color='blue',alpha=0.5)  
        
        
        #********************* SHUFFLE DISTRIUBTIONS *************
        
        #maxloc_array_shuffled = np.array(maxloc_array)
        #np.random.shuffle(maxloc_array_shuffled)
        loc_distribution = []
        for k in range(len(maxloc_array_shuffled)-1):
            if (maxloc_array_shuffled[k]!=123) and ((maxloc_array_shuffled[k+1]!=123)): 
                
                loc_distribution.append(maxloc_array_shuffled[k+1]-maxloc_array_shuffled[k])

        #print loc_distribution
        ax = plt.subplot(1,2,1)
        
        y = np.histogram(loc_distribution, bins = np.arange(0,100,bin_width))
        plt.bar(y[1][:-1], y[0], bin_width, color='red',alpha=0.5)    
        
        
        KS, p_val = stats.ks_2samp(maxloc_array, maxloc_array_shuffled)
       
        print KS, p_val
        
    
    
    
    #************************** PLOT AUTOCORRELATION FUNCTIONS ******************
    
    f3 =  plt.figure()


    # generate some data
    #x = np.arange(len(maxloc_array))
    x = np.int32(Sort_lfp.units[lfp_cluster]*1E-6*50)
    print len(x), x[-1]
    
    x_zeros = np.zeros(int(x[-1]+1), dtype=np.float32)
    print len(x_zeros)
    
    for k in range(len(x)):
        x_zeros[x[k]] = maxloc_array[k]
    
    plt.scatter(np.arange(len(x_zeros)), x_zeros)
    plt.show()
    quit()

    y = x_zeros

    # y = np.random.uniform(size=300)
    yunbiased = y-np.mean(y)
    ynorm = np.sum(yunbiased**2)
    acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
    # use only second half
    acor = acor[len(acor)/2:]

    plt.plot(acor)
        

    plt.show()



def plot_msl_continuous_single(self):

    lfp_cluster = int(self.parent.lfp_cluster.text())

    sig = int(self.sigma_width.text())


    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+str(sig)+'sigma'
    jitter_time = 50 #Time to jitter spiketrian
    #file_out_jittered = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+str(jitter_time)+"ms_jitter"
    #file_out_poisson = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+str(jitter_time)+"ms_window_poisson"
    #file_out_poisson_singles = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_poisson_singles"
    
    #Load lock times
    lock_time = np.load(file_out+'.npy')
    #lock_time_jittered = np.load(file_out_jittered+'.npy')
    #lock_time_poisson = np.load(file_out_poisson+'.npy')
    #lock_time_poisson_singles = np.load(file_out_poisson_singles+'.npy')

    compress_factor = 50

    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data
    lfp_cluster = int(self.parent.lfp_cluster.text())
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    print "... # of lfp events; ", len(pop_spikes)

    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    #tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))

    #Make plots for single cell raster
    unit = int(self.starting_cell.text())
    f1 = plt.figure()
    ax = plt.subplot(1, 1, 1)

    #************** Plot locked times rasters *********************
    even_locks=[]; odd_locks=[]
    for k in range(len(lock_time[unit])): 
        even_locks.append(lock_time[unit][k][0])
        odd_locks.append(lock_time[unit][k][1])
    
    #Plot event times scatter
    even_times = []
    for k in range(len(even_locks)):
        if even_locks[k]!=0.0:
            even_times.append([k, even_locks[k]])
    
    even_times = np.array(even_times).T
    #if len(even_times)>0: 
    #    plt.scatter(even_times[0], even_times[1], s = 10, color='blue', alpha=.6)
    
    #Plot odd times scatter
    odd_times = []
    for k in range(len(odd_locks)):
        if odd_locks[k]!=0.0:
            odd_times.append([k, odd_locks[k]])
    odd_times = np.array(odd_times).T
    #if len(odd_times)>0: 
    #    plt.scatter(odd_times[0], odd_times[1], s = 10, color='red', alpha=.6)
    
    #Plot ave times plot; search for locking times; if only an odd or even time exists, set that time; otherwise use average;
    ave_times = []
    error_array = []
    for k in range(len(even_locks)):
        temp = 0
        if even_locks[k]!=0: 
            temp+=even_locks[k]
            if odd_locks[k]!=0: 
                temp+=odd_locks[k]
                ave_times.append([k, temp/2.])
                error_array.append([min(even_locks[k],odd_locks[k]), max(even_locks[k],odd_locks[k])])
            else:
                ave_times.append([k, temp])
                error_array.append([even_locks[k], even_locks[k]])
        else:
            if odd_locks[k]!=0:
                ave_times.append([k, odd_locks[k]])
                error_array.append([odd_locks[k], odd_locks[k]])

    if len(ave_times)>0: 
        for k in range(len(ave_times)-1):
            
            #Plot magenta line average;
            if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                plt.plot([ave_times[k][0], ave_times[k+1][0]] , [ave_times[k][1], ave_times[k+1][1]], color='blue', linewidth = 8, alpha=.85)
                
                #Plot error bar at each point 
                x = np.arange(ave_times[k][0], ave_times[k+1][0],0.1)-0.5
                y1 = error_array[k][0]
                y2 = error_array[k][1]
                plt.fill_between(x, y1, y2, color='blue', alpha=0.3)

           
    #Track stats on original data
    diffs = []
    for k in range(len(even_locks)):
        if (even_locks[k]!=0) and (odd_locks[k]!=0):
            diff = abs(even_locks[k]-odd_locks[k])
            #if diff < 50: 
            diffs.append(diff)    #Only 
            
    
    if False:
        
        #************** Plot poisson times rasters *********************
        lock_time = lock_time_poisson
        even_locks=[]; odd_locks=[]
        for k in range(len(lock_time[unit])): 
            even_locks.append(lock_time[unit][k][0])
            odd_locks.append(lock_time[unit][k][1])
        
        #Plot event times scatter
        even_times = []
        for k in range(len(even_locks)):
            if even_locks[k]!=0.0:
                even_times.append([k, even_locks[k]])
        
        even_times = np.array(even_times).T
        #if len(even_times)>0: 
        #    plt.scatter(even_times[0], even_times[1], s = 10, color='green', alpha=.6)
        
        #Plot odd times scatter
        odd_times = []
        for k in range(len(odd_locks)):
            if odd_locks[k]!=0.0:
                odd_times.append([k, odd_locks[k]])
        odd_times = np.array(odd_times).T
        #if len(odd_times)>0: 
        #    plt.scatter(odd_times[0], odd_times[1], s = 10, color='cyan', alpha=.6)
        
        #Plot ave times plot
        ave_times = []
        error_array = []
        for k in range(len(even_locks)):
            temp = 0
            if even_locks[k]!=0: 
                temp+=even_locks[k]
                if odd_locks[k]!=0: 
                    temp+=odd_locks[k]
                    ave_times.append([k, temp/2.])
                    error_array.append([min(even_locks[k],odd_locks[k]), max(even_locks[k],odd_locks[k])])
                else:
                    ave_times.append([k, temp])
                    error_array.append([even_locks[k], even_locks[k]])

            else:
                if odd_locks[k]!=0:
                    ave_times.append([k, odd_locks[k]])
                    error_array.append([odd_locks[k], odd_locks[k]])

        ave_times =np.array(ave_times)
        ave_times[:,1]=ave_times[:,1]+60.0
        error_array =np.array(error_array)
        error_array=error_array+60.0
        if len(ave_times)>0: 
            for k in range(len(ave_times)-1):
                if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                    plt.plot([ave_times[k][0], ave_times[k+1][0]] , [ave_times[k][1], ave_times[k+1][1]], color='green', linewidth = 8, alpha=.85)
        
                    #Plot error bar at each point 
                    x = np.arange(ave_times[k][0], ave_times[k+1][0],0.1)-0.5
                    y1 = error_array[k][0]
                    y2 = error_array[k][1]
                    plt.fill_between(x, y1, y2, color='green', alpha=0.3)


        #Track stats on poisson data
        diffs_poisson = []
        for k in range(len(even_locks)):
            if (even_locks[k]!=0) and (odd_locks[k]!=0):
                diff = abs(even_locks[k]-odd_locks[k])
                #if diff < 50: 
                diffs_poisson.append(diff)        
                
                
        ##************** Plot poisson times rasters *********************
        #lock_time = lock_time_poisson_singles
        #even_locks=[]; odd_locks=[]
        #for k in range(len(lock_time[unit])): 
            #even_locks.append(lock_time[unit][k][0])
            #odd_locks.append(lock_time[unit][k][1])
        
        ##Plot event times scatter
        #even_times = []
        #for k in range(len(even_locks)):
            #if even_locks[k]!=0.0:
                #even_times.append([k, even_locks[k]])
        
        #even_times = np.array(even_times).T
        ##if len(even_times)>0: 
            ##plt.scatter(even_times[0], even_times[1], s = 10, color='green', alpha=.6)
        
        ##Plot odd times scatter
        #odd_times = []
        #for k in range(len(odd_locks)):
            #if odd_locks[k]!=0.0:
                #odd_times.append([k, odd_locks[k]])
        #odd_times = np.array(odd_times).T
        ##if len(odd_times)>0: 
            ##plt.scatter(odd_times[0], odd_times[1], s = 10, color='cyan', alpha=.6)
        
        ##Plot ave times plot
        #ave_times = []
        #for k in range(len(even_locks)):
            #temp = 0
            #if even_locks[k]!=0: 
                #temp+=even_locks[k]
                #if odd_locks[k]!=0: 
                    #temp+=odd_locks[k]
                    #ave_times.append([k, temp/2.])
                #else:
                    #ave_times.append([k, temp])
            #else:
                #if odd_locks[k]!=0:
                    #ave_times.append([k, odd_locks[k]])

        #if len(ave_times)>0: 
            #for k in range(len(ave_times)-1):
                #if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                    #plt.plot([ave_times[k][0],ave_times[k+1][0]] , [ave_times[k][1], ave_times[k+1][1]], color='brown', linewidth = 5, alpha=.3)
        
        
    
    #**********Plot labels and additional info
    old_label = np.arange(50, 170, 10)
    new_label = [-50,-40,-30,-20,-10,0,-50, -40,-30,-20,-10,0]
    plt.yticks(old_label,new_label, fontsize=25)
    
    #plt.plot([0, tsf.n_vd_samples/tsf.SampleFrequency/60.], [100,100], 'r--', color='black', linewidth=4, alpha=.5)
    #plt.plot([0, tsf.n_vd_samples/tsf.SampleFrequency/60.], [160,160], 'r--', color='black', linewidth=4, alpha=.5)
    
    win_len = self.sliding_window_length.text()      #Work in ms
    plt.title("Unit: "+str(unit) + ",   #spks: "+str(len(Sort_sua.units[unit]))+",    sliding window: "+win_len+" mins.\n Real data: ave(diff): " + str(round(np.mean(diffs),2))+ "   std(diff): " + str(round(np.std(diffs),2))
                ,      fontsize=25)
    ax.tick_params(axis='both', which='both', labelsize=25)
    plt.ylim(40,170)
    #plt.xlim(0, tsf.n_vd_samples/tsf.SampleFrequency/60.)

    plt.ylabel("Locking Latency (ms)\n Real Data           Poisson Data", fontsize = 30)
    plt.xlabel("Recording time (mins)", fontsize = 25)
    
    #Plot LFP Events
    spikes = pop_spikes*1E-6/60.
    ymin=np.zeros(len(spikes))
    ymax=np.zeros(len(spikes))
    ymin+=40
    ymax+=44
    plt.vlines(spikes, ymin, ymax, linewidth=1, color='brown',alpha=.5) #colors[mod(counter,7)])

    #Plot Single Spike Events
    spikes = Sort_sua.units[unit]*1E-6/60.
    ymin=np.zeros(len(spikes))
    ymax=np.zeros(len(spikes))
    ymin+=45    
    ymax+=49
    plt.vlines(spikes, ymin, ymax, linewidth=1, color='blue',alpha=.5) #colors[mod(counter,7)])

    
    #ax.set_xticklabels([])
    #ax.set_yticklabels([])
        
    plt.show()
    

def plot_msl_continuous_multi_unit(self):

    colors = ['blue', 'red', 'green', 'magenta', 'brown', 'orange', 'cyan', 'pink', 'grey', 'indigo']

   
    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    
    sig = int(self.sigma_width.text())

    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+str(sig)+'sigma'
    
    #Load lock times
    lock_time = np.load(file_out+'.npy')
    

    #file_out_locked_spikes = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_locked_spikes"
    #locked_spike_array = np.load(file_out_locked_spikes+'.npy') #This is al ist of all spikes locking to each 30min window; used for non-Luczak type studies;


    file_spikerates = file_out+"_spikerates"

    spikerates = np.load(file_spikerates+".npy")

    #Load lfp event rasters
    compress_factor = 50
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data
    lfp_cluster = int(self.parent.lfp_cluster.text())
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    print "... # of lfp events; ", len(pop_spikes)


    #Load single unit rasters
    Sort_sua = PTCS.PTCS(self.parent.sua_file) #Auto load flag for Nick's data
    #Select Singel Units to process
    if self.multiple_units.text()=='00':
        units = np.arange(Sort_sua.n_units)
    else:
        units = np.int16(self.multiple_units.text().split(","))


    ##Compute periods of synchrony from si index
    lfp = TSF.TSF(self.parent.lfp_tsf_file)
    #lfp = TSF.TSF(self.parent.sua_file.replace('_hp_butter_alltrack.ptcs','.tsf').replace('hp','lp'))
    lfp.read_ec_traces()    #Ok to read all LFP, smaller files


    #*************************** PLOT RASTERS ***********************
    ax = plt.subplot(1, 1, 1)

    # define the colormap
    cmap = plt.cm.get_cmap('viridis', Sort_sua.n_units)
    plt.set_cmap('viridis')

    #Plot Unit Rasters
    raster_offset = -40
    clr_ctr=0
    #for clr_ctr, unit in enumerate(units):
        #spikes = Sort_sua.units[unit]*1E-6/60.
        #ymin=np.zeros(len(spikes))
        #ymax=np.zeros(len(spikes))
        #ymin+=raster_offset
        #ymax+=raster_offset-5
        #plt.vlines(spikes, ymin, ymax, linewidth=1, color = cm.viridis(int(float(unit)/Sort_sua.n_units*256)), alpha=0.5) #colors[mod(counter,7)])
        
        #raster_offset = raster_offset-clr_ctr*5
    
    #Plot LFP Events
    raster_offset-=3
    spikes = pop_spikes*1E-6/60.
    ymin=np.zeros(len(spikes))
    ymax=np.zeros(len(spikes))
    ymin+=raster_offset-4-clr_ctr*3
    ymax+=raster_offset-14-clr_ctr*3
    plt.vlines(spikes, ymin, ymax, linewidth=1, color='black',alpha=.5) #colors[mod(counter,7)])


    #*************************** PLOT MSL DRIFT DATA ***********************
    Luczak_method = True
    std_method = False
    poisson_method = False
    
    
    #Make plots for single cell raster
    if Luczak_method: 
        #clr_ctr=0

        for clr_ctr, unit in enumerate(units):
            print "... unit: ", unit
            ax = plt.subplot(1, 1, 1)

            #************** Plot locked times rasters *********************
            even_locks=[]; odd_locks=[]
            for k in range(len(lock_time[unit])): 
                even_locks.append(lock_time[unit][k][0])
                odd_locks.append(lock_time[unit][k][1])
            
            #Plot event times scatter
            even_times = []
            for k in range(len(even_locks)):
                if even_locks[k]!=0.0:
                    even_times.append([k, even_locks[k]])
            
            even_times = np.array(even_times).T
            #if len(even_times)>0: 
            #    plt.scatter(even_times[0], even_times[1], s = 10, color='blue', alpha=.6)
            
            #Plot odd times scatter
            odd_times = []
            for k in range(len(odd_locks)):
                if odd_locks[k]!=0.0:
                    odd_times.append([k, odd_locks[k]])
            odd_times = np.array(odd_times).T
            #if len(odd_times)>0: 
            #    plt.scatter(odd_times[0], odd_times[1], s = 10, color='red', alpha=.6)
            
            #Plot ave times plot; search for locking times; if only an odd or even time exists, set that time; otherwise use average;
            ave_times = []
            error_array = []
            skip_unit = False
            for k in range(len(even_locks)):
                temp = 0
                if even_locks[k]!=0: 
                    temp+=even_locks[k]
                    if odd_locks[k]!=0: 
                        temp+=odd_locks[k]
                        ave_times.append([k, temp/2.])
                        error_array.append([min(even_locks[k],odd_locks[k]), max(even_locks[k],odd_locks[k])])
                    else:
                        ave_times.append([k, temp])
                        error_array.append([even_locks[k], even_locks[k]])
                else:
                    if odd_locks[k]!=0:
                        ave_times.append([k, odd_locks[k]])
                        error_array.append([odd_locks[k], odd_locks[k]])
                    else:
                        error_array.append([odd_locks[k], odd_locks[k]])
                   
                if abs(error_array[k][0]-error_array[k][1])>30: 
                    skip_unit = True
            
            if skip_unit: continue
                
            offset_temp =int(self.sliding_window_length.text())/2   #OFFSET DATA BY 1/2 WINDOW LENGTH TO REFLECT AVERAGING TO MIDDLE OF WINDOW;
               
            if len(ave_times)>0: 
                for k in range(len(ave_times)-1):
                    
                    #Plot line average;
                    if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                        plt.plot([ave_times[k][0]*int(self.sliding_window_step.text())+offset_temp, ave_times[k+1][0]*int(self.sliding_window_step.text())+offset_temp],    #X coordinates
                                 [ave_times[k][1]-100, ave_times[k+1][1]-100],                                                                                              #Y coordinates
                                 color = cm.viridis(int(float(unit)/Sort_sua.n_units*256)), linewidth = 5, alpha=.75)
                       
                        #Plot error; LUCZAK METHOD: odd vs. even location times.
                        if True: 
                            x = np.arange(ave_times[k][0]*int(self.sliding_window_step.text()), ave_times[k+1][0]*int(self.sliding_window_step.text()),0.1)-0.5+offset_temp
                            y1 = error_array[k][0]-100
                            y2 = error_array[k][1]-100
                            plt.fill_between(x, y1, y2, color = cm.viridis(int(float(unit)/Sort_sua.n_units*256)), alpha=0.4)

                   

        #if len(ave_times)>0: 
            #for k in range(len(ave_times)-1):
                
                ##Plot magenta line average;
                #if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                    #plt.plot([ave_times[k][0], ave_times[k+1][0]] , [ave_times[k][1], ave_times[k+1][1]], color=colors[unit%8], linewidth = 8, alpha=.85)
                    
                    ##Plot error bar at each point 
                    #x = np.arange(ave_times[k][0], ave_times[k+1][0],0.1)-0.5
                    #y1 = error_array[k][0]
                    #y2 = error_array[k][1]
                    #plt.fill_between(x, y1, y2, color='blue', alpha=0.3)

               
               
                           ##Track stats on original data
            #if False:
                #diffs = []
                #for k in range(len(even_locks)):
                    #if (even_locks[k]!=0) and (odd_locks[k]!=0):
                        #diff = abs(even_locks[k]-odd_locks[k])
                        ##if diff < 50: 
                        #diffs.append(diff)    #Only 
        
            ##if clr_ctr==10: break
            ##clr_ctr+=1


    #COMPUTE MEAN AND STD OF DISTRIBUTION DIRECTLY  
    if std_method:
        for k in range(len(locked_spike_array)):
            #print "... computing epoch: ", k
            
            for clr_ctr, unit in enumerate(units):
                if len(locked_spike_array[k][unit])>0: 
    
                    spikes = np.hstack(locked_spike_array[k][unit])*1E-3
                    indexes = np.where(np.logical_and(spikes>=-50, spikes<=50))[0]   #Work in milliseconds
                    if len(indexes)>0: 
                        print indexes
                        spikes = spikes[indexes]
                        temp_ave = np.mean(spikes, axis=0)
                        temp_std = np.std(spikes, axis=0)
                        plt.scatter(k, temp_ave, color=colors[clr_ctr%10], s = 30, alpha=.85)
                        plt.plot([k,k], [temp_ave-temp_std,temp_ave+temp_std], color=colors[clr_ctr%10], linewidth=2, alpha=0.3)

    
    if poisson_method:
        
        #************** Plot poisson times rasters *********************
        lock_time = lock_time_poisson
        even_locks=[]; odd_locks=[]
        for k in range(len(lock_time[unit])): 
            even_locks.append(lock_time[unit][k][0])
            odd_locks.append(lock_time[unit][k][1])
        
        #Plot event times scatter
        even_times = []
        for k in range(len(even_locks)):
            if even_locks[k]!=0.0:
                even_times.append([k, even_locks[k]])
        
        even_times = np.array(even_times).T
        #if len(even_times)>0: 
        #    plt.scatter(even_times[0], even_times[1], s = 10, color='green', alpha=.6)
        
        #Plot odd times scatter
        odd_times = []
        for k in range(len(odd_locks)):
            if odd_locks[k]!=0.0:
                odd_times.append([k, odd_locks[k]])
        odd_times = np.array(odd_times).T
        #if len(odd_times)>0: 
        #    plt.scatter(odd_times[0], odd_times[1], s = 10, color='cyan', alpha=.6)
        
        #Plot ave times plot
        ave_times = []
        error_array = []
        for k in range(len(even_locks)):
            temp = 0
            if even_locks[k]!=0: 
                temp+=even_locks[k]
                if odd_locks[k]!=0: 
                    temp+=odd_locks[k]
                    ave_times.append([k, temp/2.])
                    error_array.append([min(even_locks[k],odd_locks[k]), max(even_locks[k],odd_locks[k])])
                else:
                    ave_times.append([k, temp])
                    error_array.append([even_locks[k], even_locks[k]])

            else:
                if odd_locks[k]!=0:
                    ave_times.append([k, odd_locks[k]])
                    error_array.append([odd_locks[k], odd_locks[k]])

        ave_times =np.array(ave_times)
        ave_times[:,1]=ave_times[:,1]+60.0
        error_array =np.array(error_array)
        error_array=error_array+60.0
        if len(ave_times)>0: 
            for k in range(len(ave_times)-1):
                if (ave_times[k+1][0] - ave_times[k][0])==1:  #Only plot lines between consecutive times
                    plt.plot([ave_times[k][0], ave_times[k+1][0]] , [ave_times[k][1], ave_times[k+1][1]], color='green', linewidth = 8, alpha=.85)
        
                    #Plot error bar at each point 
                    x = np.arange(ave_times[k][0], ave_times[k+1][0],0.1)-0.5
                    y1 = error_array[k][0]
                    y2 = error_array[k][1]
                    plt.fill_between(x, y1, y2, color='green', alpha=0.3)


        #Track stats on poisson data
        diffs_poisson = []
        for k in range(len(even_locks)):
            if (even_locks[k]!=0) and (odd_locks[k]!=0):
                diff = abs(even_locks[k]-odd_locks[k])
                #if diff < 50: 
                diffs_poisson.append(diff)        
                
                
    #**********Plot labels and additional info
    #old_label = np.arange(50, 170, 10)
    #new_label = [-50,-40,-30,-20,-10,0,-50, -40,-30,-20,-10,0]
    #plt.yticks(old_label,new_label, fontsize=25)
    
    #plt.plot([0, tsf.n_vd_samples/tsf.SampleFrequency/60.], [100,100], 'r--', color='black', linewidth=4, alpha=.5)
    #plt.plot([0, tsf.n_vd_samples/tsf.SampleFrequency/60.], [160,160], 'r--', color='black', linewidth=4, alpha=.5)
    
    win_len = self.sliding_window_length.text()      #Work in ms
    ax.tick_params(axis='both', which='both', labelsize=25)
    
    #plt.ylim(raster_offset-7-clr_ctr*3,20)
    #plt.ylim(raster_offset-7-clr_ctr*3,100)
    plt.ylim(-30, 50)
    
    #plt.xlim(0, lfp.n_vd_samples/lfp.SampleFrequency/60.)
    plt.xlim(0, 320.)

    plt.ylabel("Mean-Spike-Latency (ms)\n Single Units   LFP Events", fontsize = 30)
    plt.xlabel("Recording time (mins)", fontsize = 25)

    
    #ax.set_xticklabels([])
    #ax.set_yticklabels([])
        
        
    plt.suptitle(os.path.split(self.parent.sua_file)[1], fontsize=25)
    plt.show()
    

def msl_continuous_parallel(chunk_ctr, time_chunks, pop_spikes, locked_spike_array, Sort_sua, win_len, lock_window_start, lock_window_end, cell_rasters, cell_rasters_poisson, sig):
    
    print "*******************************", chunk_ctr
    time_chunk = time_chunks[chunk_ctr]
    temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
    print time_chunk
    
    lock_time = []
    lock_time_poisson = []
    locked_spike_array = []
    
    for k in range(len(Sort_sua.units)):
        lock_time.append([])
        lock_time_poisson.append([])
        locked_spike_array.append([])
    
    #locked_spike_array.append([])
    for unit in range(len(Sort_sua.units)):
        print chunk_ctr, unit
        if (len(temp3)/(win_len*1E-3))<0.01:            #Exclude periods with LFP rates < 0.01 Hz
            lock_time[unit].append([0,0])
            lock_time_poisson[unit].append([0,0])
            locked_spike_array[chunk_ctr].append([])

            continue

        locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])

        if (len(locked_spikes)/(win_len*1E-3))<0.01:    #Exclude cells during periods with firing rates < 0.01 Hz
            lock_time[unit].append([0,0])
            lock_time_poisson[unit].append([0,0])
            locked_spike_array[chunk_ctr].append([])

            continue
    

        locked_spike_array[chunk_ctr].append(locked_spikes)
        #locked_spike_array[ctr].append(np.array(cell_rasters[unit])[temp3])

        #***********************************************************
        #*********** Compute Lock Times ****************************
        #***********************************************************

        fit_even = np.zeros(2000, dtype=np.float32)
        fit_odd = np.zeros(2000, dtype=np.float32)
        x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
        sig_gaussian = np.float32(gaussian(x, 0, sig))
        for g in range(len(locked_spikes)):
            mu = int(locked_spikes[g]*1E-3)
            if g%2==0: fit_even += np.roll(sig_gaussian, mu)
            else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

        #Search for peaks in PETH within the locking window
        if np.max(fit_even[1000+lock_window_start:1000+lock_window_end])>0:
            even_lock = np.argmax(fit_even[1000+lock_window_start:1000+lock_window_end])
        else:
            even_lock = 0

        if np.max(fit_odd[1000+lock_window_start:1000+lock_window_end])>0:
            odd_lock = np.argmax(fit_odd[1000+lock_window_start:1000+lock_window_end])
        else:
            odd_lock = 0

        lock_time[unit].append([even_lock, odd_lock])


        #***********************************************************
        #************ Compute bursty poisson process Lock Times *******************
        #***********************************************************
        #locked_spikes_poisson = np.random.poisson(np.random.randint(200), len(locked_spikes))-100       #NB: THIS IS ALREADY IN MS
        #locked_spikes_poisson = np.sort(np.random.poisson(10, len(locked_spikes))+(np.random.randint(jitter_time)-jitter_time/2.))  #Make sure spikes are time sorted
        locked_spikes_poisson = np.hstack(np.array(cell_rasters_poisson[unit])[temp3])
        #locked_spikes_poisson = np.unique(np.sort(locked_spikes_poisson))
        
        fit_even = np.zeros(2000, dtype=np.float32)
        fit_odd = np.zeros(2000, dtype=np.float32)
        x = np.linspace(-1000,1000, 2000)                   #Make an array from -1000ms .. +1000ms with microsecond precision
        sig_gaussian = np.float32(gaussian(x, 0, sig))
        for g in range(len(locked_spikes_poisson)):
            mu = int(locked_spikes_poisson[g])                  #Already in ms
            if g%2==0: fit_even += np.roll(sig_gaussian, mu)
            else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

        #Search for peaks in PETH within the locking window
        if np.max(fit_even[1000+lock_window_start:1000+lock_window_end])>0:
            even_lock = np.argmax(fit_even[1000+lock_window_start:1000+lock_window_end])
        else:
            even_lock = 0

        if np.max(fit_odd[1000+lock_window_start:1000+lock_window_end])>0:
            odd_lock = np.argmax(fit_odd[1000+lock_window_start:1000+lock_window_end])
        else:
            odd_lock = 0

        lock_time_poisson[unit].append([even_lock, odd_lock])
        
    #return zip(lock_time, lock_time_poisson, locked_spike_array)
    return locked_spike_array

def compute_msl_spikerates(self):
    
    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1000000  # window width in usec
    
    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())

    lfp_cluster = int(self.parent.lfp_cluster.text())
    #tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))

    sig = float(self.sigma_width.text())
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())

    #Load SUA Sort
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    #Load LFP Sort
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data
    
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardcoded to: ", compress_factor

    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())
      
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    pop_spikes=np.sort(np.unique(pop_spikes))*1E-3   #Exclude duplicates; convert to milisecond time
    #pop_spikes=(pop_spikes)*1E-3   #Exclude duplicates; convert to milisecond time
    
    print " ... # LFP events: ", len(pop_spikes)
    
    #Load saved cell rasters
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")
    
    
    ##Compute periods of synchrony from si index                    #***********************************REIMPLEMENT ASAP
    #lfp = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf').replace('hp','lp'))
    #lfp.read_ec_traces()    #Ok to read all LFP, smaller files
    lfp = TSF.TSF(self.parent.lfp_tsf_file)
    lfp.read_ec_traces()
    sync_ch=9

    sync_periods_file = self.parent.sua_file.replace('.ptcs','')+"_sync_periods_ch"+str(sync_ch)+'.txt'
    if os.path.exists(sync_periods_file):
        sync_periods = np.loadtxt(sync_periods_file)
    else:
        si, t, sync_periods = synchrony_index(lfp.ec_traces[sync_ch], lfp.SampleFrequency, si_limit)
        np.savetxt(sync_periods_file, sync_periods)

    temp_list = []
    for p in range(len(sync_periods)):
        indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0]*1E3, pop_spikes<=sync_periods[p][1]*1E3))[0]   #Work in milliseconds
        temp_list.extend(pop_spikes[indexes])
    pop_spikes=np.array(temp_list)
    
    print "... # LFP events during sync states: ", len(pop_spikes)
    
    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    
    #OPTION 1: Divide into chunks of recording length
    win_len = float(self.sliding_window_length.text())*60000.0                          #Work in ms
    step_len = float(self.sliding_window_step.text())*60000.0                           #Work in ms
    #rec_len = tsf.n_vd_samples/float(tsf.SampleFrequency)*1E3   #Work in ms
    rec_len = lfp.n_vd_samples/float(lfp.SampleFrequency)*1E3   #Work in ms         #USE LFP RECORDING IF High-pass not available
    print "... lfp rec_len: ", rec_len
    
    
    temp_chunks = np.arange(0, rec_len-win_len, step_len)   #set begginings of each chunk up to end of recording (less the window length)
    time_chunks = []
    for k in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[k], temp_chunks[k]+win_len])
    time_chunks.append([temp_chunks[-1], rec_len])

    #**************************************************************************
    #*************************** LOOP OVER TIME CHUNKS ************************
    #**************************************************************************

    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+ self.sigma_width.text()+"sigma"

    chunk_index = np.arange(0,len(time_chunks),1).tolist()
    
    spiking_rate_array = []
    for chunk_ctr in chunk_index:
        time_chunk = time_chunks[chunk_ctr]
        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]       #indexes of spikes in sliding window, for example 60mins; 
        print len(temp3), self.min_fire_rate.text(), len(temp3)/(win_len*1.0E-3)
        
        spiking_rate_array.append([])
        for unit in range(len(Sort_sua.units)):
            if (len(temp3)/(win_len*1.0E-3))< float(self.min_fire_rate.text()):            #Exclude periods with LFP spike rates < 0.01 Hz; NB: win_len is in milliseconds;
                spiking_rate_array[chunk_ctr].append(0)
                #print "...excluding epoch: ", chunk_ctr,
                continue

            locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])

            spiking_rate_array[chunk_ctr].append(len(locked_spikes)/(win_len*1E-3)) #append spiking rate for each epoch for each cell;
        print ''
    np.save(file_out+"_spikerates", spiking_rate_array)

    print "... done saving spikerates matrix..."
    
def compute_msl_continuous(self):
    
    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1000000  # window width in usec
    
    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())

    lfp_cluster = int(self.parent.lfp_cluster.text())
    #tsf = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf'))

    sig = float(self.sigma_width.text())
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())

    #Load SUA Sort
    Sort_sua = PTCS.PTCS(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    #Load LFP Sort
    Sort_lfp = PTCS.PTCS(self.parent.lfp_event_file) #Auto load flag for Nick's data
    
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardcoded to: ", compress_factor

    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())
      
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    pop_spikes=np.sort(np.unique(pop_spikes))*1E-3   #Exclude duplicates; convert to milisecond time
    #pop_spikes=(pop_spikes)*1E-3   #Exclude duplicates; convert to milisecond time
    
    print " ... # LFP events: ", len(pop_spikes)
    
    #Load saved cell rasters
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")
    
    
    ##Compute periods of synchrony from si index                    #***********************************REIMPLEMENT ASAP
    #lfp = Tsf_file(self.parent.sua_file.replace('.ptcs','.tsf').replace('hp','lp'))
    #lfp.read_ec_traces()    #Ok to read all LFP, smaller files
    lfp_filename = self.parent.lfp_tsf_file
    lfp = TSF.TSF(lfp_filename)
    lfp.read_ec_traces()
    sync_ch=9




    sync_periods_file = self.parent.sua_file.replace('.ptcs','')+"_sync_periods_ch"+str(sync_ch)+'.txt'
    if os.path.exists(sync_periods_file):
        sync_periods = np.loadtxt(sync_periods_file)
    else:
        si, t, sync_periods = synchrony_index(lfp.ec_traces[sync_ch], lfp.SampleFrequency, si_limit)
        np.savetxt(sync_periods_file, sync_periods)

    temp_list = []
    for p in range(len(sync_periods)):
        indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0]*1E3, pop_spikes<=sync_periods[p][1]*1E3))[0]   #Work in milliseconds
        temp_list.extend(pop_spikes[indexes])
    pop_spikes=np.array(temp_list)
    
    print "... # LFP events during sync states: ", len(pop_spikes)
    
    #**************************************************************************
    #********* CHUNK UP TIME - 3 OPTIONS: TIME, # SPIKES, # EVENTS ************
    #**************************************************************************
    
    #OPTION 1: Divide into chunks of recording length
    win_len = float(self.sliding_window_length.text())*60000.0                          #Work in ms
    step_len = float(self.sliding_window_step.text())*60000.0                           #Work in ms
    #rec_len = tsf.n_vd_samples/float(tsf.SampleFrequency)*1E3   #Work in ms
    rec_len = lfp.n_vd_samples/float(lfp.SampleFrequency)*1E3   #Work in ms         #USE LFP RECORDING IF High-pass not available
    print "... lfp rec_len: ", rec_len
    
    
    temp_chunks = np.arange(0, rec_len-win_len, step_len)   #set begginings of each chunk up to end of recording (less the window length)
    time_chunks = []
    for k in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[k], temp_chunks[k]+win_len])
    time_chunks.append([temp_chunks[-1], rec_len])

    #**************************************************************************
    #*************************** LOOP OVER TIME CHUNKS ************************
    #**************************************************************************

    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+ self.sigma_width.text()+"sigma"

    jitter_time = 50 #Time to jitter spiketrian
    file_out_poisson = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_"+str(jitter_time)+"ms_window_poisson_"+ self.sigma_width.text()+"sigma"

    file_out_locked_spikes = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_locked_spikes_"+ self.sigma_width.text()+"sigma"

    locked_spike_array = [] #This is al ist of all spikes locking to each window; used for non-Luczak type studies;


    cell_rasters_poisson = []
    for unit in range(len(Sort_sua.units)):
        cell_rasters_poisson.append([])
        for p in range(len(cell_rasters[unit])):
            #poisson_spikes = np.random.poisson(10, len(cell_rasters[unit][p]))+np.random.randint(jitter_time)-jitter_time/2.
            #poisson_spikes = np.random.poisson(np.random.randint(jitter_time), len(cell_rasters[unit][p])) -jitter_time
            poisson_spikes = np.random.poisson(10, len(cell_rasters[unit][p])) - np.random.randint(jitter_time)
            cell_rasters_poisson[unit].append(poisson_spikes)
        
            #if len(cell_rasters[unit][p])>0: 
                #print cell_rasters[unit][p]
                #print cell_rasters_poisson[unit][p]
                #print ""
        #return
    
    #if os.path.exists(file_out+'.npy'): 
    
    if False: 
        lock_time = np.load(file_out+'.npy')
    else:
        img=[]; lock_time=[]; lock_time_jittered=[]; lock_time_shifted=[]; lock_time_poisson=[]; lock_time_poisson_singles=[]
        for k in range(len(Sort_sua.units)): 
            lock_time.append([]); lock_time_jittered.append([]); lock_time_shifted.append([]); lock_time_poisson.append([]); lock_time_poisson_singles.append([])
        chunk_index = np.arange(0,len(time_chunks),1)
    
        for ctr, chunk_ctr in enumerate(chunk_index):
            time_chunk = time_chunks[chunk_ctr]
            temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
            print time_chunk
            
            locked_spike_array.append([])
            for unit in range(len(Sort_sua.units)):
                if (len(temp3)/(win_len*1E-3))<0.01:            #Exclude periods with LFP rates < 0.01 Hz
                    lock_time[unit].append([0,0])
                    lock_time_jittered[unit].append([0,0])
                    lock_time_shifted[unit].append([0,0])
                    lock_time_poisson[unit].append([0,0])
                    lock_time_poisson_singles[unit].append([0,0])
        
                    locked_spike_array[ctr].append([])

                    continue

                locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])

                if (len(locked_spikes)/(win_len*1E-3))<0.01:    #Exclude periods with firing rates < 0.01 Hz
                    lock_time[unit].append([0,0])
                    lock_time_jittered[unit].append([0,0])
                    lock_time_shifted[unit].append([0,0])
                    lock_time_poisson[unit].append([0,0])
                    lock_time_poisson_singles[unit].append([0,0])

                    locked_spike_array[ctr].append([])

                    continue
            

                locked_spike_array[ctr].append(locked_spikes)
                #locked_spike_array[ctr].append(np.array(cell_rasters[unit])[temp3])

                #***********************************************************
                #*********** Compute Lock Times ****************************
                #***********************************************************

                fit_even = np.zeros(2000, dtype=np.float32)
                fit_odd = np.zeros(2000, dtype=np.float32)
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(locked_spikes)):
                    mu = int(locked_spikes[g]*1E-3)
                    if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                    else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

                #Search for peaks in PETH within the locking window
                if np.max(fit_even[1000+lock_window_start:1000+lock_window_end])>0:
                    even_lock = np.argmax(fit_even[1000+lock_window_start:1000+lock_window_end])
                else:
                    even_lock = 0

                if np.max(fit_odd[1000+lock_window_start:1000+lock_window_end])>0:
                    odd_lock = np.argmax(fit_odd[1000+lock_window_start:1000+lock_window_end])
                else:
                    odd_lock = 0

                lock_time[unit].append([even_lock, odd_lock])


                #***********************************************************
                #************ Compute bursty poisson process Lock Times *******************
                #***********************************************************
                #locked_spikes_poisson = np.random.poisson(np.random.randint(200), len(locked_spikes))-100       #NB: THIS IS ALREADY IN MS
                #locked_spikes_poisson = np.sort(np.random.poisson(10, len(locked_spikes))+(np.random.randint(jitter_time)-jitter_time/2.))  #Make sure spikes are time sorted
                locked_spikes_poisson = np.hstack(np.array(cell_rasters_poisson[unit])[temp3])
                #locked_spikes_poisson = np.unique(np.sort(locked_spikes_poisson))
                
                fit_even = np.zeros(2000, dtype=np.float32)
                fit_odd = np.zeros(2000, dtype=np.float32)
                x = np.linspace(-1000,1000, 2000)                   #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(locked_spikes_poisson)):
                    mu = int(locked_spikes_poisson[g])                  #Already in ms
                    if g%2==0: fit_even += np.roll(sig_gaussian, mu)
                    else: fit_odd += np.roll(sig_gaussian, mu , axis=0)

                #Search for peaks in PETH within the locking window
                if np.max(fit_even[1000+lock_window_start:1000+lock_window_end])>0:
                    even_lock = np.argmax(fit_even[1000+lock_window_start:1000+lock_window_end])
                else:
                    even_lock = 0

                if np.max(fit_odd[1000+lock_window_start:1000+lock_window_end])>0:
                    odd_lock = np.argmax(fit_odd[1000+lock_window_start:1000+lock_window_end])
                else:
                    odd_lock = 0

                lock_time_poisson[unit].append([even_lock, odd_lock])
                
        np.save(file_out, lock_time)
        #np.save(file_out_jittered, lock_time_jittered)
        np.save(file_out_poisson, lock_time_poisson)
        #np.save(file_out_poisson_singles, lock_time_poisson_singles)
        
        np.save(file_out_locked_spikes, locked_spike_array)
    
        
    for unit in range(len(Sort_sua.units)):
        ax = plt.subplot(11, 11, unit+1)
    #for unit in range(20):
    #    ax = plt.subplot(4, 5, unit+1)

        even_locks=[]; odd_locks=[]
        for k in range(len(lock_time[unit])): 
            even_locks.append(lock_time[unit][k][0])
            odd_locks.append(lock_time[unit][k][1])

        times = []
        for k in range(len(even_locks)):
            if even_locks[k]!=0.0:
                times.append([k, even_locks[k]])
        
        times = np.array(times).T
        if len(times)>0: 
            plt.scatter(times[0], times[1], s = 10, color='blue', alpha=.6)
        
        times = []
        for k in range(len(odd_locks)):
            if odd_locks[k]!=0.0:
                times.append([k, odd_locks[k]])
        times = np.array(times).T
        if len(times)>0: 
            plt.scatter(times[0], times[1], s = 10, color='red', alpha=.6)
        

        diffs = []
        for k in range(len(even_locks)):
            if (even_locks[k]!=0) and (odd_locks[k]!=0):
                diff = abs(even_locks[k]-odd_locks[k])
                if diff < 50: diffs.append(diff)
        #print diffs
        
        plt.title(str(unit) + " "+str(len(Sort_sua.units[unit]))+" " + str(round(np.mean(diffs),2))+ " " + str(round(np.std(diffs),2)), 
        fontsize=8, y=.95)
        ax.tick_params(axis='both', which='both', labelsize=8)
        plt.ylim(50,110)
        plt.xlim(0, len(time_chunks))
        
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
    plt.show()

def parallel_pval(p, spk1, k, lfp_period_indexes, locked_spike_array, lock_window_start, lock_window_end):
             #epochs, spk1, lfp_period_indexes, locked_spike_array, unit, self,
    
    #print "...epoch #1: ", k, " epoch #2: ", p


    if len(lfp_period_indexes[p])==0:
        return 1.0

    spk2 = np.hstack(locked_spike_array[lfp_period_indexes[p]])*1E-3

    indexes = np.where(np.logical_and(spk2>=lock_window_start, spk2<=lock_window_end))[0]   #Work in milliseconds            #********************* ONLY USING SPIKES WITHIN 100ms of t=0 for KS TEST

    if len(indexes)==0:
        return 1.0 #DESYNCH STATES; APPEND DUMMY VALUE

    spk2 = np.sort(spk2[indexes])

    KS, p_val = stats.ks_2samp(spk1, spk2)
    return p_val




def compute_window_pval(self):
    
    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    sort_lfp = PTCS.PTCS(self.parent.lfp_event_file)
    lfp_spikes = sort_lfp.units[lfp_cluster]*1E-6*50

    #Load cell rasters from scratch
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)+"_"+self.sigma_width.text()+'ms_histograms'
    locked_spike_array = np.load(cell_rasters_filename+'.npy', mmap_mode='c')
    print locked_spike_array.shape


    unit = int(self.starting_cell.text())

    #******************* PVALUE MATRIX ***************
    lfp_period_indexes = []
    for k in range(320):
        lfp_period_indexes.append(np.where(np.logical_and(lfp_spikes>=k*60, lfp_spikes<=(k+30)*60))[0])         #Find indexes of lfp events in sliding 30 minute period 

    #COMPUTE KS P-VALUES FOR CURRENT LFP CLUSTER
    n_epochs = 320
    cell_pval_filename = cell_rasters_filename+"_unit"+self.starting_cell.text()+"_"+str(n_epochs)+"epochs_pvals.npy"

    empty_row = np.zeros(len(locked_spike_array),dtype=np.float32)+1.0
    if os.path.exists(cell_pval_filename)==False:
        kstest_array = []
        #for k in range(len(locked_spike_array[unit])/30):
        #for k in range(10):
        for k in range(n_epochs):
            print "...epoch: ", k
            kstest_array.append([])    
            
            #Find spikes in base epoch
            #lfp_period_indexes = np.where(np.logical_and(lfp_spikes>=k*60, lfp_spikes<=(k+30)*60))[0]         #Find indexes of lfp events in sliding 30 minute period 
            if len(lfp_period_indexes[k])==0:
                kstest_array[k] = empty_row     #If no more events
                continue
            
            spk1 = np.hstack(locked_spike_array[unit][lfp_period_indexes[k]])*1E-3
            indexes = np.where(np.logical_and(spk1>=int(self.parent.lock_window_start.text()), spk1<=int(self.parent.lock_window_end.text())))[0]   #Work in milliseconds            #********************* ONLY USING SPIKES WITHIN 100ms of t=0 for KS TEST
            if len(indexes)==0:
                kstest_array[k] = empty_row             #DESYNCH STATES; APPEND DUMMY VALUE
                continue

            spk1 = np.sort(spk1[indexes])
            
            #Find spikes in all other epochs
            #for p in range(len(locked_spike_array[unit])/30):
            #for p in range(10):
            import parmap
            if True: 
                epochs=np.arange(n_epochs).tolist()
                #print epochs
                kstest_array[k] = parmap.map(parallel_pval, epochs, spk1, k, lfp_period_indexes, locked_spike_array[unit], 
                                  int(self.parent.lock_window_start.text()), int(self.parent.lock_window_end.text()), processes=20)

            else:

                for p in range(n_epochs):

                    print "... unit: ", unit, " epoch: ", k, " vs epoch: ", p 
                    #lfp_period_indexes = np.where(np.logical_and(lfp_spikes>=p*60, lfp_spikes<=(p+30)*60))[0]     #Find indexes of lfp events in sliding 30 minute period 
                    if len(lfp_period_indexes[p])==0:
                        kstest_array[k].append(1.0)
                        continue
                    spk2 = np.hstack(locked_spike_array[unit][lfp_period_indexes[p]])*1E-3

                    indexes = np.where(np.logical_and(spk2>=int(self.parent.lock_window_start.text()), spk2<=int(self.parent.lock_window_end.text())))[0]   #Work in milliseconds            #********************* ONLY USING SPIKES WITHIN 100ms of t=0 for KS TEST

                    if len(indexes)==0:
                        kstest_array[k].append(1.0) #DESYNCH STATES; APPEND DUMMY VALUE
                        continue 

                    spk2 = np.sort(spk2[indexes])
                    print "...in..."
                    KS, p_val = stats.ks_2samp(spk1, spk2)
                    print "...out...\n"
                    kstest_array[k].append(p_val)
       
            kstest_array[k]=np.hstack(kstest_array[k])
       
        np.save(cell_pval_filename, kstest_array)

    else:
        kstest_array = np.load(cell_pval_filename)

    print kstest_array
    print kstest_array[0]
    
    kstest_array += 1E-20
    #Compute log of P-value array
    kstest_array = np.vstack(kstest_array)[::-1]        #Invert array 
    kstest_array = np.log10(kstest_array)
    vmin_value = np.min(kstest_array)
    if vmin_value < -6: vmin_value = -6
    vmax_value = max(0,np.max(kstest_array))

    f2 = plt.figure()

    cax = plt.imshow(kstest_array, vmin=vmin_value, vmax=vmax_value, cmap=self.cmap.text(), interpolation='none')


    #cax = ax.imshow(ks_img, vmax=float(self.vmin_value.text()), vmin=float(self.vmax_value.text()), cmap=cm.jet_r, interpolation='none')
    #plt.ylim(len(ks_img),0)

    #ax = plt.subplot(2,3,4)
    
    ##Plot SUA Events
    #spikes = Sort_sua.units[unit]*1E-6/60.
    #ymin=np.zeros(len(spikes))
    #ymax=np.zeros(len(spikes))
    #ymin+=0
    #ymax+=10
    #plt.vlines(spikes, ymin, ymax, linewidth=1, color='black',alpha=.5) #colors[mod(counter,7)])
    
    plt.tick_params(axis='both', which='both', labelsize=30)
    plt.xlabel("Epochs", fontsize=30)
    plt.ylabel("Epochs", fontsize=30)
    
    tstep= 1
    old_xlabel = np.arange(0, len(kstest_array), tstep)
    new_xlabel = np.arange(1,len(kstest_array)+1,tstep)
    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    
    
    old_ylabel = np.arange(len(kstest_array)-1, -1, -tstep)
    new_ylabel = np.arange(1,len(kstest_array)+1,tstep)
    plt.yticks(old_ylabel, new_ylabel, fontsize=30) #, rotation='vertical')    
    
    #vmin_value = np.min(kstest_array)
    #vmax_value = 1.0
    print vmin_value, vmax_value
    #x_ticks = np.arange(float(self.vmin_value.text()), float(self.vmax_value.text())+1, 1)
    x_ticks = np.arange(vmax_value, vmin_value, -1)
    x_ticks = np.append(x_ticks, vmin_value)
    print x_ticks
    #x_ticks = [float(vmax_value), float(vmin_value)]
    cbar = f2.colorbar(cax, ticks=x_ticks)
    
    x_tick_labels = []
    for k in range(0,len(x_ticks)-1,1):
        x_tick_labels.append('$10^{'+str(int(x_ticks[k]))+"}$")
    
    x_tick_labels.append('$10^{'+str(round(x_ticks[-1],1))+"}$")
    

    print x_tick_labels
    cbar.ax.set_yticklabels(x_tick_labels)  # vertically oriented colorbar
    #cbar.ax.set_xticklabels([str(vmax_value), str(vmin_value)])  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=30) 
    
    #plt.suptitle('$log_{10}(P-Value)$', fontsize=30)
    plt.suptitle('P-Value', fontsize=30, fontweight='bold')

    
    plt.show()
   

def compute_msl_state_space(self):
    
    colors = ['blue', 'red', 'green', 'magenta', 'brown', 'orange', 'cyan', 'pink', 'grey', 'indigo']

    lfp_cluster = int(self.parent.lfp_cluster.text())
        

    #Load lock times
    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step"
    lock_time = np.load(file_out+'.npy')

    jitter_time = 50 #Time to jitter spiketrian
    file_out_poisson = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step"
    lock_time_poisson = np.load(file_out_poisson+'.npy')

    file_out_locked_spikes = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_locked_spikes"
    locked_spike_array = np.load(file_out_locked_spikes+'.npy') #This is al ist of all spikes locking to each 30min window; used for non-Luczak type studies;


    file_spikerates = file_out+"_spikerates"

    spikerates = np.load(file_spikerates+".npy")

    #Load lfp event rasters
    compress_factor = 50
    Sort_lfp = Ptcs(self.parent.lfp_event_file) #Auto load flag for Nick's data
    lfp_cluster = int(self.parent.lfp_cluster.text())
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor

    #Load single unit rasters
    Sort_sua = PTCS.PTCS(self.parent.sua_file) #Auto load flag for Nick's data
    units = np.arange(0,Sort_sua.n_units,1)
    #units = np.int16(self.multiple_units.text().split(","))

    ##Compute periods of synchrony from si index
    lfp = TSF.TSF(self.parent.sua_file.replace('_hp_butter_alltrack.ptcs','_lfp_250hz_alltrack.tsf'))
    lfp.read_ec_traces()    #Ok to read all LFP data as they are smaller files


    ax = plt.subplot(1, 1, 1)

  
    mean_locks=[]
    for k in range(200): mean_locks.append([])
    for unit in units:
        print "...unit: ", unit, "  overall fire rate: ", len(Sort_sua.units[unit])/float(lfp.n_vd_samples/lfp.SampleFrequency), " color: ", int(float(unit)/Sort_sua.n_units*256)

        even_locks=[]; odd_locks=[]; temp_locks = []
        for k in range(len(lock_time[unit])): 
            if spikerates[k][unit] > float(self.min_fire_rate.text()):       #*** Check if fire rate during epoch was sufficiently high
                
                #if ((lock_time[unit][k][0]-lock_time[unit][k][1])<10.0) and (lock_time[unit][k][0]!=0):
                if (lock_time[unit][k][0]==0) and (lock_time[unit][k][1]==0):    #If either odd or even lock times are zero
                    temp_locks.append(0.0)
                else:
                    if (lock_time[unit][k][0]==0):    #If either odd or even lock times are zero
                        temp_locks.append(lock_time[unit][k][1])
                    elif (lock_time[unit][k][1]==0):
                        temp_locks.append(lock_time[unit][k][0])
                    else:
                        temp_locks.append((lock_time[unit][k][0]+lock_time[unit][k][1])/2.)
            else:
                temp_locks.append(0.0)
        
        mean_locks[unit] = temp_locks
            
    #mean_locks = np.float32(mean_locks)
    
    lowest_index = 1000
    min_index = 100
    cell_msl_array = []
    for unit in units:
        indexes = np.where(np.array(mean_locks[unit])==0.0)[0]
        if len(indexes)>0:
            print "... unit: ", unit, " is stable to time point: ", indexes[0]
            if indexes[0]>=min_index:
                cell_msl_array.append(mean_locks[unit][:min_index])
        else:
            print "... unit: ", unit, " is stable to end of recording..."
            cell_msl_array.append(mean_locks[unit][:min_index])
        
       
    print "... # of cells with at least: ", min_index, " epochs: ", len(cell_msl_array)

    methods = ['MDS', 'tSNE', 'PCA', 'BH_tSNE']
    method = 2
    filename = file_out

    data_in = np.array(cell_msl_array).T
    print data_in.shape
    recompute = True
    data_out = dim_reduction_general(data_in, method, filename, recompute).T
    print "...data_out: ", data_out.shape
    
    #Generate control data:
    #shuffle_indexes = np.random.choice(150,150)
    #print shuffle_indexes
    #data_in_control = data_in[shuffle_indexes]
    data_in_control = np.random.rand(min_index,len(cell_msl_array))*50+75
    print data_in_control.shape
    filename = file_out+'_control'
    data_out_control = dim_reduction_general(data_in_control, method, filename).T
    

    
    #Plot data: 
    fig = plt.figure(1, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, 1, 1], elev=48, azim=134)

    colors = []
    for k in range(min_index):
        colors.append(int(float(k)/min_index*256))

    cmap = plt.cm.get_cmap('viridis', min_index)
    plt.set_cmap('viridis')
    for k in range(len(data_out[0])-1):
        ax.scatter(data_out[0][k], data_out[1][k], data_out[2][k], s=200, c=cm.viridis(int(float(k)/min_index*256)), cmap=cmap)
        ax.plot([data_out[0][k],data_out[0][k+1]], [data_out[1][k],data_out[1][k+1]], [data_out[2][k],data_out[2][k+1]], color='black')

    ax.scatter(data_out[0][-1], data_out[1][-1], data_out[2][-1], s=200, c=cm.viridis(int(float(k+1)/min_index*256)), cmap=cmap)
        
    
    fig = plt.figure(2, figsize=(4, 3))
    plt.clf()
    ax = Axes3D(fig, rect=[0, 0, 1, 1], elev=48, azim=134)

    #cmap = plt.cm.get_cmap('magma', Sort_sua.n_units)
    #plt.set_cmap('magma')
    for k in range(len(data_out[0])):
        ax.scatter(data_out_control[0][k], data_out_control[1][k], data_out_control[2][k], s=200, c=cm.viridis(int(float(k)/min_index*256)), cmap=cmap)
        #ax.scatter(data_out_control[0], data_out_control[1], data_out_control[2], s=200, c=colors, cmap=cmap)


    #ax.w_xaxis.set_ticklabels([])
    #ax.w_yaxis.set_ticklabels([])
    #ax.w_zaxis.set_ticklabels([])
    
    #ax.set_xlim(np.min(X[:, 0])-1, np.max(X[:, 0])+1)
    #ax.set_ylim(np.min(X[:, 1])-1, np.max(X[:, 1])+1)
    #ax.set_zlim(np.min(X[:, 2])-1, np.max(X[:, 2])+1)

    plt.show()

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)
    
    
    
def compute_msl_single_lfpevent(self):
    ''' Forgot what this does... :(
        I think it plots MSL for average (normal measure) vs. 1st spike vs ?
    '''
    
    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1000000  # window width in usec
    
    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())

    lfp_cluster = int(self.parent.lfp_cluster.text())

    tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))


    sig = float(self.sigma_width.text())
    lock_window_start = int(self.parent.lock_window_start.text())
    lock_window_end = int(self.parent.lock_window_end.text())

    #Load SUA Sort
    Sort_sua = PTCS.PTCS(self.parent.sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    #Load LFP Sort
    Sort_lfp = PTCS.PTCS(self.parent.lfp_event_file) #Auto load flag for Nick's data
    
    #Load pop events during synch periods (secs)
    compress_factor = 50
    print "... compress_factor is hardcoded to: ", compress_factor

    starting_cell = int(self.starting_cell.text()); ending_cell = int(self.ending_cell.text())
      
    #Load LFP Cluster events                                     #*******************ENSURE THAT NO DUPLICATE POP SPIKES MAKE IT THROUGH
    pop_spikes = np.uint64(Sort_lfp.units[lfp_cluster])*compress_factor
    pop_spikes=np.sort(np.unique(pop_spikes))*1E-3   #Exclude duplicates; convert to milisecond time
    #pop_spikes=(pop_spikes)*1E-3   #Exclude duplicates; convert to milisecond time
    
    print " ... # LFP events: ", len(pop_spikes)


    
    #REDO THE DISTANCE DIFF
    min_diff = 0.050    #50ms
    pop_spikes=pop_spikes*1E-3
    t = [0,31000]

    if False:
        for unit in range(20):
            ax = plt.subplot(4,5,unit+1)

            spikes = Sort_sua.units[unit]*1E-6

            recomputed_array1 = []
            diff_array1 = []
            spikes_perlfp_array = []
            for k in range(len(pop_spikes)): spikes_perlfp_array.append([])
            for k in range(len(spikes)):
                nearest_lfp_index = find_nearest(pop_spikes, spikes[k])
                temp_diff = spikes[k]- pop_spikes[nearest_lfp_index]
                if np.abs(temp_diff)<min_diff:
                    recomputed_array1.append(spikes[k])
                    diff_array1.append(temp_diff)
                    spikes_perlfp_array[nearest_lfp_index].append(temp_diff)
            
            plt.scatter(recomputed_array1, np.array(diff_array1)*1E3, color='blue')

            #Compute average lfp lock time for each LFP event:
            ave_lfp_lock = []
            x = []
            for k in range(len(spikes_perlfp_array)):
                if len(spikes_perlfp_array[k])>0:
                    ave_lfp_lock.append(np.mean(spikes_perlfp_array[k]))
                    x.append(pop_spikes[k])
            
            plt.scatter(x, np.array(ave_lfp_lock)*1E3-50, color='magenta')
                
            ax.ticklabel_format(useOffset=False)
            plt.ylim(-100,50)
            plt.xlim(0,t[-1])
            plt.title("Unit: " + str(unit), fontsize=20)
                
            ymin=np.zeros(len(spikes)); ymax=ymin+2.0
            plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=.8) #colors[mod(counter,7)])

            ymin=np.zeros(len(pop_spikes))+2; ymax=ymin+2.0
            plt.vlines(pop_spikes, ymin, ymax, linewidth=3, color='green',alpha=.8) #colors[mod(counter,7)])

            ymin=np.zeros(len(recomputed_array1))+4; ymax=ymin+2.0
            plt.vlines(recomputed_array1, ymin, ymax, linewidth=3, color='red',alpha=.8) #colors[mod(counter,7)])
            
            ax.tick_params(axis='both', which='both', labelsize=8)
        
    def runningMeanFast(x, N):
        return np.convolve(x, np.ones((N,))/N)[(N-1):]
    
    
    unit = starting_cell
    ax = plt.subplot(1,1,1)

    spikes = Sort_sua.units[unit]*1E-6

    #Mean ? 
    recomputed_array1 = []
    diff_array1 = []
    spikes_perlfp_array = []
    for k in range(len(pop_spikes)): spikes_perlfp_array.append([])
    for k in range(len(spikes)):
        nearest_lfp_index = find_nearest(pop_spikes, spikes[k])
        temp_diff = spikes[k]- pop_spikes[nearest_lfp_index]
        if np.abs(temp_diff)<min_diff:
            recomputed_array1.append(spikes[k])
            diff_array1.append(temp_diff)
            spikes_perlfp_array[nearest_lfp_index].append(temp_diff)
    
    plt.scatter(recomputed_array1, np.array(diff_array1)*1E3, color='blue', alpha=0.6)

    running_mean = runningMeanFast(np.array(diff_array1)*1E3, 10)
    plt.plot(recomputed_array1, running_mean, color='blue', linewidth=5, alpha=0.5)
    

    #Compute average lfp lock time for each LFP event:
    ave_lfp_lock = []
    x = []
    for k in range(len(spikes_perlfp_array)):
        if len(spikes_perlfp_array[k])>0:
            ave_lfp_lock.append(np.mean(spikes_perlfp_array[k]))
            x.append(pop_spikes[k])
    
    plt.scatter(x, np.array(ave_lfp_lock)*1E3-50, color='magenta', alpha=0.6)
    running_mean = runningMeanFast(np.array(ave_lfp_lock)*1E3-50, 10)
    plt.plot(x, running_mean, color='magenta', linewidth=5, alpha=0.5)


    first_spike = []
    x = []
    for k in range(len(spikes_perlfp_array)):
        if len(spikes_perlfp_array[k])>0:
            first_spike.append(spikes_perlfp_array[k][0])
            x.append(pop_spikes[k])
    
    plt.scatter(x, np.array(first_spike)*1E3-100, color='brown', alpha=0.6)
    running_mean = runningMeanFast(np.array(first_spike)*1E3-100, 10)
    plt.plot(x, running_mean, color='brown', linewidth=5, alpha=0.5)
    
    ax.ticklabel_format(useOffset=False)
    plt.ylim(-150,50)
    plt.xlim(0,t[-1])
    plt.title("Unit: " + str(unit), fontsize=20)
        
    #Plot rasters
    ymin=np.zeros(len(spikes))+20; ymax=ymin+10.0
    plt.vlines(spikes, ymin, ymax, linewidth=2, color='black',alpha=.8) #colors[mod(counter,7)])

    ymin=np.zeros(len(pop_spikes))+30; ymax=ymin+10.0
    plt.vlines(pop_spikes, ymin, ymax, linewidth=2, color='green',alpha=.8) #colors[mod(counter,7)])

    ymin=np.zeros(len(recomputed_array1))+40; ymax=ymin+10.0
    plt.vlines(recomputed_array1, ymin, ymax, linewidth=2, color='red',alpha=.8) #colors[mod(counter,7)])
    
    ax.tick_params(axis='both', which='both', labelsize=15)


    old_ylabel = [-150,-125,-100,-75, -50,-25,  0, 50]
    new_ylabel = [-50,  -25,   0,-25,   0,-25,  0, 50]
    plt.yticks(old_ylabel, new_ylabel, fontsize=15)
    plt.xlabel("Time (sec)", fontsize=25)
    plt.plot([0,t[-1]], [0,0], 'r--', color='black', alpha=0.7)
    plt.plot([0,t[-1]], [-50,-50], 'r--', color='black', alpha=0.7)
    plt.plot([0,t[-1]], [-100,-100], 'r--', color='black', alpha=0.7)
            
    plt.show()
    
    
def Compute_MSL_depth(self):
    
    '''Align msl latencies and also show by depth 
    '''
    
    #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'

    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file


    top_channel = np.loadtxt(os.path.split(os.path.split(self.parent.sua_file)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.

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
    n_units_incortex = len(np.where(Sort_sua.maxchan>=top_channel)[0])                                          #Number of units that are in tissue, i.e. below the top_channel.txt (see file)


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
    
    #Load cell rasters from file
    
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")

    #********************************************* START OVER LOADING DATA FROM DISK *****************************************

    chunk_index = np.arange(0,len(time_chunks),1).tolist()

    #chunk_index = []
    #for chk in self.chunks_to_plot.text().split(','):
        #chunk_index.append(int(chk))
    
    print "...chunk_index: ", chunk_index
    
    f1 = plt.figure()
    img_means = []
    sig = float(self.sigma_width.text())
    #for ctr, chunk_ctr in enumerate(chunk_index):                      #LOOP OVER ALL CHUNKS
    for ctr in range(1):                                    #COMPUTE ONLY FIRST CHUNK
        chunk_ctr = chunk_index[0]
        print len(time_chunks), chunk_ctr
        time_chunk = time_chunks[chunk_ctr]
        
        print time_chunk
        fit_sum = np.zeros((total_units, 2000), dtype=np.float32)

        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        
        for unit in range(total_units):
        #for unit in range(1):
            if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it
            if len(Sort_sua.units[unit])<min_spikes: continue
            locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])
            
            n_spikes_pre = len(locked_spikes)
            locked_spikes = np.unique(np.sort(locked_spikes))
               
            print "... chunk: ", time_chunk, " ...unit: ", unit, " #spikes locked: ", len(locked_spikes), " / ", len(Sort_sua.units[unit]), \
            "   duplicates: ", n_spikes_pre - len(locked_spikes)

            #NEW METHOD: JUST COMPUTE GAUSSIAN ONCE THEN TRANSLATE THE ARRAY
            if len(locked_spikes)>0: 
                x = np.linspace(-1000, 1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(locked_spikes)):
                    mu = int(locked_spikes[g]*1E-3)
                    fit_sum[unit] += np.roll(sig_gaussian, mu)
        
        img=[]
        lock_time=[]
        #Add each cell histogram normalized to itself.
        for unit in range(total_units):
            if np.max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end])>0:
                lock_time.append(np.argmax(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
                img.append(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]/max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
            else:
                lock_time.append(lock_window_end*4)
                img.append(np.zeros(lock_window_end-lock_window_start, dtype=np.float32))
                

        #**************** PLOT MSL - BY ORDER **********************
        #ax = plt.subplot(len(chunk_index), 2, ctr*2+1)
        ax = plt.subplot(1,2,1)

        #ORDER MSL IMG BY LOCK TIME OF FIRST EPOCH
        if (ctr ==0): inds = np.array(lock_time).argsort()
        print "Order: ", inds

        img_ordered=np.array(img)[inds]
        temp_img = []
        for i in range(len(img_ordered)):
            #if np.max(img[i])!=0:
                temp_img.append(img_ordered[i])
        img_ordered=np.array(temp_img)
        
        im = ax.imshow(img_ordered, origin='upper', extent=[0,(lock_window_end-lock_window_start), len(img),0], aspect='auto', interpolation='none')

        #SET LABELS
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if ctr ==0:
            xx = np.arange(0,(lock_window_end-lock_window_start)+1,(lock_window_end-lock_window_start)/2.)
            x_label = np.arange(lock_window_start,lock_window_end+1,(lock_window_end-lock_window_start)/2.)
            plt.xticks(xx,x_label, fontsize=25)
            plt.xlabel("Time (ms)", fontsize=30, fontweight='bold')

            yy = np.arange(0,len(img),20)
            y_label = np.arange(0,len(img),20)
            plt.yticks(yy, y_label, fontsize=30)

            plt.ylabel("Cell #", fontsize=30, fontweight='bold')

        #plt.title(str(int(time_chunk[0]/(60.*1E6)))+".."+str(int(time_chunk[1]/(60.*1E6)))+" mins", fontsize=15)
        plt.title("Latency Order", fontsize=25)
        plt.ylim(len(inds),0)
    
        plt.plot([-lock_window_start,-lock_window_start],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.xaxis.labelpad = 0

        img_means.append(np.mean(img, axis=0))
    
        #**************** PLOT MSL - BY DEPTH **********************
        #ax = plt.subplot(len(chunk_index), 2, ctr*2+2)
        ax = plt.subplot(1,2,2)

        im = ax.imshow(img, origin='upper', extent=[0,(lock_window_end-lock_window_start), len(img),0], aspect='auto', interpolation='none')

        #SET LABELS
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if ctr ==0:
            xx = np.arange(0,(lock_window_end-lock_window_start)+1,(lock_window_end-lock_window_start)/2.)
            x_label = np.arange(lock_window_start,lock_window_end+1,(lock_window_end-lock_window_start)/2.)
            plt.xticks(xx,x_label, fontsize=20)
            plt.xlabel("Time (ms)", fontsize=30, fontweight='bold')

            yy = np.arange(0,len(img),20)
            y_label = np.arange(0,len(img),20)
            plt.yticks(yy, y_label, fontsize=30)

            #plt.ylabel("Cell #", fontsize=30, fontweight='bold')

        #plt.title(str(int(time_chunk[0]/(60.*1E6)))+".."+str(int(time_chunk[1]/(60.*1E6)))+" mins", fontsize=15)
        plt.title("Depth Order", fontsize=25)
        plt.yticks([])
        plt.ylim(len(inds),0)
        
        plt.plot([-lock_window_start,-lock_window_start],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.xaxis.labelpad = 0

        img_means.append(np.mean(img, axis=0))

    plt.suptitle("Group: "+self.parent.lfp_cluster.text(), fontsize=25, fontweight='bold')
    #plt.suptitle(self.parent.sua_file.replace('.ptcs','')+",  sigma: " + str(sig) +"(ms)", fontsize=20)


    ##*******************PLOT MSL ALL CELL DISTRIBUTIONS ************************
    ##cmap = plt.cm.get_cmap('viridis', Sort_sua.n_units)
    ##plt.set_cmap('viridis')
    
    #f2 = plt.figure()
    #for ctr,img_mean in enumerate(img_means):
        #print img_mean
        #plt.plot(img_mean, linewidth = 5, color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
    
    #old_xlabel = np.linspace(0, window*2, 3)
    #new_xlabel = np.linspace(-window, window, 3)
    
    
    #plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    
    #plt.xlabel("Time (ms)", fontsize = 30, fontweight='bold')
    #plt.xlim(0,lock_window_end*2)
    #plt.ylim(bottom=0)
     
    
    ##*******************PLOT PERCENT SPIKING HISTOGRAMS ************************
    ##cmap = plt.cm.get_cmap('viridis', Sort_sua.n_units)
    ##plt.set_cmap('viridis')
    
    #f3 = plt.figure()
    #img_means = []
    ##for ctr, chunk_ctr in enumerate(chunk_index):
    ##    time_chunk = time_chunks[chunk_ctr]
    ##    temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]

    #x = np.arange(0,n_units_incortex,1)

    #cumulative_bars = np.zeros(n_units_incortex, dtype=np.float32)
    #for k in range(10):
        #cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(k)+".npy"
        #if os.path.exists(cell_rasters_filename): 
            
            #cell_rasters = np.load(cell_rasters_filename)

            #percent_array = []
            #for unit in range(total_units):
                #if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it

                #locked_spikes = np.hstack(np.array(cell_rasters[unit]))*1E-3        #Convert from usec to msec
                
                #indexes = np.where(np.logical_and(locked_spikes>=-50, locked_spikes<=50))[0]   #Look for spikes between -50 to +50msec around LFP time

                #percent_array.append(float(len(indexes))/len(Sort_sua.units[unit])*1E2) #Convert to %
                
            #print percent_array

            ##plt.bar(x, percent_array, 1, color='blue')
            #print len(x), len(percent_array)
            #p2 = plt.bar(x, percent_array, 0.95, bottom=cumulative_bars, color = colors[k])
            
            #cumulative_bars=cumulative_bars + np.float32(percent_array)

    #plt.ylim(0,100)

    plt.show()
        
    
        
def Compute_MSL_chunks(self):
    '''Align msl latencies by lfp cluster 
    '''
    
    #self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'

    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file


    top_channel = np.loadtxt(os.path.split(os.path.split(self.parent.sua_file)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.

    #colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

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
    n_units_incortex = len(np.where(Sort_sua.maxchan>=top_channel)[0])                                          #Number of units that are in tissue, i.e. below the top_channel.txt (see file)


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
    #OPTION 1: Divide into equal chunks of recording length
    
    if False: 
        self.parent.tsf = TSF.TSF(self.parent.sua_file.replace('.ptcs','.tsf'))
        temp_chunks = np.linspace(0,self.parent.tsf.n_vd_samples*1E6/float(self.parent.tsf.SampleFrequency), int(self.parent.time_chunks.text()))
    
    else: 
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
    
    #Load cell rasters from file
    
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")

    #********************************************* START OVER LOADING DATA FROM DISK *****************************************

    chunk_index = np.arange(0,len(time_chunks),1).tolist()

    #chunk_index = []
    #for chk in self.chunks_to_plot.text().split(','):
    #    chunk_index.append(int(chk))
    
    
    print "...chunk_index: ", chunk_index
    
    f1 = plt.figure()
    img_means = []
    img_sums = []
    sig = float(self.sigma_width.text())
    unstable_cells = []
    img_array = []
    for ctr, chunk_ctr in enumerate(chunk_index):
        print len(time_chunks), chunk_ctr
        time_chunk = time_chunks[chunk_ctr]
        
        print time_chunk
        fit_sum = np.zeros((total_units,2000), dtype=np.float32)

        temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]
        
        for unit in range(total_units):
        #for unit in range(1):
            if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it
            if len(Sort_sua.units[unit])<min_spikes: continue
            locked_spikes = np.hstack(np.array(cell_rasters[unit])[temp3])
            
            n_spikes_pre = len(locked_spikes)
            locked_spikes = np.unique(np.sort(locked_spikes))
               
            print "... chunk: ", time_chunk, " ...unit: ", unit, " #spikes locked: ", len(locked_spikes), " / ", len(Sort_sua.units[unit]), \
            "   duplicates: ", n_spikes_pre - len(locked_spikes)

            #NEW METHOD: JUST COMPUTE GAUSSIAN ONCE THEN TRANSLATE THE ARRAY
            if len(locked_spikes)>0 : 
                x = np.linspace(-1000,1000, 2000)    #Make an array from -1000ms .. +1000ms with microsecond precision
                sig_gaussian = np.float32(gaussian(x, 0, sig))
                for g in range(len(locked_spikes)):
                    mu = int(locked_spikes[g]*1E-3)
                    fit_sum[unit] += np.roll(sig_gaussian, mu)
                                
        img=[]
        lock_time=[]
        #Add each cell histogram normalized to itself.
        for unit in range(total_units):
            #if np.max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end])>0:
            if (np.max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end])==np.max(fit_sum[unit])) and (np.max(fit_sum[unit])!=0):
                lock_time.append(np.argmax(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
                img.append(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]/max(fit_sum[unit][1000+lock_window_start:1000+lock_window_end]))
            else:
                #CELL DOESN"T HAVE A LOCK DURING THIS PERIOD
                unstable_cells.append(unit)
                lock_time.append(lock_window_end*4)
                img.append(np.zeros(lock_window_end-lock_window_start, dtype=np.float32))
                
        #ORDER MSL IMG BY LOCK TIME OF FIRST EPOCH
        if (ctr ==0): 
            inds = np.array(lock_time).argsort()
        print "Order: ", inds

        
        img=np.array(img)[inds]
        temp_img = []
        for i in range(len(img)):
            #if np.max(img[i])!=0:
                temp_img.append(img[i])
        img=np.array(temp_img)
        img_array.append(img)
    
        if ctr==2: break        #Exit after 3 epochs
    
    unstable_list = np.unique(unstable_cells)
    for k in range(len(unstable_list)):
        unstable_list[k] = np.argwhere(inds== unstable_list[k])
    print "..unstable_list: ", unstable_list
    #inds = np.array(lock_time).argsort()
    #unstable_cells = np.array(unstable_cells)
    
    stable_lindexes = np.arange(len(inds))
    if True:
        stable_list = []
        for k in range(total_units):
            if k not in unstable_list:
                stable_list.append(k)
        print stable_list
        stable_lindexes=np.int16(stable_list)
    else:
        print "... not excluding any cells..."
    
    #******************** PLOT IMAGES ONLY FOR STABLE CELLS **************************    
    for ctr, chunk_ctr in enumerate(chunk_index):

        #********** PLOTING ROUTINES **************
        #ax=plt.subplot(1,int(self.chunks_to_plot.text()),chunk_ctr+1)
        #ax = plt.subplot(1, len(chunk_index), ctr+1)
        ax = plt.subplot(1, 3, ctr+1)
        
        img = img_array[ctr][stable_lindexes]
        im = ax.imshow(img, origin='upper', extent=[0,(lock_window_end-lock_window_start), len(img),0], aspect='auto', interpolation='none')

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if ctr ==0:
            xx = np.arange(0,(lock_window_end-lock_window_start)+1,(lock_window_end-lock_window_start)/2.)
            x_label = np.arange(lock_window_start,lock_window_end+1,(lock_window_end-lock_window_start)/2.)
            plt.xticks(xx,x_label, fontsize=25)
            plt.xlabel("Time from LFP event (ms)", fontsize=30, fontweight='bold')

            yy = np.arange(0,len(img),20)
            y_label = np.arange(0,len(img),20)
            plt.yticks(yy, y_label, fontsize=30)

            plt.ylabel("Cell #", fontsize=30, fontweight='bold')
        
        plt.title(str(int(time_chunks[ctr][0]/(60.*1E6)))+" mins", fontsize=15)

        plt.ylim(len(img),0)

        
        #plt.ylabel("Neuron (lock order)", fontsize=30)
        plt.plot([-lock_window_start,-lock_window_start],[0,total_units], 'r--', linewidth=4, color='white')
        ax.tick_params(axis='both', which='both', labelsize=25)
        ax.xaxis.labelpad = 0

        img_means.append(np.mean(img, axis=0))
        img_sums.append(np.sum(img, axis=0))

        if ctr==2: break        #Exit after 3 epochs

    #plt.suptitle("Group: "+self.parent.lfp_cluster.text(), fontsize=25, fontweight='bold')
    plt.suptitle(self.parent.sua_file.replace('.ptcs','')+",  sigma: " + str(sig) +"(ms)", fontsize=20)

    plt.show()
    
    return

    #*******************PLOT MSL ALL CELL DISTRIBUTIONS ************************
    #cmap = plt.cm.get_cmap('viridis', Sort_sua.n_units)
    #plt.set_cmap('viridis')
    
    #colors = ['green', 'blue', 'red']
    f2 = plt.figure()
    for ctr in range(len(img_means)):
        #img_sum=img_sums[ctr]/np.max(np.array(img_sums))
        img_mean=img_means[ctr]/np.max(np.array(img_means))
        
        plt.plot(img_mean, linewidth = 15, color = colors[int(self.parent.lfp_cluster.text())])#color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
        #plt.plot(img_mean, linewidth = 15, color = colors[int(self.parent.lfp_cluster.text())-1])#color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
        #plt.plot(img_sum, linewidth = 15, color = 'red')#color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
    
    old_xlabel = np.linspace(0, window*2, 3)
    new_xlabel = np.linspace(-window, window, 3)
    
    
    plt.plot([len(img_mean)/2,len(img_mean)/2], [0,1.1], 'r--', color='black', linewidth=5, alpha=0.8)
    
    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    
    plt.xlabel("Time (ms)", fontsize = 30, fontweight='bold')
    plt.xlim(0,lock_window_end*2)
    plt.ylim(0, 1.1)
    
    f2_1 = plt.figure()
    ctr = int(self.starting_cell.text())
    #img_sum=img_sums[ctr]/np.max(np.array(img_sums))
    img_mean=img[ctr]/np.max(np.array(img[ctr]))
    
    plt.plot(img_mean, linewidth = 15, color = colors[int(self.parent.lfp_cluster.text())])#color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
    #plt.plot(img_sum, linewidth = 15, color = 'red')#color = cm.viridis(int(float(ctr)/len(chunk_index)*256)))
    
    old_xlabel = np.linspace(0, window*2, 3)
    new_xlabel = np.linspace(-window, window, 3)


    plt.plot([len(img_mean)/2,len(img_mean)/2], [0,1.1], 'r--', color='black', linewidth=5, alpha=0.8)

    plt.xticks(old_xlabel, new_xlabel, fontsize=30) #, rotation='vertical')    
    plt.xlabel("Time (ms)", fontsize = 30, fontweight='bold')
    plt.xlim(0,lock_window_end*2)
    plt.ylim(0, 1.1)



    #*******************PLOT PERCENT SPIKING HISTOGRAMS ************************
    #cmap = plt.cm.get_cmap('viridis', Sort_sua.n_units)
    #plt.set_cmap('viridis')
    
    f3 = plt.figure()
    img_means = []
    #for ctr, chunk_ctr in enumerate(chunk_index):
    #    time_chunk = time_chunks[chunk_ctr]
    #    temp3 = np.where(np.logical_and(pop_spikes>=time_chunk[0], pop_spikes<=time_chunk[1]))[0]

    x = np.arange(0,n_units_incortex,1)

    cumulative_bars = np.zeros(n_units_incortex, dtype=np.float32)
    for k in range(10):
        cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(k)+".npy"
        if os.path.exists(cell_rasters_filename): 
            
            cell_rasters = np.load(cell_rasters_filename)

            percent_array = []
            for unit in range(total_units):
                if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it

                locked_spikes = np.hstack(np.array(cell_rasters[unit]))*1E-3        #Convert from usec to msec
                
                indexes = np.where(np.logical_and(locked_spikes>=-50, locked_spikes<=50))[0]   #Look for spikes between -50 to +50msec around LFP time

                percent_array.append(float(len(indexes))/len(Sort_sua.units[unit])*1E2) #Convert to %
                
            #print percent_array

            #plt.bar(x, percent_array, 1, color='blue')
            #print len(x), len(percent_array)
            p2 = plt.bar(x, percent_array, 0.95, bottom=cumulative_bars, color = colors[k])
            
            cumulative_bars=cumulative_bars + np.float32(percent_array)

    plt.ylim(0,100)

    plt.show()
    

def compute_msl_pvals(self):
    #COMPUTE KS P-VALUES FOR CURRENT LFP CLUSTER
    
    colors = ['blue', 'red', 'green', 'magenta', 'brown', 'orange', 'cyan', 'pink', 'grey', 'indigo']

    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    #units = np.int16(self.multiple_units.text().split(","))
    #units = np.arange(int(self.starting_cell.text()), int(self.ending_cell.text()), 1)
    
    unit = int(self.starting_cell.text())
    
    #Load lock times
    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step"
    lock_time = np.load(file_out+'.npy')

    file_out_locked_spikes = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_locked_spikes"
    locked_spike_array = np.load(file_out_locked_spikes+'.npy') #This is al ist of all spikes locking to each XXmin window; used for non-Luczak type studies;

    print len(locked_spike_array)
    print len(locked_spike_array[0])
    
    units = np.arange(0,118,1)
    #unit = [56,99,100]
    empty_row = np.zeros(len(locked_spike_array),dtype=np.float32)+123.0
    for unit in units: 
        kstest_array = []
        for k in range(len(locked_spike_array)):
            print "... unit: ", unit, " epoch: ", k
            kstest_array.append([])    
            
            temp = locked_spike_array[k][unit][:10]
            if len(temp)==0:
                kstest_array[k] = empty_row             #DESYNCH STATES; APPEND DUMMY VALUE
                continue
    
            spk1 = np.hstack(locked_spike_array[k][unit])*1E-3
            indexes = np.where(np.logical_and(spk1>=-50, spk1<=50))[0]   #Work in milliseconds
            spk1 = spk1[indexes]
              
            for p in range(len(locked_spike_array)):
            #for p in range(200):
                temp = locked_spike_array[p][unit][:10]
                if len(temp)==0:
                    kstest_array[k].append(123.0)       #DESYNCH STATES; APPEND DUMMY VALUE     
                    continue
                
                spk2 = np.hstack(locked_spike_array[p][unit])*1E-3
                indexes = np.where(np.logical_and(spk2>=-50, spk2<=50))[0]   #Work in milliseconds
                spk2 = spk2[indexes]
                 
                KS, p_val = stats.ks_2samp(spk1, spk2)
                
                
                #print "... p_val: ", p_val
                kstest_array[k].append(p_val)
       
            kstest_array[k]=np.hstack(kstest_array[k])
       
        kstest_array = np.vstack(kstest_array)

        np.save(file_out_locked_spikes+"_unit_"+str(unit), kstest_array)
    
    plt.imshow(kstest_array, vmin=0.05, vmax=1.0, interpolation='none')
    plt.show()
    
    
    
def view_msl_pval(self):
    #VIEW KS P-VALUES FOR SINGLE CELL
    
    colors = ['blue', 'red', 'green', 'magenta', 'brown', 'orange', 'cyan', 'pink', 'grey', 'indigo']

    lfp_cluster = int(self.parent.lfp_cluster.text())
    
    #units = np.int16(self.multiple_units.text().split(","))
    #units = np.arange(int(self.starting_cell.text()), int(self.ending_cell.text()), 1)
    
    unit = int(self.starting_cell.text())
    
    Sort_sua = Ptcs(self.parent.sua_file) #Auto load flag for Nick's data
    Sort_lfp = Ptcs(self.parent.lfp_event_file)
    
    #ax = plt.subplot(1,1,1)
    fig, ax = plt.subplots()

    #Load lock times
    file_out = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step"
    lock_time = np.load(file_out+'.npy')

    file_out_locked_spikes = self.parent.sua_file.replace('.ptcs','')+"_"+str(lfp_cluster)+"lfpcluster_"+self.sliding_window_length.text()+"window_"+self.sliding_window_step.text()+"step_locked_spikes"

    ks_img = np.load(file_out_locked_spikes+"_unit_"+str(unit)+'.npy')
    ks_img = np.flipud(ks_img)
    
    #plt.imshow(ks_img, vmax=float(self.vmin_value.text()), vmin=float(self.vmax_value.text()))
    #plt.show()
    #return
    
    print ks_img.shape
    
    #Mask desynch periods
    if True: 
        mask_indexes = np.where(ks_img.ravel()==123.0)
        mask_array = np.zeros(len(ks_img.ravel()), dtype=np.int8)
        mask_array[mask_indexes]=1
        ks_img = np.ma.masked_array(ks_img, mask = mask_array)

        #Recombine to make img
        print ks_img.shape
        ks_img = np.ma.array(np.split(ks_img, 453))
        print ks_img.shape
        ks_img = np.ma.vstack(ks_img)
        print ks_img.shape
    
    
    #Remove desynch periods
    if False:
        new_array = []
        for k in range(len(ks_img)):
           if np.min(ks_img[k])<100:
                new_array.append(ks_img[k])
        new_array = np.vstack(new_array).T
        
        ks_img = []
        for k in range(len(new_array)):
            if np.min(new_array[k])<100:
                ks_img.append(new_array[k])
        ks_img=np.vstack(ks_img).T        
    
        
    #print ks_img
    #ks_img = 1-ks_img    

    #ks_img = np.ma.log10(np.clip(1-ks_img,0,1))
    ks_img = np.ma.log10(ks_img)
    
    #print ks_img
    
    cax = ax.imshow(ks_img, vmax=float(self.vmin_value.text()), vmin=float(self.vmax_value.text()), cmap=cm.jet_r, interpolation='none')
    #plt.ylim(len(ks_img),0)

    #ax = plt.subplot(2,3,4)
    
    ##Plot SUA Events
    #spikes = Sort_sua.units[unit]*1E-6/60.
    #ymin=np.zeros(len(spikes))
    #ymax=np.zeros(len(spikes))
    #ymin+=0
    #ymax+=10
    #plt.vlines(spikes, ymin, ymax, linewidth=1, color='black',alpha=.5) #colors[mod(counter,7)])
    
    plt.tick_params(axis='both', which='both', labelsize=30)
    plt.xlabel("Time (min)", fontsize=30)
    plt.ylabel("Time (min)", fontsize=30)
    
    tstep= 50
    #old_xlabel = np.arange(len(ks_img), 0, -tstep)
    new_xlabel = np.arange(0,len(ks_img),tstep)
    plt.xticks(new_xlabel, fontsize=30) #, rotation='vertical')    
    
    old_ylabel = np.arange(len(ks_img), 0, -tstep)
    new_ylabel = np.arange(0,len(ks_img),tstep)
    plt.yticks(old_ylabel, new_ylabel, fontsize=30) #, rotation='vertical')    
        
    x_ticks = np.arange(float(self.vmin_value.text()), float(self.vmax_value.text())+1, 1)
    cbar = fig.colorbar(cax, ticks=[x_ticks])
    #cbar.ax.set_yticklabels(['1^'+self.vmin_value.text(), '1^'+self.vmax_value.text()])  # vertically oriented colorbar
    cbar.ax.tick_params(labelsize=30) 

    plt.suptitle("2-Sample KS-Test: P-values "+ " Unit: " + str(unit), fontsize=30)
    plt.show()
    
    

    #quit()

def load_lfp_length(file_name):     #Nick/Martin data has different LFP structure to their data.
    data_in = np.load(file_name)
    t1 = data_in['t1']  #Save length of lfp recording

    return t1


def compute_natscene_pairisi(self):
    


    #Load time interval
    seq_isi_start = float(self.seqisi_start.text())
    seq_isi_end = float(self.seqisi_end.text())
    print seq_isi_start, seq_isi_end
    
    unit1 = int(self.starting_cell.text())
    unit2 = int(self.ending_cell.text())

    #Find offset for stimulus processing
    recordings = sorted(glob.glob(self.rec_path+"/*"))
    offset = 0
    for recording in recordings:
        print recording
        temp= recording.replace(self.rec_path+'/','')
        if temp in self.rec_name: 
            din_file = self.rec_path+'/'+temp+'/'+temp+'.din'       #Save .din filename before exiting
            break
        t1 = load_lfp_length(glob.glob(recording+"/*.lfp.zip")[0])
        offset += t1

    offset = offset*1E-6
    print "...offset: ", offset



    #******************************************************
    #*************** LOAD DIN STIMULUS INFO ***************
    #******************************************************
    print din_file
    f = open(din_file, 'rb')
    din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols containing [time_stamp in usec : frame number]
    f.close()
    print "No of frames in din file: ", len(din)
    print "Length of recording: ", float(din[-1][0])*1E-6/60., " mins." 

    start_indexes = np.where(din[:,1]==0)[0][::3]   #each frame is on screen for 3 refreshes (NOT ALWAYS TRUE*********)
    stimuli = din[:,0][start_indexes]*1E-6
      


    #Load LFP
    Sort_lfp = Ptcs(self.parent.lfp_event_file) 
    lfp_cluster = int(self.parent.lfp_cluster.text())
    lfp_events = Sort_lfp.units[lfp_cluster]*1E-6
    compress_factor = 50.   #Needed to uncompress the LFP compressed sorts

    print lfp_events[-10:]
    print lfp_events[-10:]*compress_factor
    print len(lfp_events)

    #Load saved cell rasters
    cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(lfp_cluster)
    cell_rasters = np.load(cell_rasters_filename+".npy")
    print len(cell_rasters)
    print len(cell_rasters[0])
        

    #Load LFP rasters falling within
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    lfp_spikes_array = []
    all_spikes=[]
    lfp_events = lfp_events*compress_factor
    lfp_indexes = np.where(np.logical_and(lfp_events>=stimuli[0]+offset, lfp_events<=stimuli[-1]+offset+1.0))[0]
    
    print stimuli[0]+offset, stimuli[-1]+offset
    print lfp_indexes
    print len(lfp_indexes)
    print lfp_events[lfp_indexes]
    
    print cell_rasters[unit1][lfp_indexes]
    print cell_rasters[unit2][lfp_indexes]
    
    return

    

    #Load spike rasters
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    spikes_array = []
    all_spikes=[]
    for unit in range(len(Sort_sua.units)):
        spikes = Sort_sua.units[unit]*1E-6
        spikes = spikes[np.where(np.logical_and(spikes>=stimuli[0]+offset, spikes<=stimuli[-1]+offset+1.0))[0]]     #Save all spikes that fall within experiment time
        spikes_array.append(spikes-offset)
        all_spikes.extend(spikes-offset)
           
           
    #*********** Plot unit 1 rasters rasters

    ax = plt.subplot(3,2,1)
    all_spikes1 = []
    print "...unit: ", unit1
    plt.ylim(0,400)
    plt.xlim(0,5.5)
    depth=0
    for stimulus in stimuli:
        ymin=np.zeros(len(spikes_array[unit1]))
        ymax=np.zeros(len(spikes_array[unit1]))
        ymin+=depth
        ymax+=depth+.9
        depth+=1
        
        spikes = spikes_array[unit1][np.where(np.logical_and(spikes_array[unit1]>=stimulus+seq_isi_start, spikes_array[unit1]<=stimulus+seq_isi_end))[0]]-stimulus
        all_spikes1.append(spikes)
        
        plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=.8) #colors[mod(counter,7)])

    #Plot distribution of spikes

    ax = plt.subplot(3,2,2)
    all_spikes2 = []
    print "...unit: ", unit2
    plt.ylim(0,400)
    plt.xlim(0,5.5)
    depth=0
    for stimulus in stimuli:
        ymin=np.zeros(len(spikes_array[unit2]))
        ymax=np.zeros(len(spikes_array[unit2]))
        ymin+=depth
        ymax+=depth+.9
        depth+=1
        
        spikes = spikes_array[unit2][np.where(np.logical_and(spikes_array[unit2]>=stimulus+seq_isi_start, spikes_array[unit2]<=stimulus+seq_isi_end))[0]]-stimulus
        all_spikes2.append(spikes)
        
        plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=.8) #colors[mod(counter,7)])


    #Plot averge ISI during each stimulus 
    ax = plt.subplot(3,2,3)
    sequential_isi = []
    for k in range(len(all_spikes1)):   #Loop over each stimulus presentation
        temp_array=[]
        for p in range(len(all_spikes2[k])):
            if (len(all_spikes1[k])>0) and (len(all_spikes2[k])>0):
                temp_array.extend(all_spikes1[k]-all_spikes2[k][p])
        print k, temp_array
        sequential_isi.append(np.mean(temp_array))

    print len(sequential_isi), len(np.hstack(sequential_isi))
    plt.plot(np.arange(len(sequential_isi)), sequential_isi, linewidth=3)
    plt.scatter(np.arange(len(sequential_isi)), sequential_isi, s=10)




    plt.show()



def compute_natscene_seqisi(self):
    
    #Load time interval
    seq_isi_start = float(self.seqisi_start.text())
    seq_isi_end = float(self.seqisi_end.text())
    print seq_isi_start, seq_isi_end
    
    
    #Find offset for stimulus processing
    recordings = sorted(glob.glob(self.rec_path+"/*"))
    offset = 0
    for recording in recordings:
        print recording
        temp= recording.replace(self.rec_path+'/','')
        if temp in self.rec_name: 
            din_file = self.rec_path+'/'+temp+'/'+temp+'.din'       #Save .din filename before exiting
            break
        t1 = load_lfp_length(glob.glob(recording+"/*.lfp.zip")[0])
        offset += t1

    offset = offset*1E-6
    print "...offset: ", offset



    #**** Load .din stimulus timestamps
    print din_file
    f = open(din_file, 'rb')
    din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols containing [time_stamp in usec : frame number]
    f.close()
    print "No of frames in din file: ", len(din)
    print "Length of recording: ", float(din[-1][0])*1E-6/60., " mins." 

    start_indexes = np.where(din[:,1]==0)[0][::3]   #each frame is on screen for 3 refreshes (NOT ALWAYS TRUE*********)
    stimuli = din[:,0][start_indexes]*1E-6
      
    

    #Load spike rasters
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    spikes_array = []
    all_spikes=[]
    for unit in range(len(Sort_sua.units)):
        spikes = Sort_sua.units[unit]*1E-6
        spikes = spikes[np.where(np.logical_and(spikes>=stimuli[0]+offset, spikes<=stimuli[-1]+offset+1.0))[0]]     #Save all spikes that fall within experiment time
        spikes_array.append(spikes-offset)
        all_spikes.extend(spikes-offset)
           
           
    #plot spike rasters
    #for unit in range(len(Sort_sua.units)):
    unit = int(self.starting_cell.text())

    ax = plt.subplot(2,2,1)
    all_spikes = []
    print "...unit: ", unit
    plt.ylim(0,400)
    plt.xlim(0,5.5)
    depth=0
    for stimulus in stimuli:
        ymin=np.zeros(len(spikes_array[unit]))
        ymax=np.zeros(len(spikes_array[unit]))
        ymin+=depth
        ymax+=depth+.9
        depth+=1
        
        spikes = spikes_array[unit][np.where(np.logical_and(spikes_array[unit]>=stimulus+seq_isi_start, spikes_array[unit]<=stimulus+seq_isi_end))[0]]-stimulus
        all_spikes.append(spikes)
        
        plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=.8) #colors[mod(counter,7)])

    #Plot distribution of spikes
    ax = plt.subplot(2,2,3)
    plt.xlim(0,5.5)

    bin_width = 0.001   #1ms bins
    y = np.histogram(np.hstack(all_spikes), bins = np.arange(0, 5.5, bin_width))
    plt.plot(y[1][:-1], y[0], linewidth=3, color='blue')    
    
    #Plot distribution of sequential pairwise distances
    ax = plt.subplot(2,2,2)
    plt.xlim(0,0.01)
    sequential_isi = []
    for k in range(len(all_spikes)-1):
        for p in range(len(all_spikes[k+1])):
            if (len(all_spikes[k])>0) and (len(all_spikes[k+1])>0):
                sequential_isi.extend(np.abs(all_spikes[k]-all_spikes[k+1][p]))
    
    print np.sort(sequential_isi)[:100]
    if len(sequential_isi)>0:
        bin_width = float(self.bin_width.text())
        y = np.histogram(np.hstack(sequential_isi), bins = np.arange(0,0.150,bin_width))
        plt.plot(y[1][:-1], y[0], linewidth=3, color='blue')    

    #Plot scrambled order pairwise distances
    shuffled_indexes = np.random.choice(len(all_spikes)-1,len(all_spikes)) 
    
    sequential_isi = []
    for k in range(len(shuffled_indexes)-1):
        print shuffled_indexes[k], shuffled_indexes[k+1]
        for p in range(len(all_spikes[shuffled_indexes[k+1]])):
            if (len(all_spikes[shuffled_indexes[k]])>0) and (len(all_spikes[shuffled_indexes[k+1]])>0):
                sequential_isi.extend(np.abs(all_spikes[shuffled_indexes[k]]-all_spikes[shuffled_indexes[k+1]][p]))
    
    print np.sort(sequential_isi)[:100]
    if len(sequential_isi)>0:
        y = np.histogram(np.hstack(sequential_isi), bins = np.arange(0,0.150,bin_width))
        plt.plot(y[1][:-1], y[0], linewidth=3, color='red')    

    #PLOT SCRAMBLED RASTERS
    ax = plt.subplot(2,2,1)
    print "...unit: ", unit
    plt.ylim(0,400)
    plt.xlim(0,5.5)
    depth=0
    for k in range(len(all_spikes)):
        ymin=np.zeros(len(all_spikes[shuffled_indexes[k]]))
        ymax=np.zeros(len(all_spikes[shuffled_indexes[k]]))
        ymin+=depth
        ymax+=depth+.9
        depth+=1
        
        spikes = all_spikes[shuffled_indexes[k]]
        
        plt.vlines(spikes, ymin, ymax, linewidth=3, color='red',alpha=.8) #colors[mod(counter,7)])

    print all_spikes
    



    plt.show()



def compute_natscene_rasters(self):

    #Find offset for stimulus processing
    recordings = sorted(glob.glob(self.rec_path+"/*"))
    offset = 0
    for recording in recordings:
        print recording
        temp= recording.replace(self.rec_path+'/','')
        if temp in self.rec_name: 
            din_file = self.rec_path+'/'+temp+'/'+temp+'.din'       #Save .din filename before exiting
            break
        t1 = load_lfp_length(glob.glob(recording+"/*.lfp.zip")[0])
        offset += t1

    offset = offset*1E-6
    print "...offset: ", offset



    #**** Load .din stimulus timestamps
    print din_file
    f = open(din_file, 'rb')
    din = np.fromfile(f, dtype=np.int64).reshape(-1, 2) # reshape to 2 cols containing [time_stamp in usec : frame number]
    f.close()
    print "No of frames in din file: ", len(din)
    print "Length of recording: ", float(din[-1][0])*1E-6/60., " mins." 

    start_indexes = np.where(din[:,1]==0)[0][::3]   #each frame is on screen for 3 refreshes (NOT ALWAYS TRUE*********)
    stimuli = din[:,0][start_indexes]*1E-6
    #print stimuli
       
    

    #Load spike rasters
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    spikes_array = []
    all_spikes=[]
    for unit in range(len(Sort_sua.units)):
        spikes = Sort_sua.units[unit]*1E-6
        spikes = spikes[np.where(np.logical_and(spikes>=stimuli[0]+offset, spikes<=stimuli[-1]+offset+1.0))[0]]
        spikes_array.append(spikes-offset)
        all_spikes.extend(spikes-offset)
           
           
    #plot spike rasters
    #for unit in range(len(Sort_sua.units)):
    unit = int(self.starting_cell.text())

    for unit in [unit]:
        print "...unit: ", unit
        #ax = plt.subplot(6,7,unit+1)
        plt.ylim(0,400)
        plt.xlim(0,5.5)
        depth=0
        for stimulus in stimuli:
            ymin=np.zeros(len(spikes_array[unit]))
            ymax=np.zeros(len(spikes_array[unit]))
            ymin+=depth
            ymax+=depth+.9
            depth+=1
            
            spikes = spikes_array[unit][np.where(np.logical_and(spikes_array[unit]>=stimulus, spikes_array[unit]<=stimulus+5.5))[0]]-stimulus

            
            plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=1) #colors[mod(counter,7)])
        
        
    plt.show()


def compute_triplet_sequences(self):
    #np.savez(out_file, sequences=out_x, controls=out_color, lfp_events = pop_spikes)    

    #Load SUA
    Sort_sua = Ptcs(self.parent.sua_file)
    for k in range(len(Sort_sua.units)): print k, len(Sort_sua.units[k])

    #Load LFP
    Sort_lfp = Ptcs(self.parent.lfp_event_file) 
    lfp_cluster = int(self.parent.lfp_cluster.text())
    lfp_events = Sort_lfp.units[lfp_cluster]*1E-6
    print lfp_events[:10]

    #Load saved sequences
    lfp_cluster = int(self.parent.lfp_cluster.text())
    filename = self.parent.sua_file.replace('.ptcs','')+'_sequences_LPF_'+str(lfp_cluster)
    data = np.load(filename+'.npz')
    
    sequences = data['sequences']
    controls = data['controls']
    print len(sequences)
    print len(sequences[0])

    #Compute relative locations of units
    units = [11, 49, 65]
    #units = [72, 99, 100]
    #units = [49, 60, 108]

    #Plot LFP drift data
    ax=plt.subplot(3,1,1)
    for k in range(len(sequences)):
        for p in units:
            plt.scatter(k, sequences[k][p], color=controls[k][p]*float(p)/len(sequences[0]))
    
    plt.xlim(0, len(sequences)+100)
    plt.tick_params(axis='both', which='both', labelsize=18)

    #Plot LFP drift data - relative first 
    ax=plt.subplot(3,1,2)
    for k in range(len(sequences)):
        #print k, 
        for p in units[1:]:
            #print sequences[k][0]-sequences[k][p],
            #print sequences[k][p],
            plt.scatter(k, sequences[k][units[0]]-sequences[k][p], color=controls[k][p]*float(p)/len(sequences[0]))
        #print ''
    plt.xlim(0, len(sequences)+100)


    #***************************************************************************
    #Compute pairwise distances
    spikes = []
    for p in units:
        spikes.append(Sort_sua.units[p]*1E-6)


    locations = np.zeros ((len(sequences)+100,2), dtype=np.float32)
    alphas = np.zeros((len(sequences)+100,2), dtype=np.float32)
    min_dist = 0.020
    for s in range(len(spikes[0])):
        t0 = spikes[0][s]

        #Check nearest LFP event
        temp_lfp = lfp_events[find_nearest(lfp_events,t0)]
        #print "...lfp distance: ", temp_lfp-t0
        if abs(temp_lfp-t0)<min_dist: continue

        #Find nearest spikes from other two units
        temp_nearest1 = spikes[1][find_nearest(spikes[1],t0)]
        temp_nearest2 = spikes[2][find_nearest(spikes[2],t0)]
        
        #print spikes[0][s], temp_nearest1, temp_nearest2
        dist1 = t0 - temp_nearest1
        dist2 = t0 - temp_nearest2

        #Check for spikes within 50ms of each other
        if (abs(dist1))<min_dist: 
            locations[int(t0/60.),0] = dist1*1E3
            alphas[int(t0/60.),0] = 1.0

        if (abs(dist2))<min_dist: 
            locations[int(t0/60.),1] = dist2*1E3
            alphas[int(t0/60.),1] = 1.0

        #if s > 10000: break
    plt.ylabel("LFP Time Lock (ms)", fontsize = 20)
    plt.tick_params(axis='both', which='both', labelsize=18)


    #Compute least squares fit
    ax = plt.subplot(3,1,3)

    # Fit a polynomial 
    x = np.arange(0,len(locations),1)
    y = locations[:, 0]

    trialX = np.linspace(x[0],x[-1],1000)
    fitted = np.polyfit(x, y, 12)[::-1]
    y = np.zeros(len(trialX))
    for i in range(len(fitted)):
       y += fitted[i]*trialX**i
    plt.plot(trialX,   y, 'r--', color='blue', linewidth=3)

    y = locations[:, 1]

    trialX = np.linspace(x[0],x[-1],1000)
    fitted = np.polyfit(x, y, 12)[::-1]
    y = np.zeros(len(trialX))
    for i in range(len(fitted)):
       y += fitted[i]*trialX**i
    plt.plot(trialX,   y, 'r--', color='magenta', linewidth=3)



    print len(locations)
    for k in range(len(locations)):
        plt.scatter(k, locations[k,0], color='green', alpha=alphas[k,0])
        plt.scatter(k, locations[k,1], color='red', alpha=alphas[k,1])

    plt.xlim(0, len(sequences)+100)
    plt.ylim(-5, +5)
    plt.tick_params(axis='both', which='both', labelsize=18)
    plt.ylabel("Pairwise ISI\n(ms)", fontsize = 20)
    plt.xlabel("Recording Time\n(mins)", fontsize = 20)
    plt.show()


def compute_csd_rasters(self):
    
    #Find offset for stimulus processing
    recordings = sorted(glob.glob(self.rec_path+"/*"))
    offset = 0
    for recording in recordings:
        temp= recording.replace(self.rec_path+'/','')
        if temp in self.rec_name: break
        t1 = load_lfp_length(glob.glob(recording+"/*.lfp.zip")[0])
        offset += t1

    offset = offset*1E-6
    print "...offset: ", offset


    #Load ontimes of stimulus; use existing value in 'recording' variable
    #print recording+temp+"/*ontimes.txt"
    stimuli = np.loadtxt(glob.glob(recording+"/*ontimes.txt")[0])
    #print stimuli
    stimuli=stimuli*1E-6

    ax = plt.subplot(1,1,1)
    for stimulus in stimuli:
        plt.axvspan(stimulus, stimulus+0.100, facecolor='black', alpha=0.2)

    #Load spike rasters
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    spikes_array = []
    all_spikes=[]
    for unit in range(len(Sort_sua.units)):
        spikes = Sort_sua.units[unit]*1E-6
        spikes = spikes[np.where(np.logical_and(spikes>=stimuli[0]+offset, spikes<=stimuli[-1]+offset+1.0))[0]]
        spikes_array.append(spikes-offset)
        all_spikes.extend(spikes-offset)

    #plot spike rasters
    for spikes in sorted(spikes_array, key=len):
        ymin=np.zeros(len(spikes))
        ymax=np.zeros(len(spikes))
        ymin+=depth
        ymax+=depth+.9
        depth+=1

        plt.vlines(spikes, ymin, ymax, linewidth=3, color='black',alpha=1) #colors[mod(counter,7)])

    #plot MUA histogram
    bin_width = 0.100   #100ms bins
    y = np.histogram(all_spikes, bins = np.arange(0,np.max(all_spikes),bin_width))
    #plt.bar(y[1][:-1], y[0], bin_width, color='blue')
    plt.plot(y[1][:-1], y[0]-40, linewidth=3, color='blue')

    #plot lfp data
    tsf = load_lfp_all(recording+'/'+temp+'.lfp.zip')
    
    t = np.linspace(0, tsf.n_vd_samples/float(tsf.SampleFrequency), tsf.n_vd_samples)
    traces = tsf.ec_traces[9]/1000.
    plt.plot(t, traces-15, color='red', linewidth=3)
    
    #Labels
    plt.plot([0,tsf.n_vd_samples/float(tsf.SampleFrequency)], [0, 0], 'r--', linewidth = 3, alpha=0.5)
    
    plt.ylabel(" MUA  LFP   Units", fontsize=25)
    plt.xlabel("Time (sec)", fontsize= 25)
    plt.tick_params(axis='both', which='both', labelsize=20)

    plt.show()

def compute_csd_histogram(self):

    #Find offset for stimulus processing
    recordings = sorted(glob.glob(self.rec_path+"/*"))
    offset = 0
    for recording in recordings:
        temp= recording.replace(self.rec_path+'/','')
        if temp in self.rec_name: break
        t1 = load_lfp_length(glob.glob(recording+"/*.lfp.zip")[0])
        offset += t1

    offset = offset*1E-6
    print "...offset: ", offset


    #Load ontimes of stimulus; use existing value in 'recording' variable
    #print recording+temp+"/*ontimes.txt"
    stimuli = np.loadtxt(glob.glob(recording+"/*ontimes.txt")[0])
    stimuli=stimuli*1E-6

    #Load spike rasters
    Sort_sua = Ptcs(self.parent.sua_file)
    depth = 0
    spikes_array = []
    all_spikes=[]
    for unit in range(len(Sort_sua.units)):
        spikes_array.append([])
        spikes = Sort_sua.units[unit]*1E-6
        for stimulus in stimuli: 
            temp_spikes = spikes[np.where(np.logical_and(spikes>=stimulus+offset, spikes<=stimulus+offset+1.0))[0]]
            spikes_array[unit].extend(temp_spikes-stimulus-offset)
            #print temp_spikes
    
    
    #Plot histograms
    #plot MUA histogram
    bin_width = 0.025   #100ms bins

    for unit in range(len(spikes_array)):
        y = np.histogram(spikes_array[unit], bins = np.arange(0,1,bin_width))
        plt.plot(y[1][:-1], y[0]+10*unit, linewidth=1, color='blue')

    #Labels
    
    plt.ylabel("Unit Fire Rate (au)", fontsize=25)
    plt.xlabel("Time from CSD (sec)", fontsize= 25)
    plt.tick_params(axis='both', which='both', labelsize=20)
    plt.title(temp)

    plt.show()


def tsf_subsample(self):
    
    tsf = TSF.TSF(self.selected_recording)
    tsf.read_ec_traces()
    
    tsf.SampleFrequency = tsf.SampleFrequency/2
    
    temp_traces=[]
    for k in range(tsf.n_electrodes):
        temp_traces.append(tsf.ec_traces[k][::2])

    tsf.ec_traces = np.int16(temp_traces)
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    
    tsf.file_names=[]; tsf.n_samples=[]; tsf.n_digital_chs=[]; tsf.digital_chs=[]
    tsf.layout = np.arange(tsf.n_electrodes)    
    
    if tsf.n_cell_spikes!=0:
        tsf.fake_spike_times = np.int32(tsf.fake_spike_times/2.)
        #np.int32(tsf.fake_spike_times).tofile(fout)
        
    save_tsf_single(tsf, self.selected_recording[:-4]+"_subsampled.tsf")

def load_lfp_all(file_name):     #Nick/Martin data has different LFP structure to their data.
    
    class Object_empty(object):
        def __init__(self):
            pass
    
    tsf = Object_empty()
    tsf.file_name = file_name
    tsf.iformat = 1002
    tsf.header = 'Test spike file '
    
    data_in = np.load(file_name)
    tsf.SampleFrequency = data_in['tres']
    tsf.chans = data_in['chans']
    tsf.n_electrodes = len(tsf.chans)
    tsf.Siteloc = data_in['chanpos']
    tsf.vscale_HP = data_in['uVperAD']
    tsf.ec_traces = data_in['data']
    tsf.ec_traces = tsf.ec_traces

    #Home made notch filter; filter.notch doesn't seem to work...
    offset = np.zeros(int(data_in['t0']*1E-3), dtype=np.int16)      #Convert microsecond offset to miliseconds;
    temp_traces = []
    for k in range(tsf.n_electrodes):
        tsf.ec_traces[k] = Notch_Filter(tsf.ec_traces[k])
        temp_traces.append(np.append(offset, tsf.ec_traces[k]))
    
    tsf.ec_traces = np.int16(temp_traces)
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    
    return tsf

def sua_lock_percentage(self):

    colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']


    #******************************* OLD CODE *************************
    if False:
        
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
        tsf = TSF.TSF(self.parent.sua_file.replace('_hp.ptcs','_lp.tsf')) 
        tsf.read_trace(int(self.specgram_ch.text()))
           
            
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
        start_window = self.start_window*1E-3                            #*********************************** CODE THESE INTO THE GUI AS TEXT BOXES
        end_window = self.end_window*1E-3
        min_lfp_isi = 0.100   #Lockout lfp events more 
        
        si_limit = 0.7                                  #*********************************** MAYBE ALSO THIS!?
        
        

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





#self.parent.animal.ptcsName = self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'

    min_spikes = float(self.min_spikes.text())

    print self.parent.sua_file 
    print self.parent.lfp_event_file

    font_size = 25

    top_channel = np.loadtxt(os.path.split(os.path.split(self.parent.sua_file)[0])[0]+"/top_channel.txt") - 1      #Load top channel for track; convert to 0-based ichannel values.


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
    n_units_incortex = len(np.where(Sort_sua.maxchan>=top_channel)[0])                                          #Number of units that are in tissue, i.e. below the top_channel.txt (see file)


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



    #************************** COMPUTE EXPECTED LOCK *******************
    all_pop_spikes = 0
    for k in range(Sort_lfp.n_units):
        all_pop_spikes+=len(Sort_lfp.units[k])
    
    sync_period_length = 0
    for k in range(Sort_sua.n_units):
        if np.max(Sort_sua.units[k])>sync_period_length: 
            sync_period_length = np.max(Sort_sua.units[k])

    sync_period_length = sync_period_length * 1E-6
    start_window = float(self.start_window.text())                            #*********************************** CODE THESE INTO THE GUI AS TEXT BOXES
    end_window = float(self.end_window.text())
    
    print "# pop spikes: ", all_pop_spikes
    print "len sync_period: ", sync_period_length
    exp_lock = all_pop_spikes*(end_window-start_window)*1E-3/sync_period_length*1E2
    print "expected lock: ", exp_lock
    
    #************************* FIND PEAK ACROSS ALL LFP 


    f3 = plt.figure()
    
    #************************* PLOT % PLOTS BY DEPTH *********************
    ax = plt.subplot(1,2,1)
    img_means = []

    x = np.arange(0,n_units_incortex,1)


    cumulative_bars = np.zeros(n_units_incortex, dtype=np.float32)
    no_clusters = 10
    for k in range(no_clusters):

        
        #if k==0: continue

        cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(k)+".npy"
        if os.path.exists(cell_rasters_filename): 
            
            cell_rasters = np.load(cell_rasters_filename)

            percent_array = []
            for unit in range(total_units):
                if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it

                locked_spikes = np.hstack(np.array(cell_rasters[unit]))*1E-3        #Convert from usec to msec
                print "lfp: ", k, " unit: ", unit, "...# locked spikes to LFP event: ", len(locked_spikes), " / ", len(Sort_sua.units[unit]), 
                
                #y = np.histogram(locked_spikes, bins = np.arange(-100,100,5))
                #window_offset = np.argmax(y[0], axis=0)*5-100
                window_offset = 0
                
                indexes = np.where(np.logical_and(locked_spikes>=start_window+window_offset, locked_spikes<=end_window+window_offset))[0]   #Look for spikes between -50 to +50msec around LFP time

                print "... % within window: ", round(float(len(indexes))/len(Sort_sua.units[unit])*1E2, 2), "%"
                percent_array.append(float(len(indexes))/len(Sort_sua.units[unit])*1E2) #Convert to %
                
            #print percent_array

            #plt.bar(x, percent_array, 1, color='blue')
            #print len(x), len(percent_array)
            #p2 = plt.barh(x, percent_array, 0.95, bottom=cumulative_bars, color = colors[k])
            p2 = plt.barh(x, percent_array, 0.95, left=cumulative_bars, color = colors[k])

            
            cumulative_bars=cumulative_bars + np.float32(percent_array)

    plt.plot([exp_lock, exp_lock], [Sort_sua.n_units-0.5,-0.5], 'r--', color='cyan', linewidth = 3, alpha=0.8)

    plt.xlim(0,100)
    plt.xlabel("Percent Locking", fontsize = font_size, fontweight = 'bold')
    plt.ylabel("Deep <---- Depth ----> Superficial", fontsize = font_size, fontweight = 'bold')
    plt.ylim(Sort_sua.n_units-0.5,-0.5)

    plt.tick_params(axis='both', which='both', labelsize=font_size)


    #*********************** PLOT % PLOTS BY FIRING RATE ********************
    ax = plt.subplot(1,2,2)

    #Order units by firing rate
    spike_arrays = []
    for k in range(Sort_sua.n_units):
        spike_arrays.append(len(Sort_sua.units[k]))
    indexes_frate = np.array(spike_arrays).argsort()

    x = np.arange(0,n_units_incortex,1)

    cumulative_bars = np.zeros(n_units_incortex, dtype=np.float32)
    for k in range(no_clusters):
        cell_rasters_filename = self.parent.sua_file.replace('.ptcs','')+"_cell_rasters_lfp"+str(k)+".npy"
        if os.path.exists(cell_rasters_filename): 
            
            cell_rasters = np.load(cell_rasters_filename)

            percent_array = []
            for unit in indexes_frate:
                if Sort_sua.maxchan[unit]<top_channel: continue         #If unit is above cortex exclude it

                locked_spikes = np.hstack(np.array(cell_rasters[unit]))*1E-3        #Convert from usec to msec
                print "...# locked spikes to LFP event: ", len(locked_spikes),
                
                y = np.histogram(locked_spikes, bins = np.arange(-100,100,5))
                window_offset = np.argmax(y[0], axis=0)*5-100
                #print window_offset
                #plt.plot(y[1][:-1], y[0], linewidth=1, color='blue')
                #print window_offset
                #plt.show()
                
                indexes = np.where(np.logical_and(locked_spikes>=start_window+window_offset, locked_spikes<=end_window+window_offset))[0]   #Look for spikes between -50 to +50msec around LFP time

                print round(float(len(indexes))/len(Sort_sua.units[unit])*1E2, 2), "%"
                percent_array.append(float(len(indexes))/len(Sort_sua.units[unit])*1E2) #Convert to %
                
            #print percent_array

            #plt.bar(x, percent_array, 1, color='blue')
            #print len(x), len(percent_array)
            #p2 = plt.barh(x, percent_array, 0.95, bottom=cumulative_bars, color = colors[k])
            p2 = plt.barh(x, percent_array, 0.95, left=cumulative_bars, color = colors[k])

            
            cumulative_bars=cumulative_bars + np.float32(percent_array)
   
    plt.plot([exp_lock, exp_lock], [Sort_sua.n_units-0.5,-0.5], 'r--', color='cyan', linewidth = 3, alpha=0.8)

    plt.xlim(0,100)
    plt.yticks([])
    plt.xlabel("Percent Locking", fontsize = font_size, fontweight = 'bold')
    plt.ylabel("Higher <--- Firing Rate ---> Lower", fontsize = font_size, fontweight = 'bold')#, labelpad=-10)
    plt.ylim(Sort_sua.n_units-0.5,-0.5)

    plt.tick_params(axis='both', which='both', labelsize=font_size)

    plt.suptitle("Window width: "+str(end_window-start_window)+"(ms)", fontsize=font_size)

    print "....expected % lock: ", exp_lock
    
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


    
def dim_reduction_general(matrix_in, method, filename, recompute=False):
    
    import sklearn
    from sklearn import metrics, manifold
    
    methods = ['MDS', 'tSNE', 'PCA', 'BH_tSNE']    
    
    print "Computing dim reduction, size of array: ", matrix_in.shape
    
    if method==0:
        #MDS Method - SMACOF implementation Nelle Varoquaux
        
        if os.path.exists(filename+'_'+methods[method]+'.npy')==False or recompute:

            print "... MDS-SMACOF..."
            print "... pairwise dist ..."
            dists = metrics.pairwise.pairwise_distances(matrix_in)
            adist = np.array(dists)
            amax = np.amax(adist)
            adist /= amax
            
            print "... computing MDS ..."
            mds_clf = manifold.MDS(n_components=3, metric=True, n_jobs=-1, dissimilarity="precomputed", random_state=6)
            results = mds_clf.fit(adist)
            Y = results.embedding_         
         
            np.save(filename+'_'+methods[method], Y)
        
        else:
            Y = np.load(filename+'_'+methods[method]+'.npy')


    elif method==1:
        ##t-Distributed Stochastic Neighbor Embedding; Laurens van der Maaten
        if os.path.exists(filename+'_'+methods[method]+'.npy')==False or recompute:


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
        
            np.save(filename+'_'+methods[method], Y)
        
        else:
            Y = np.load(filename+'_'+methods[method]+'.npy')

    elif method==2:

        if os.path.exists(filename+'_'+methods[method]+'.npy')==False or recompute:
            Y = PCA_reduction(matrix_in, 3)
            np.save(filename+'_'+methods[method], Y)
        else:
            Y = np.load(filename+'_'+methods[method]+'.npy')

                
    elif method==3:
        
        from tsne import bh_sne
        if os.path.exists(filename+'_'+methods[method]+'.npy')==False or recompute:

        #if os.path.exists(mouse.home_dir+mouse.name+'/tSNE_barnes_hut.npy')==False:
            print "... computing Barnes-Hut tSNE..."
            Y = bh_sne(np.float64(matrix_in))
        
            np.save(filename+'_'+methods[method], Y)
        else:
            Y = np.load(filename+'_'+methods[method]+'.npy')

    return Y




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
        


def on_click_OLD(event):
    
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



      
def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return b, a


      
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff/nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a
    
    
def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y





def find_nearest(array, value):
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
    
    
def MCD_read_imagingtimes(MCDFilePath):
    np.set_printoptions(suppress=True)      #Supress scientific notation printing

    import neuroshare as ns
 
    #open file using the neuroshare bindings
    fd = ns.File(MCDFilePath)
 
    #create index
    indx = 0
 
    #create empty dictionary
    data = dict()
 
    #loop through data and find all analog entities
 
    for entity in fd.list_entities():
        #print "looping over entities: ", entity

        if entity.entity_type == 1:
            data["extra"] = fd.entities[indx]

    print "... searching trigger on channel 17th .mcd file... looking for epochs..."
    temp_data = []
    for i in range(data['extra'].item_count):
        temp_data.append(data['extra'].get_data(i)) 
    temp_data = np.array(temp_data)[:,0]    #Select time column only from 17th channel

    start_array = []
    end_array = []
    start_array.append(temp_data[0])
    for i in range(1,len(temp_data)-1,1):
        if temp_data[i+1]-temp_data[i]>1.0:
            end_array.append(temp_data[i])
            start_array.append(temp_data[i+1])
    end_array.append(temp_data[-1])
    
    out_data = []
    for i in range(len(start_array)):
        out_data.append([start_array[i], end_array[i]])

    path_dir = os.path.dirname(MCDFilePath)
    np.savetxt(path_dir+'/epochs.txt', out_data, fmt='%5.5f')
    
    filename = os.path.dirname(MCDFilePath)+"/rec_index.txt"
    if os.path.exists(filename): 
        rec_index = int(np.loadtxt(filename))-1
    else:
        rec_index =0

    out_array = np.hstack((round(start_array[rec_index],5), round(end_array[rec_index],5)))     #NEED TO SAVE THE CORRECT PERIOD FOR MULTI EPOCH RECORDINGS
    
    np.savetxt(MCDFilePath[:-4]+'_imagingonoff.txt', out_array, fmt='%5.5f')
    
    #Visualize different epoch starts/ends
    #print temp_data
    #for i in range(len(start_array)):
    #    plt.plot([start_array[i],start_array[i]],[0,len(temp_data)])
    #    plt.plot([end_array[i],end_array[i]],[0,len(temp_data)])
    #    plt.axvspan(start_array[i], end_array[i], color='red', alpha=0.35)

    #plt.show()
    #quit()
    
    
    
