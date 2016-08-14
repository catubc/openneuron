#Animal class and subclasses 
#Code: Cat 

from analysis import *

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
import cPickle as pickle
import gc
from skimage.measure import block_reduce
import shutil
from load_intan_rhd_format import *
from scipy.signal import butter, filtfilt


class Mouse(object):      
    
    def __init__(self, animal_name, home_dir):
        
        self.name = animal_name
        self.home_dir = home_dir

        self.move_files()

        self.load_filenames()
        
        self.probe = Probe()   #Load intan probe map;
        try: 
            self.img_rate = float(np.loadtxt(self.home_dir+self.name+'/img_rate.txt'))
        except:
            self.img_rate = 0
            print "No img rate file..."

    def move_files(self):
        print "... making default directories and moving files..."
        dirs = ['camera_files','movie_files','rhd_files','stm_files','tif_files','tsf_files']
        for dir_ in dirs:
            new_dir = self.home_dir+self.name+'/'+dir_
            if not os.path.exists(new_dir):
                print "...making dir: ",  new_dir;     os.makedirs(new_dir)
        
        #Move .bin and .rhd files to correct foloders
        temp_names = glob.glob(self.home_dir+self.name+'/*.rhd')  
        for r in range(len(temp_names)):
            new_file = self.home_dir+self.name+'/rhd_files'+temp_names[r].replace(self.home_dir+self.name,'')
            os.rename(temp_names[r], new_file)

        #Move .bin and .rhd files to correct foloders
        temp_names = glob.glob(self.home_dir+self.name+'/*.bin')  
        for r in range(len(temp_names)):
            new_file = self.home_dir+self.name+'/tif_files'+temp_names[r].replace(self.home_dir+self.name,'')
            os.rename(temp_names[r], new_file)
            
        
    
    def load_filenames(self):
        print "...loading filenames..."
        self.filenames = glob.glob(self.home_dir+self.name+'/rhd_files/*.rhd')   #use .tif extension otherwise will load .npy files
        for f in range(len(self.filenames)):
            self.filenames[f] = self.filenames[f][:-4]
        
        #self.tsf_filenames = glob.glob( self.home_dir+self.name+'/tsf_files/*.tsf')   #use .tif extension otherwise will load .npy files
    
    
    def load_tsf_header(self, file_name):
        self.tsf = Tsf_file(file_name)


    def load_tsf(self, file_name):
        self.tsf = Tsf_file(file_name)
        self.tsf.read_ec_traces()
       
        
    def load_channel(self, file_name, channel):
        
        self.tsf = Tsf_file(file_name)
        self.tsf.file_name = file_name
        self.tsf.read_trace(channel)
    
    
    def rhd_digital_save(self):
        '''Read .rhd files, and save digital channels.
        NB: there can be 2, 4 or 6 digital channels inside Intan file
        chs 1 and 2 are laser meta data and laser pulse times (these are off for other experiments)
        chs 3 and 4 are camera pulse times and on/off times from clampx computer (these are chs 1 and 2 usually as laser chs are off)
        chs 5 and 6 are vis stim pulse times and meta data (these are usually 3 and 4 as not recorded w. laser on)
        '''
        
        print "...reading digital amp data..."

        for file_name in self.filenames:

            camera_frames_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:]+'_camera_pulses'
            camera_onoff_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:]+'_camera_onoff'
            
            if os.path.exists(camera_onoff_filename+'.npy')==True: continue

            data = read_data(file_name+'.rhd')

            SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
            print "SampleFrequency: ", SampleFrequency

            counter=0
            #laser_pulses = data['board_dig_in_data'][counter]; np.save(laser_filename, laser_pulses); counter+=1
            #meta_data = data['board_dig_in_data'][counter]; np.save(meta_filename, meta_data); counter+=1
            
            camera_frames = data['board_dig_in_data'][counter]; np.save(camera_frames_filename, camera_frames); counter+=1
            camera_onoff = data['board_dig_in_data'][counter]; np.save(camera_onoff_filename, camera_onoff); counter+=1
            
            #Save stim info as well
            if len(data['board_dig_in_data'])>2:
                stim_pulses_filename = file_name[0:file_name.find('rhd_files')]+'stm_files/'+file_name[file_name.find('rhd_files')+10:]+'_stim_pulses'
                stim_meta_filename = file_name[0:file_name.find('rhd_files')]+'stm_files/'+file_name[file_name.find('rhd_files')+10:]+'_stim_meta'

                stim_pulses = data['board_dig_in_data'][counter]; np.save(stim_pulses_filename, stim_pulses); counter+=1
                stim_meta = data['board_dig_in_data'][counter]; np.save(stim_meta_filename, stim_meta); counter+=1
                
    

    def rhd_to_tsf(self):
        '''Read .rhd files, convert to correct electrode mapping and save to .tsf file
        NB: There are 2 possible mapping depending on the insertion of the AD converter 
        TODO: implement a wavelet high pass filter directly to avoid SpikeSorter Butterworth filter artifacts
        '''
        
        print "...reading amp data..."

        for file_name in self.filenames:
            ec_traces = 0.; ec_traces_hp = 0.; data=0. #Delete previous large arrays;
            
            file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_hp.tsf'
            if os.path.exists(file_out)==True: continue

            print "Processing: \n", file_name+'.rhd'

            data = read_data(file_name+'.rhd')
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

            #SAVE RAW DATA - ******NB:  SHOULD CLEAN THIS UP: the write function should be shared by all, just data is changing so no need to repeat;
            if True:
                print "Writing raw data ..."
                #print "CHANGE THIS TO WORK THROUGH FUNCTION WITHOUT REPEATING"
                file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_raw.tsf'
                fout = open(file_out, 'wb')
                fout.write(header)
                fout.write(struct.pack('i', 1002))
                fout.write(struct.pack('i', SampleFrequency))
                fout.write(struct.pack('i', self.probe.n_electrodes))
                fout.write(struct.pack('i', n_vd_samples))
                fout.write(struct.pack('f', vscale_HP))
                
                for i in range (self.probe.n_electrodes):
                    fout.write(struct.pack('h', self.probe.Siteloc[i][0]))
                    fout.write(struct.pack('h', self.probe.Siteloc[i][1]))
                    fout.write(struct.pack('i', i+1))

                for i in range(self.probe.n_electrodes):
                    print i,
                    ec_traces[self.probe.layout[i]].tofile(fout)  #Frontside

                fout.write(struct.pack('i', n_cell_spikes))
                fout.close()
                
            #SAVE HIGH PASS WAVELET FILTERED DATA
            if True:
                print "Writing hp data ..."
                file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_hp.tsf'
                fout = open(file_out, 'wb')
                fout.write(header)
                fout.write(struct.pack('i', 1002))
                fout.write(struct.pack('i', SampleFrequency))
                fout.write(struct.pack('i', self.probe.n_electrodes))
                fout.write(struct.pack('i', n_vd_samples))
                fout.write(struct.pack('f', vscale_HP))
                
                for i in range (self.probe.n_electrodes):
                    fout.write(struct.pack('h', self.probe.Siteloc[i][0]))
                    fout.write(struct.pack('h', self.probe.Siteloc[i][1]))
                    fout.write(struct.pack('i', i+1))

                print "Wavelet filtering..."
                ec_traces_hp = wavelet(ec_traces, wname="db4", maxlevel=6)
                print ec_traces_hp.shape

                for i in range(self.probe.n_electrodes):
                    print i,
                    ec_traces_hp[self.probe.layout[i]].tofile(fout)  #Frontside

                fout.write(struct.pack('i', n_cell_spikes))
                fout.close()
            
            #OPTIONAL: SAVE LOW PASS DATA
            if False:
                file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp.tsf'
                fout = open(file_out, 'wb')
                fout.write(header)
                fout.write(struct.pack('i', 1002))
                fout.write(struct.pack('i', SampleFrequency))
                fout.write(struct.pack('i', self.probe.n_electrodes))
                fout.write(struct.pack('i', n_vd_samples))
                fout.write(struct.pack('f', vscale_HP))
                
                for i in range (self.probe.n_electrodes):
                    fout.write(struct.pack('h', self.probe.Siteloc[i][0]))
                    fout.write(struct.pack('h', self.probe.Siteloc[i][1]))
                    fout.write(struct.pack('i', i+1))
               
                for i in range(self.probe.n_electrodes):
                    (ec_traces[self.probe.layout[i]]-ec_traces_hp[self.probe.layout[i]]).tofile(fout)  #Frontside

                fout.write(struct.pack('i', n_cell_spikes))
                fout.close()

    def tsf_to_lfp(self):
        '''Read .tsf files - subsample to 1Khz, save as *_lp.tsf
        '''
        
        print "...making low-pass tsf files (1Khz sample rates)..."

        for file_name in self.filenames:
            
            file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp.tsf'
            if os.path.exists(file_out)==True: continue

            file_in = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_raw.tsf'

            print "Processing: \n", file_in

            self.load_tsf(file_in)
            #self.load_channel(file_in, 0)
            print self.tsf.Siteloc.shape
            
            n_vd_samples = len(self.tsf.ec_traces[0]); print "Number of samples: ", n_vd_samples
            
            print "...converting raw to .lfp (1Khz) sample rate tsf files ..."
            temp_array=[]
            lowcut = 0.1; highcut=110; fs=1000
            #import matplotlib.pyplot as plt
            #lowcut = 5; highcut=240; fs=1000    #Use 5Hz low cutoff to reduce slower oscillations for spike sorting;
            for k in range(len(self.tsf.ec_traces)):      #Subsample to 1Khz and notch filter
                print "ch: ", k
                temp = np.array(butter_bandpass_filter(self.tsf.ec_traces[k][::int(self.tsf.SampleFrequency/1000)], lowcut, highcut, fs, order = 2), dtype=np.int16)
                #Home made notch filter; filter.notch doesn't seem to work...
                notch = np.array(butter_bandpass_filter(temp, 59.9, 60.1, fs, order = 2), dtype=np.int16)
                temp = temp-notch
                temp_array.append(temp)
            
            ec_traces = np.int16(temp_array) 
            print ec_traces.shape

            #SAVE RAW DATA ******NB:  SHOULD CLEAN THIS UP: the write function should be shared by all, just data is changing so no need to repeat;
            print "Writing raw data ..."

            fout = open(file_out, 'wb')
            fout.write(self.tsf.header)
            fout.write(struct.pack('i', 1002))
            fout.write(struct.pack('i', 1000))
            fout.write(struct.pack('i', self.tsf.n_electrodes))
            fout.write(struct.pack('i', len(ec_traces[0])))
            fout.write(struct.pack('f', self.tsf.vscale_HP))
            
            for i in range (self.tsf.n_electrodes):
                fout.write(struct.pack('h', self.tsf.Siteloc[i*2]))
                fout.write(struct.pack('h', self.tsf.Siteloc[i*2+1]))
                fout.write(struct.pack('i', i+1))

            for i in range(self.tsf.n_electrodes):
                print i,
                ec_traces[i].tofile(fout)  #Frontside

            fout.write(struct.pack('i', self.tsf.n_cell_spikes))
            fout.close()
            #quit()


    def bin_to_npy(self):
        print "... converting .bin to .npy ..."

        for file_name in self.filenames:
            file_out = file_name[:file_name.find('rhd_files')]+'tif_files/'+ file_name[file_name.find('rhd_files')+10:]+'.npy'
            if os.path.exists(file_out)==True: continue

            print "...reading raw bin: ", file_out
                                    
            #file_name = '/media/cat/12TB/in_vivo/tim/cat/2016_06_09_test/test_laser_imaging'

            #IS THIS DOUBLING DATA SIZE UNECESSARILY?  ORIGINAL DATA MAY BE UINT8 SO NO NEED TO MAKE IT INT16
            data = np.fromfile(file_out[:-4]+'.bin', dtype=np.int16)
            print "...reshaping array..."
            data = data.reshape((-1, 128, 128))
            print "...saving .npy array..."
            np.save(file_out, data)
            
            #import matplotlib.pyplot as plt
            #print data[0]
            #plt.imshow(data[0], cmap=plt.get_cmap('gray'))
            #plt.show()


    #def save_dig_input(self):
        ##Read and save digital input data
        #laser_filename = file_name[0:file_name.find('rhd_files')]+'laser_files/'+file_name[file_name.find('rhd_files')+10:]+'_laser_times'
        #meta_filename = file_name[0:file_name.find('rhd_files')]+'laser_files/'+file_name[file_name.find('rhd_files')+10:]+'_meta_data'
        #camera_pulses_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:]+'_camera_pulses'
        #camera_onoff_filename = file_name[0:file_name.find('rhd_files')]+'camera_files/'+file_name[file_name.find('rhd_files')+10:]+'_camera_onoff'

        #SampleFrequency = data['frequency_parameters']['board_adc_sample_rate']
        #print "SampleFrequency: ", SampleFrequency

        #counter=0
        ##laser_pulses = data['board_dig_in_data'][counter]; np.save(laser_filename, laser_pulses); counter+=1
        ##meta_data = data['board_dig_in_data'][counter]; np.save(meta_filename, meta_data); counter+=1
        ##camera_frames = data['board_dig_in_data'][counter]; np.save(camera_frames_filename, camera_frames); counter+=1
        #camera_onoff = data['board_dig_in_data'][counter]; np.save(camera_onoff_filename, camera_onoff); counter+=1
        ##plt.plot(camera_onoff)
        ##plt.show()

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
        #print " ... reading data.....chs: ", self.n_electrodes, " ... nsamples: ", self.n_vd_samples
        #self.ec_traces =  np.fromfile(self.fin, dtype=np.int16, count=self.n_electrodes*self.n_vd_samples)
        #self.ec_traces.shape = self.n_electrodes, self.n_vd_samples

        #self.n_cell_spikes = struct.unpack('i',self.fin.read(4))[0] 
        
        #print "No. ground truth cell spikes: ", self.n_cell_spikes
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
 


def convert_tif_DUPLICATE(self):

    if (os.path.exists(self.tif_file[:-4] +'.npy')==False):
        print "...read: ", self.tif_file
        images_raw = tiff.imread(self.tif_file)

        print "... saving .npy"
        np.save(self.tif_file[:-4], images_raw)



def find_previous(self, array, value):
    temp = (np.abs(array-value)).argmin()
    if array[temp]>value: return temp-1
    else: return temp
