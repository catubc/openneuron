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


class Object_empty(object):
    def __init__(self):
        pass
            
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
        
def save_tsf(tsf,file_name):
    
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
        fout.write(struct.pack('i', i))

    tsf.ec_traces.tofile(fout)

    fout.write(struct.pack('i', tsf.n_cell_spikes))
    #print "No. cell spikes: ", n_cell_spikes
    #if (n_cell_spikes>0):
    #    if (iformat==1001):
    #        fout.write(struct.pack('i', vertical_site_spacing)) # = struct.unpack('i',fin.read(4))[0] 
    #        fout.write(struct.pack('i', n_cell_spikes)) #Write # of fake spikes
    #    fake_spike_times.tofile(fout)
    #    fake_spike_assignment.tofile(fout) 
    #    fake_spike_channels.tofile(fout) 
    fout.close()

        

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
    
    self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'_blue_light_frames.npy'
    
    if os.path.exists(self.blue_light_filename)==True: 
        print "...Blue Light Boundaries already found... returning..."
        return

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
    
    movie_data = np.load(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'.npy')
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

    if plotting: 
        ax = plt.subplot(2,1,2)
        plt.plot(blue_light_roi)
        plt.ylim(bottom=0)
        plt.show()
        

def event_triggered_movies(self):
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


def event_triggered_movies_Ca(self):
    """ Load [Ca] imaging and behavioural camera data and align to selected trial"""

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
    self.blue_light_filename = self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'_blue_light_frames.npy'
    self.movie_data = movie_data[np.load(self.blue_light_filename)]
    
    #Find movie frame corresponding to lever pull trigger
    movie_times = np.linspace(0, self.abstimes[-1], self.movie_data.shape[0])
    self.movie_04frame_locations = find_nearest(movie_times, self.selected_locs_44threshold)
    print "... frame event triggers: ", self.movie_04frame_locations

    #Make movie stack
    temp_img_rate = 15
    self.movie_stack = self.movie_data[self.movie_04frame_locations-3*temp_img_rate: self.movie_04frame_locations+3*temp_img_rate]

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
        print j
        plt.suptitle(self.selected_dff_filter+'  ' +self.dff_method + "\nFrame: "+str(j)+"  " +str(format(float(j)/30-3.,'.2f'))+"sec", fontsize = 15)

        # set the data in the axesimage object
        im[0].set_array(self.ca_stack[j])
        im[1].set_array(self.movie_stack[j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(self.movie_stack)), interval=100, blit=False, repeat=True)

    if True:
        ani.save(self.parent.root_dir+self.parent.animal.name+"/video_files/"+self.selected_session+'_'+str(len(self.movie_stack))+'_'+str(self.selected_trial)+'trial.mp4', writer=writer)
    plt.show()



def filter_data(self):
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
        print "...data already filtered..."
        return

    #Load aligned images
    print "... loading aligned imgs..."
    images_aligned = np.load(images_file)
    
    #Save mean of images_aligned if not already done
    if os.path.exists(images_file[:-4]+'_mean.npy')==False: 
        np.save(images_file[:-4]+'_mean', np.mean(images_aligned, axis=0))
            
    #Load mask for display:
    n_pixels = len(images_aligned[0])
    generic_coords = np.loadtxt(self.parent.animal.home_dir + self.parent.animal.name+'/genericmask.txt')
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)): generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True
        
    if plotting: f = plt.figure(); ax = f.gca(); f.show()
    
    #Butter
    if filter_type == 'butterworth':
        data_out = images_aligned.copy()*0.0
        for row in range(len(images_aligned[0])):
            print "...filtering row: ", row
            for col in range(len(images_aligned[0,row])):
                data_out[:,row,col] = butter_bandpass_filter(images_aligned[:,row,col], lowcut, highcut, fs, order = 2)
            
            if plotting: 
                ax.clear(); ax.imshow(np.ma.masked_array(data_out[1000], mask=generic_mask_indexes))
                ax.set_xticks([]); ax.set_yticks([]) 
                ax.set_title("Lowcut: "+str(lowcut)+"hz, highcut: "+ str(highcut)+ "hz" + #, maxDFF: " + str(round(np.nanmax(data_out[1000]),1)) + " minDFF: " + str(round(np.nanmin(data_out[1000]),1))+ 
                            "\nFrame: 1000 / "+str(len(images_aligned)))
                plt.pause(0.000001) 
        
        print "... saving filtered data..."
        temp_out = np.float32(data_out+np.mean(images_aligned, axis=0))
        np.save(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz', temp_out)
        print "...DONE..."

    #Cheby
    elif filter_type == 'chebyshev':

        nyq = fs / 2.0
        order = 4
        rp = 0.1
        Wn = [lowcut / nyq, highcut / nyq]
        b, a = cheby1(order, rp, Wn, 'bandpass', analog=False)
        
        data_out = images_aligned.copy()*0.0
        for row in range(len(images_aligned[0])):
            print "...filtering row: ", row
            for col in range(len(images_aligned[0,row])):
                data_out[:,row,col] = filtfilt(b, a, images_aligned[:,row,col], axis=0)
            
            if plotting: 
                ax.clear(); ax.imshow(np.ma.masked_array(data_out[1000], mask=generic_mask_indexes))
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_title("Lowcut: "+str(lowcut)+"hz, highcut: "+ str(highcut)+ "hz" +#, maxDFF: " + str(round(np.nanmax(data_out[1000]),1)) + " minDFF: " + str(round(np.nanmin(data_out[1000]),1))+ 
                            "\nFrame: 1000 / "+str(len(images_aligned)))
                plt.pause(0.000001) 
                
        print "... saving filtered data...",
        temp_out = np.float32(data_out+np.mean(images_aligned, axis=0))
        np.save(images_file[:-4]+'_'+filter_type+'_'+str(lowcut)+'hz_'+str(highcut)+'hz', temp_out)
        print "...DONE..."
        
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

    #Load aligned/filtered data and find ON/OFF light
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
        
        lever_trace = self.abspositions[lever_position_index-lever_window:lever_position_index+lever_window]

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



def view_static_stm(self):
    
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
        ani.save(filename+'_.mp4', writer=writer)

    plt.show()

    #quit() 
    


def quick_mask(self, data):
        
    generic_mask_file = self.parent.animal.home_dir+self.parent.animal.name + '/genericmask.txt'
    if (os.path.exists(generic_mask_file)==True):
        generic_coords = np.int32(np.loadtxt(generic_mask_file))
    else:
        print "...generic mask not found..."
        return
        #fname = glob.glob(animal.filenames[0].replace('rhd_files/','tif_files/')[:animal.filenames[0].find('rhd_files/')+10]+"*std*")[0]
        #images_temp = np.load(fname)
        #Define_generic_mask_single_frame(images_temp, animal)
        #generic_coords = np.int32(np.loadtxt(generic_mask_file))
                        
    generic_mask_indexes=np.zeros((128,128))
    for i in range(len(generic_coords)):
        generic_mask_indexes[int(generic_coords[i][0])][int(generic_coords[i][1])] = True

    for i in range(int(self.midline_mask.text())):
        generic_mask_indexes[:,64+int(int(self.midline_mask.text())/2)-i]=True

    n_pixels = 128
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
        
    #Mask all frames; NB: PROBABLY FASTER METHOD
    print "... masking data: ", data.shape
    for i in range(0, len(data),1):
        #if len(data)!=128:
        temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes, fill=np.nan)
        #else:
        #    temp_array[i] = np.ma.masked_array(data, mask=generic_mask_indexes)
        #    print "single frame"
        #    return temp_array[0]
    
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
    
    lowband = [0.1, 4]; highband = [8,100]     #Martin's suggestions; Saleem 2010?
    #lowband = [0.1, 4]; highband = [15,100]     #variations
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
        tsf.Siteloc[i*2+1]=i*50 #GUESSING each tetrode is 50um apart
            
    tsf.ec_traces = []
    min_samples = 1E14
    for ctr, file_name in enumerate(ncs_files):
        
        data = loadNcs(file_name)
        tsf.vscale_HP = float(data[2][14].replace('-ADBitVolts ',''))*1E6
        
        tsf.SampleFrequency =float(data[2][12].replace('-SamplingFrequency ' ,''))
        print "SampleFrequency = ", tsf.SampleFrequency,
        
        tsf.n_vd_samples = len(data[0])
        print "... #samples: ", tsf.n_vd_samples
        if tsf.n_vd_samples<min_samples: min_samples = tsf.n_vd_samples

        tsf.ec_traces.append(data[0]) 
    
    #Trunkate extra voltage values (some chs in Neuralynx recs have more/less values than others)
    tsf.n_vd_samples = min_samples 
    for k in range(len(tsf.ec_traces)):
        tsf.ec_traces[k]=tsf.ec_traces[k][:min_samples]
        
    tsf.ec_traces = np.array(tsf.ec_traces, dtype=np.int16)
    
    tsf.ec_traces_hp = tsf.ec_traces
    tsf.ec_traces_lp = tsf.ec_traces.copy()
    
    #******************SAVE HIGH PASS RECORD******************
    #Wavelet filter record first
    tsf.ec_traces = wavelet(tsf.ec_traces)

    print ''; print "...saving alltrack _hp.tsf..."
    file_name = ncs_files[0][:-4]+"_alltrack_hp.tsf"
    save_tsf(tsf, file_name)
    
    #*************SAVE LOW PASS RECORD @ 1KHZ***************
    tsf.ec_traces = tsf.ec_traces_lp-tsf.ec_traces #Subtract wavelet filtered hp record from original record
    
    temp_array=[]
    lowcut = 0.1; highcut=110; fs=1000
    for k in range(len(tsf.ec_traces)):      #Subsample to 1Khz and notch filter
        print "...low pass ch: ", k
        temp = np.array(butter_bandpass_filter(tsf.ec_traces[k][::int(tsf.SampleFrequency/1000)], lowcut, highcut, fs, order = 2), dtype=np.int16)
        notch = np.array(butter_bandpass_filter(temp, 59.9, 60.1, fs, order = 2), dtype=np.int16)
        temp = temp-notch
        temp_array.append(temp)
    
    tsf.ec_traces = np.int16(temp_array) 
    
    #saving low pass record
    tsf.SampleFrequency = 1000
    print "#samples: ", tsf.n_vd_samples
    tsf.n_vd_samples = len(tsf.ec_traces[0])
    print "#samples: ", tsf.n_vd_samples
    
    print ''; print "...saving alltrack _lp.tsf..."
    file_name = ncs_files[0][:-4]+"_alltrack_lp.tsf"
    save_tsf(tsf, file_name)


    #************ SAVE LOW PASS COMPRESSED RECORD *************
    tsf.SampleFrequency = 50000
    print ''; print "...saving alltrack _lp_compress.tsf..."
    file_name = ncs_files[0][:-4]+"_alltrack_lp_compress.tsf"
    save_tsf(tsf, file_name)



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
    save_tsf(tsf, file_name)

    

def concatenate_tsf(self):
    """ Function doc """
    
    print "...concatenate multiple .tsf..."
    
    for k in range(len(self.parent.animal.tsf_files)):  print self.parent.animal.tsf_files[k]
        
    #mouse intan data: select .tsf files directly, load and concatenate them
    if '.tsf' in self.parent.animal.tsf_files[0]:

        for ctr, file_name in enumerate(self.parent.animal.tsf_files):
            if ctr==0:
                #Read first tsf file
                tsf = Tsf_file(file_name)
                tsf.read_ec_traces()
                print len(tsf.ec_traces)
            else:
                #Read additional tsf files
                temp_tsf = Tsf_file(file_name)
                temp_tsf.read_ec_traces()
                
                temp_ec_traces=[]
                print "...updating ch: ",
                for ch in range(len(tsf.ec_traces)):
                    print ch,
                    temp_ec_traces.append(np.append(tsf.ec_traces[ch],temp_tsf.ec_traces[ch]))
                
                tsf.ec_traces=np.int16(temp_ec_traces)
                tsf.n_vd_samples += temp_tsf.n_vd_samples
        
        print ''; print "...saving alltrack .tsf..."
        file_name = self.parent.animal.tsf_files[0][:-4]+"_alltrack.tsf"
        save_tsf(tsf, file_name)
        
        
    #cat data: select outer directory and then find .tsf inside each directory to load and concatenate
    else:
        
        for ctr, dir_name in enumerate(self.parent.animal.tsf_files):
            
            file_name = dir_name+dir_name[dir_name.rfind('/'):]+'.tsf'
            if ctr==0:
                tsf = Tsf_file(file_name)
                tsf.read_ec_traces()
                print len(tsf.ec_traces)
            else:
                temp_tsf = Tsf_file(file_name)
                temp_tsf.read_ec_traces()
                print len(temp_tsf.ec_traces)
                
                temp_ec_traces=[]
                print "...updating chs: ",
                for ch in range(len(tsf.ec_traces)):
                    print ch,
                    temp_ec_traces.append(np.append(tsf.ec_traces[ch],temp_tsf.ec_traces[ch]))
                
                tsf.ec_traces=np.int16(temp_ec_traces)
                tsf.n_vd_samples += temp_tsf.n_vd_samples
        
        print ''; print "...saving alltrack .tsf..."
        file_name = self.parent.animal.tsf_files[0] + "_alltrack.tsf"
        save_tsf(tsf, file_name)    


def concatenate_lfp_zip(self):
    """ Function doc """

    print "...concatenate lfp.zip files..."
    
    for k in range(len(self.parent.animal.tsf_files)):  print self.parent.animal.tsf_files[k]
        
    #mouse intan data: select .tsf files directly, load and concatenate them
    for ctr, dir_name in enumerate(self.parent.animal.tsf_files):
        
        file_name = dir_name+dir_name[dir_name.rfind('/'):]+'.lfp.zip'
        if ctr==0:
            self.parent.animal.load_lfp_all(file_name)
            tsf = self.parent.animal.tsf
            tsf.iformat = '1002'
            tsf.header = 'Test spike file '
            tsf.Siteloc = "LOAD FROM NCHANS" 
            tsf.vscale_HP = 1.0
            tsf.n_electrodes = 10 #OR LOAD FROM .lfp.zip file...
            
            #******************************NB: NEED TO LOAD LENGTH OF RECORDING FROM HIGHPASS FILE AND APPEND ZEROS TO END OF PADDED LFP****
            print "TODO ***************** NEED TO APPEND ZEROS HERE ***********************"
        else:
            self.parent.animal.load_lfp_all(file_name)
            tsf_temp = self.parent.animal.tsf
            
            temp_ec_traces=[]
            print "...updating chs: ",
            for ch in range(len(tsf.ec_traces)):
                print ch,
                temp_ec_traces.append(np.append(tsf.ec_traces[ch],temp_tsf.ec_traces[ch]))
            
            tsf.ec_traces=np.int16(temp_ec_traces)
            tsf.n_vd_samples += temp_tsf.n_vd_samples
    
    print ''; print "...saving alltrack .tsf..."
    file_name = self.parent.animal.tsf_files[0] + "_alltrack.tsf"
    save_tsf(tsf, file_name)    

def compress_lfp(self):
    
    print "...making compressed lfp files ..."

    for file_name in self.filenames:
        
        file_out = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp_compressed.tsf'
        if os.path.exists(file_out)==True: continue

        #Load low-pass .tsf file
        file_in = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'_lp.tsf'
        print "Processing: \n", file_in
        self.load_tsf(file_in)
        print self.tsf.Siteloc.shape

        #Save .lfp.zip file; Martin format; still used by some routines (may wish to remove eventually)
        out_file = file_name[:file_name.find('rhd_files')]+'tsf_files/'+ file_name[file_name.find('rhd_files')+10:]+'.lfp.zip'
        print "Saving LFP : ", out_file
        t0 = 0
        t1 = len(self.tsf.ec_traces[0])*1E6/self.tsf.SampleFrequency  #time of end of file in usec 
        tres = 1000     #1Khz sample rate
        uVperAD = 1.0
        chans = np.arange(0, self.tsf.n_electrodes, 1)
        chanpos = self.tsf.Siteloc

        np.savez(out_file, chans=chans, chanpos=chanpos, data=self.tsf.ec_traces, t0=t0, t1=t1, tres=tres, uVperAD=uVperAD)
        os.rename(out_file+'.npz', out_file)
        
        #COMPRESS LFP 
        #PAD DATA - Required for Nick/Martin recs as LFP starts bit after highpass data;
        #print "Loading LFP (.lfp.zip) file: ", file_lfp
        #data_in  = np.load(file_lfp+'.lfp.zip')
        #t_start = data_in['t0']*1E-6      #Convert to seconds
        #t_end = data_in['t1']*1E-6        #Convert to seconds
        #start_offset = np.zeros((len(data_in['data']), int(t_start*1E3)), dtype=np.int16)  #make padding at front of data
        #data_out = np.concatenate((start_offset, data_out), axis=1)

        #*********** SAVE COMPRESSED LOW PASS .TSF FILE *********
        header = 'Test spike file '
        iformat = 1002
        n_electrodes = self.tsf.n_electrodes
        SampleFrequency = 50000
        vscale_HP = self.tsf.vscale_HP #use the same as above
        #n_vd_samples = len(self.tsf.ec_traces[0][::Compress_factor])
        n_vd_samples = len(self.tsf.ec_traces[0])

        f1 = open(file_out, 'wb')
        f1.write(header)
        f1.write(struct.pack('i', iformat))
        f1.write(struct.pack('i', SampleFrequency))
        f1.write(struct.pack('i', n_electrodes))
        f1.write(struct.pack('i', n_vd_samples))
        f1.write(struct.pack('f', vscale_HP))
        for i in range (n_electrodes):
            f1.write(struct.pack('h', self.tsf.Siteloc[i*2]))
            f1.write(struct.pack('h', self.tsf.Siteloc[i*2+1]))
            f1.write(struct.pack('i', i+1)) #Need to add extra value for Fortran arrays

        print "Writing data"
        for i in range(n_electrodes):
            self.tsf.ec_traces[i].tofile(f1)

        f1.write(struct.pack('i', 0)) #Write # of fake spikes
        #text = "Compression:" + str(overall_compression)
        #n_bytes = len(text)
        #f1.write(struct.pack('i', n_bytes))
        #f1.write(text)
        f1.close()


def Specgram_syncindex(self):
    
    channel=int(self.parent.specgram_ch.text())
    
    if self.parent.exp_type=="mouse":
        self.parent.animal.load_channel(self.parent.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+
                        '_lp.tsf', channel) #Loads single channel as animal.tsf
    elif self.parent.exp_type=="cat":
        temp_name = self.parent.animal.recName.replace('.ptcs','.lfp.zip').replace('.tsf','.lfp.zip')
        self.parent.animal.load_lfp_channel(temp_name, channel) #Loads single channel as animal.tsf
    elif self.parent.exp_type=='rat':
        temp_name = self.parent.animal.recName.replace('_hp.tsf','_lp.tsf')
        self.parent.animal.load_channel(temp_name, channel) #Loads single channel as animal.tsf
        
        #print "exp_type unknown"
        #return
    
    colors=['blue','green','violet','lightseagreen','lightsalmon','dodgerblue','mediumvioletred','indianred','lightsalmon','pink','darkolivegreen']

    tsf = self.parent.animal.tsf
    samp_freq = tsf.SampleFrequency
    print "rec length: ", len(tsf.ec_traces)/float(tsf.SampleFrequency), " sec."

    ax = plt.subplot(1,1,1)
    font_size = 30
    height = 25
    width_plots = 35 #int(max(20, int(math.ceil(len(SUA_sort.sec_len)*3)))*1.16)

    #Compute Specgram
    print "computing specgram..."
    data_in = tsf.ec_traces
    P, extent = Compute_specgram_signal(data_in, samp_freq)
    plt.imshow(P, extent=extent, aspect='auto')

    #Compute sync index
    si_limit=0.7
    si, t, sync_periods = synchrony_index(data_in, samp_freq, si_limit)
    
   # print "t: ", t
    
    plt.plot(t, si*50-60, linewidth=3, color='blue')
    plt.plot([0,max(t)],[-50*.3-10,-50*.3-10], 'r--', color='red', linewidth = 3, alpha=0.8)
    plt.plot([0,max(t)],[-10,-10], color='black', linewidth = 2, alpha=0.8)
    plt.plot([0,max(t)],[-60,-60], color='black', linewidth = 2, alpha=0.8)
    
    xx = np.linspace(0,max(t),5)
    x_label = np.round(np.linspace(0, max(t)/60.,5))
    plt.xticks(xx, x_label, fontsize=20)
       
    ax.set_xlim((0,P.shape[1]/2))    
    ax.set_ylim((-60,110))    
    
    old_ylabel = [-60, -25, -10, 0, 50, 100]
    new_ylabel = [0, 0.7, 1, 0, 50, 100]
    plt.yticks(old_ylabel, new_ylabel, fontsize=font_size)
    
    ax.tick_params(axis='both', which='both', labelsize=font_size)

    plt.ylabel("Synchrony Index             Specgram Frequency (Hz)      ", fontsize=font_size-5)           
    plt.xlabel("Time (mins)", fontsize = font_size)
    plt.title(tsf.file_name, fontsize=font_size-10)
    plt.show()

def Compute_specgram_signal(data, SampleFrequency):

    t0=None
    t1=None
    f0=0.1
    f1=110
    p0=-60
    p1=None
    chanis=-1
    width=2
    tres=.5
    cm=None
    colorbar=False
    title=True
    figsize=(20, 6.5)
    
    F0, F1 = 0.2, 110 # Hz #THESE ARE UNUSED AT THIS TIME;
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
    end_offset = float(indexes[-1])/25000       #NB: compressed data samplerate is 50Khz or other so not reliable.
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


def msl_plots(self):
    #(sim_dir, Sorts_sua, Sorts_lfp, lfp, recs, track_name):
    '''Align msl latencies by lfp cluster 
    '''

    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']
    #colors=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','pink','darkolivegreen','cyan']

    si_limit = 0.7
    window=1.0  #NB: ************* SET THIS VALUE TO ALLOW ARBITRARY ZOOM IN AND OUT ALONG WITH 2000 sized arrays below
    

    lock_window = int(self.lock_window.text())
    
    #Loading length of .tsf: NB: LOADING CH NOT NECESSARY; SHOULD BE ABLE TO READ JUST HEADER
    self.animal.load_channel(self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp.tsf', channel=55) #Loads single channel as animal.tsf
    
    temp_chunks = np.linspace(0,self.animal.tsf.n_vd_samples/float(self.animal.tsf.SampleFrequency), int(self.time_chunks.text())+1)
    time_chunks=[]
    for t in range(len(temp_chunks)-1):
        time_chunks.append([temp_chunks[t],temp_chunks[t+1]])
    
    #UPDATE PARAMS FROM CURRENT WIDGET TEXTBOXES
    self.animal.name = self.animal_name.text()
    self.animal.recName = self.root_dir+self.animal.name+'/rhd_files/'+self.rec_name.text()
        
    #Load SUA Sort
    sua_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_hp.ptcs'
    Sort_sua = Ptcs(sua_file) #Auto load flag for Nick's data
    total_units = len(Sort_sua.units)

    #Load LFP Sort
    lfp_file = self.animal.recName.replace('rhd_files','tsf_files').replace('.rhd','')+'_lp_compressed.ptcs'
    Sort_lfp = Ptcs(lfp_file) #Auto load flag for Nick's data

    start_lfp = min(int(self.start_lfp.text()),len(Sort_lfp.units)-1)
    end_lfp = min(int(self.end_lfp.text()),len(Sort_lfp.units)-1)
    lfp_ctr=0
    for lfp_cluster in range(start_lfp,end_lfp, 1):
        #Pick a particular LFP event to lock to
        fit_sum = np.zeros((total_units,2000), dtype=np.float32)

        print "LFP Cluster # : ", lfp_cluster, " / ", len(Sort_lfp.units),

        #Load pop events during synch periods (secs)
        compress_factor = 50.
        if '2016_07_11' in self.animal.recName: compress_factor=40.  #July 11 recording uses different factors
        pop_spikes = np.array(Sort_lfp.units[lfp_cluster])/Sort_lfp.samplerate*compress_factor#*1E-3  
        print " ... # events: ", len(pop_spikes)

        ##Compute periods of synchrony from si index #**********SKIP FOR MSL PLOT ONLY ********
        #data_in = lfp.data[rec_index][9]
        #si, t, sync_periods = synchrony_index(data_in, lfp, rec_index, si_limit)
        #temp_list = []
        #for p in range(len(sync_periods)):
        #    indexes = np.where(np.logical_and(pop_spikes>=sync_periods[p][0], pop_spikes<=sync_periods[p][1]))[0]
        #    temp_list.extend(pop_spikes[indexes])
        #pop_spikes=np.array(temp_list)

        #CELL RASTER LISTS TO HOLD SPIKES FOR EACH LFP EPOCH
        cell_rasters = []
        for chunk_ctr, time_chunk in enumerate(time_chunks):
            print time_chunk
            #Loop over single units, get LFP triggered rasters, build histograms
            
            #latencies=np.zeros((total_units, 2), dtype=np.float32)

            cell_rasters.append([])            
            for unit in range(len(Sort_sua.units)):
            #for unit in range(10):
                #Load unique track-wide unit id 
                unique_unit = Sort_sua.uid[unit]
                #Load sua spikes during synch periods (secs); use default unit
                spike_array = np.array(Sort_sua.units[unit],dtype=np.float32)/float(Sort_sua.samplerate)  #This converts to seconds

                temp3 = np.where(np.logical_and(spike_array>=time_chunk[0], spike_array<=time_chunk[1]))[0]
                spike_array= spike_array[temp3]
                
                print "...unit: ", unit, " uid: ", unique_unit,  " #spikes: ", len(spike_array), " / ", len(Sort_sua.units[unit])
                

                #SKIP SYNC PERIOD FOR NOW
                #temp_list = []
                #for p in range(len(sync_periods)):
                #    indexes = np.where(np.logical_and(spike_array>=sync_periods[p][0], spike_array<=sync_periods[p][1]))[0]
                #    temp_list.extend(spike_array[indexes])
                #spike_array=np.array(temp_list)

                xx_even=[]          #collect even spikes
                xx_odd=[]           #collect odd spikes
                xx1=[]              #collect all spikes for KS stat test
                for j in range(len(pop_spikes)):
                    #Skip pop spikes that occur w/in 100ms of each other
                    if j<(len(pop_spikes)-1):
                        if (pop_spikes[j+1]-pop_spikes[j])<0.100: continue

                    #find spikes that fall w/in +/- window of pop event
                    temp2 = np.where(np.logical_and(spike_array>=pop_spikes[j]-window, spike_array<=pop_spikes[j]+window))[0]

                    #NB: NOT excluding duplicate spikes from broader window; make sure to take into account for analysis
                    x=(spike_array[temp2]-pop_spikes[j])*1E3 #Offset time of spike to t_0; convert to ms
                    xx1.append(x)
                    if j % 2 == 0:
                        xx_even.append(x)
                    else:
                        xx_odd.append(x)
                
                if len(xx_even)>0: xx_even = np.hstack(xx_even)
                if len(xx_odd)>0: xx_odd = np.hstack(xx_odd)
                if len(xx1)>0: xx1 = np.hstack(xx1)

                #Convolve all spike train data w. 20ms gaussian  - Martin's suggestion; also in Luczak 2007, 2009;
                sig = 20
                fit_even = np.zeros(2000, dtype=np.float32)
                def gaussian(x, mu, sig):
                    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
                if len(xx_even)>0 : 
                    for g in range(len(xx_even)):
                        mu = np.array(xx_even[g]) 
                        fit_even += gaussian(np.linspace(-1E3, 1E3, 2000), mu, sig)

                fit_odd = np.zeros(2000, dtype=np.float32)
                if len(xx_odd)>0 : 
                    for g in range(len(xx_odd)):
                        mu = xx_odd[g] 
                        fit_odd += gaussian(np.linspace(-1E3,1E3, 2000), mu, sig)

                ##Plot control plots
                #if False:
                    #plt.plot(fit_even[960:1040], color='blue')
                    #plt.plot(fit_odd[960:1040], color='red')
                    #plt.plot([40,40],[0,max(np.max(fit_even),np.max(fit_odd))*1.2], 'r--', linewidth=3, color='black')
                    #plt.title(track_name)
                    #plt.show()

                fit_sum[unit] += (fit_even + fit_odd)/2.

                cell_rasters[chunk_ctr].append(xx1)

            #Pick an lfp cluster to order rest of data: 
            img=[]
            lock_time=[]
            for unit in range(total_units):
                if np.max(fit_sum[unit][1000-lock_window:1000+lock_window])>0:
                    lock_time.append(np.argmax(fit_sum[unit][1000-lock_window:1000+lock_window]))
                    img.append(fit_sum[unit][1000-lock_window:1000+lock_window]/max(fit_sum[unit][1000-lock_window:1000+lock_window]))
                else:
                    lock_time.append(200)
                    img.append(np.zeros(lock_window*2, dtype=np.float32))
            
            #ORDER MSL IMG BY LOCK TIME
            if (chunk_ctr ==0): inds = np.array(lock_time).argsort()
            print "Order: ", inds
            #if (chunk_ctr ==0) and (start_lfp==lfp_cluster): inds = np.array(lock_time).argsort()

            img=np.array(img)[inds]
            temp_img = []
            for i in range(len(img)):
                #if np.max(img[i])!=0:
                    temp_img.append(img[i])
            img=np.array(temp_img)
            
            

            #********** PLOTING ROUTINES **************
            ax=plt.subplot(end_lfp-start_lfp,len(time_chunks),lfp_ctr+1)
            im = ax.imshow(img, origin='upper', extent=[0,lock_window*2, len(img),0], aspect='auto', interpolation='none')
            #plt.plot([100,100],[1,total_units-1], 'r--', linewidth=2, color='white')

            ax.set_xticklabels([])
            if lfp_ctr>((end_lfp-start_lfp-1)*(len(time_chunks)-1)):
                xx = np.arange(0,lock_window*2+1,lock_window)
                x_label = np.arange(-lock_window,lock_window+1,lock_window)
                plt.xticks(xx,x_label, fontsize=25)
                plt.xlabel("Time (ms)", fontsize=30)

            yy = np.arange(0,len(img),10)
            y_label = np.arange(0,len(img),10)
            plt.yticks(yy, y_label, fontsize=25)
            
            if (lfp_ctr%len(time_chunks))==0:
                plt.ylabel("lfp: "+str(lfp_cluster)+"\n#"+str(len(pop_spikes)), fontsize=25)

            if (lfp_ctr<len(time_chunks)):
                plt.title(str(int(time_chunk[0]/60.))+".."+str(int(time_chunk[1]/60.))+" mins", fontsize=25)

            plt.ylim(len(inds),0)
            
            #plt.ylabel("Neuron (lock order)", fontsize=30)
            plt.plot([lock_window,lock_window],[0,total_units], 'r--', linewidth=4, color='white')
            ax.tick_params(axis='both', which='both', labelsize=25)
            ax.xaxis.labelpad = 0
            lfp_ctr+=1
    
    self.cell_rasters = cell_rasters
    self.Sort_sua = Sort_sua
    self.t_chunks = time_chunks
    plt.show()

def compute_msl_pvals(self):
    
    colors=['blue','green','cyan','magenta','red','pink','orange', 'brown', 'yellow']


    #COMPUTE P-VALUES FOR CURRENT LFP CLUSTER
    from scipy import stats

    n_spikes = int(self.n_spikes.text()) #Min number of spikes in clusters considered
    cell_rasters = self.cell_rasters
    Sort_sua = self.Sort_sua
    time_chunks = self.t_chunks
    
    p_val_array = []
    for p in range(len(Sort_sua.units)):
        p_val_array.append([])
        if len(Sort_sua.units[p])<n_spikes: continue
        for c in range(len(time_chunks)):
            if (len(cell_rasters[0][p])>0) and (len(cell_rasters[c][p])>0):
                KS, p_val = stats.ks_2samp(cell_rasters[0][p], cell_rasters[c][p])
                p_val_array[p].append(p_val)
    
    for p in range(len(Sort_sua.units)):
        if len(Sort_sua.units[p])<250: continue
        plt.scatter([p]*len(p_val_array[p]), p_val_array[p], color=colors[p%9])
        
    plt.ylim(0,1)
    plt.xlim(-1,len(Sort_sua.units))
    plt.show()

    #quit()

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
        
        for k in range(self.n_units):
            self.readUnit(f)
            self.units[k]= self.spikes

            if 'martin' in self.full_path:
                self.uid[k]= self.nid
            else: #All other sorts are from Nick's SS so should be the same
                self.uid[k]= self.nid-1
               
            #print "SAMPLERATE: ", self.samplerate
            #if ptcs_flag: #Martin's data has wrong flag for saves
            self.units[k]=[x*self.samplerate/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps
            #else:
            #    self.units[k]=[x*self.samplerate/2/1E+6 for x in self.units[k]] #Converts spiketimes from usec to timesteps

            self.n_sorted_spikes[k] = len(self.units[k])
            self.size.append(self.nspikes)
            self.maxchan.append(self.maxchanu)
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
        self.xpos = float(np.fromfile(f, dtype=np.float64, count=1)) # xpos (um)
        self.ypos = float(np.fromfile(f, dtype=np.float64, count=1)) # ypos (um)
        self.zpos = float(np.fromfile(f, dtype=np.float64, count=1)) # zpos (um)
        self.nchans = int(np.fromfile(f, dtype=np.uint64, count=1)) # nchans
        self.chans = np.fromfile(f, dtype=np.uint64, count=self.nchans) #NB: Some errors here from older .ptcs formats
        self.maxchanu = int(np.fromfile(f, dtype=np.uint64, count=1)) # maxchanid

        self.nt = int(np.fromfile(f, dtype=np.uint64, count=1)) # nt: number of time points in template

        self.nwavedatabytes, self.wavedata = self.read_wave(f) #TEMPLATE

        self.nwavestdbytes, self.wavestd = self.read_wave(f) #STANDARD DEVIATION
        self.nspikes = int(np.fromfile(f, dtype=np.uint64, count=1)) # nspikes

        self.spikes = np.fromfile(f, dtype=np.uint64, count=self.nspikes) # spike timestamps (us):

        # convert from unsigned to signed int for calculating intervals:
        self.spikes = np.asarray(self.spikes, dtype=np.float64)

            
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
    

