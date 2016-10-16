#Mouse class and subclasses for processing lever pull data projec Greg Silasi and Tim Murphy lab
#Code: Cat Mitelut

import os
import glob
import numpy as np
import struct
import string, re
import scipy
import tifffile as tiff
import cPickle as pickle
import gc
#from skimage.measure import block_reduce
import matplotlib.pyplot as plt
import shutil


class Mouse_lever(object):      
    
    def __init__(self, mouse_name, home_dir, n_sec):
        
        self.name = mouse_name
        self.home_dir = home_dir
        self.n_sec = n_sec

        #Read sampling frequency from
        img_rate_file = home_dir+mouse_name + '/img_rate.txt'
        self.img_rate = int(np.loadtxt(img_rate_file))
        self.window = int(self.img_rate*n_sec)

        #Fix pixels to 128 for now... 
        self.n_pixels = 128

        #Make lists for trial chunking
        self.traces = []
        self.DFF_files = []
        self.codes = []
        self.trial_times = []
        self.days_counter = []

        print "...img rate: ", self.img_rate, "hz"

        #self.preprocess_mouse_lever()      #***************MAY WISH TO RUN THIS AT LEAST ONCE PER FILE AS PREPROCESSING STEP TO INITIALIZE REALTIME FILES **********


    def preprocess_mouse_lever(self):

        self.load_filenames()   #Loads tif files, event_files, lever_files
        
        self.process_sessions()    #Loads session data

        self.stroke = Stroke(self.sessions, self.home_dir, self.name)      #Load stroke information; need to have processed session info first to place location of stroke

        #self.save_mouse()       #Save mouse file names and trace data to pickle object (NB: no DFF data saved)
        


    def load_traces(self):
        
        counter=0
        for session in self.sessions:
            print "...traces session: ", counter
            session.traces=np.load(session.tif_file[:-4]+"_"+str(self.window/self.img_rate)+"sec_traces.npy")
            counter+=1
            
            
    def load_filenames(self):
        '''Load event files, tiffs and lever position files'''
        
        #LOAD EVENT FILE NAMES
        event_files = os.listdir(self.home_dir + self.name + '/event_files')

        #ORDER EVENT FILES BY DATE
        file_order = []
        for event_file in event_files:
            if '2014' in event_file: 
                temp = event_file.replace(self.name+'_','')
                temp = temp[:temp.find('2014')-1]
                month = temp[:temp.find('-')]
                day = temp[temp.find('-')+1:]
                #print month,day
                file_order.append('2014'+month.zfill(2)+day.zfill(2))
            
            if '2015' in event_file: 
                temp = event_file.replace(self.name+'_','')
                temp = temp[:temp.find('2015')-1]
                month = temp[:temp.find('-')]
                day = temp[temp.find('-')+1:]
                #print month,day
                file_order.append('2015'+month.zfill(2)+day.zfill(2))
            
            if '2016' in event_file: 
                temp = event_file.replace(self.name+'_','')
                temp = temp[:temp.find('2016')-1]
                month = temp[:temp.find('-')]
                day = temp[temp.find('-')+1:]
                #print month,day
                file_order.append('2016'+month.zfill(2)+day.zfill(2))


        indexes= np.argsort(np.int32(file_order)) #Sorted indexes for input files.
        file_order = list(np.array(file_order)[indexes])
        event_files = list(np.array(event_files)[indexes])
            
        #LOAD LEVER FILES NAMES; REMOVE ONES THAT DONT" HAVE CORRESPONDING LEVER POS FILES
        lever_files = []
        del_file_counter = []
        counter=0
        for event_file in event_files:
            temp_index = event_file.find('201')
            temp_filename = event_file.replace('data','')
            load_filename = self.home_dir+self.name+'/leverPos/'+temp_filename[:temp_index+5]+ 'leverPos' + temp_filename[temp_index+5:]

            if (os.path.exists(load_filename)==True):
                lever_files.append(load_filename)
            else:
                print "Missing lever pos file: ", load_filename
                del_file_counter.append(counter)
            counter+=1
            
        #del_file_counter = np.array(del_file_counter)
        for i in reversed(del_file_counter):
            del event_files[i]; del file_order[i]
            
        #LOAD TIF FILE NAMES; need to convert to letter dates;
        tif_files = []
        months = {'1':'Jan', '2':'Feb', '3':'Mar','4':'Apr','5':'May','6':'Jun',
                  '7':'Jul','8':'Aug','9':'Sep','10':'Oct','11':'Nov','12':'Dec'}
        counter=0
        while True:
            event_file=event_files[counter].replace(self.name+'_','')
            month = event_file[:2].replace('-','')
            day = event_file[len(month):len(month)+3].replace('-','')
            #print event_file, "  month: ", month,  "day: ", day
            
            if 'AM' in event_file: am_pm = 'am'
            elif 'PM' in event_file: am_pm = 'pm'
            else: am_pm=''
            
            if self.img_rate == 30: #Newer files use 30Hz and have am/pm and moth/day inverted
                temp_name = self.home_dir+self.name+'/tif_files/'+self.name+am_pm+'_'+months[month]+day
            if self.img_rate == 15:
                temp_name = self.home_dir+self.name+'/tif_files/'+self.name+'_'+months[month]+day+am_pm 
                
            dir_name = glob.glob(temp_name+"_*")   #use .tif extension otherwise will load .npy files

            
            if len(dir_name)==1: #Found tif file
                make_dir = dir_name[0]
                fname = glob.glob(make_dir+"*")[0]
                new_name = make_dir +'/'+ fname.replace(self.home_dir+self.name+'/tif_files/','')+'.tif'

                tif_files.append(new_name)
                counter+=1
                
            elif len(dir_name)>1:
                print dir_name
                print "---------TOO MANY TIF FILES-------"
                quit()
                
            elif len(dir_name)==0:
                print "Missing tif file: ", temp_name+"*"
                del lever_files[counter]
                del event_files[counter]
                del file_order[counter]
                
            if counter==len(event_files): break

        for k in range(len(event_files)):
            event_files[k] = self.home_dir+self.name+'/event_files/'+event_files[k]
        
        self.event_files = event_files
        self.lever_files = lever_files
        self.tif_files = tif_files
        self.file_order = file_order
        
    
        np.save(self.home_dir+self.name+ '/event_files', self.event_files)      #NB: Some of the sessions may still be incorrect; see DFF ca
        np.save(self.home_dir+self.name+ '/lever_files', self.lever_files)
        np.save(self.home_dir+self.name+ '/tif_files', self.tif_files)
        np.save(self.home_dir+self.name+ '/session_realtime', self.file_order)

    def load_reclength(self, filename):
        """ Load realtime length of a single session. Probably should be in session, but was quicker to dump here"""

        text_file = open(filename, "r")
        lines = text_file.read().splitlines()
        event_text = []
        for line in lines:
            event_text.append(re.split(r'\t+',line))
        
        #Delete false starts from event file
        for k in range(len(event_text)-1,-1,-1):        #Search backwards for the 1st occurence of "date" indicating last imaging start
                                                        #NB: There can be multiple false starts WHICH DON"T LINE UP - NEED TO IGNORE SUCH SESSIONS
            if event_text[k][0]=='date': 
                event_text = event_text[k+2:]         #Remove first 2 lines
                break
        
        if len(event_text)==0:
            print "empty session"
            self.reclength = 0
        else:
            if event_text[-1][2] != "None": 
                print "Missing event data end ...skipping session"
                self.reclength = 0
            else: 
                self.reclength = float(event_text[-1][3])
        
        return self.reclength
                
    
    def process_sessions(self):

        self.sessions = []

        #if os.path.exists(home_dir+mouse+"/"+mouse+"_imaging.npy")==False:
        counter=-1
        for tif_file, event_file, lever_file in zip(self.tif_files, self.event_files, self.lever_files):
            print ""
            print "******************************************************"
            print "Session: ", len(self.sessions), " / ", len(self.tif_files)
            print ".......: ", tif_file
            print ".......: ", event_file
            print ".......: ", lever_file
    
            #counter+=1; print counter
            #if counter!=6: continue
            
            session = Session(tif_file, event_file, lever_file, self.window, len(self.sessions),self.home_dir, self.name, self.img_rate)
            session.process_session()
            self.sessions.append(session)
            #if len(self.sessions)>15: return
            
            
    def save_mouse(self):
        print "Saving mouse to disk... NOT CURRENTLY USED"

    def load(self):

        print "Loading mouse from disk..."

        self.load_filenames()

        #Load file names; LOAD ONLY DATA THAT HAS BEEN DOUBLECHECKED FOR ERRORS...
        self.event_files = np.load(self.home_dir+self.name+ '/event_files.npy')
        self.tif_files = np.load(self.home_dir+self.name+ '/tif_files.npy')
        self.lever_files = np.load(self.home_dir+self.name+ '/lever_files.npy')
        self.session_realtimes = np.load(self.home_dir+self.name+ '/session_realtime.npy')
        self.img_rate = int(np.loadtxt(self.home_dir+self.name + '/img_rate.txt'))
        
        #Loop over processed sessions:
        counter=-1
        self.sessions = []
        for tif_file, event_file, lever_file, session_realtime in zip(self.tif_files, self.event_files, self.lever_files, self.session_realtimes):
            print len(self.sessions)
            #if len(self.sessions)>35: break
            
            #Check to see if traces data was generated; NB: the incorrect img_rate files aren't saved anyways, so img_rate check below is redundant
            if (os.path.exists(tif_file[:-4]+"_"+str(int(self.window/self.img_rate))+"sec_traces.npy")==False): 
                print "********** SESSION NOT PROCESSED....SKIPPING********"
                continue
            
            #Check to see if img_rate was ~30.00Hz; otherwise skip
            img_rate = np.load(tif_file[:-4]+'_img_rate.npy') #LOAD IMG_RATE
            if abs(img_rate-float(self.img_rate))>0.01: 
                print "******* img_rate mismatch *************"
                continue
            
            session = Session(tif_file, event_file, lever_file, self.window, len(self.sessions),self.home_dir, self.name, self.img_rate)

            session.abstimes = np.load(tif_file[:-4]+'_abstimes.npy')
            session.abspositions = np.load(tif_file[:-4]+'_abspositions.npy')
            session.abscodes = np.load(tif_file[:-4]+'_abscodes.npy')
            session.locs_44threshold = np.load(tif_file[:-4]+'_locs44threshold.npy')
            session.code_44threshold = np.load(tif_file[:-4]+'_code44threshold.npy')
            
            #Load traces
            session.traces=np.load(tif_file[:-4]+"_"+str(int(self.window/self.img_rate))+"sec_traces.npy")
            
            #Generate absolute time for each trace:
            session.realtimes= [session_realtime]*len(session.traces)
            
            #Save DFF_File names to list; IS THIS REDUNDANT? Isn't self.event_files already containing this info
            session.DFF_files=[]
            for p in range(len(session.traces)):
                session.DFF_files.append(tif_file[:-4]+'_'+str(p).zfill(4)+'.npy')

            self.sessions.append(session)

        print len(self.sessions)
        print ''
        self.stroke = Stroke(self.sessions, self.home_dir, self.name)      #Load stroke information; need to have processed session info first to place location of stroke

    def move_tifs(self):
        import shutil
        
        #THIS WORKS AS STAND ALONE... MODIFY TO WORK WITH OBJECT
        home_dir = '/media/cat/12TB/in_vivo/tim/yuki/I1/tif_files/'
        fnames = glob.glob(home_dir+"*.tif")   #use .tif extension otherwise will load .npy files

        for fname in fnames: 
            if '.tif' in fname: 
                dir_name = fname.replace('.tif','/')
                os.makedirs(dir_name)

                new_name = dir_name + fname.replace(home_dir,'')
                print fname
                print new_name
                shutil.move(fname, new_name)    

    #def save_DFF(self):
        #''' NOT USED CURRENTLY '''
        
        #print "Loading all DFF - saving to single file (unit8)"
        
        #DFF_array = []
        #if (os.path.exists(self.home_dir+self.name+'/'+self.name+'_DFF.npy')==False):
            #for session in self.sessions:
                #print session.tif_file
                #DFF_array.extend(np.load(session.tif_file[:-4]+'_3sec.npy'))
        
            ##Convert list to np array
            #DFF_array = np.array(DFF_array)
            #print DFF_array.shape

            #print "Writing DFF_array to disk..."
            #np.save(self.home_dir+self.name+'/'+self.name+'_DFF.npy', DFF_array)

            ##CODE FOR CONVERTING TO 64pixel by averaging... not great, but may be ok for debugging.
            ##DFF_array_64 = block_reduce(DFF_array, block_size=(1,1,2,2), func=np.mean)
            ##np.save(self.home_dir+self.name+'/'+self.name+'_DFF_64.npy', DFF_array_64)
     
        #else:
            #print "DFF already saved..."
            
            #DFF_array = np.load(self.home_dir+self.name+'/'+self.name+'_DFF.npy')

        #print DFF_array.shape
        #quit()
        

    def chunk_data(self,pre_stroke_days, post_stroke_days, post_post_stroke_days):
        
        print pre_stroke_days, post_stroke_days, post_post_stroke_days
        return
        
        all_traces=[]
        all_codes=[]
        all_DFF_files = []
        all_realtimes = []
        for session in self.sessions:
            all_traces.extend(session.traces)
            all_codes.extend(session.code_44threshold)
            all_DFF_files.extend(session.DFF_files)
            all_realtimes.extend(session.realtimes)

        #Mark each trial w. relative time;
        self.n_trials = len(all_traces)
        trial_times = np.linspace(0, 1, self.n_trials)     #Distribute list of trials uniformly
        all_traces = np.abs(all_traces)

        #Define first and last sessions
        self.realtimes_first = all_realtimes[0]
        self.realtimes_last = all_realtimes[-1]

        #Convert all_realtimes from string to sequential day activity;
        self.stroke.day=0
        days_counter = []
        days_counter.append(0)
        self.month_counter = []
        self.month_text = []
        for k in range(1,len(all_realtimes),1):
            if all_realtimes[k][4:6] == all_realtimes[k-1][4:6]:  #Same month; just subtract days
                days_counter.append(days_counter[k-1]+int(all_realtimes[k][6:8]) - int(all_realtimes[k-1][6:8]))
            else:   #Following month; add remaining days in previous month; plus days from beginning of month; asume 31 days/month
                temp1 = 31 - int(all_realtimes[k-1][6:8])
                days_counter.append(days_counter[k-1]+temp1 + int(all_realtimes[k][6:8]))
                
                self.month_counter.append(days_counter[k-1]+temp1 + int(all_realtimes[k][6:8]))
                self.month_text.append(all_realtimes[k])
            
            if (self.stroke.trial<=k) and (self.stroke.day==0): 
                self.stroke.day = days_counter[k]-3

        print len(days_counter), len(all_traces)
        self.n_days = days_counter[-1]

        #Define data chunks
        epochs = []
        epochs.append(np.where(np.logical_and(np.array(days_counter)>=(self.stroke.day-pre_stroke_days), np.array(days_counter)<(self.stroke.day)))[0])
        epochs.append(np.where(np.logical_and(np.array(days_counter)>=(self.stroke.day), np.array(days_counter)<(self.stroke.day+post_stroke_days)))[0])
        epochs.append(np.where(np.logical_and(np.array(days_counter)>=(self.stroke.day+post_stroke_days), np.array(days_counter)<(self.stroke.day+post_post_stroke_days)))[0])

        #Load 
        for epoch in epochs:                  
            print epoch
            if len(epoch)==0: continue
            start_trial = epoch[0]
            end_trial = epoch[-1]

            chunk_traces = np.array(all_traces[start_trial:end_trial])
            chunk_DFF_files = np.array(all_DFF_files[start_trial:end_trial])
            chunk_codes = np.array(all_codes[start_trial:end_trial])
            chunk_days_counter = np.array(days_counter[start_trial:end_trial])
            chunk_times = trial_times[start_trial:end_trial]

            codes_02 = np.where(chunk_codes=='02')[0]
            codes_04 = np.where(chunk_codes=='04')[0]
            codes_07 = np.where(chunk_codes=='07')[0]
            print "# 04 trials: ", len(codes_04)

            #Select only traces and DFF files for code
            self.traces.append(chunk_traces[codes_04])
            self.DFF_files.append(chunk_DFF_files[codes_04])
            self.trial_times.append(chunk_times[codes_04])
            self.days_counter.append(chunk_days_counter[codes_04])


class Session(object):
    
    def __init__(self, tif_file, event_file, lever_file, window, index, home_dir, name, img_rate):

        self.home_dir = home_dir
        self.name = name
        self.tif_file = tif_file
        self.event_file = event_file
        self.lever_file = lever_file
        self.window = window
        self.index = index                  #Save session index for access inside session object
        self.aligned_images = []            #Make empty list - NOT SURE IF REQUIRED, POSSIBLY EASIER WAY TO DO THIS
        self.img_rate = img_rate
        
    def process_session(self):

        self.load_lever()                   #Load lever positions

        if self.reclength == 0: return
    
        self.convert_tif()                  #Convert .tif -> .npy format; save to disk

        self.align_images()                 #Align raw_images to first session frame 1000

        #REMOVED OUT OF PRE-PROCESSING
        #self.compute_DFF()                  #Compute DFF on aligned images; ################THIS IS OLD VERSION OF DFF PROCESSING, SAVED INDIVIDUALLY; DO NOT USE

        #self.empty_session()                #Remove data already saved to disk; otherwise object too large
        
        #self.make_trials()                  #Make and populate individual reward trials for each session

    def load_lever(self):

        print " ... loading lever positions..."
        #Make lever object
        self.trials = []       #List for each trial 
        self.positions = []    #All lever positions 
        self.threshold = []    #Holds '44' values which indicates threshold crossing; '77' otherwise
        self.times = []        #Corresponding times for each lever position
        self.code = []         #Code for trial
        
        #****************LOAD EVENT FILE****************
        text_file = open(self.event_file, "r")
        lines = text_file.read().splitlines()
        event_text = []
        for line in lines:
            event_text.append(re.split(r'\t+',line))
        
        #Delete false starts from event file
        for k in range(len(event_text)-1,-1,-1):        #Search backwards for the 1st occurence of "date" indicating last imaging start
                                                        #NB: There can be multiple false starts WHICH DON"T LINE UP - NEED TO IGNORE SUCH SESSIONS
            if event_text[k][0]=='date': 
                event_text = event_text[k+2:]         #Remove first 2 lines
                break
        
        if len(event_text)==0:
            print "empty session"
            self.reclength = 0
            return 
        else:
            if event_text[-1][2] != "None": 
                print "Missing event data end ...skipping session"
                self.reclength = 0
                return
                
            else: 
                self.reclength = float(event_text[-1][3])


        #**************** LOAD LEVER POSTIION FILE **********
        text_file = open(self.lever_file, "r")
        lines = text_file.read().splitlines()
        lever_text = []
        for line in lines:
            lever_text.append(re.split(r'\t+',line))
        lever_text = np.array(lever_text)

        #Delete false starts from lever position file
        for k in range(0,len(lever_text),1):     #Search for 1st occurence of first correct value in event file
            if lever_text[k][0]==event_text[0][0]: 
                lever_text = lever_text[k:]        
                break
        
        #Convert lever_text array into data for the lever object
        trial_counter=-1
        lever_counter=0
        event_counter=-1
        for event in event_text[:-1]:       #SHITY CODE USES 3 COUNTERS... NEED TO REDO THESE LOOPS
            event_counter+=1
            if lever_counter==len(lever_text): break
            if (event[0]==lever_text[lever_counter][0]) and (event[1]==lever_text[lever_counter][2]):
                trial_counter+=1
                #Save code from file
                self.code.append(lever_text[lever_counter][2])
                
                #Save relative time
                self.trials.append(float(event_text[event_counter][3]))

                #Save threshold value, positions and times
                self.threshold.append([])
                self.positions.append([])
                self.times.append([])

                while True: 
                    lever_counter+=1
                    if lever_counter==len(lever_text): break
                    #if (lever_text[k][0][:5]=='2016-') or (lever_text[k][0][:5]=='2015-'):   #CAREFUL THIS MAY BE FOOLED BY OTHER  DATES
                    if (lever_text[lever_counter][1]=='None'):   break 
                    
                    #Save lever positions
                    self.times[trial_counter].append(float(lever_text[lever_counter][0]))       #Time values from leverpos file; e.g. 2.91, 2.99, 3.07...5.02
                    self.positions[trial_counter].append(float(lever_text[lever_counter][1]))   #Position vals from levpos file; e.g. 0, 0, 1, 11, 34, 44, ...
                    self.threshold[trial_counter].append(float(lever_text[lever_counter][2]))   #Threshold vals from levps file; e.g. 77, 77, 77, 44, 77, 77...
          
       
        #*************** GENERATE REALTIME LEVER LOCATIONS AND POSITIONS *******************
        ''' NB: This data is saved in chunks and not in continuous time;  i.e. sequential indexes of lever position data (abspositions) 
            not temporally continuous; values (abstimes) needed to be reliably used to match imaging data - so we interpolate and match all times
            TODO: Convert data to complete stream
        '''
        plotting=True
        self.abstimes = []; self.abstimes.extend([0])                    #Generate continuous time stream array
        self.abspositions = []; self.abspositions.extend([0])               #Absolute positions
        self.abscodes = []; self.abscodes.extend([0])                       #Trial code for each recorded position
        self.locs_44threshold = []                                            #Successful codes: '44' value;
        self.code_44threshold = []                                          #Double save of code value
        times = [0]     #Set default in case not successful trials
        
        #Loop over each trial's threshold markers; i.e. 77.0, 77.0, 44.0, ... values to look for threshold code '44' 
        for k in range(len(self.threshold)):
            if 44 in self.threshold[k]:  #Search trial for 44 value to match to lever.trials times
                #FIND ABSOLUTE TIME LOCATION OF TIMES BY OFFSETTING USING EVENT CODE FILE DTSTART
                thresh_index = self.threshold[k].index(44)  #Find location of the 44 in order to align times and leverpos to this value
                times = np.array(self.times[k]) - self.times[k][thresh_index] + self.trials[k]      #Save absolute times to nearest milisecond
                
                #INSERT MISSING TIMES
                missing_times = np.arange(self.abstimes[-1], times[0], 1./120.) #INTERPLOATE MISSING DATA AT @120HZ
                self.abstimes.extend(missing_times)
                self.abstimes.extend(times)
                
                #Save absolute lever positions
                self.abspositions.extend(np.zeros(len(missing_times), dtype=np.float32))
                self.abspositions.extend(self.positions[k])

                #Save code for each position: 02, 04, 07...
                self.abscodes.extend(np.zeros(len(missing_times), dtype=np.int16))
                self.abscodes.extend([self.code[k]]*(len(self.positions[k])))

                #Save threshold location in seconds; and value; BIT REDUNDANT AS ALREADY STORED IN abscodes  
                temp_44 = self.trials[k]
                self.locs_44threshold.append(temp_44)
                self.code_44threshold.append(self.code[k])

        #FILL IN LEVER POS BETWEEN LAST LEVER DATA AND END OF RECORDING
        missing_times = np.arange(times[-1],self.reclength, 1./120.)
        self.abstimes.extend(missing_times)
        self.abstimes = np.array(self.abstimes)

        missing_positions = np.zeros(len(missing_times), dtype=np.int16)
        self.abspositions.extend(missing_positions)
        self.abspositions = np.array(self.abspositions)
        
        missing_codes = np.zeros(len(missing_times), dtype=np.int16)
        self.abscodes.extend(missing_codes)
        self.abscodes = np.array(self.abscodes)

        plotting = False
        if plotting:
            print "reclength: ", self.reclength, len(self.abstimes), len(self.abspositions), len(self.abscodes)

            for p in range(len(self.locs_44threshold)):
                plt.plot([self.locs_44threshold[p],self.locs_44threshold[p]], [0,90], 'r--', color='black', alpha=.5, linewidth=3)

            plt.plot(self.abstimes, self.abspositions, color='blue', linewidth = 3)
            plt.plot(self.abstimes, self.abscodes, color='red', linewidth = 2)

            plt.plot([self.abstimes[0],self.abstimes[-1]], [2,2], 'r--', color='red', alpha=.5, linewidth=3)
            plt.plot([self.abstimes[0],self.abstimes[-1]], [4,4], 'r--', color='blue', alpha=.5, linewidth=3)
            plt.plot([self.abstimes[0],self.abstimes[-1]], [7,7], 'r--', color='black', alpha=.5, linewidth=3)
            
            plt.title(self.tif_file, fontsize=20)
            plt.ylim(0,80)
            plt.show()
        
        #THIS LINE GIVES WARNINGS FOR EMPTY ARRAYS
        print "Event Codes:  # 02: ", len(np.where(np.array(self.code)=='02')[0]), ",   # 04: ", len(np.where(np.array(self.code)=='04')[0]), \
        ",   # 07: ", len(np.where(np.array(self.code)=='07')[0]), ",   # 08: ", len(np.where(np.array(self.code)=='08')[0]), \
        ",   # 00: ", len(np.where(np.array(self.code)=='00')[0])

        self.trigger_indexes = self.locs_44threshold
    
        #SAVE DATA TO FILE:
        np.save(self.tif_file[:-4]+'_abstimes', self.abstimes)
        np.save(self.tif_file[:-4]+'_abspositions', self.abspositions)
        np.save(self.tif_file[:-4]+'_abscodes', self.abscodes)
        np.save(self.tif_file[:-4]+'_locs44threshold', self.locs_44threshold)  #***** NEED TO SAVE LOCS and CODES FOR ANALYSIS BELOW ****** 
        np.save(self.tif_file[:-4]+'_code44threshold', self.code_44threshold)  #Analysis of DFF will automatically exclude triggers that are out of bounds.

    
    def convert_tif(self):
    
        if (os.path.exists(self.tif_file[:-4] +'.npy')==False) and (os.path.exists(self.tif_file[:-4] +'_aligned.npy')==False):
            print "...read: ", self.tif_file
            images_raw = tiff.imread(self.tif_file)

            print "... saving .npy"
            np.save(self.tif_file[:-4], images_raw)
        else:
            print "... *_aligned.npy file exists...skipping conversion..."


    def compute_DFF(self):

        #IF TRACES SAVED, DFF already computed; skip
        if os.path.exists(self.tif_file[:-4]+"_"+str(int(self.window/self.img_rate))+"sec_traces.npy"): 
            print "3sec DFF already computed ...skiping load..."
            #self.traces=np.load(self.tif_file[:-4]+"_"+str(self.window/30)+"sec_traces.npy")
            return

        #Load rotated imaging data
        if len(self.aligned_images)==0: 
            print "...loading aligned img data..."
            self.aligned_images = np.load(self.tif_file[:-4]+"_aligned.npy")
        
        #FIND BEGINNING AND END OF IMAGING
        blue_light_threshold = 400
        #Case #1: imaging starts with light on; only need to remove end
        if np.average(self.aligned_images[0])> blue_light_threshold:
            for k in range(len(self.aligned_images)):
                if np.average(self.aligned_images[k])< blue_light_threshold:
                    self.aligned_images = self.aligned_images[k:]
                    break
        #Case #2: imaging starts with light off; need to remove starting and end chunk;
        else:
            #Find first light on
            for k in range(len(self.aligned_images)):
                if np.average(self.aligned_images[k])> blue_light_threshold:
                    self.aligned_images = self.aligned_images[k:]
                    break

            #Find light off - count backwards from end of imaging data
            for k in range(len(self.aligned_images)-1,0,-1):
                if np.average(self.aligned_images[k])> blue_light_threshold:
                    self.aligned_images = self.aligned_images[:k]
                    #print k
                    break

        self.n_images=len(self.aligned_images)

        #CHECK img_rate 
        session_img_rate = self.n_images/self.reclength
        print "# img frames: ", self.n_images, " rec length: ", self.reclength, " img_rate: ", session_img_rate
        
        if abs(session_img_rate-float(self.img_rate))<0.01:         #Compare computed session img_rate w. experimentally set img_rate
            print "Correct img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
            np.save(self.tif_file[:-4]+'_img_rate', session_img_rate)
        else:
            print "***Incorrect img rate: ", session_img_rate, ",  # img frames: ", self.n_images, ",  rec length: ", self.reclength
            np.save(self.tif_file[:-4]+'_img_rate', session_img_rate)
            return

        #COMPUTE absolute times of triggers from lever pull data
        trigger_times = self.locs_44threshold

        #FIND NEAREST PREVIOUS IMAGING FRAME TIME FOR EACH LEVER POS THRESHOLD
        frame_times = np.linspace(0, self.reclength, self.n_images)             #Set time for each frame in imaging data
        img_frame_triggers = []
        for i in range(len(trigger_times)):
            #img_frame_triggers.append(self.find_previous(frame_times, trigger_times[i]))
            img_frame_triggers.append(self.find_nearest(frame_times, trigger_times[i]))

        
        #BASELINE FOR GLOBAL BASLINE REMOVAL
        global_DFF_filename = self.tif_file[:-4]+"_global_DFF"
        if os.path.exists(global_DFF_filename+'.npy'): 
            global_DFF = np.load(global_DFF_filename +'.npy')
        else:
            global_DFF = np.average(self.aligned_images[1000:-1000], axis=0)
            np.save(global_DFF_filename, global_DFF)


        print "...computing DF/F..."
        data_3sec = []; data_global=[]; traces = []; locs = []; codes = []
        counter=-1
        plotting=False
        print "Trigger frame: ", 
        for trigger in img_frame_triggers:
            counter+=1
            print trigger,
            #NB: Skip first 100 frames or so as the camera is not on and DF/F averaging will not work
            if trigger < (2*self.window) or trigger>(self.n_images-self.window): 
                #print "..skip: too close to start/end ..."
                continue  #Skip if too close to start/end

            #add locs and codes
            locs.append(self.locs_44threshold[counter])
            codes.append(self.code_44threshold[counter])

            #***PROCESS IMAGING; load data before and after pull
            data_chunk = self.aligned_images[trigger-self.window:trigger+self.window]
            
            #DF/F computation; uses baseline -2*window .. -window; e.g. for 3sec windows, baseline: average(-6sec..-3sec)
            baseline = np.average(self.aligned_images[trigger-2*self.window:trigger-self.window], axis=0)
            data_3sec.append((data_chunk-baseline)/baseline)
            data_global.append((data_chunk-global_DFF)/global_DFF)
            
             
            #***PROCESS TRACES - WORKING IN DIFFERENT TIME SCALE *************
            lever_window = 120*3    #NB: Lever window is computing in real time steps @120Hz; and the trigger_indexes are discontinuous;
            t = np.linspace(-lever_window*0.0082,lever_window*0.0082, lever_window*2)
            lever_position_index = self.find_nearest(np.array(self.abstimes), self.locs_44threshold[counter])
            
            lever_trace = self.abspositions[lever_position_index-lever_window:lever_position_index+lever_window]

            if len(lever_trace)!=len(t):    #Extraplote missing data
                print "...missing lever trace data ... extrapolating..."
                lever_trace = np.zeros(lever_window*2,dtype=np.float32)
                for k in range(-lever_window,lever_window,1):
                    lever_trace[k+lever_window] = self.abspositions[k+lever_window]     #Double check this...

            traces.append(lever_trace)

            plotting=False
            if plotting: self.plot_traces()

        print '\n'
        #Save traces, and 44 threshold locations and codes for trials within boundaries
        np.save(self.tif_file[:-4]+'_'+str(int(self.window/self.img_rate))+"sec_traces", traces)
        np.save(self.tif_file[:-4]+'_locs44threshold', locs)
        np.save(self.tif_file[:-4]+'_code44threshold', codes)

        #Save individual trial time dynamics
        data_3sec = np.float16(data_3sec)
        data_global = np.float16(data_global)
        print "Saving trial DFF...",
        for k in range(len(data_3sec)):
            print k,
            np.save(self.tif_file[:-4]+'_3sec_'+str(k).zfill(4), data_3sec[k])
            np.save(self.tif_file[:-4]+'_global_'+str(k).zfill(4), data_global[k])

        print ''
        #Set aligned_images to empty 
        self.aligned_images = []

    def plot_traces(self):
        import matplotlib.gridspec as gridspec
        generic_mask_indexes = load_generic_mask(self.home_dir+self.name, data_chunk)
        
        gs = gridspec.GridSpec(10,1)
        ax = plt.subplot(gs[:2,:])

        for p in range(len(self.locs_44threshold)):
            plt.plot([self.locs_44threshold[p],self.locs_44threshold[p]], [0,90], 'r--', color='black', alpha=.5, linewidth=3)

        plt.plot(self.abstimes, self.abspositions, color='blue', linewidth = 3)
        plt.plot(self.abstimes, self.abscodes, color='red', linewidth = 2)

        plt.plot([self.abstimes[0],self.abstimes[-1]], [2,2], 'r--', color='red', alpha=.5, linewidth=3)
        plt.plot([self.abstimes[0],self.abstimes[-1]], [4,4], 'r--', color='blue', alpha=.5, linewidth=3)
        plt.plot([self.abstimes[0],self.abstimes[-1]], [7,7], 'r--', color='black', alpha=.5, linewidth=3)
        
        plt.title(self.tif_file, fontsize=20)
        plt.ylim(0,80)
        
        
        ax = plt.subplot(gs[2:4,:])
        t = np.linspace(0,len(lever_trace),len(lever_trace))/120.-3.0
        plt.plot(t, lever_trace)
        
        block=6
        ax = plt.subplot(gs[4:6,:]) 
        data_temp = data_chunk
        data_temp = fast_mask(data_temp, generic_mask_indexes)

        v_max = np.max(np.abs(data_temp))
        v_min = -v_max
        plot_array = []
        for k in range(0,len(data_temp), block):                          
            plot_array.append(np.ma.average(data_temp[k:k+block], axis=0))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.imshow(np.ma.hstack((plot_array)), vmin=v_min, vmax=v_max)
        plt.title("Raw Data", fontsize=12)
        
        
        ax = plt.subplot(gs[6:8,:]) 
        data_temp = (data_chunk-baseline)/baseline
        data_temp = fast_mask(data_temp, generic_mask_indexes)

        v_max = np.max(np.abs(data_temp))
        v_min = -v_max
        plot_array = []
        for k in range(0,len(data_temp), block):                          
            plot_array.append(np.ma.average(data_temp[k:k+block], axis=0))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.imshow(np.ma.hstack((plot_array)), vmin=v_min, vmax=v_max)
        plt.title("DFF: 3sec baseline ", fontsize=12)
                        
        ax=plt.subplot(gs[8:10,:])
        data_temp = (data_chunk-global_DFF)/global_DFF
        data_temp = fast_mask(data_temp, generic_mask_indexes)
        v_max = np.max(np.abs(data_temp))
        v_min = -v_max
        
        plot_array = []
        for k in range(0,len(data_temp), block):                          
            plot_array.append(np.ma.average(data_temp[k:k+block], axis=0))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.imshow(np.ma.hstack((plot_array)), vmin=v_min, vmax=v_max)
        plt.title("DFF: average signal removal ", fontsize=12)

        plt.show()


    def align_images(self):

        if os.path.exists(self.tif_file[:-4]+"_aligned.npy")==False:
            
            raw_images=np.load(self.tif_file[:-4]+'.npy')
            np.save(self.tif_file[:-4]+'_frame1000', raw_images[1000])
            
            #Save 1st session frame 1000 to be used by all other code
            print self.index
            if self.index==0: np.save(self.home_dir+self.name+'/'+self.name+'_align_frame', raw_images[1000])
            
            #Load 1st session frame1000:
            first_img = np.load(self.home_dir+self.name+'/'+self.name+"_align_frame.npy")

            #Load current session frame1000
            current_img = raw_images[1000]
            
            #Find translation between frame1000 current session and 1st session
            #im2, scale, angle, t = similarity(im0, plot_img)
            t = translation(first_img, current_img)
            
            #Load raw .npy imaging data
            aligned =np.load(self.tif_file[:-4]+".npy")
            for k in range(len(aligned)):
                aligned[k] = np.roll(aligned[k], t[0], axis=0)
                aligned[k] = np.roll(aligned[k], t[1], axis=1)
            
            #Save rotated imaging data
            np.save(self.tif_file[:-4]+"_aligned", aligned)
            self.aligned_images = aligned

    def empty_session(self):
        
        self.locs_44threshold=[]

        self.abstimes = []

        self.abspositions = []
        
        self.abscodes = []
        

    def make_trials(self):
        ''' Saves rewarded events from each session into individual trial objects for later analysis
            NB: Maybe more efficient to skip this object
        '''
        
        print "...making trials..."
        
        self.trials = []
        for k in range(len(self.traces)):
            trial = Trial(self.traces[k])
            self.trials.append(trial)

        return temp
            
    def find_nearest(self, array, value):
        return (np.abs(array-value)).argmin()

    def find_previous(self, array,value):
        temp = (np.abs(array-value)).argmin()
        if array[temp]>value: return temp-1
        else: return temp
        
class Trial(object):
    ''' Make individual trial swithin each session representing individual DF/F stacks and lever position traces'''
    def __init__(self, traces):

        #self.DFF = DFF
        self.traces = traces


class Stroke(object):
    ''' Make individual trial swithin each session representing individual DF/F stacks and lever position traces'''
    
    def __init__(self, sessions, home_dir, name):

        self.stroke_time(sessions, home_dir, name)
        
        print name
        if name in ['IA1', 'IA2', 'IA3', 'IJ1', 'IJ2']: 
            self.stroke_mask(sessions, home_dir, name)
        else:
            #NO stroke-day single img files
            contour_save = [0,0]
            whole_save = [0,0]

            np.save(sessions[0].home_dir+sessions[0].name+'/stroke_contour', contour_save)
            np.save(sessions[0].home_dir+sessions[0].name+'/stroke_mask', whole_save)

    def stroke_time(self, sessions, home_dir, name):
        
        stroke_data = np.loadtxt(home_dir+name+'/stroke.txt', dtype=str)
        
        if len(stroke_data)==3: 
            self.ampm = stroke_data[1]
            self.kind = stroke_data[2]
            self.day = stroke_data[0][3:5]
            self.month = stroke_data[0][0:2]
            self.year = '20'+stroke_data[0][6:8]

        else:
            print "-----Stroke file problem"
            quit()
        
        self.fulldate = self.year+self.month+self.day

        #DEFAULT VALUES FOR TESTING
        self.trial = 0
        self.session = 0

        #Find relative location of stroke in rewarded trials for later referencing/indexing
        trial_counter = 0
        session_counter = 0
        for session in sessions:
            date = session.event_file.replace(home_dir,'').replace(name+'/event_files/','').replace('.txt','').replace(name+'_','')
            month = date[0:date.index('-')];    date = date[date.index('-')+1:]
            day = date[0:date.index('-')];    date = date[date.index('-')+1:]
            year = date[0:date.index('_')]; 
            
            if year == self.year:
                if (int(month) >= int(self.month)) and (int(day)>=int(self.day)):
                    #print self.month, self.day, self.year
                    #print session.event_file
                    self.trial = trial_counter
                    self.session = session_counter
                    print "Stroke trial: ", self.trial, " session: ", self.session
                    break

            trial_counter+=len(session.locs_44threshold)
            session_counter+=1
    
       
    def stroke_mask(self, sessions, home_dir, name):
        
        print sessions[0].home_dir+sessions[0].name+'/tif_files/*stroke*'
        file_name = glob.glob(sessions[0].home_dir+sessions[0].name+'/tif_files/*stroke*')

        if len(file_name)==0:
            print "MISSING STROKE IMG FOR SHAM MOUSE"
            
        else:
            file_name = file_name[0]
            data = tiff.imread(file_name)
            data = np.rot90(data)

            plotting = False
            #if os.path.exists(sessions[0].home_dir+sessions[0].name+'/stroke_mask.npy')==False:
            if True:
                #Load 1st session frame1000:
                first_img = np.load(sessions[0].home_dir+sessions[0].name+'/'+sessions[0].name+"_align_frame.npy")
                first_img = first_img - ndimage.gaussian_filter(first_img, sigma=2)
                #first_img[:,:64]=0
                
                #Load current session frame1000
                current_img = data
                current_img = data - ndimage.gaussian_filter(current_img, sigma=2)
                #current_img[:,:64]=0
                
                #Find translation between frame1000 current session and 1st session
                #im2, scale, angle, t = similarity(im0, plot_img)
                t = translation(first_img, current_img)
                
                #Load raw .npy imaging data
                aligned = np.roll(current_img, t[0], axis=0)
                aligned = np.roll(aligned, t[1], axis=1)

                data_out = np.roll(data, t[0], axis=0)
                data_out = np.roll(data_out, t[1], axis=1)
            
                if plotting: 
                    ax1 = plt.subplot(321)
                    plt.imshow(data)
                    
                    ax1 = plt.subplot(322)
                    plt.imshow(aligned)
                    
                    ax1 = plt.subplot(323)
                    plt.imshow(first_img)

                    ax1 = plt.subplot(324)
                    plt.imshow(data_out)

            print np.max(data_out), np.min(data_out)
            temp_plot = []
            data_out = ndimage.gaussian_filter(data_out, sigma=2)
            
            data_1d = data_out.reshape(data_out.shape[0]*data_out.shape[1])
            val999 = np.percentile(data_1d, 97.5)               #Mark stroke as the 97.5 percentile and higher values; 
            thresh = val999                                     #Due to rotation of data.
            #thresh = np.max(data_out)*0.9                          #Should theoretically be the max val; but some smoothing likely occuring
            thresh_delta = thresh*.175

            sum_temp0 = []
            sum_temp1 = []
            contour_save = []
            whole_save = []
            for k in range(len(data_out)):
                temp0 = np.where(np.logical_and(data_out[k]>=thresh, data_out[k]<=thresh+thresh_delta))[0]  #Use 1.0 for coverage maps
                contour_save.append(temp0)
                line = np.zeros(128, dtype=np.int32)
                line[temp0]=1.0
                sum_temp0.append(line)
                
                temp0 = np.where(data_out[k]>=thresh)[0]  #Use 1.0 for coverage maps
                whole_save.append(temp0)

                line = np.zeros(128, dtype=np.int32)
                line[temp0]=1.0
                sum_temp1.append(line)
            
            if plotting: 
                sum_temp0 = np.vstack((sum_temp0))
                ax1 = plt.subplot(325)
                plt.imshow(sum_temp0)

                sum_temp1 = np.vstack((sum_temp1))
                ax1 = plt.subplot(326)
                plt.imshow(sum_temp1)
                plt.show()
            
            np.save(sessions[0].home_dir+sessions[0].name+'/stroke_contour', contour_save)
            np.save(sessions[0].home_dir+sessions[0].name+'/stroke_mask', whole_save)
        
        #quit()




