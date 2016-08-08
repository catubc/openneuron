import os
import glob
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import struct
import string
import scipy
import time

from threading import *
import time

import sklearn  
from sklearn import manifold
#from tsne import bh_sne
import scipy.ndimage as ndimage

import matplotlib.gridspec as gridspec
from matplotlib import animation
import matplotlib.animation as animation
from matplotlib.pyplot import figure
from mpl_toolkits.axes_grid1 import make_axes_locatable #Used for colorbar; allocates a bit of space for it


colors_names=['blue','red', 'green','violet','lightseagreen','lightsalmon','indianred','lightsalmon','pink','darkolivegreen']

colors = {
'#E24A33': (0.8862745098039215, 0.2901960784313726, 0.2), '#C4AD66': (0.7686274509803922, 0.6784313725490196, 0.4), '#fa8174': (0.9803921568627451, 0.5058823529411764, 0.4549019607843137), 
'#E8000B': (0.9098039215686274, 0.0, 0.043137254901960784), '#B0E0E6': (0.6901960784313725, 0.8784313725490196, 0.9019607843137255), '#7A68A6': (0.47843137254901963, 0.40784313725490196, 0.6509803921568628), 
'#ccebc4': (0.8, 0.9215686274509803, 0.7686274509803922), '.8': (0.8, 0.8, 0.8), 'w': (1.0, 1.0, 1.0), 'r': (1.0, 0.0, 0.0), 
'm': (0.75, 0, 0.75), '#4878CF': (0.2823529411764706, 0.47058823529411764, 0.8117647058823529), 'b': (0.0, 0.0, 1.0), 
'#03ED3A': (0.011764705882352941, 0.9294117647058824, 0.22745098039215686), '#bc82bd': (0.7372549019607844, 0.5098039215686274, 0.7411764705882353), '0.70': (0.7, 0.7, 0.7), 
'0.00': (0.0, 0.0, 0.0), '#56B4E9': (0.33725490196078434, 0.7058823529411765, 0.9137254901960784), '#6d904f': (0.42745098039215684, 0.5647058823529412, 0.30980392156862746), 
'#eeeeee': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), '#77BEDB': (0.4666666666666667, 0.7450980392156863, 0.8588235294117647), '#A60628': (0.6509803921568628, 0.023529411764705882, 0.1568627450980392), 
'#bcbcbc': (0.7372549019607844, 0.7372549019607844, 0.7372549019607844), '#55A868': (0.3333333333333333, 0.6588235294117647, 0.40784313725490196), '#8dd3c7': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451), 
'#64B5CD': (0.39215686274509803, 0.7098039215686275, 0.803921568627451), '#F0E442': (0.9411764705882353, 0.8941176470588236, 0.25882352941176473), '#81b1d2': (0.5058823529411764, 0.6941176470588235, 0.8235294117647058), 
'k': (0.0, 0.0, 0.0), '#bfbbd9': (0.7490196078431373, 0.7333333333333333, 0.8509803921568627), '#8b8b8b': (0.5450980392156862, 0.5450980392156862, 0.5450980392156862), 
'#777777': (0.4666666666666667, 0.4666666666666667, 0.4666666666666667), 'c': (0.0, 0.75, 0.75), '#348ABD': (0.20392156862745098, 0.5411764705882353, 0.7411764705882353), 
'#CC79A7': (0.8, 0.4745098039215686, 0.6549019607843137), '0.50': (0.5, 0.5, 0.5), '#fc4f30': (0.9882352941176471, 0.30980392156862746, 0.18823529411764706), 
'#FF9F9A': (1.0, 0.6235294117647059, 0.6039215686274509), 'blue': (0.0, 0.0, 1.0), '#afeeee': (0.6862745098039216, 0.9333333333333333, 0.9333333333333333), '#97F0AA': (0.592156862745098, 0.9411764705882353, 0.6666666666666666), '.15': (0.15, 0.15, 0.15), '0.8': (0.8, 0.8, 0.8), 'cyan': (0.0, 1.0, 1.0), '#ffed6f': (1.0, 0.9294117647058824, 0.43529411764705883), '#00D7FF': (0.0, 0.8431372549019608, 1.0), '#FFFEA3': (1.0, 0.996078431372549, 0.6392156862745098), '#cbcbcb': (0.796078431372549, 0.796078431372549, 0.796078431372549), '#D0BBFF': (0.8156862745098039, 0.7333333333333333, 1.0), '#B47CC7': (0.7058823529411765, 0.48627450980392156, 0.7803921568627451), '#92C6FF': (0.5725490196078431, 0.7764705882352941, 1.0), '0.5': (0.5, 0.5, 0.5), '#B8860B': (0.7215686274509804, 0.5254901960784314, 0.043137254901960784), '#6ACC65': (0.41568627450980394, 0.8, 0.396078431372549), '#FFB5B8': (1.0, 0.7098039215686275, 0.7215686274509804), '#988ED5': (0.596078431372549, 0.5568627450980392, 0.8352941176470589), '#467821': (0.27450980392156865, 0.47058823529411764, 0.12941176470588237), '#e5ae38': (0.8980392156862745, 0.6823529411764706, 0.2196078431372549), 'y': (0.75, 0.75, 0), '#8EBA42': (0.5568627450980392, 0.7294117647058823, 0.25882352941176473), '#b3de69': (0.7019607843137254, 0.8705882352941177, 0.4117647058823529), '#fdb462': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961), 'purple': (0.5019607843137255, 0.0, 0.5019607843137255), '#003FFF': (0.0, 0.24705882352941178, 1.0), 'magenta': (1.0, 0.0, 1.0), '0.75': (0.75, 0.75, 0.75), '#0072B2': (0.0, 0.4470588235294118, 0.6980392156862745), '#8A2BE2': (0.5411764705882353, 0.16862745098039217, 0.8862745098039215), '#30a2da': (0.18823529411764706, 0.6352941176470588, 0.8549019607843137), '0.40': (0.4, 0.4, 0.4), '#feffb3': (0.996078431372549, 1.0, 0.7019607843137254), '#8172B2': (0.5058823529411764, 0.4470588235294118, 0.6980392156862745), '#7600A1': (0.4627450980392157, 0.0, 0.6313725490196078), '#EAEAF2': (0.9176470588235294, 0.9176470588235294, 0.9490196078431372), '#EEEEEE': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), '#E5E5E5': (0.8980392156862745, 0.8980392156862745, 0.8980392156862745), '#f0f0f0': (0.9411764705882353, 0.9411764705882353, 0.9411764705882353), '#001C7F': (0.0, 0.10980392156862745, 0.4980392156862745), '#C44E52': (0.7686274509803922, 0.3058823529411765, 0.3215686274509804), '#017517': (0.00392156862745098, 0.4588235294117647, 0.09019607843137255), 'g': (0.0, 0.5, 0.0), '#D65F5F': (0.8392156862745098, 0.37254901960784315, 0.37254901960784315), '#8C0900': (0.5490196078431373, 0.03529411764705882, 0.0), '#006374': (0.0, 0.38823529411764707, 0.4549019607843137), '#555555': (0.3333333333333333, 0.3333333333333333, 0.3333333333333333), '#CCB974': (0.8, 0.7254901960784313, 0.4549019607843137), 'darkgoldenrod': (0.7215686274509804, 0.5254901960784314, 0.043137254901960784), '0.60': (0.6, 0.6, 0.6), 'gray': (0.5019607843137255, 0.5019607843137255, 0.5019607843137255), 'green': (0.0, 0.5019607843137255, 0.0), '#00FFCC': (0.0, 1.0, 0.8), '#FBC15E': (0.984313725490196, 0.7568627450980392, 0.3686274509803922), 'black': (0.0, 0.0, 0.0), '#4C72B0': (0.2980392156862745, 0.4470588235294118, 0.6901960784313725), 'firebrick': (0.6980392156862745, 0.13333333333333333, 0.13333333333333333), '#FFC400': (1.0, 0.7686274509803922, 0.0), 'red': (1.0, 0.0, 0.0), 'white': (1.0, 1.0, 1.0), '0.6': (0.6, 0.6, 0.6), 'yellow': (1.0, 1.0, 0.0), '#D55E00': (0.8352941176470589, 0.3686274509803922, 0.0), '#009E73': (0.0, 0.6196078431372549, 0.45098039215686275)}

#**************************************************************        

def PCA(X, n_components):
    from sklearn import decomposition

    pca = decomposition.PCA(n_components)
    pca.fit(X)
    X=pca.transform(X)

    coords = []
    for i in range(len(X)):
         coords.append([X[i][0], X[i][1], X[i][2]])
    
    return X, np.array(coords).T #THIS IS REDUNDANT... REDUCE IT

def Meanshift(data):
    from sklearn.cluster import MeanShift, estimate_bandwidth
    from sklearn.datasets.samples_generator import make_blobs
    
    # The following bandwidth can be automatically detected using
    bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=1000)

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(data)
    labels = ms.labels_
    cluster_centers = ms.cluster_centers_

    labels_unique = np.unique(labels)
    n_clusters_ = len(labels_unique)

    print("number of estimated clusters : %d" % n_clusters_)

    return labels
    
def KMEANS(data, n_clusters):

    from sklearn import cluster, datasets
    clusters = cluster.KMeans(n_clusters, max_iter=1000, n_jobs=-1, random_state = 121)
    clusters.fit(data)
    
    return clusters.labels_

 

def average_parallel((data_in)):
    
    print "in... "
    ave_data = np.ma.average(data_in, axis=0)
    print "...done"

    return ave_data

def fast_mask(mouse, data, generic_mask_indexes):
                
    n_pixels = mouse.n_pixels
    temp_array = np.ma.array(np.zeros((len(data),n_pixels,n_pixels),dtype=np.float32), mask=True)
        
    #Mask all frames; NB: PROBABLY FASTER METHOD
    for i in range(0, len(data),1):
        if len(data)!=128:
            temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes, fill=np.nan)
        else:
            temp_array[i] = np.ma.masked_array(data, mask=generic_mask_indexes)
            print "single frame"
            return temp_array[0]
    
    return temp_array

def load_generic_mask(mouse):
    
    n_pixels=mouse.n_pixels
    generic_mask_file = mouse.home_dir + mouse.name+'/genericmask.txt'
    generic_coords = np.loadtxt(generic_mask_file)
        
    generic_mask_indexes=np.zeros((n_pixels,n_pixels))
    for i in range(len(generic_coords)):
        generic_mask_indexes[generic_coords[i][0]][generic_coords[i][1]] = True

    #Mask centreline also
    for x in range(n_pixels/2-10, n_pixels/2+10, 1):
        for k in range(n_pixels):
            generic_mask_indexes[k][x]=True
    
    return generic_mask_indexes


def plot_pyqt(app, coords, color_index):
    
    w = gl.GLViewWidget()
    w.setBackgroundColor([0,0,0])
    w.opts['distance'] = 20
    w.show()
    w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

    g = gl.GLGridItem()
    w.addItem(g)

    pos = []
    
    size = []
    color = []
    for i in range(len(coords)):
        if coords.shape[1] == 3:         pos.append([coords[i][0], coords[i][1], coords[i][2]])
        else:        pos.append([coords[i][0], coords[i][1], 0])

        size.append(10)
        color.append(colors[str(colors.keys()[color_index[i]])])  #May wish to randomize color loading
        
    pos = np.array(pos)
    size = np.array(size)
    color = np.array(color)
    
    sp1 = gl.GLScatterPlotItem(pos=pos, size=size, color=color, pxMode=False)
    w.addItem(sp1)


    ## Start Qt event loop unless running in interactive mode.
    QtGui.QApplication.instance().exec_()
    app.closeAllWindows()

def pre_post_stroke_traces(mouse): 
    
    #OLD CODE
    mouse.temp_std = []
    mouse.temp_ave = []
    mouse.trace_time = []
    mouse.prestroke = []
    mouse.poststroke = []
    n_traces = []
    for k in range(n_clusters):
        temp_event_times = []
        ave_trace = []
        ave_trace_prestroke = []
        ave_trace_poststroke = []
        
        for j in range(len(mouse.all_traces)):
            if cluster_labels[j]==k:
                ave_trace.append(mouse.all_traces[j])
                temp_event_times.append(j)

                if j < mouse.stroke_index:
                    ave_trace_prestroke.append(mouse.all_traces[j])
                else:
                    ave_trace_poststroke.append(mouse.all_traces[j])
        
        if len(ave_trace)==0: continue

        mouse.trace_time.append(((np.array(temp_event_times)/float(len(mouse.all_traces)))*mouse.window*2/mouse.img_rate)-mouse.window/mouse.img_rate)

        n_traces.append(len(ave_trace))
        mouse.temp_std.append(np.std(np.array(ave_trace), axis=0))
        mouse.temp_ave.append(np.average(np.array(ave_trace),axis=0))

        mouse.prestroke.append(ave_trace_prestroke)
        mouse.poststroke.append(ave_trace_poststroke)


def plot_traces_pyqt(app, mouse):
    
    print "PLOTING TRACES..."
    
    #global app, selected_window, n_traces
    #app = QtGui.QApplication([])    #Starts up the QtGui; 

    print "gui initiatied..."
    
    
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    window = pg.GraphicsWindow(title="Clustered traces")
    window.resize(1000,600)
    window.setWindowTitle(mouse.name+" #: "+str(len(mouse.traces)))

    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)

    all_traces = []
    for k in range(len(mouse.traces)):
        all_traces.extend(mouse.traces[k])
    all_traces = np.array(all_traces)
    
    all_trial_times = []
    for k in range(len(mouse.trial_times)):
        all_trial_times.extend(mouse.trial_times[k])
    all_trial_times = np.array(all_trial_times)
    
    
    print len(all_traces)
    print len(mouse.cluster_labels)


    #OLD CODE
    mouse.temp_std = []
    mouse.temp_ave = []
    mouse.trace_time = []
    mouse.prestroke = []
    mouse.poststroke = []
    n_traces = []
    for k in range(mouse.n_clusters):
        ave_trace = []
        temp_event_times = []
        for j in range(len(all_traces)):
            if mouse.cluster_labels[j]==k:
                ave_trace.append(all_traces[j])
                temp_event_times.append(all_trial_times[j])

        if len(ave_trace)==0: continue

        mouse.trace_time.append((np.array(temp_event_times)*mouse.window*2/mouse.img_rate)-mouse.window/mouse.img_rate)

        n_traces.append(len(ave_trace))
        mouse.temp_std.append(np.std(np.array(ave_trace), axis=0))
        mouse.temp_ave.append(np.average(np.array(ave_trace),axis=0))


    #PLOT CURVES
    print "Plotting curves..."
    xx = np.linspace(-mouse.window/mouse.img_rate,mouse.window/mouse.img_rate,len(all_traces[0]))
    print len(xx), len(mouse.temp_ave[0])
    for k in range(len(n_traces)):
        print "Cluster : ", k
        win = window.addPlot(title="Cluster: "+str(k)+", "+str(n_traces[k]) )
        win.setXRange(-mouse.window/mouse.img_rate, mouse.window/mouse.img_rate, padding=0)
        win.setYRange(0, 90, padding=0)

        win.plot(pen = {'color':255*np.array(colors[str(colors.keys()[k])]), 'width':5}, x=xx, y=mouse.temp_ave[k])

        p1 = win.plot(xx, mouse.temp_ave[k]+mouse.temp_std[k])
        p2 = win.plot(xx, mouse.temp_ave[k]-mouse.temp_std[k])
        temp_color = list(255*np.array(colors[str(colors.keys()[k])]))
        fill = pg.FillBetweenItem(p1, p2, temp_color)
        win.addItem(fill)

        #Add dashed lines
        win.plot(pen = {'color':(0,0,0), 'width':2, 'style':QtCore.Qt.DashLine}, x=[xx[0],xx[-1]], y = [11,11])
        win.plot(pen = {'color':(0,0,0), 'width':2, 'style':QtCore.Qt.DashLine}, x=[xx[0],xx[-1]], y = [59,59])
        win.plot(pen = {'color':(0,0,0), 'width':2, 'style':QtCore.Qt.DashLine}, x=[0,0], y = [0,80])
        

        #Add locations of pulls
        x = []
        y = []
        for p in range(len(mouse.trace_time[k])):
            x.extend([mouse.trace_time[k][p], mouse.trace_time[k][p]])
            y.extend([11, 59])

        x = np.array(x)
        y = np.array(y)
        connect = np.ones(len(y), dtype=np.ubyte)
        connect[1::2] = 0  #  disconnect segment between lines
        path = pg.arrayToQPath(x, y, connect)   #CAN ALSO COLAPSE A 2D ARRAY TO A 1D ARRAY USING: pg.arrayToQPath(x.reshape(lines*points), y.reshape(lines*points), connect.reshape(lines*points))
        item = pg.QtGui.QGraphicsPathItem(path)
        item.setPen(pg.mkPen((0,0,0), width=.2))
        win.addItem(item)

        #Add location of stroke
        temp_color = (255,0,0)
        if mouse.stroke.kind == 'sham': temp_color = (0,255,0)
        
        #This plots stroke_index time relative indexes for all lever pulls; NB: DFF plotter needs stroke_index relative specific cluster lever pulls
        stroke_loc = (mouse.stroke.trial/float(mouse.n_trials)-0.5)*6
        win.plot(pen = {'color':temp_color, 'width':5}, x=[stroke_loc,stroke_loc], y = [5,65])
                
        #Loop for every 
        if (k%np.sqrt(len(n_traces)))==(np.sqrt(len(n_traces))-1): window.nextRow()

    #mouse.stroke_mark = mouse.stroke_index/float(len(mouse.all_traces))*mouse.window*2/30-mouse.window/30  #Save relative time for use in later plots

    #selected_window = 0
    def onClick(event):
        if event.button() == 1: #Mouse left button; right click can still save data...
            items = window.scene().items(event.scenePos())
            mouse_xy = [x for x in items if isinstance(x, pg.PlotItem)]
            size = window.size(); width = size.width(); height = size.height()

            p_width = width/np.sqrt(mouse.n_clusters)
            p_height = height/np.sqrt(mouse.n_clusters)
                        
            mouse.selected_cluster = int(int(mouse_xy[0].x()*1.1/p_width)+np.sqrt(mouse.n_clusters)*int(mouse_xy[0].y()*1.1/p_height))
            #print "CLUSTER SELECTED: ", mouse.selected_cluster, "  n_clusters:", n_clusters
            
            app.closeAllWindows()

    window.scene().sigMouseClicked.connect(onClick)
    
    #Send to qtgui
    QtGui.QApplication.instance().exec_()


def plot_array(mouse, data_array, ax1, generic_mask_indexes):

    masked = np.ma.array(np.zeros((mouse.n_pixels,mouse.n_pixels),dtype=np.float32), mask=True)
    temp_plot1 = []
    for k in range(0, len(data_array[1]), 128):
        masked = data_array[:, k:k+128] #Chunk the data back to individual frames
        temp_plot1.append(masked)
    
    temp_plot1 = fast_mask(mouse, temp_plot1, generic_mask_indexes)
    data_array = np.ma.hstack((temp_plot1))
    ax1.yaxis.set_ticks([])
    ax1.xaxis.set_ticks([])
    
    v_max = np.max(np.abs(data_array))
    v_min = -v_max
    
    return data_array, v_min, v_max

#def plot_array_allframes(mouse, index, chunk_size, ax1, generic_mask_indexes):

    #data = mouse.non_norm_plot
    #temp_plot = []
    #for p in range(len(data)):
        #temp_plot.append(fast_mask(mouse, data[p], generic_mask_indexes))
    
    #stacks = []
    #for p in range(len(temp_plot)):
        #stacks.append(np.ma.hstack((temp_plot[p][index*len(temp_plot[p])/chunk_size:(index+1)*len(temp_plot[p])/chunk_size])))
    
    #stacked_plot = np.ma.vstack((stacks))
           
    #ax1.yaxis.set_ticks([])
    #ax1.xaxis.set_ticks([])
    
    #v_max = np.max(np.abs(stacks))
    #v_min = -v_max
    
    #return stacked_plot, v_min, v_max


def plot_array_aveframes(mouse, data, block, generic_mask_indexes):
    
    #Load masks
    contour = np.load(mouse.home_dir+mouse.name+'/stroke_contour.npy')
    
    #data = mouse.non_norm_plot
    stacks = []
    for p in range(len(data)):  #Mask each epoch of 180 frames
        temp_img = fast_mask(mouse, data[p], generic_mask_indexes)
        if mouse.name in ['IA1', 'IA2', 'IA3', 'IJ1', 'IJ2']: 
            for f in range(len(temp_img)):
                for l in range(len(temp_img[f])):
                    temp_img[f][l][contour[l]]=np.nan
            
        stacks.append(temp_img)
    
    stacks_block = []
    for p in range(len(stacks)):
        temp_stack = []
        for q in range(0,len(stacks[p]), block):
            temp_stack.append(np.ma.average(stacks[p][q:q+block], axis=0))
        stacks_block.append(np.ma.hstack((temp_stack)))
   
    return stacks_block, stacks

def average_sum(mouse, data_array, k):

    print " averaging epoch: ", k, " ... "

    #non_norm += np.array(data_array[p])
    mouse.non_norm_plot[k] = np.average(data_array, axis=0)

def plot_1D(mouse, generic_mask_indexes):
    
    epoch_sums = np.zeros((3,2,mouse.img_rate*6), dtype=np.float32)
    for p in range(len(mouse.non_norm_plot)):
        temp = fast_mask(mouse, mouse.non_norm_plot[p], generic_mask_indexes)
        for k in range(len(temp)):
            epoch_sums[p][0][k] = np.ma.average(temp[k][:70,:64])
            epoch_sums[p][1][k] = np.ma.average(temp[k][:70,64:])
    
    labels = ['Pre-stroke', 'Post-stroke (2 weeks)', 'Post-stroke (4 weeks)']
    t = np.linspace(-3,3,mouse.img_rate*6)
    for p in range(3):
        ax=plt.subplot(2,3,p+1)
        plt.title(labels[p])
        plt.plot(t, epoch_sums[p][0]*100, linewidth=5, color='red')
        plt.plot(t, epoch_sums[p][1]*100, linewidth=5, color='blue')
        plt.xlabel("Time (sec)", fontsize=15)
        plt.ylabel("Average pixel DFF (%DFF)", fontsize=15, labelpad=-4)
        plt.ylim(-5,15)
        plt.plot([0,0],[-5,15], 'r--', color='black', linewidth=2, alpha=0.5)
        plt.plot([-3,3],[0,0], 'r--', color = 'black', linewidth=2, alpha=0.5)
        

    plt.suptitle(mouse.name+"\nRight (blue) vs. Left (red) Hemisphere Average Pixel Activity", fontsize=20)
    plt.show()
    
def pixel_plot(mouse, generic_mask_indexes):
    
    #epoch_sums = np.zeros((3,2,180), dtype=np.float32)
    n_pixels_left=[]; n_pixels_right = []
    for p in range(len(mouse.non_norm_plot)): #LOOP OVER ALL 3 EPOCHS
        n_pixels_left.append([]); n_pixels_right.append([])
        
        #Mask all 180 frames for each epoch
        temp = fast_mask(mouse, mouse.non_norm_plot[p], generic_mask_indexes)

        #Find max DFF value over entire stack
        v_max_left = np.ma.max(temp[:, :70, :])
        v_max_right = np.ma.max(temp[:, :70, :])
        
        print "v_max: ", v_max_left, v_max_right
        
        for k in range(len(temp)):
            n_pixels_left[p].append(len(np.where(temp[k][:70,:64]>(v_max_left/2.))[0]))
            n_pixels_right[p].append(len(np.where(temp[k][:70,64:]>(v_max_right/2.))[0]))
        
        np.savetxt(mouse.home_dir+mouse.name+'/left_epoch_'+str(p), n_pixels_left[p])
        np.savetxt(mouse.home_dir+mouse.name+'/right_epoch_'+str(p), n_pixels_right[p])

    labels = ['Pre-stroke', 'Post-stroke (2 weeks)', 'Post-stroke (4 weeks)']
    t = np.linspace(-3,3,mouse.img_rate*6)
    h_max =max(np.ma.max(n_pixels_left), np.ma.max(n_pixels_right))
    for p in range(len(mouse.non_norm_plot)):
        ax=plt.subplot(2,3,p+1)
        plt.title(labels[p])
        plt.plot(t, n_pixels_left[p], linewidth=5, color='red')
        plt.plot(t, n_pixels_right[p], linewidth=5, color='blue')
        plt.xlabel("Time (sec)", fontsize=15)
        plt.ylabel("No. Pixels > 50% Max(DFF)", fontsize=15, labelpad=-2)
        
        plt.ylim(-5,h_max+100)
        plt.plot([0,0],[-5,h_max+100], 'r--', color='black', linewidth=2, alpha=0.5)
        plt.plot([-3,3],[0,0], 'r--', color = 'black', linewidth=2, alpha=0.5)
        

    plt.suptitle(mouse.name+"\nRight (blue) vs. Left (red) Hemisphere Average Pixel Activity", fontsize=20)
    plt.show()
    
    #quit()

def plot_selected_cluster_DFF(mouse):
    ''' Function loads DFF for individual trials based on clustered ids
        Plots [Ca] dynamics
     '''

    data_chunks = []
    #n_trials = 5000
    #Loop over chunks/epochs of traces and
    for p in range(len(mouse.DFF_files)): 
        data_chunks.append([])
        for k in range(len(mouse.DFF_files[p])):
            trial_DFF = np.load(mouse.DFF_files[p][k][:-9]+'_3sec'+mouse.DFF_files[p][k][-9:])

            #OPTION MASK TRIAL BEFORE CHECKING FOR MAX VALS! OR DO np.array max search
            if np.isnan(trial_DFF).any() or (np.max(trial_DFF)>0.5) or (np.min(trial_DFF)<-0.5):
                print "---------------------SKIPPING NAN OR TOO LARGE/SMALL VALS----------------"
                continue
                
            print "Chunk: ", p, " trial #: ", len(data_chunks[p]), " /  ~", str(len(mouse.DFF_files[p]))
            data_chunks[p].append(trial_DFF)
    
            #if len(data_chunks[p])>n_trials: break
    
    #************** MAKE SINGLE TRIAL MOVIES
    if False: make_movies_singletrials(mouse, data_chunks)  #NOT WORKING YET
    
    
    #***************SUM DFF TRIALS IN PARALLEL
    functions = []
    for k in range(len(data_chunks)):
        functions.append(average_sum) #*len(data_chunks)        #Function that averages data is called "average_sum"
    
    mouse.non_norm_plot=[[]]*len(data_chunks) #Make lists to hold the averaged [ca] data for each chunk in data_array
    
    par_funcs = []
    for k in range(len(functions)):
        par_funcs.append(Thread(target=functions[k], args=(mouse, data_chunks[k], k)))
    
    for k in range(len(par_funcs)):
        par_funcs[k].start()
    
    for k in range(len(par_funcs)):
        par_funcs[k].join()
    
  
    #*********************PLOTING PARAMETERS 
    clr = 'red'
    if mouse.stroke.kind == 'sham': clr='green'
    left_border = -3    #Second where to start data plot
    right_border = 3    #Second where to end data plot

    #Load generic mask for 128 pixel data
    generic_mask_indexes = load_generic_mask(mouse)


    
    #**************1D PLOTS
    if True: plot_1D(mouse, generic_mask_indexes)   

    #************** PIXEL PLOTS
    if True: pixel_plot(mouse, generic_mask_indexes)

    #Generate epochs in list of stacks; both raw and averaged every "block" frames
    block = 6               #pick the # frame block to average together.
    stacks_block, stacks  = plot_array_aveframes(mouse, mouse.non_norm_plot, block, generic_mask_indexes)
    temp_plot = np.ma.vstack((stacks_block));  v_max = np.max(np.abs(temp_plot));  v_min = -v_max

    #Make movies of activity
    if False: make_movies(mouse, stacks, v_max, v_min)
    
    
    #BEGIN PLOT NON-NORMALIZED DATA
    fig = plt.figure()
    gs = gridspec.GridSpec(18,9)        #Make plotting size array to accomodate all chunks

    
    #*****************PLOT TRIAL LOCATIONS
    for q in range(len(data_chunks)):
        ax1 = fig.add_subplot(gs[0:1,q*3:(q+1)*3])                    #NB: data is already in frame format, i.e. not hstacked!
        ax1.yaxis.set_ticks([])
        plt.title("# trials: " + str(len(mouse.traces[q])), fontsize=15)
        plt.xlim(0, mouse.n_days)
        plt.ylim(0,1)
        
        #Plot trial locations
        for k in range(len(mouse.days_counter[q])):
            plt.plot([mouse.days_counter[q][k],mouse.days_counter[q][k]], [0,1], linewidth = 3, color = 'black', alpha=.15)
                
        #Plot stroke location
        clr='red'
        if mouse.stroke.kind == 'sham': clr='green'
        plt.plot([mouse.stroke.day,mouse.stroke.day], [0,1], linewidth=4, color=clr)
        
        #Plot locations of each month
        for k in range(len(mouse.month_counter)):
            plt.plot([mouse.month_counter[k],mouse.month_counter[k]], [0,1], 'r--', linewidth=1, color='black', alpha=0.15)
        
        #Generate labels
        old_xlabel = []; old_xlabel.append(0);
        for k in range(len(mouse.month_counter)): old_xlabel.append(mouse.month_counter[k])
        old_xlabel.append(mouse.stroke.day); old_xlabel.append(mouse.n_days)
        
        new_xlabel = []; new_xlabel.append(mouse.realtimes_first)
        for k in range(len(mouse.month_text)): new_xlabel.append(mouse.month_text[k])
        new_xlabel.append(mouse.stroke.fulldate); new_xlabel.append(mouse.realtimes_last)
        plt.xticks(old_xlabel, new_xlabel, rotation=90, fontsize=11)

       
    #******************PLOT LEVERPOS
    t = np.linspace(-mouse.window/mouse.img_rate,mouse.window/mouse.img_rate, len(mouse.traces[0][0]))
    new_ylabel = np.linspace(0,80,5)
    for q in range(len(data_chunks)):
        ax1 = fig.add_subplot(gs[3:5,q*3:(q+1)*3])                    #NB: data is already in frame format, i.e. not hstacked!
        trace_ave = np.average(mouse.traces[q], axis=0)
        trace_std = np.std(mouse.traces[q], axis=0)
        plt.plot(t, trace_ave, color='black')
        ax1.fill_between(t, trace_ave+trace_std, trace_ave-trace_std, color='blue', alpha=0.4)
        plt.ylim(0,90)
        if q==0: plt.ylabel("Lever Position", fontsize=12)
        plt.yticks(new_ylabel, fontsize=10)

        plt.plot([0,0],[0,90], 'r--', color='black', linewidth=3)
        
                
    
    #********************PLOT AVERAGED BLOCKS
    ax1 = fig.add_subplot(gs[6:12,:])                   
    #Set ticks to nada
    ax1.yaxis.set_ticks([])
    ax1.xaxis.set_ticks([])
    
    #Make temp_plot for raw data
    temp_plot = np.ma.vstack((stacks_block));  v_max = np.max(np.abs(temp_plot));  v_min = -v_max
    temp_plot[:,temp_plot.shape[1]/2-3:temp_plot.shape[1]/2+3] = v_min #Set t=0 line
    
    #Cut data to borders intended
    left_index = (3+left_border)*30/block*128
    right_index = (3+right_border)*30/block*128
    temp_plot = temp_plot[:,left_index:right_index]
    
    im = plt.imshow(temp_plot, vmin=v_min, vmax=v_max)

    plt.ylabel("Epochs:\n3rd - 2nd - 1st", fontsize=16)
    
    old_xlabel = np.linspace(0,temp_plot.shape[1], (right_border - left_border)*15/block+1)
    new_xlabel = np.around(np.linspace(left_border,right_border, (right_border - left_border)*15/block+1), decimals=2)
    plt.xticks(old_xlabel, new_xlabel, fontsize=10)
    
    #Colorbar
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "2%", pad="1%")
    cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
    cbar.ax.set_yticklabels([str(int(v_min*100))+'%', '0%', str(int(v_max*1000)/10)+'%'])  # horizontal colorbar


    #*****************PLOT AVERAGED BLOCKS - FIRST EPOCHS
    ax1 = fig.add_subplot(gs[12:15,:])                   
    ax1.yaxis.set_ticks([])
    ax1.xaxis.set_ticks([])
    
    #Make temp_plot for diff data
    temp_plot = stacks_block[1]-stacks_block[0];  v_max = np.max(np.abs(temp_plot));  v_min = -v_max
    temp_plot[:,temp_plot.shape[1]/2-3:temp_plot.shape[1]/2+3] = v_min #Set t=0 line
    
    #Cut data to borders intended
    left_index = (3+left_border)*30/block*128
    right_index = (3+right_border)*30/block*128
    temp_plot = temp_plot[:,left_index:right_index]
    
    im = plt.imshow(temp_plot, vmin=v_min, vmax=v_max)

    #plt.xlabel("Time (sec)", fontsize=20)
    plt.ylabel("Diff:\n2nd - 1st", fontsize=10)
    plt.xticks(old_xlabel, new_xlabel, fontsize=10)
    
    #Colorbar
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "2%", pad="1%")
    cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
    cbar.ax.set_yticklabels([str(int(v_min*100))+'%', '0%', str(int(v_max*1000)/10)+'%'])  # horizontal colorbar


    if len(stacks_block)==3: 
        ax1 = fig.add_subplot(gs[15:18,:])                   
        ax1.yaxis.set_ticks([])
        ax1.xaxis.set_ticks([])
        
        #Make temp_plot for diff data
        temp_plot = stacks_block[2]-stacks_block[0]; v_max = np.max(np.abs(temp_plot));  v_min = -v_max
        temp_plot[:,temp_plot.shape[1]/2-3:temp_plot.shape[1]/2+3] = v_min #Set t=0 line
        
        #Cut data to borders intended
        left_index = (3+left_border)*30/block*128
        right_index = (3+right_border)*30/block*128
        temp_plot = temp_plot[:,left_index:right_index]
        
        im = plt.imshow(temp_plot, vmin=v_min, vmax=v_max)

        plt.xlabel("Time (sec)", fontsize=20)
        plt.xticks(old_xlabel, new_xlabel, fontsize=10)
        plt.ylabel("Didff:\n3rd - 1st", fontsize=10)

        #Colorbar
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "2%", pad="1%")
        cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
        cbar.ax.set_yticklabels([str(int(v_min*100))+'%', '0%', str(int(v_max*1000)/10)+'%'])  # horizontal colorbar


    plt.suptitle(mouse.name + " Absolute DFF scale ", fontsize=20)
    plt.show()

def plot_single_trial(mouse):
    ''' Function loads DFF for individual trials based on clustered ids
        Plots [Ca] dynamics
     '''

    n_trials=20
    epoch_selected = 2
    data_chunks = []
    data_chunks_large = []
    ctr=0
    img_ctr=[]
    #Loop over chunks/epochs of traces and
    #for p in range(len(mouse.DFF_files)): 
    for k in range(len(mouse.DFF_files[epoch_selected])):
        trial_DFF = np.load(mouse.DFF_files[epoch_selected][k][:-9]+'_3sec'+mouse.DFF_files[epoch_selected][k][-9:])

        #OPTION MASK TRIAL BEFORE CHECKING FOR MAX VALS! OR DO np.array max search
        if (np.max(trial_DFF)>0.5) or (np.min(trial_DFF)<-0.5):
            img_ctr.append(ctr)
            print "max threshold reached: ", ctr
        if np.isnan(trial_DFF).any():
            print "---------------------SKIPPING NAN----------------"
            continue
            
        print "Chunk: ", epoch_selected, " trial #: ", len(data_chunks), " /  ~", str(len(mouse.DFF_files[epoch_selected]))
        data_chunks.append(trial_DFF)
        
        ctr+=1
    #*********************PLOTING PARAMETERS 
    left_border = -3    #Second where to start data plot
    right_border = 3    #Second where to end data plot

    #Load generic mask for 128 pixel data
    generic_mask_indexes = load_generic_mask(mouse)
    
    #Generate epochs in list of stacks; both raw and averaged every "block" frames
    block = 6               #pick the # frame block to average together.
    data = data_chunks[:n_trials]
    print len(data)
    stacks_block, stacks  = plot_array_aveframes(mouse, data, block, generic_mask_indexes)
    temp_plot = np.ma.vstack((stacks_block));  v_max = np.max(np.abs(temp_plot));  v_min = -v_max

    #BEGIN PLOT NON-NORMALIZED DATA
    fig = plt.figure()
    #gs = gridspec.GridSpec(18,9)        #Make plotting size array to accomodate all chunks

       
    ##******************PLOT LEVERPOS
    #t = np.linspace(-mouse.window/mouse.img_rate,mouse.window/mouse.img_rate, len(mouse.traces[0][0]))
    #new_ylabel = np.linspace(0,80,5)
    #for q in range(len(data_chunks)):
        #ax1 = fig.add_subplot(gs[3:5,q*3:(q+1)*3])                    #NB: data is already in frame format, i.e. not hstacked!
        #trace_ave = np.average(mouse.traces[q], axis=0)
        #trace_std = np.std(mouse.traces[q], axis=0)
        #plt.plot(t, trace_ave, color='black')
        #ax1.fill_between(t, trace_ave+trace_std, trace_ave-trace_std, color='blue', alpha=0.4)
        #plt.ylim(0,90)
        #if q==0: plt.ylabel("Lever Position", fontsize=12)
        #plt.yticks(new_ylabel, fontsize=10)

        #plt.plot([0,0],[0,90], 'r--', color='black', linewidth=3)
        
    
    #********************PLOT AVERAGED BLOCKS
    ctr=0
    #for stack in stacks_block:
    for stack in stacks_block[:n_trials]:
        print "...plotting trial: ", ctr
        #ax1 = fig.add_subplot(gs[ctr+1,:])                   
        ax1 = plt.subplot(len(stacks_block), 1, ctr+1)
        #Set ticks to nada
        ax1.yaxis.set_ticks([])
        ax1.xaxis.set_ticks([])

        #Make temp_plot for raw data
        temp_plot = stack;  v_max = np.max(np.abs(temp_plot));  v_min = -v_max
        temp_plot[:,temp_plot.shape[1]/2-3:temp_plot.shape[1]/2+3] = v_min #Set t=0 line
        
        #Cut data to borders intended
        left_index = (3+left_border)*30/block*128
        right_index = (3+right_border)*30/block*128
        temp_plot = temp_plot[:,left_index:right_index]
        
        im = plt.imshow(temp_plot, vmin=v_min, vmax=v_max)

        label = plt.ylabel(str(ctr), fontsize=16)
        if ctr in img_ctr: label.set_color("red")
        
        if ctr==(len(stacks_block)-1):
            old_xlabel = np.linspace(0,temp_plot.shape[1], (right_border - left_border)*15/block+1)
            new_xlabel = np.around(np.linspace(left_border,right_border, (right_border - left_border)*15/block+1), decimals=2)
            plt.xticks(old_xlabel, new_xlabel, fontsize=10)
        
        #Colorbar
        divider = make_axes_locatable(plt.gca())
        cax = divider.append_axes("right", "2%", pad="1%")
        cbar = plt.colorbar(im, cax=cax, ticks=[v_min, 0, v_max], orientation='vertical')
        cbar.ax.set_yticklabels([str(int(v_min*100))+'%', '0%', str(int(v_max*1000)/10)+'%'], fontsize=6)  # horizontal colorbar
        
        ctr+=1

    plt.suptitle(mouse.name + " Absolute DFF scale, epoch: "+str(epoch_selected)+", #trials: "+str(len(stacks_block))+"\n("+str(block)+"-block frame averages)", fontsize=20)
    plt.show()


def make_movies_singletrials(mouse, stacks):

    #******* LOAD MASKS
    #Load generic mask for 128 pixel data
    generic_mask_indexes = load_generic_mask(mouse)
    
    #Load masks
    contour = np.load(mouse.home_dir+mouse.name+'/stroke_contour.npy')

    #Mask individual frames of data
    data = stacks
    
    #Loop over 3 epochs:
    stacks_out = []
    n_trials = 49
    for k in range(len(data)):
        print "...masking epoch: ", k
        temp_data = []
        for p in range(len(data[k][:n_trials])):  #Mask each epoch of 180 frames
            print "...trial: ", p
            temp_img = fast_mask(mouse, data[k][p], generic_mask_indexes)
            #for f in range(len(temp_img)):
                #for l in range(len(temp_img[f])):
                    #temp_img[f][l][contour[l]]=np.nan
            #temp_array[i] = np.ma.masked_array(data[i], mask=generic_mask_indexes, fill=np.nan)

            temp_data.append(temp_img)
            
        stacks_out.append(temp_data)

    stacks = stacks_out

    left_border = -3    #Second where to start data plot
    right_border = 3    #Second where to end data plot

    #***********GENERATE ANIMATIONS
    for e in range(len(stacks)):
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

        #Convert list of images to masked arrays required
        vid_array = np.ma.array(np.zeros((n_trials,90,128,128),dtype=np.float32), mask=True) #Not sure if necessary or overkill
        for k in range(len(stacks[e])):
            vid_array[k] = np.ma.array(stacks[e][k][60:150])
        
        titles = ["Pre-stroke (3 weeks)", "Post-stroke (0-2 weeks)", "Post-stroke (2-4 weeks)"]
        fig = plt.figure()
        #fig.tight_layout()

        im=[]
        for k in range(len(stacks[e])):
            
            im.append([])
            ax = plt.subplot(np.sqrt(n_trials), np.sqrt(n_trials),k+1) 
            
            v_max = np.ma.max(np.ma.abs(vid_array[k])); v_min = -v_max
            print v_max, v_min
            print vid_array[k].shape
            
            ax.get_xaxis().set_visible(False)
            ax.yaxis.set_ticks([])
            plt.ylabel(str(int(v_max*100))+ ".."+str(int(v_min*100)), labelpad=-1, fontsize=6)
            im[k] = plt.imshow(vid_array[k][0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
            
        #function to update figure
        def updatefig(j):
            print j,
            plt.suptitle(mouse.name + "  "+titles[e]+"  Frame: "+str(j)+"  " +str(format(float(j)/mouse.img_rate-1.,'.2f'))+"sec", fontsize = 15)

            # set the data in the axesimage object
            for k in range(len(vid_array)):
                im[k].set_array(vid_array[k][j])

            # return the artists set
            return im
            
        # kick off the animation
        ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array[0])), interval=100, blit=False, repeat=True)

        if True:
        #if save_animation:
            ani.save(mouse.home_dir+mouse.name+'/'+mouse.name+'_epoch_'+str(e)+'.mp4', writer=writer)

        plt.show()

    #quit()


    
    
    
def make_movies(mouse, stacks, v_max, v_min):

    #***********GENERATE ANIMATIONS
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)

    #Convert list of images to masked arrays required
    vid_array = np.ma.array(np.zeros((3,90,128,128),dtype=np.float32), mask=True) #Not sure if necessary or overkill
    for k in range(len(stacks)):
        vid_array[k] = np.ma.array(stacks[k][60:150])
    
    titles = ["Pre-stroke (3 weeks)", "Post-stroke (0-2 weeks)", "Post-stroke (2-4 weeks)"]
    fig = plt.figure()
    gs = gridspec.GridSpec(16,150)        #Make plotting size array to accomodate all chunks
    locs=[[0,10,0,45],[0,10,50,95],[0,10,100,150]]
    im=[]
    for k in range(len(stacks)+1):
        
        im.append([])

        ##Colorbar
        if k == 3: 
            #divider = make_axes_locatable(plt.gca())
            #cax = divider.append_axes("right", "2%", pad="1%")
            im[k] = plt.colorbar(im[0], ticks=[v_min, 0, v_max], orientation='vertical', fraction=0.046)
            im[k].ax.set_yticklabels([str(int(v_min*100))+'%', '0%', str(int(v_max*1000)/10)+'%'])  # horizontal colorbar
            continue
    
        ax = fig.add_subplot(gs[locs[k][0]:locs[k][1], locs[k][2]:locs[k][3]]) #q*3:(q+1)*3])                    #NB: data is already in frame format, i.e. not hstacked!
        
        #v_max = np.ma.max(np.ma.abs(np.ma.vstack((stacks[k]))));  v_min = -v_max
        
        print v_max, v_min
        print vid_array[k].shape
        
        ax.get_xaxis().set_visible(False)
        ax.yaxis.set_ticks([])
        ax.yaxis.labelpad = 0
        #ax.set_ylabel("",fontsize=6)
        ax.set_title(titles[k], fontsize=11)
        #if k ==0: 
        #    plt.ylabel("DFF: "+str(int(v_min*100))+"%.."+str(int(v_max*100))+"%", fontsize=12)
        
        im[k] = plt.imshow(vid_array[k][0], cmap=plt.get_cmap('jet'), vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)
        #im[k] = plt.imshow(vid_array[k][0], cmap=blue_red1, vmin=v_min, vmax=v_max, interpolation='none')#, vmin=0, vmax=v_max)

        
    #function to update figure
    def updatefig(j):
        print j,
        plt.suptitle(mouse.name + "\nFrame: "+str(j)+"  " +str(format(float(j)/mouse.img_rate-1.,'.2f'))+"sec", fontsize = 20)

        # set the data in the axesimage object
        for k in range(len(vid_array)):
            im[k].set_array(vid_array[k][j])

        # return the artists set
        return im
        
    # kick off the animation
    ani = animation.FuncAnimation(fig, updatefig, frames=range(len(vid_array[0])), interval=100, blit=False, repeat=True)

    if True:
    #if save_animation:
        ani.save(mouse.home_dir+mouse.name+'/'+mouse.name+'.mp4', writer=writer)

    plt.show()

    #quit()



def multi_dim_scaling(mouse, method):
    
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
    
    
    





def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

def find_previous(array,value):
    temp = (np.abs(array-value)).argmin()
    if array[temp]>value: return temp-1
    else: return temp
    



#def plot_traces(mouse, all_traces, cluster_labels, n_clusters):

    #temp_std = []
    #temp_ave = []
    #trace_time = []
    #n_traces = []
    #for k in range(n_clusters):
        #temp_event_times = []
        #ave_trace = []
        #for j in range(len(all_traces)):
            #if cluster_labels[j]==k:
                #ave_trace.append(all_traces[j])
                #temp_event_times.append(j)
        
        #if len(ave_trace)==0: continue

        #trace_time.append(((np.array(temp_event_times)/float(len(all_traces)))*.84)-.42)

        #n_traces.append(len(ave_trace))
        #temp_std.append(np.std(np.array(ave_trace), axis=0))
        #temp_ave.append(np.average(np.array(ave_trace),axis=0))

    ##SORT CURVES
    ##cc = []
    ##for k in range(len(temp_ave)):
        ##cc.append(np.corrcoef(temp_ave[0], temp_ave[k])[0][1])

    ##indexes = np.argsort(cc[::-1])
    ##temp_std = np.array(temp_std)[indexes]
    ##temp_ave = np.array(temp_ave)[indexes]
    ##n_traces = np.array(n_traces)[indexes]


    ##PLOT CURVES
    #print "Plotting curves..."
    #xx = np.linspace(-.42,.42,100)
    #for k in range(len(n_traces)):
        #ax1=plt.subplot(np.sqrt(n_clusters),np.sqrt(n_clusters),k+1)
        #plt.plot([0,0],[-1,1], color='black')
        #plt.plot([-3,3],[0,0], color='black')
        #plt.ylim(0,80)
        #plt.xlim(-.42,.42)        
            
        #for p in range(len(trace_time[k])):
            #ax1.axvline(trace_time[k][p], ymin=0.25, ymax=0.75, color = 'black', alpha=0.3)
        
        #ax1.plot([0,0],[0,80], 'r--', linewidth=2, color='black', alpha=0.5)
        #ax1.plot([-.42,.42],[11,11], 'r--', linewidth=2, color='black', alpha=0.5)
        #ax1.plot([-.42,.42],[59,59], 'r--', linewidth=2, color='black', alpha=0.5)

        #ax1.fill_between(xx, temp_ave[k]+temp_std[k], temp_ave[k]-temp_std[k], facecolor=colors[k%10], alpha=0.55)
        #ax1.plot(xx, temp_ave[k], color=colors[k%10], linewidth=4, alpha=1.0)
        #plt.title("Cluster: "+str(k+1)+ "  #traces: "+str(n_traces[k]), fontsize=10)
        #ax1.tick_params(axis='x', which='both', labelsize=10)

        #if k<(n_clusters-np.sqrt(n_clusters)): ax1.xaxis.set_ticks([])

    #plt.suptitle(mouse.name + " tot reward pulls: " + str(len(all_traces)), fontsize = 25)

    #plt.show()
    
