import numpy as np
import matplotlib.pyplot as plt
import os

from ncs import *

file_names= [
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC29.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC21.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC14.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC13.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC9.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC6.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC5.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC3.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC1.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC32.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC31.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC30.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC28.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC27.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC26.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC25.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC24.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC23.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC22.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC20.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC19.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC18.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC17.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC16.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC15.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC12.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC11.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC10.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC8.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC7.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC4.ncs',
'/media/cat/12TB/in_vivo/barak/Front_Beh_2014-04-26_12-27-45/CSC2.ncs',
]

data = loadNcs(file_names[0])

SampleFrequency=float(data[2][12].replace('-SamplingFrequency ' ,''))
print "SampleFrequency = ", SampleFrequency

counter=0
for file_name in file_names:

    print "Loading: ", file_name
    if (os.path.exists(file_name[:-4]+'.npy')==False):
        data = loadNcs(file_name)
        np.save(file_name[:-4], data[0])

    else:
        temp = np.load(file_name[:-4]+'.npy')
        data = temp[0:10000000]
        t = np.linspace(0,len(data),len(data))/SampleFrequency
        plt.plot(t,data+counter*10000)
        
        print "Length of rec: ", len(temp)/SampleFrequency, " sec."
    
    counter+=1

plt.show()
