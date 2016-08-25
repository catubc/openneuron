#Sergey Gratiy's CSD function

import numpy as np
import scipy.linalg as la
from math import *
import pylab as plt

def inverse_tikhonov(A,SNR):

    '''
        calcuate matrix W which where CSD_est=W*LFP
        R = resolution matrix R=W*A
    '''
    At=A.T
    AAt = np.dot(A,At)
    nz = A.shape[0]
    np.mean(np.diag(AAt),axis=0)
    AAt_reg = AAt+1/SNR**2*np.mean(np.diag(AAt),axis=0)*np.eye(nz)
    #print AAt_reg[0]
    invAAt_reg = la.inv(AAt_reg)
    W = np.dot(At,invAAt_reg)
    R = np.dot(W,A)

    return [W,R]


def forward_operator(zs,b,sigma):
    
    '''
        caclulate matrix A where LFP=A*CSD
    '''    
    K = np.zeros((len(zs),len(zs)))

    dz = zs[1]-zs[0]
    for iz,z in enumerate(zs):
        for izp, zp in enumerate(zs):
            K[iz,izp] = sqrt(((z-zp)/b)**2 + 1) - abs((z-zp)/b)
    
    
    A = b*dz/(2*sigma)*K
      
    return [A,K]  
 
 
def show_operators(A,A0,W,W0,R,R0):
    
    
    fig100 = plt.figure(100)
    plt.subplot(2,3,1);
    VspanA = np.array([-1, 1])*abs(A).max()
    plt.imshow(A, vmin = VspanA[0], vmax = VspanA[1],cmap='jet'); plt.colorbar();
    plt.title('forward operator')
    
    plt.subplot(2,3,4);
    VspanA0 = np.array([-1, 1])*abs(A0).max()
    plt.imshow(A0, vmin = VspanA0[0], vmax = VspanA0[1],cmap='jet'); plt.colorbar();
    plt.title('forward operator 0')
    
    
    plt.subplot(2,3,2);
    VspanW = np.array([-1, 1])*abs(W).max()
    plt.imshow(W, vmin = VspanW[0], vmax = VspanW[1],cmap='jet'); plt.colorbar();
    plt.title('Inverse operator')
    
    plt.subplot(2,3,3);
    VspanR = np.array([-1, 1])*abs(R).max()
    plt.imshow(R, vmin = VspanR[0], vmax = VspanR[1],cmap='jet'); plt.colorbar();
    plt.title('Resolution matrix')
    
    
    plt.subplot(2,3,5);
    VspanW0 = np.array([-1, 1])*abs(W0).max()
    plt.imshow(W0, vmin = VspanW0[0], vmax = VspanW0[1],cmap='jet'); plt.colorbar();
    plt.title('Inverse operator 0')
    
    plt.subplot(2,3,6);
    VspanR0 = np.array([-1, 1])*abs(R0).max()
    plt.imshow(R0, vmin = VspanR0[0], vmax = VspanR0[1],cmap='jet'); plt.colorbar();
    plt.title('Resolution matrix 0')
    


    fig101 = plt.figure(101)
    nw=W.shape[0]

    plt.subplot(2,1,1); 
    plt.plot(W[:,[0,nw/2,nw-1]],'.-')
    plt.title('Inverse operator for different channels')

    plt.subplot(2,1,2); 
    plt.plot(W0[:,[0,nw/2,nw-1]],'.-')
        
    return [fig100,fig101]
     
    
