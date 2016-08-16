#Visualization toolbox 
#Created for Allen Institute - Modelling, Analysis and Theory biophysically detailed network models
#Authors: Catalin Mitelut, Sergey Gratiy
#MIT License

import numpy as np

from PyQt4 import QtGui, QtCore, QtOpenGL   #QtOpenGL installs via "sudo apt-get install python-qt4-gl"
from OpenGL.GL import *                     #OpenGL installs in ubuntu via "sudo pip install PyOpenGL PyOpenGL_accelerate"
from OpenGL.GLU import *                    #only used in one location... not clear if necessary


RED = 255, 0, 0
ORANGE = 255, 127, 0
YELLOW = 255, 255, 0
GREEN = 0, 255, 0
CYAN = 0, 255, 255
LIGHTBLUE = 0, 127, 255
BLUE = 0, 0, 255
VIOLET = 127, 0, 255
MAGENTA = 255, 0, 255  
GREY = 85, 85, 85
WHITE = 255, 255, 255
DARK_GREY = 30, 30, 30
PURPLE = 154, 44, 209
CMAP = np.array([RED, ORANGE, YELLOW, GREEN, CYAN, LIGHTBLUE, BLUE, VIOLET, MAGENTA,
                 GREY, WHITE, PURPLE], dtype=np.uint8)

COLOUR_DICT = {'pyramid1': 1, 'pyramid2': 2, 'pyramid3': 3, 'pyramid4': 4, 'pyramid21': 5,
               'pyramid26': 6, 'pyramid28': 7, 'bc1': 0, 'bc2': 8, 'bc3': 1, 'bc4': 9, 'bc5': 10, 'bc6': 11,
	       'pyramid23': 1, 'pyramid5': 2, 'pyramid6': 3, 
	       'basket23': 5, 'basket4': 6, 'basket5': 7, 'basket6': 8}



class GLWindow(QtGui.QWidget):
    
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.glWidget = GLWidget(parent=self)
        
        self.glWidget.bkgr=0
        
        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.glWidget)
        
        self.setLayout(mainLayout)
        self.setWindowTitle(self.tr("OpenGL test"))

        #Default: don't plot segments or somas
        self.glWidget.plot_seg=0
        self.glWidget.plot_soma=0
        self.glWidget.plot_frame=1
        
        self.show()

class GLWidget(QtOpenGL.QGLWidget):
    
    def __init__(self, parent=None):
        QtOpenGL.QGLWidget.__init__(self, parent)
        self.lastPos = QtCore.QPoint()
        self.focus = np.float32([0, 0, 0]) # init camera focus
        self.axes = 'both' # display both mini and focal xyz axes by default
        

        format = QtOpenGL.QGLFormat()
        
        #format.setVersion(3, 0) # not available in PyQt 4.7.4
        # set to color index mode, unsupported in OpenGL >= 3.1, don't know how to load
        # GL_ARB_compatibility extension, and for now, can't force OpenGL 3.0 mode.
        # Gives "QGLContext::makeCurrent(): Cannot make invalid context current." error:
        #format.setRgba(False)
        
        format.setDoubleBuffer(True) # works fine
        self.setFormat(format)
        #QtOpenGL.QGLFormat.setDefaultFormat(format)
        
        '''
        c = QtGui.qRgb
        cmap = [c(255, 0, 0), c(0, 255, 0), c(0, 0, 255), c(255, 255, 0), c(255, 0, 255)]
        colormap = QtOpenGL.QGLColormap()
        colormap.setEntries(cmap)
        self.setColormap(colormap)
        '''
        
        

    def minimumSizeHint(self):
        return QtCore.QSize(50, 50)

    def sizeHint(self):
        return QtCore.QSize(1000, 800)

    def initializeGL(self):
        if self.bkgr==1:
            glClearColor(0.0, 0.0, 0.0, 1.0) # Black / White toggle switch
        else:
            glClearColor(1.0, 1.0, 1.0, 1.0)

        glClearDepth(10.0) # same as default
        glEnable(GL_DEPTH_TEST) # display points according to occlusion, not order of plotting
        #GL.glEnable(GL.GL_POINT_SMOOTH) # doesn't seem to work right, proper way to antialiase?
        #GL.glEnable(GL.GL_LINE_SMOOTH) # works better
        #GL.glPointSize(1.5) # truncs to the nearest pixel if antialiasing is off
        #GL.glShadeModel(GL.GL_FLAT)
        #GL.glEnable(GL.GL_CULL_FACE) # only useful for solids
        glTranslate(0, 750, -3000) # init camera distance from origin



        #GL.glEnable(GL_LIGHTING)
        #GL.glEnable(GL_LIGHTING)
        #GL.glEnable(GL_LIGHT0)
        
        #lightpos = [1.,1.,1., 0.]
        #GL.glLightfv(GL_LIGHT0, GL_POSITION, lightpos)

        #white = [0.8, 0.8, 0.8, 1.0]
        #cyan = [0., .8, .8, 1.]
        #GL.glMaterialfv(GL.GL_FRONT, GL_DIFFUSE, cyan)
        #GL.glMaterialfv(GL.GL_FRONT, GL_SPECULAR, white)
        #shininess = [50]
        #GL.glMaterialfv(GL.GL_FRONT, GL_SHININESS, shininess)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        # Don't load identity matrix. Do all transforms in place against current matrix
        # and take advantage of OpenGL's state-machineness.
        # Sure, you might get round-off error over time, but who cares? If you wanna return
        # to a specific focal point on 'f', that's when you really need to first load the
        # identity matrix before doing the transforms
        #GL.glLoadIdentity() # loads identity matrix into top of matrix stack

        # viewing transform for camera: where placed, where pointed, which way is up:
        #GLU.gluLookAt()
        #GL.glScale() # modelling transformation, lets you stretch your objects

        glEnableClientState(GL_COLOR_ARRAY);
        glEnableClientState(GL_VERTEX_ARRAY);
        #GL.glEnable(GL_LINE_SMOOTH);
        #GL.glHint(GL_LINE_SMOOTH_HINT,  GL_NICEST);
        #GL.glColorPointerub(cell.colors_3dpts) # unsigned byte, ie uint8
        #GL.glVertexPointerf(cell.points_3dpts) # float32
        #GL.glDrawArrays(GL.GL_LINES, 0, cell.npoints_3dpts)
        
        #for cid in range(0, window.n_cells, window.skip):
        #GL.glVertexPointerf(cells.points_seg) # float32

        #************** Plots somas ************
        #if self.plot_soma==1:
        glColorPointerub(self.soma_colours) # unsigned byte, ie uint8
        glVertexPointerf(self.soma_points) # float32
        glDrawArrays(GL_TRIANGLES, 0, len(self.soma_points))

        if self.axes: # paint xyz axes
            glClear(GL_DEPTH_BUFFER_BIT) # make axes paint on top of data points
            #if self.axes in ['both', 'mini']:
            self.paint_mini_axes()
            #if self.axes in ['both', 'focal']:
            self.paint_focal_axes()

        # might consider using buffer objects for even more speed (less unnecessary vertex
        # data from ram to vram, I think). Apparently, buffer objects don't work with
        # color arrays?

        #GL.glFlush() # forces drawing to begin, only makes difference for client-server?
        self.swapBuffers() # doesn't seem to be necessary, even though I'm in double-buffered
                           # mode with the back buffer for RGB sid encoding, but do it anyway
                           # for completeness

        # print the modelview matrix
        #print(self.MV)

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        # fov (deg) controls amount of perspective, and as a side effect initial apparent size
        gluPerspective(45, width/height, 0.01, 1000000) # fov, aspect, nearz & farz
                                                           # clip planes
        glMatrixMode(GL_MODELVIEW)
    

    def paint_mini_axes(self):
        """Paint mini xyz axes in bottom left of widget"""
        w, h = self.width(), self.height()
        vt = self.getTranslation() # this is in eye coordinates
        glViewport(0, 0, w//4, h//4) # mini viewport at bottom left of widget
        self.setTranslation((-1, -.5, -3)) # draw in center of this mini viewport
        self.paint_axes()
        self.setTranslation(vt) # restore translation vector to MV matrix
        glViewport(0, 0, w, h) # restore full viewport

    def paint_focal_axes(self):
        """Paint xyz axes proportional in size to sigma, at focus"""
        glTranslate(*self.focus) # translate to focus
        #self.paint_axes(self.sigma)
        glTranslate(*-self.focus) # translate back

    def update_focal_axes(self):
        """Called every time sigma is changed in main spyke window"""
        #self.update_sigma()
        self.updateGL()

    def paint_axes(self, l=1):
        """Paint axes at origin, with lines of length l"""
        glBegin(GL_LINES)
        glColor3f(1, 0, 0) # red x axis
        glVertex3f(0, 0, 0)
        glVertex3f(l, 0, 0)
        glColor3f(0, 1, 0) # green y axis
        glVertex3f(0, 0, 0)
        glVertex3f(0, l, 0)
        glColor3f(0, 0, 1) # blue z axis
        glVertex3f(0, 0, 0)
        glVertex3f(0, 0, l)
        glEnd()

    def get_MV(self):
        """Return modelview matrix"""
        return glGetDoublev(GL_MODELVIEW_MATRIX) # I think this acts like a copy

    def set_MV(self, MV):
        glLoadMatrixd(MV)

    MV = property(get_MV, set_MV)

    # modelview matrix is column major, so we work on columns instead of rows
    def getViewRight(self):
        """View right vector: 1st col of modelview matrix"""
        return self.MV[:3, 0]

    def getViewUp(self):
        """View up vector: 2nd col of modelview matrix"""
        return self.MV[:3, 1]

    def getViewNormal(self):
        """View normal vector: 3rd col of modelview matrix"""
        return self.MV[:3, 2]

    def getTranslation(self):
        """Translation vector: 4th row of modelview matrix"""
        return self.MV[3, :3]

    def setTranslation(self, vt):
        """Translation vector: 4th row of modelview matrix"""
        MV = self.MV
        MV[3, :3] = vt
        self.MV = MV

    def getDistance(self):
        v = self.getTranslation()
        #return np.sqrt((v**2).sum()) # from data origin
        return np.sqrt(((v-self.focus)**2).sum()) # from focus

    def pan(self, dx, dy):
        """Translate along view right and view up vectors"""
        d = self.getDistance()
        vr = self.getViewRight()
        vr *= dx*d
        glTranslate(vr[0], vr[1], vr[2])
        vu = self.getViewUp()
        vu *= dy*d
        glTranslate(vu[0], vu[1], vu[2])

    def zoom(self, dr):
        """Translate along view normal vector"""
        d = self.getDistance()
        vn = self.getViewNormal()
        vn *= dr*d
        glTranslate(vn[0], vn[1], vn[2])

    def pitch(self, dangle): # aka elevation
        """Rotate around view right vector"""
        vr = self.getViewRight()
        glTranslate(*self.focus)
        glRotate(dangle, *vr)
        glTranslate(*-self.focus)

    def yaw(self, dangle): # aka azimuth
        """Rotate around view up vector"""
        vu = self.getViewUp()
        glTranslate(*self.focus)
        glRotate(dangle, *vu)
        glTranslate(*-self.focus)

    def roll(self, dangle):
        """Rotate around view normal vector"""
        vn = self.getViewNormal()
        glTranslate(*self.focus)
        glRotate(dangle, *vn)
        glTranslate(*-self.focus)

    def panTo(self, p=None):
        """Translate along view right and view up vectors such that data point p is
        centered in the viewport. Not entirely sure why or how this works, figured
        it out using guess and test"""
        if p == None:
            p = self.focus
        MV = self.MV
        vr = self.getViewRight()
        vu = self.getViewUp()
        p = -p
        x = np.dot(p, vr) # dot product
        y = np.dot(p, vu)
        MV[3, :2] = x, y # set first two entries of 4th row to x, y
        self.MV = MV

    def mousePressEvent(self, event):
        self.lastPos = QtCore.QPoint(event.pos())

    def mouseMoveEvent(self, event):
        buttons = event.buttons()
        modifiers = event.modifiers()
        dx = event.x() - self.lastPos.x()
        dy = event.y() - self.lastPos.y()

        if buttons == QtCore.Qt.LeftButton:
            if modifiers == QtCore.Qt.ControlModifier:
                self.roll(-0.5*dx - 0.5*dy)
            elif modifiers == QtCore.Qt.ShiftModifier:
                self.pan(dx/600., -dy/600.) # qt viewport y axis points down
            else:
                self.yaw(0.5*dx)
                self.pitch(0.5*dy)
        elif buttons == QtCore.Qt.RightButton:
            self.zoom(-dy/500.) # qt viewport y axis points down

        self.updateGL()
        self.lastPos = QtCore.QPoint(event.pos())

    def wheelEvent(self, event):
        self.zoom(event.delta() / 1000.)
        self.updateGL()

    #PICKER FUNCTIONS
    
    def get_sids(self):
        return self._sids

    def set_sids(self, sids):
        """Set up rgbsids array for later use in self.pick()"""
        self._sids = sids
        print sids
        # encode sids in RGB
        r = sids // 256**2
        rem = sids % 256**2 # remainder
        g = rem // 256
        b = rem % 256
        self.rgbsids = np.zeros((self.npoints, 3), dtype=np.uint8)
        self.rgbsids[:, 0] = r
        self.rgbsids[:, 1] = g
        self.rgbsids[:, 2] = b
    
    sids = property(get_sids, set_sids) #Odd way of updating a functions'

    def pick(self, x, y, pb=10, multiple=False):
        """Return sid of point at window coords x, y (bottom left origin),
        or first or multiple sids that fall within a square 2*pb+1 pix on a side,
        centered on x, y. pb is the pixel border to include around x, y"""
        width = self.size().width()
        height = self.size().height()
        #print('coords: %d, %d' % (x, y))
        # constrain to within border 1 pix smaller than widget, for glReadPixels call
        if not (pb <= x < width-pb and pb <= y < height-pb): # cursor out of range
            return
        if self.npoints > 2**24-2: # the last one is the full white background used as a no hit
            raise OverflowError("Can't pick from more than 2**24-2 sids")
        # draw encoded RGB values to back buffer
        #GL.glDrawBuffer(GL_BACK) # defaults to back
        GL.glClearColor(1.0, 1.0, 1.0, 1.0) # highest possible RGB means no hit
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glEnableClientState(GL.GL_COLOR_ARRAY)
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glColorPointerub(self.rgbsids) # unsigned byte, ie uint8
        GL.glVertexPointerf(self.points) # float32
        GL.glDrawArrays(GL.GL_POINTS, 0, self.npoints) # to back buffer
        GL.glClearColor(0.0, 0.0, 0.0, 1.0) # restore to default black
        # grab back buffer:
        #GL.glReadBuffer(GL.GL_BACK) # defaults to back
        # find rgb at or around cursor coords, decode sid:
        backbuffer = GL.glReadPixels(x=x-pb, y=y-pb, width=2*pb+1, height=2*pb+1,
                                     format=GL.GL_RGB, type=GL.GL_UNSIGNED_BYTE,
                                     array=None, outputType=None)
        # NOTE: outputType kwarg above must be set to something other than str to ensure
        # that an array is returned, instead of a string of bytes
        if (backbuffer == 255).all(): # no hit
            return
        if not multiple:
            sid = self.decodeRGB(backbuffer[pb, pb]) # check center of backbuffer
            if sid != None:
                #print('hit at exact cursor pos')
                return sid # hit at exact cursor position
        # 2D array with nonzero entries at hits:
        hitpix = (backbuffer != [255, 255, 255]).sum(axis=2)
        if not multiple:
            ri = np.where(hitpix.ravel())[0][0] # get ravelled index of first hit
            i, j = np.unravel_index(ri, dims=hitpix.shape) # unravel to 2D index
            #print('hit at %d, %d' % (i, j))
            return self.decodeRGB(backbuffer[i, j]) # should be a valid sid
        ijs = zip(*np.where(hitpix)) # list of ij tuples
        
        out_array = np.asarray([ self.decodeRGB(backbuffer[i, j]) for i, j in ijs ])
        
        print out_array
        
        return out_array


