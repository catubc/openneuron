import glob, os

from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class FileDialog(QtGui.QFileDialog):
    def __init__(self, *args, **kwargs):
        super(FileDialog, self).__init__(*args, **kwargs)
        self.setOption(QtGui.QFileDialog.DontUseNativeDialog, True)
        self.setFileMode(QtGui.QFileDialog.ExistingFiles)
        #self.setFileMode(QtGui.QFileDialog.DirectoryOnly)

        self.out_files=[]

    def accept(self):
        self.out_files = self.selectedFiles()   #not sure why can't index into this directly, but seems ok to use "out_files"
        self.deleteLater()
        #super(FileDialog, self).accept()
