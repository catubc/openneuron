import sip
sip.setapi('QString', 2) #Sets the qt string to native python strings so can be read without weird stuff

import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from mouse import *
from mouse_lever import *
from cat import *
from rat import *

from analysis import *
from mouse_lever_analysis import *

np.set_printoptions(suppress=True)      #Supress scientific notation printing

from window import Window


def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())


run()
