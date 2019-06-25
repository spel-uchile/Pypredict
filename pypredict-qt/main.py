from PyQt5 import QtWidgets
from ui.main_window import Ui_MainWindow
import sys

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime

import numpy as np


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        # Configure window
        super(ApplicationWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Add the map widget to the empty 'frame_plot' frame
        self.map = MapCanvas(self, self.ui.frame_plot.width()/100, self.ui.frame_plot.height()/100)
        self.ui.frame_plot.layout().addWidget(self.map)

        # Connect signals and slots
        self.ui.button_text.clicked.connect(self.update_text)
        self.ui.button_plot.clicked.connect(self.map.plot_update)
        self.ui.datetime.dateTimeChanged.connect(self.update_nightshade)

        # Other variables
        self.i = 0

    def update_text(self):
        """Slot to update labels text"""
        # Just and example to update all labels
        labels = [self.ui.label_1, self.ui.label_2, self.ui.label_3,
                  self.ui.label_4, self.ui.label_5, self.ui.label_6,
                  self.ui.label_7, self.ui.label_8, self.ui.label_9]
        self.i += 1
        for i, label in enumerate(labels):
            label.setText("Label {}-{}".format(i, self.i))

    def update_nightshade(self):
        """ Slot to update time"""
        # Read datetime widget as python datetime object
        dt = self.ui.datetime.dateTime().toPyDateTime()
        self.map.plot_nightshade(dt)


class MapCanvas(FigureCanvas):
    """Matplotlib map-plot widget"""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        # Init the plot widget
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        # Save variable to update the map later
        self.axes = self.fig.add_subplot(111)
        self.plt = None
        self.map = None
        self.cs = None

        # Initial plot
        self.plot_map()

    def plot_map(self, date=None):
        """Plot the initial map"""
        self.map = Basemap(ax=self.axes)
        # Plot coastlines, draw label meridians and parallels.
        self.map.drawcoastlines()
        self.map.drawparallels(np.arange(-90, 90, 30), labels=[1, 0, 0, 0])
        self.map.drawmeridians(np.arange(self.map.lonmin, self.map.lonmax + 30, 60), labels=[0, 0, 0, 1])
        # Set background
        self.map.bluemarble(scale=0.25)
        # Set initial nightshade
        date = date if date else datetime.utcnow()
        self.cs = self.map.nightshade(date)
        # Init the x,y plot, empty at this moment @see plot_update
        self.plt, = self.axes.plot([], [], 'g--')
        self.draw()

    def plot_nightshade(self, dt):
        """Update the map nightshade"""
        # Remove previous nightshade (kind of a hack)
        for c in self.cs.collections:
            self.axes.collections.remove(c)
        # Update the nightshade
        self.cs = self.map.nightshade(dt)
        # Update the GUI
        self.draw()
        self.flush_events()

    def plot_update(self):
        """Update the line plot"""
        # New data
        x = np.linspace(-180, 180, 1000)
        a = np.random.randint(0, 90)
        phi = np.random.random(1)*0.5*np.pi
        data = a*np.cos(2*np.pi*np.deg2rad(x)+phi)
        # Update the plot
        self.plt.set_ydata(data)
        self.plt.set_xdata(x)
        # Update the GUI
        self.draw()
        self.flush_events()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    application = ApplicationWindow()
    application.show()
    sys.exit(app.exec_())
