"""specview.py: A spectrum viewing tool.

To run, `python specview.py`

"""

from pathlib import Path
import pickle
import tkinter as tk
import tkinter.ttk as ttk

import cartopy.crs as ccrs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# Implement the default Matplotlib key bindings.
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
from PIL import Image


# data

data_dir = Path(__file__).parent / "data"
"Path to pickle and image data."

europa_img = Image.open(data_dir / "europa_basemap_with_chaos.png")
with open(data_dir / "fullwaves.pkl", 'rb') as f:
    waves = pickle.load(f)
with open(data_dir / "corr_cubes_101520.pkl", 'rb') as f:
    data = pickle.load(f)


# main window

window = tk.Tk()
window.title("Europa Specview")

class MapPicker(ttk.Frame):
    """Map Picker widget.
    
    This widget displays a map and lets users select positions within the map, using
    Matplotlib under the hood. Users may select regions by clicking and dragging
    their mouse over the map. When a region has been selected, 'select' callback
    functions are invoked with the latitude and longitude data of the selection.

    Args:
      parent: The parent widget to attach this to.
      img: The map to display, as an Image.
      figsize: The size of the matplotlib figure.
    """
    def __init__(self, parent, img, figsize=(10,5)):
        super().__init__(parent)

        fig = Figure(figsize=figsize)
        plate = ccrs.PlateCarree(central_longitude=-180)
        ax = fig.add_subplot(projection=plate)
        ax.imshow(europa_img, transform=ccrs.PlateCarree(central_longitude=0), extent=[0,360,-90,90])
        gl = ax.gridlines(draw_labels=True,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': 16}
        gl.ylabel_style = {'size': 16}

        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0)

        self.status = tk.StringVar()
        self.statusbar = ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.grid(row=1,column=0, sticky='we')

        self.canvas.callbacks.connect('button_press_event', self.on_click)
        self.canvas.callbacks.connect('button_release_event', self.on_release)
        self.canvas.callbacks.connect('motion_notify_event', self.on_motion)

        self.selection = None
        self.clicked = False
        self.fig = fig
        self.ax = ax

        self.select_fns = []

    
    def register_select_callback(self, fn):
        """Register a callback function to invoke when the user selects a region.

        The function will be invoked with the coordinates of the upper left and lower
        right corners of the selected region.
        
        Args:
          fn: A function of four arguments: lat1,lon1 (upper left) and lat2,lon2 (lower right).
        """
        self.select_fns.append(fn)


    def on_motion(self, event):
        if event.inaxes is not None:
            # we are inside the 'map' part of the plot
            self.status.set(f"{event.xdata:.3f}, {event.ydata:.3f}")
            if self.clicked and self.selection is not None:
                # we have a selection -- update it!
                x,y = self.selection.get_xy()
                self.selection.set(width=event.xdata-x, height=event.ydata-y)
                self.canvas.draw()
        else:
            # we're outside of the 'map' part of the plot
            self.status.set("")

    def on_click(self, event):
        if event.inaxes is not None:
            # we are inside of the 'map' part of the plot
            # set 'clicked' to True to indicate that we have an active selection
            self.clicked = True
            if self.selection is not None:
                self.selection.remove()
                self.selection = None
            self.selection = mpatches.Rectangle(
                (event.xdata,event.ydata), 1, 1,
                fill=False,
                edgecolor='red'
            )
            self.ax.add_artist(self.selection)
            self.ax.draw_artist(self.selection)
            self.canvas.draw()
        
    def on_release(self, event):
        if self.clicked:
            # we have completed an active selection -- invoke callbacks
            lat1, lon1 = self.selection.get_xy()
            lat2 = lat1 + self.selection.get_width()
            lon2 = lon1 + self.selection.get_height()
            for fn in self.select_fns:
                fn(lat1, lon1, lat2, lon2)
        self.clicked = False


class SpecViewer(ttk.Frame):
    """Display a spectrum corresponding to some spatial region.
    
    This uses Matplotlib, and the actual plotting is controlled
    by the `plot_spectrum` method.

    Args:
      parent: The parent tkinter widget.
      wave: The array of wavelengths.
      data: The data cube -- a dictionary of them!
      figsize: The size of the Matplotlib figure to produce."""
    def __init__(self, parent, wave, data_cube, figsize=(6,5)):
        super().__init__(parent)
        self.wave = wave
        self.data_cube = data_cube
    
        self.fig = Figure(figsize=figsize)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0)

        self.status = tk.StringVar()
        self.statusbar = ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.E)
        self.statusbar.grid(row=1,column=0, sticky='we')


    def plot_spectrum(self, lat1, lon1, lat2, lon2):
        """Plot the spectrum associated to a specified lat/lon region.
        
        This assumes that the region is rectangular (in the projection), and that
        we are provided the upper left and lower right coordinates of the rectangle.

        Args:
            lat1: The latitude of the upper-left corner.
            lon1: The longitude of the upper-left corner.
            lat2: The latitude of the lower-right corner.
            lon2: The longitude of the lower-right corner.
        """
        self.fig.clear()
        ax = self.fig.add_subplot()
        ax.set_ylabel("Albedo", fontsize=22)
        ax.set_xlabel("Wavelength (microns)", fontsize=22)
        # TODO: Fix this!!!
        ax.plot(self.wave[90:-80], self.data_cube['225'][:,64,64][90:-80], label=f'{lat1:.3f}$^\circ E$, {lon1:.3f}$^\circ$ N')
        ax.legend(fontsize=16)
        # helpful: print the region in the statusbar
        self.status.set(f"selection ({lat1:.3f},{lon1:.3f}) to ({lat2:.3f},{lon2:.3f})")
        self.canvas.draw()

mp = MapPicker(window, europa_img)
sv = SpecViewer(window, waves, data)
mp.register_select_callback(sv.plot_spectrum)
mp.grid(row=0, column=0)
sv.grid(row=0, column=1)

window.mainloop()