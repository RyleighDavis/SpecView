"""specview.py: A spectrum viewing tool.

To run, `python specview.py`

"""

from pathlib import Path
import pickle
import tkinter as tk
import tkinter.ttk as ttk

import cartopy.crs as ccrs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# Implement the default Matplotlib key bindings.
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from dataclasses import dataclass
import time
from itertools import cycle


# data
@dataclass
class Pixel:
    lats: np.array #latitude range (lower left, upper right)
    lons: np.array #longitude range (lower left, upper right)
    wave: np.array # avelengths
    intensity: np.array # spectra (intensity)
    cen_dist: float #distance from center of disk when observed
        
    def __init__(self, lats: np.array, lons: np.array, wave:np.array,
                intensity: np.array, cen_dist: float):
        self.lats = lats
        self.lats = lats
        self.lons = lons
        self.wave = wave
        self.intensity = intensity
        self.cen_dist = cen_dist

# Load data
data_dir = Path(__file__).parent / "data"
"Path to pickle and image data."

europa_img = Image.open(data_dir / "europa_basemap_with_chaos.png")
with open(data_dir / "fullwaves.pkl", 'rb') as f:
    waves = pickle.load(f)
with open(data_dir / "corr_cubes_101520.pkl", 'rb') as f:
    data = pickle.load(f)
pixels = np.load(str(data_dir)+"/pixels.npy", allow_pickle=True)


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
    def __init__(self, parent, img, figsize=(12,3.5)):
        super().__init__(parent)

        self.colors = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        self.img = img

        fig = Figure(figsize=figsize)

        self.fig = fig
        self.ax = self.draw_map()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0)

        self.status = tk.StringVar()
        self.statusbar = ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.W)
        self.statusbar.grid(row=1,column=0, sticky='we')

        self.canvas.callbacks.connect('button_press_event', self.on_click)
        self.canvas.callbacks.connect('button_release_event', self.on_release)
        self.canvas.callbacks.connect('motion_notify_event', self.on_motion)

        self.clearselection = tk.Button(parent, text='Clear ROIs', command=self.on_remove)
        self.clearselection.grid(row=0, column=1, sticky='nw')

        self.selection = None
        self.clicked = False

        self.select_fns = []

    def draw_map(self):
            fig = self.fig
            plate = ccrs.PlateCarree(central_longitude=0)
            ax = fig.add_subplot(projection=plate)
            ax.imshow(self.img, transform=ccrs.PlateCarree(central_longitude=0), extent=[-180,180,-90,90])
            gl = ax.gridlines(draw_labels=True,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')
            gl.xlabel_style = {'size': 16}
            gl.ylabel_style = {'size': 16}

            return ax

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
                #self.selection.remove()
                self.selection = None
            self.selection = mpatches.Rectangle(
                (event.xdata,event.ydata), 1, 1,
                fill=False,
                edgecolor=next(self.colors),
                linewidth=2
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
                if (lat1 < lat2) & (lon1 < lon2):
                    fn(lat1, lon1, lat2, lon2)
                elif (lat2 < lat1) & (lon1 < lon2):
                    fn(lat2, lon1, lat1, lon2)
                elif (lat1 < lat2) & (lon2 < lon1):
                    fn(lat1, lon2, lat2, lon1)
                else:
                    fn(lat2, lon2, lat1, lon1)

        self.clicked = False

    def on_remove(self):
        # Remove all past selections
        self.clicked = False

        self.selection.remove()
        self.selection = None

        self.fig.clear()
        self.ax = self.draw_map()
        self.canvas.draw()

        self.colors = cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
        



class SpecViewer(ttk.Frame):
    """Display a spectrum corresponding to some spatial region.
    
    This uses Matplotlib, and the actual plotting is controlled
    by the `plot_spectrum` method.

    Args:
      parent: The parent tkinter widget.
      wave: The array of wavelengths.
      data: The data cube -- a dictionary of them!
      figsize: The size of the Matplotlib figure to produce."""
    def __init__(self, parent, wave, data_cube, pixels, figsize=(12,4)):
        super().__init__(parent)
        self.wave = wave
        self.data_cube = data_cube
        self.pixels = pixels
    
        self.fig = Figure(figsize=figsize)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0)

        toolbarFrame = ttk.Frame(master=self)
        toolbarFrame.grid(row=0,column=0, sticky='sw')
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbarFrame)
        self.toolbar.update()
        #self.toolbar.grid(row=0, column=0)

        self.status = tk.StringVar()
        self.statusbar = ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.E)
        self.statusbar.grid(row=1,column=0, sticky='we')

        self.clearbutton = tk.Button(parent, text='Clear Plot', command=self.clear_plot)
        self.clearbutton.grid(row=1, column=1, sticky='nw')

        self.savebutton = tk.Button(parent, text='Save Plot', command=self.save_plot)
        self.savebutton.grid(row=1, column=1, sticky='nw', pady=40)

    def clear_plot(self):
        self.fig.clear()
        self.canvas.draw()

    def save_plot(self):
        t = time.time() 
        filenm = str(data_dir)+'/plots/plot_'+str(t)+'.png'
        self.fig.savefig(filenm)


    def match_pixels(self, lat1, lon1, lat2, lon2):
        matches = []
        for p in self.pixels:
            if (((p.lats[0] >= lat1) & (lat2 >= p.lats[1])) & ((p.lons[0] >= lon1) & (p.lons[1] <= lon2))):
                if np.sign(p.lats[0]) == np.sign(p.lats[1]):
                    #then not wrapping around 180/-180 boundary
                    matches.append(p)
                else:
                    # TO DO: Add logic to deal with +/- 180 boundary
                    print('No matches found!')
            
        return matches


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
        ax = self.fig.add_subplot()
        ax.set_ylabel("Albedo", fontsize=16)
        ax.set_xlabel("Wavelength (microns)", fontsize=16)

        #Plot relevant spectra!
        matches = self.match_pixels(lat1, lon1, lat2, lon2)
        
        if len(matches) == 1:
            ax.plot(self.wave[100:-80], matches[0].intensity[100:-80], 
                    label=f'{lat1:.2f}:{lat2:.2f}$^\circ E$, {lon1:.2f}:{lon2:.2f}$^\circ$ N')
        elif len(matches) > 1:
            # Median combine if multiple data pixels
            med = np.nanmedian([m.intensity[100:-80] for m in matches], axis=0)
            ax.plot(self.wave[100:-80], med, 
                    label=f'{lat1:.2f}:{lat2:.2f}$^\circ E$, {lon1:.2f}:{lon2:.2f}$^\circ$ N')

        self.fig.subplots_adjust(bottom=0.15, right=0.7)
        ax.legend(bbox_to_anchor=(1.5, 1), fontsize=10)

        if ax.get_ylim()[0] < 0:
            ax.set_ylim(0, ax.get_ylim()[1])
        if ax.get_ylim()[1] > 0.55:
            ax.set_ylim(ax.get_ylim()[0], 0.55)
        #ax.set_ylim(0, 0.6)
        # helpful: print the region in the statusbar
        self.status.set(f"selection ({lat1:.3f},{lon1:.3f}) to ({lat2:.3f},{lon2:.3f})")
        self.canvas.draw()

mp = MapPicker(window, europa_img)
sv = SpecViewer(window, waves, data, pixels)
mp.register_select_callback(sv.plot_spectrum)
mp.grid(row=0, column=0)
sv.grid(row=1, column=0)

window.mainloop()