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
import matplotlib.ticker as mticker

from PIL import Image
import numpy as np
from dataclasses import dataclass
import time
from itertools import cycle
from copy import copy

@dataclass
class Coordinates:
    lat1: float
    lon1I: float
    lat2: float
    lon2I: float
    lon1W: float
    lon2W: float
    lon1EW: float
    lon2EW: float
        
    def __init__(self,  lon1I: float, lon2I:float, lat1: float, lat2: float):
        self.lat1 = lat1
        self.lat2 = lat2
        
        self.lon1I = lon1I
        self.lon2I = lon2I
        
        self.posW_transform()
        self.posEW_transform()
        
    def posW_transform(self):
        """ Translate lats from image coords to [360, 0] (W)"""
        def transW(l):
            return 180. - l
        
        if transW(self.lon1I) <= transW(self.lon2I):
            self.lon1W = transW(self.lon1I)
            self.lon2W = transW(self.lon2I)
        else:
            self.lon1W = transW(self.lon2I)
            self.lon2W = transW(self.lon1I)
       
    def posEW_transform(self):
        """ Translate lats from image coords to [-180, 180] (EW)"""
        def transEW(l):
            return -1*np.sign(l)*(180. - np.abs(l))
        
        if transEW(self.lon1I) <= transEW(self.lon2I):
            self.lon1EW = transEW(self.lon1I)
            self.lon2EW = transEW(self.lon2I)
        else:
            self.lon1EW = transEW(self.lon2I)
            self.lon2EW = transEW(self.lon1I)
            
    def within_region(self, lonEW, latEW) -> bool:
        if not ((self.lon1EW <= lonEW) and (lonEW <= self.lon2EW)):
            return False
        if (self.lat1 <= latEW) and (latEW <= self.lat2):
            return True
        elif (self.lat2 <= latEW) and (latEW <= self.lat1):
            return True
        return False


# data
@dataclass
class Pixel:
    lats: np.array #latitude range (lower left, upper right)
    lons: np.array #longitude range (lower left, upper right)
    wave: np.array # avelengths
    intensity: np.array # spectra (intensity)
    cen_dist: float #distance from center of disk when observed
    scaled_intensity: np.array #allow arbitrary 2 wavelength scaling to be plotted
        
    def __init__(self, lats: np.array, lons: np.array, wave:np.array,
                intensity: np.array, cen_dist: float):
        self.lats = lats
        self.lons = lons
        self.wave = wave
        self.intensity = intensity
        self.cen_dist = cen_dist
        self.scaled_intensity = intensity

# Load data:
# TO DO: Lat/lons are swapped, need to fix this!
data_dir = Path(__file__).parent / "data"
"Path to pickle and image data."

europa_img = Image.open(data_dir / "europa_basemap_with_chaos.png")
with open(data_dir / "fullwaves.pkl", 'rb') as f:
    waves = pickle.load(f)
with open(data_dir / "corr_cubes_101520.pkl", 'rb') as f:
    data = pickle.load(f)
ps = np.load(str(data_dir)+"/pixels.npy", allow_pickle=True) 

# To Do: more elegant solution for intensity scaling?
pixels = [Pixel(p.lats, p.lons, p.wave, p.intensity, p.cen_dist) for p in ps]


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
        def west_positive_formatted(longitude, num_format='g'):
            fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    
            if longitude >= 0:
                longitude = 360 - longitude
            else:
                -1*longitude
            #Transform 180W to 180E to 0-360 W
            return fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere='W',
                             degree=u'\u00B0')
        
        fig = self.fig
        plate = ccrs.PlateCarree(central_longitude=180)
        ax = fig.add_subplot(projection=plate)
        ax.imshow(self.img, transform=ccrs.PlateCarree(central_longitude=0))#, extent=[-360, -0, -90,90])
        gl = ax.gridlines(draw_labels=True,
                          linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xformatter = mticker.FuncFormatter(lambda v, pos:
                                             west_positive_formatted(v))

        gl.xlabel_style = {'size': 16}
        gl.ylabel_style = {'size': 16}
        ax.set_global()

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
            pos = Coordinates(event.xdata, event.xdata, event.ydata, event.ydata)
            self.status.set(f"{pos.lon1W:.3f}, {pos.lat1:.3f}")
            #self.status.set(f"{event.xdata:.3f}, {event.ydata:.3f}")
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
            lon1, lat1 = self.selection.get_xy()
            lon2 = lon1 + self.selection.get_width()
            lat2 = lat1 + self.selection.get_height()
            self.coords = Coordinates(lon1, lon2, lat1, lat2)
            for fn in self.select_fns:
                fn(self.coords)

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
        

def valid_coord_region(lon1,lon2,lat1,lat2):
    "Check that a coordinate region is valid."
    return (lon1 <= lon2) and (lat1 <= lat2)


# def correct_pixels(pixels):
#     "Last minute corrections to pixel formatting."
    
#     corrected = []
    
#     for p in pixels:
#         if p.intensity[100:-80].count() == 0:
#             continue
#         # todo: confirm that what we fix here (lat/lon flipped) matches with how the data was generated
#         corrected.append(Pixel(np.array(p.lons), np.array(p.lats), p.wave, p.intensity, p.cen_dist))
        
#     return corrected


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
        # TO DO: fix underlying data in pixels and remove correct function
        #self.pixels = correct_pixels(pixels)
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
        
        # Widget to scale the intensity
        #tk.Label(parent, text="Lambda 1:").grid(column=1, row=0)
        #tk.Label(parent, text="Lambda 2:").grid(column=1, row=0)

        w1 = tk.Entry(parent, text="Lambda 1: ")
        w2 = tk.Entry(parent, text="Lambda 2: ")
        
        w1.grid(row=0, column=1, sticky='nw', pady=40)
        w2.grid(row=0, column=1, sticky='nw', pady=80)

        self.scaleintbutton = tk.Button(parent, text='Scale Spectra',
                                       command = lambda: self.scale_pixels(w1.get(), w2.get()))
        self.scaleintbutton.grid(row=0, column=1, sticky='nw', pady=120)

    def clear_plot(self):
        self.fig.clear()
        self.canvas.draw()

    def save_plot(self):
        t = time.time() 
        filenm = str(data_dir)+'/plots/plot_'+str(t)+'.png'
        self.fig.savefig(filenm)

    def match_pixels(self, coords):
        matches = []
        for p in self.pixels:
            if coords.within_region(p.lons[0], p.lats[0]) and coords.within_region(p.lons[1], p.lats[1]):
                matches.append(p)
                
        if not matches:
            print("No matches found :(")
            
        return matches
    
    def scale_pixels(self, wave1, wave2):
        """ Scale pixel values."""
        wave1, wave2 = float(wave1), float(wave2)
        print(wave1, wave2)
        if wave1 is None:
            for p in self.pixels:
                p.scaled_intensity = copy(p.intensity)
        else: #scale based on wave values
            waves = copy(self.pixels[0].wave)

            idx1 = np.argmin(np.abs(waves - wave1))
            idx2 = np.argmin(np.abs(waves - wave2))

            w1val = np.mean([p.intensity[idx1] for p in pixels])
            w2val = np.mean([p.intensity[idx2] for p in pixels])

            for p in self.pixels:
                scale = (w2val-w1val)/(p.intensity[idx2]-p.intensity[idx1])
                p.scaled_intensity = p.intensity*(scale)

                sub = w1val - p.scaled_intensity[idx1]
                p.scaled_intensity = p.scaled_intensity + sub



    def plot_spectrum(self, coords):
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
        matches = self.match_pixels(coords)
        
        if len(matches) == 1:
            ax.plot(self.wave[100:-80], matches[0].intensity[100:-80], 
                    label=f'{coords.lon1W:.2f}:{coords.lon2W:.2f}$^\circ W$, {coords.lat1:.2f}:{coords.lat2:.2f}$^\circ$ N')
        elif len(matches) > 1:
            # Median combine if multiple data pixels
            pixdata = [m.scaled_intensity[100:-80] for m in matches]
            #np.set_printoptions(threshold=np.inf)
            #print([m.count() for m in pixdata])
            
            # TO  DO: actually bin spatial pixels
            binned = np.nanmedian(pixdata, axis=0)
            ax.plot(self.wave[100:-80], binned,
                    label=f'{coords.lon1W:.2f}:{coords.lon2W:.2f}$^\circ W$, {coords.lat1:.2f}:{coords.lat2:.2f}$^\circ$ N')

        self.fig.subplots_adjust(bottom=0.15, right=0.7)
        ax.legend(bbox_to_anchor=(1.5, 1), fontsize=10)

        if ax.get_ylim()[0] < 0:
            ax.set_ylim(0, ax.get_ylim()[1])
        if ax.get_ylim()[1] > 0.6:
            ax.set_ylim(ax.get_ylim()[0], 0.6)
        #ax.set_ylim(0, 0.6)
        # helpful: print the region in the statusbar
        self.status.set(f"selection ({coords.lat1:.3f},{coords.lon1W:.3f}) to ({coords.lat2:.3f},{coords.lon2W:.3f})")
        self.canvas.draw()

mp = MapPicker(window, europa_img)
sv = SpecViewer(window, waves, data, pixels)
mp.register_select_callback(sv.plot_spectrum)
mp.grid(row=0, column=0)
sv.grid(row=1, column=0)

window.mainloop()