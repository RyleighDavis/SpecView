from pathlib import Path
import tkinter as tk
import tkinter.ttk as ttk

import cartopy.crs as ccrs
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# Implement the default Matplotlib key bindings.
from matplotlib.figure import Figure
import matplotlib.patches as mpatches
from PIL import Image


europa_img = Image.open(Path(__file__).parent / "data" / "europa_basemap_with_chaos.png")

window = tk.Tk()
window.title("Europa Specview")

class MapPicker(ttk.Frame):
    def __init__(self, master, img):
        super().__init__(master)

        fig = Figure(figsize=(12,6))
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
        self.statusbar = ttk.Label(self, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.E)
        self.statusbar.grid(row=1,column=0, sticky='we')

        self.canvas.callbacks.connect('button_press_event', self.on_click)
        self.canvas.callbacks.connect('button_release_event', self.on_release)
        self.canvas.callbacks.connect('motion_notify_event', self.on_motion)

        self.selection = None
        self.clicked = False
        self.fig = fig
        self.ax = ax


    def on_motion(self, event):
        if event.inaxes is not None:
            self.status.set(f"{event.xdata:.3f}, {event.ydata:.3f}")
            if self.clicked and self.selection is not None:
                x,y = self.selection.get_xy()
                self.selection.set(width=event.xdata-x, height=event.ydata-y)
                self.canvas.draw()
        else:
            self.status.set("")

    def on_click(self, event):
        if event.inaxes is not None:
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
        else:
            print('Clicked ouside axes bounds but inside plot window') 

    def on_release(self, event):
        self.clicked = False


mp = MapPicker(window, europa_img)
mp.grid(row=0, column=0)

window.mainloop()