import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from ..consts import *
from ..utils import flux2mag
from matplotlib.figure import Figure

def plot_curve(curve_file, time_range=(0, 13), **kwargs) -> Figure:
    curve = pd.read_csv(curve_file)

    w1mag = flux2mag(curve.w1flux*u.Jy + 0*u.Jy, FLUX_ZERO_W1)
    w2mag = flux2mag(curve.w2flux*u.Jy + 0*u.Jy, FLUX_ZERO_W2)

    fig = plt.figure(**kwargs)
    ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2, fig=fig)
    ax2 = plt.subplot2grid((3, 3), (2, 0), colspan=2, fig=fig)
    ax3 = plt.subplot2grid((3, 3), (0, 2), rowspan=3, fig=fig)

    min_time, max_time = time_range
    mask = (curve.time >= min_time) & (curve.time < max_time)

    ax1.plot(curve.time[mask], w1mag[mask], linewidth=2, color="tab:blue")
    ax1.scatter(curve.time[mask], w1mag[mask], c=curve.time[mask], zorder=10, s=4)
    ax1.plot(curve.time[mask], w2mag[mask], linewidth=2, color="tab:red")
    ax1.yaxis.set_inverted(True)
    ax1.xaxis.set_visible(False)

    ax2.plot(curve.time[mask], (w1mag-w2mag)[mask], linewidth=2, color="k")


    ax3.plot(w1mag[mask], (w1mag-w2mag)[mask], color="k", linewidth=0.3)
    ax3.scatter(w1mag[mask], (w1mag-w2mag)[mask], c=curve.time[mask], zorder=10, s=10)
    ax3.xaxis.set_inverted(True)

    return fig

# def plot_dust
