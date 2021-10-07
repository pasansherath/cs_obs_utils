import numpy as np
from obspy.signal.polarization import polarization_analysis
from obspy import read, Trace, Stream
from datetime import datetime
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


def calc_hpol_angle(h1_h2_window):
    '''
    Function to calculate the horizontal polarisation angle
    between H1 and H2 components of an OBS of a given time 
    window using single value decomposition.
    Inputs:
        h1_h2_window - H1 and H2 amplitudes of the sliding window (numpy array)
    Returns:
        hpol_anlge - polarisation angle in the horizontal plane (float)
    '''

    ##### covariance matrix of h1_h2_window #####
    cov_mat = np.cov(h1_h2_window)

    ##### calculate eigenvectors and eigenvalues of the cov_mat using svd #####
    eigenvec, eigenval, v = np.linalg.svd(cov_mat)

    ##### calculate horizontal polarisation angle using eigenvectors #####
    hpol_angle = np.rad2deg(np.arctan(eigenvec[1][0]/eigenvec[0][0]))

    return hpol_angle

def hodogram_patch(h1_window_data, h2_window_data, offset, time):
    '''
    Function to draw patches of particle motion between H1 and H2 components 
    of an OBS.
    Inputs: 
        h1_window_data - sliding window data of H1 component (numpy array)
        h2_window_data - sliding window data of H2 component (numpy array)
        offset - offset at trace (float)
        time - central time to start plotting particle motion (float)

    Returns:
        patch - matplotlib patch of particle motion of the sliding window
    '''
    
    ##### define an array to save vertices of particle motion #####
    verts = np.zeros((len(h1_window_data)+1,2))
    verts[0,0] = offset
    verts[0,1] = time

    ##### assign particle motion coordinates to vertices array #####
    for k in range(len(h1_window_data)):
        verts [k+1,0] = offset + h2_window_data[k]/5
        verts [k+1,1] = time + h1_window_data[k]/5

    ##### assign identification values to vertices #####
    codes = np.array([1])
    for l in range(np.shape(verts)[0]-1):
        codes = np.append(codes, 2)

    ##### create matplotlib path #####
    path = Path(verts, codes)

    ##### assign patch #####
    patch = patches.PathPatch(path, facecolor='None', lw=0.5)

    return patch

def hodogram_polarisation(h1_segy_file, h2_segy_file, channel_range, time_length, window_size, step):

    '''
    Function to read in SEGY files and plot polarisation angles
    and hodograms of H1 and H2 components of an OBS
    Inputs:
        h1_segy_file, h2_segy_file - segy files for H1 and H2 components of an OBS (segy format files)
        channel_range -  channel range for first arrivign direct water wave (numpy array)
        time_length - maximum time length to calculate horizontal polarisation angle and particle motion (float)
        window_size - sliding window time length (float)
        step - sliding window time step
    
    Output:
        Plots hodogram and horizontal polarisation angle with offset and time
    '''
    ##### read segy #####
    h1 = read(h1_segy_file, unpack_trace_headers = True)
    h2 = read(h2_segy_file, unpack_trace_headers = True)   
    
    ##### define empty array to store offsets #####
    offsets = np.array([])

    ##### define a figure object to plot the figure #####
    fig = plt.figure(figsize=(10, 6)) 

    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[1, 1])
    ##### ax0 for hodogram, ax1 for eigen decomposition #####
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    # ax2 = fig.add_subplot(gs[2]) # for colorbar

    ##### define an array to store the angle between H1 and shot line from eigen decomposition of covariance matrix #####
    hpol_angles = np.zeros((int(time_length/step),len(channel_range)))

    for i in range(len(channel_range)):
        channel = int(channel_range[i])
        print ("Reading trace {:d} of {:d}".format(i+1, len(channel_range)))

        ##### Loop through identical traces in the vertical and two horizontal components #####
        for j, trace_h1 in enumerate(h1):
            ##### normalize
            h1[j].normalize()
            h2[j].normalize()

            ##### trim to desired length #####
            h1[j].trim(h1[j].stats.starttime, h1[j].stats.starttime+time_length)
            h2[j].trim(h2[j].stats.starttime, h2[j].stats.starttime+time_length)

            ##### get trace headers and assign trace stats #####
            trace_h1_hdr = trace_h1.stats.segy.trace_header
            offset = trace_h1_hdr.get('distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group')/1000
            trace_h1_ensemble_number = int(trace_h1_hdr.get('ensemble_number'))
            trace_h1.stats.channel = str(trace_h1_ensemble_number)

            ##### create overlapping time windows in the traces #####
            if trace_h1_ensemble_number == channel:
                trace_h1_windows = h1[j].slide(window_length = window_size, step = step, offset =0.0, include_partial_windows =False,  nearest_sample = True)
                trace_h2_windows = h2[j].slide(window_length = window_size, step = step, offset =0.0, include_partial_windows =False,  nearest_sample = True)

                ##### append to offsets #####
                offsets = np.append(offsets, offset)
                ##### define t=0 indicate begining of time for each trace #####
                time = 0

                ##### loop through the overlapping time windows of the two traces of horizontal components #####
                for k, (window_h1, window_h2) in enumerate(zip(trace_h1_windows, trace_h2_windows)):
                    
                    ##### copy the windowed data to a new trace object due to limitations in ObsPy #####
                    trace_h1_window_copy = window_h1.copy()
                    trace_h2_window_copy = window_h2.copy()

                    ##### find the mid point of the time window #####
                    dt = (trace_h1_window_copy.stats.endtime - trace_h1_window_copy.stats.starttime)/2

                    ##### define and assign the windowed data into an array #####
                    h1_h2_window = np.zeros((2, len(trace_h1_window_copy.data)))
                    h1_h2_window[0, :] = trace_h1_window_copy.data
                    h1_h2_window[1, :] = trace_h2_window_copy.data

                    ##### increment time #####
                    time = time + dt

                    ##### get horizontal polarisation angle #####
                    hpol_angles[k,i] = calc_hpol_angle(h1_h2_window)

                    ##### get particle motion and add patch #####
                    patch = hodogram_patch(trace_h1_window_copy.data, trace_h2_window_copy.data, offset, time)
                    ax0.add_patch(patch)


    ##### hodogram #####
    plt.sca(ax0)
    ax0.set_xlabel("OFFSET (km)")
    ax0.set_ylabel("TIME (s)")
    ax0.set_xlim(np.min(offsets), np.max(offsets))
    ax0.set_ylim(0, time_length)
    ax0.set_aspect('equal')
    ax0.set_title("Hodogram", weight = 'bold')

    ##### polarisation_angle #####
    plt.sca(ax1)
    hpol = ax1.imshow(hpol_angles, cmap='rainbow', extent = (np.min(offsets), np.max(offsets), time_length, 0), aspect='equal')
    hpol.set_clim([-90, 90])
    ax1.set_xlabel("OFFSET (km)")
    ax1.set_ylabel("TIME (s)")
    ax1.set_xlim(np.min(offsets), np.max(offsets))
    ax1.set_ylim(0, time_length)
    ax1.set_aspect('equal')
    ax1.set_title("Horizontal Polarisation Angle", weight = 'bold')

    divider = make_axes_locatable(ax1)
    cax1 = divider.append_axes("right", size="5%", pad=0.1)
    cbar = plt.colorbar(hpol, cax=cax1, orientation='vertical', ticks=np.arange(-90,100, 30), label="Polarisation Angle")
    plt.savefig("H1_shotline_angle.png", dpi=600, bbox_inches="tight")
    plt.show()

##### segy file names H1 and H2
h1_segy_file="H1.segy"
h2_segy_file="H2.segy"

##### direct water arrival channel range
direct_arrival_channel_begin = 1001
direct_arrival_channel_end = 1176

##### sliding window parameters
time_length = 5
window_size = 0.1
step = 0.05

##### read every 5th trace from the channel range
direct_arr_channel_range = np.arange(direct_arrival_channel_begin, direct_arrival_channel_end+1, 1)
direct_arr_channel_range = direct_arr_channel_range[::5]

##### do the work
hodogram_polarisation(h1_segy_file, h2_segy_file, direct_arr_channel_range, time_length, window_size, step)