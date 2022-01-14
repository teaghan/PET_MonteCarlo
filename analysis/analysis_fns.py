import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares

#plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def plot_sino(scan, vmin=0, vmax=75, slice_indx=[31,100,100],
              spacing=np.array([0.3125, 0.07, 0.07]), draw_lines=True, 
              diff_img=False, savename=None):
    
    # Extent of each axis in physical units
    mm_extent = scan.shape*spacing
    bins = np.arange(0, mm_extent[2], spacing[2])
    projections = np.arange(0, mm_extent[1], spacing[1])
    planes = np.arange(0, mm_extent[0], spacing[0])
    
    # Aspect ratio of images
    aspects = [1.5,0.5,0.2]
    
    # Current indices to plot
    iz, iy, ix = slice_indx
    ix = [ix]
    iy = [iy]
    iz = [iz]
    zmax, ymax, xmax = (np.array(scan.shape)-1)

    # Create subplots
    fig = plt.figure(figsize=(6, 9))#, dpi=300)
    gs1 = gridspec.GridSpec(7, 6)
    ax1 = plt.subplot(gs1[:2,:])
    ax2 = plt.subplot(gs1[2:4,:3])
    ax3 = plt.subplot(gs1[2:4,3:])
    ax4 = plt.subplot(gs1[4,:])
    ax5 = plt.subplot(gs1[5,:])
    ax6 = plt.subplot(gs1[6,:])
    
    for ax in [ax1, ax2, ax3]:
        ax.tick_params(axis='both', which='both',
                       bottom=False, top=False, labelbottom=False,
                       left=False, right=False, labelleft=False)
    
    # Plot scan
    def plot_slices(x_indx, y_indx, z_indx):
        #ax1.set_title('Plane: %i' % (z_indx))
        #ax2.set_title('Projection: %i' % (y_indx))
        #ax3.set_title('Bin: %i' % (x_indx))
        
        img = ax1.imshow(scan[z_indx], extent=[0, xmax, ymax, 0],
               aspect=aspects[0], vmin=vmin, vmax=vmax)
        ax2.imshow(scan[:,y_indx], origin='lower',
                   aspect=aspects[1], vmin=vmin, vmax=vmax)
        ax3.imshow(np.flip(scan[:,:,x_indx], 1), origin='lower',
                   aspect=aspects[2], vmin=vmin, vmax=vmax)

        if diff_img:
            cax = fig.add_axes([0.86, 0.5, .015, 0.4])
            cb = plt.colorbar(img, cax=cax)
            cb.set_label(r'(Simulation - Measurements)')
            #cb.ax.tick_params(labelsize=25,width=1,length=10) 

    def plot_profiles(x_indx, y_indx, z_indx):
        z_prof = scan[:,y_indx,x_indx]
        y_prof = scan[z_indx,:,x_indx]
        x_prof = scan[z_indx,y_indx,:]
        
        ax4.plot(bins, x_prof,
                 c='r')
        ax5.plot(projections, y_prof,
                 c='g')
        ax6.plot(planes, z_prof,
                c='dodgerblue')
        ax4.set_xlim(bins[0],bins[-1])
        ax5.set_xlim(projections[0],projections[-1])
        ax6.set_xlim(planes[0],planes[-1])
        
        ax4.set_ylabel('Counts')
        ax5.set_ylabel('Counts')
        ax6.set_ylabel('Counts')
        ax4.set_xlabel('Bin (mm)')
        ax5.set_xlabel('Projection (${^\circ}$)')
        ax6.set_xlabel('Plane (mm)')
    
    def plot_lines(x_indx, y_indx, z_indx):
        ax1.plot([0,xmax],[y_indx, y_indx], lw=1, c='r')
        ax1.plot([x_indx, x_indx],[0,ymax], lw=1, c='g')
        ax1.set_xlim(0,xmax)
        ax1.set_ylim(ymax,0)
        
        ax2.plot([0,xmax],[z_indx, z_indx], lw=1, c='r')
        ax2.plot([x_indx, x_indx],[0,zmax], lw=1, c='dodgerblue')
        ax2.set_xlim(0,xmax)
        ax2.set_ylim(0,zmax)
        
        ax3.plot([0,ymax],[z_indx, z_indx], lw=1, c='g')
        ax3.plot([y_indx, y_indx],[0,zmax], lw=1, c='dodgerblue')
        ax3.set_xlim(0,ymax)
        ax3.set_ylim(0,zmax)
    
    def onclick(event):
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        if event.inaxes == ax1:
            ix.append(event.xdata)
            iy.append(event.ydata)
        elif event.inaxes == ax2:
            ix.append(event.xdata)
            iz.append(event.ydata)
        elif event.inaxes == ax3:
            iy.append(event.xdata)
            iz.append(event.ydata)
        plot_slices(int(ix[-1]), int(iy[-1]), int(iz[-1]))
        if draw_lines:
            plot_lines(ix[-1], iy[-1], iz[-1])
        plot_profiles(int(ix[-1]), int(iy[-1]), int(iz[-1]))

    connection_id = fig.canvas.mpl_connect('button_press_event', onclick)

    plot_slices(int(ix[-1]), int(iy[-1]), int(iz[-1]))
    if draw_lines:
        plot_lines(ix[-1], iy[-1], iz[-1])
    plot_profiles(int(ix[-1]), int(iy[-1]), int(iz[-1]))
    
    if diff_img:
        fig.subplots_adjust(right=0.85, hspace=0.6)
    else:
        plt.tight_layout()
    
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
def plot_projections(simu_sino, scan_sino, 
                     spacing=np.array([0.3125, 0.07, 0.07]), 
                     lims=None, savename=None):

    # Extent of each axis in physical units
    mm_extent = simu_sino.shape*spacing
    bins = np.arange(0, mm_extent[2], spacing[2])
    projections = np.arange(0, mm_extent[1], spacing[1])
    planes = np.arange(0, mm_extent[0], spacing[0])      

    # Calculate projections onto each axis
    simu_plane_proj = np.sum(simu_sino, axis=(1,2))
    simu_proj_proj = np.sum(simu_sino, axis=(0,2))
    simu_bin_proj = np.sum(simu_sino, axis=(0,1))
    scan_plane_proj = np.sum(scan_sino, axis=(1,2))
    scan_proj_proj = np.sum(scan_sino, axis=(0,2))
    scan_bin_proj = np.sum(scan_sino, axis=(0,1))
    
    # Plotting limits
    if lims is None:
        lims = [(bins[0],bins[-1]), (projections[0],projections[-1]), (0,len(simu_plane_proj))]
    
    # Z-score normalize
    #simu_plane_proj = (simu_plane_proj-np.mean(simu_plane_proj))/np.std(simu_plane_proj)
    #simu_proj_proj = (simu_proj_proj-np.mean(simu_proj_proj))/np.std(simu_proj_proj)
    #simu_bin_proj = (simu_bin_proj-np.mean(simu_bin_proj))/np.std(simu_bin_proj)
    #scan_plane_proj = (scan_plane_proj-np.mean(scan_plane_proj))/np.std(scan_plane_proj)
    #scan_proj_proj = (scan_proj_proj-np.mean(scan_proj_proj))/np.std(scan_proj_proj)
    #scan_bin_proj = (scan_bin_proj-np.mean(scan_bin_proj))/np.std(scan_bin_proj)
    
    fig, axs = plt.subplots(3,1, figsize=(8,6))
    
    fig.suptitle('Sinogram Projections', x=0.43, y=0.93, fontsize=15)

    axs[0].set_xlabel('Bin (mm)', fontsize=12)
    scan_line, = axs[0].plot(bins, scan_bin_proj, lw=1, c='k')
    simu_line, = axs[0].plot(bins, simu_bin_proj, lw=1, c='r')
    axs[0].set_xlim(lims[0])
    axs[0].set_ylim(0,1.1*np.max([scan_bin_proj,simu_bin_proj]))
    
    axs[1].set_xlabel('Projection (${^\circ}$)', fontsize=12)
    axs[1].plot(projections, scan_proj_proj, lw=1, c='k')
    axs[1].plot(projections, simu_proj_proj, lw=1, c='r')
    axs[1].set_xlim(lims[1])
    axs[1].set_ylim(0,1.1*np.max([scan_proj_proj,simu_proj_proj]))
    
    #axs[2].set_xlabel('Plane (mm)', fontsize=12)
    #axs[2].plot(planes, scan_plane_proj, lw=1, c='r')
    #axs[2].plot(planes, simu_plane_proj, lw=1, c='g')
    axs[2].set_xlabel('Plane', fontsize=12)
    axs[2].plot(scan_plane_proj, lw=1, c='k')
    axs[2].plot(simu_plane_proj, lw=1, c='r')
    axs[2].set_xlim(lims[2])
    axs[2].set_ylim(0,1.1*np.max([scan_plane_proj,simu_plane_proj]))
    
    for ax in axs:
        ax.tick_params(labelsize=10)
    
    plt.subplots_adjust(right=0.7, hspace=0.5)
    fig.legend([scan_line, simu_line], ['Measurement', 'Simulation'], loc='center right', fontsize=15)
    
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    
    plt.show()
    
def sum_sino_segs(sino, seg_sizes=np.array([109,97,97,75,75,53,53,31,31])):
    
    ## Calculate number of detection within each plane
    
    # Take sum of all of the segments
    # by centering each segment in the axial direction
    sino_sum = np.zeros((seg_sizes[0], sino.shape[1], sino.shape[2]))
    for i in range(len(seg_sizes)):
        # Collect current segment
        seg = sino[np.sum(seg_sizes[:i]):np.sum(seg_sizes[:i+1])]
        # Add to centre
        sino_sum[int(sino_sum.shape[0]/2 - seg_sizes[i]/2): 
                      int(sino_sum.shape[0]/2 + seg_sizes[i]/2)] += seg  
    # Project onto plane axis
    return sino_sum
    
def calc_sens(sino, seg_sizes=np.array([109,97,97,75,75,53,53,31,31]), dz=2.027, D=200,
             act_conc=12350, t=600, central_frac=3/5):
    '''
    dz (float) = Plane thickness in mm/pixel
    D (float): Diamter of cylinder 200 in mm
    act_conc (float): Initial activity concentration in Bq/mL
    t (float): Duration of scan in seconds
    '''
    
    # Sum the sinogram segments together
    sino_sum = sum_sino_segs(sino, seg_sizes)
    # Project onto plane axis
    sino_sum = np.sum(sino_sum, axis=(1,2)) # Counts
    
    # Turn into cps
    sino_sum /= t # cps
    
    ## Calculate number of decays per plane
    
    # Volume of plane
    plane_vol = np.pi*(D/2)**2*dz # mm^3
    plane_vol *= 0.001 # mL
    # Initial activity of a single plane
    A0 = act_conc*plane_vol
    
    # Select only central fraction
    start_plane = int(np.rint((1-central_frac)/2*sino_sum.shape[0]))
    end_plane = int(np.rint((1-(1-central_frac)/2)*sino_sum.shape[0]))
        
    glob_sens = np.sum(sino_sum[start_plane:end_plane])/(A0*len(sino_sum[start_plane:end_plane]))
    print('The global sensitivity in the central %i%% is %0.2f cps/kBq' % 
          (central_frac*100, 
           glob_sens*1000))
    
    # Return in units of cps/kBq
    return sino_sum[start_plane:end_plane]/A0*1000, glob_sens*1000

def plot_sens(scan_sino_sens, scan_glob_sens, simu_sino_sens, simu_glob_sens, savename=None):

    # Difference relative to the average
    perc_diff = (scan_sino_sens - simu_sino_sens)/(scan_sino_sens + simu_sino_sens)*2
    plane_nums = np.arange((109 - len(scan_sino_sens))/2,  (len(scan_sino_sens) + (109 - len(scan_sino_sens))/2))

    fig = plt.figure()#, dpi=300)
    gs1 = gridspec.GridSpec(4, 1)
    ax1 = plt.subplot(gs1[:3])
    ax2 = plt.subplot(gs1[3], sharex=ax1)

    ax1.plot(plane_nums, scan_sino_sens, 
             label='Measured \n'+r'$S_{\mathrm{global}}=%0.2f$ cps/kBq'%(scan_glob_sens), c='k')
    ax1.plot(plane_nums, simu_sino_sens, 
             label='Simulated \n'+r'$S_{\mathrm{global}}=%0.2f$ cps/kBq'%(simu_glob_sens), c='r')
    ax1.tick_params(labelbottom=False)
    ax1.set_ylabel('Sensitivity (cps/kBq)')
    ax1.set_ylim(5.5, 14.5)
    ax1.set_xlim(plane_nums[0], plane_nums[-1])
    ax1.grid(alpha=0.5)

    ax2.plot(plane_nums, perc_diff*100, c='navy')
    ax2.set_ylim(-0.5, 8)
    ax2.set_xlabel('Plane Number')
    ax2.set_ylabel('Relative\nDifference (%)')
    ax2.grid(alpha=0.5)

    ax1.legend()
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()

def calc_primaries(p, r, D=200, mu=0.0096):
    P = np.zeros_like(r)
    inside = np.abs(r)<=D/2 
    arg = np.sqrt(D**2/4 - r[inside]**2)
    P[inside] = 2*p*np.exp(-2*mu*arg)*arg
    return P

def calc_scatter_pw(r, y_i):
    # Indices of coordinates to evaluate
    r_indices =  np.rint(np.linspace(0, len(r)-1, len(y_i))).astype(int)
    # Piecewise linear interpolation
    return np.interp(r, r[r_indices], y_i)

def calc_scatter(r, A, sigma):
    return A*np.exp(-1/2*(r/sigma)**2)

def apply_PSF(profile, fwhm):
    return gaussian_filter1d(profile, fwhm)

def calc_SF(S, T, r, D):
    sum_indices = np.where((r<=3/10*D) & (r>=-3/10*D))[0]
    return np.sum(S[sum_indices]) / np.sum(T[sum_indices])

def T_resid(params, r, D, mu, fwhm, T_data, resolution):
    # Calculate primary events
    #P = calc_primaries(1.5, r, D, mu)
    P = calc_primaries(params[0], r, D, mu)
    # Calculate scatter events
    S = calc_scatter(r, params[1:])
    # Return residual
    return apply_PSF(P+S, fwhm=fwhm/resolution) - T_data 

def T_resid(params, r, D, mu, fwhm, T_data, resolution, pw=False):
    # Calculate primary events
    #P = calc_primaries(1.5, r, D, mu)
    P = calc_primaries(params[0], r, D, mu)
    # Calculate scatter events
    if pw:
        S = calc_scatter_pw(r, params[1:])
    else:
        S = calc_scatter(r, params[1], params[2])
    # Return residual
    return apply_PSF(P+S, fwhm=fwhm/resolution) - T_data 

def calc_PST(sino, plane, D=200, fwhm=2, resolution=2.005, mu=0.0096, pw=False):
    # Select a profile along the bin direction based on averaging the projections
    T_data = np.mean(sino[plane,:,:],axis=0)
    # Bin positions
    r = (np.arange(0,len(T_data)) - len(T_data)/2 + 0.5) * resolution # mm

    # Initial guesses for the parameters
    if pw:
        p0 = np.ones((6,))
    else:
        p0 = np.ones((3,))

    # Perform fitting
    res_meas = least_squares(T_resid, p0, 
                            args=(r, D, mu, fwhm, T_data, resolution, pw), method='trf')
    # Fitted params
    p_meas = res_meas.x[0]
    if pw:
        yi_meas = res_meas.x[1:]
    else:
        A_meas = res_meas.x[1]
        sigma_meas = res_meas.x[2]

    # Calculate primary events
    P_fit = calc_primaries(p_meas, r, D, mu)

    # Calculate scatter events
    if pw:
        S_fit = calc_scatter_pw(r, yi_meas)
    else:
        S_fit = calc_scatter(r, A_meas, sigma_meas)

    # Calculate profile
    T_fit = apply_PSF(P_fit+S_fit, fwhm=fwhm/resolution)
    
    return P_fit, S_fit, T_fit, T_data, r, calc_SF(S_fit, T_fit, r, D)

def perform_SF_cal(scan_sino, simu_sino, savename=None):
    # Take sum of each segment
    scan_sino_sum = sum_sino_segs(scan_sino, 
                                  seg_sizes=np.array([109,97,97,75,75,53,53,31,31]))
    simu_sino_sum = sum_sino_segs(simu_sino, 
                                  seg_sizes=np.array([109,97,97,75,75,53,53,31,31]))
    # Compute SF for central 60% of planes
    central_frac = 6/10
    # Select only central fraction
    start_plane = int(np.rint((1-central_frac)/2*scan_sino_sum.shape[0]))
    end_plane = int(np.rint((1-(1-central_frac)/2)*scan_sino_sum.shape[0]))

    SF_meas_all = []
    SF_simu_all = []
    for plane in range(start_plane, end_plane):
        P_meas, S_meas, T_meas, T_meas_data, r, SF_meas = calc_PST(scan_sino_sum, 
                                                                   plane=plane, fwhm=4, pw=True)
        P_simu, S_simu, T_simu, T_simu_data, r, SF_simu = calc_PST(simu_sino_sum,
                                                                   plane=plane, fwhm=4, pw=True)

        SF_meas_all.append(SF_meas)
        SF_simu_all.append(SF_simu)

    print('Avg SF for measurements: (%0.2f +- %0.2f)%%' % (100*np.mean(SF_meas_all), 100*np.std(SF_meas_all)))
    print('Avg SF for simulation: (%0.2f +- %0.2f)%%' % (100*np.mean(SF_simu_all), 100*np.std(SF_simu_all)))

    ## Plot a single plane
    # Fit scatter and primaries
    P_meas, S_meas, T_meas, T_meas_data, r, SF_meas = calc_PST(scan_sino_sum, 
                                                               plane=55, fwhm=4, pw=True)
    P_simu, S_simu, T_simu, T_simu_data, r, SF_simu = calc_PST(simu_sino_sum,
                                                               plane=55, fwhm=4, pw=True)

    fontsize=13

    fig, axs = plt.subplots(1,2, figsize=(10,4), sharey=True)

    a = axs[0].scatter(r, T_meas_data, s=5, c='r')#, facecolor='None', edgecolor='k')
    d, = axs[0].plot(r, T_meas, c='k', lw=1.8)
    b, = axs[0].plot(r, P_meas, '--', c='k', lw=1.3)
    c, = axs[0].plot(r, S_meas, '-.', c='k', lw=1.5)

    axs[1].scatter(r, T_simu_data, s=5, c='r')
    axs[1].plot(r, T_simu, c='k', lw=1.8)
    axs[1].plot(r, P_simu, '--', c='k', lw=1.3)
    axs[1].plot(r, S_simu, '-.', c='k', lw=1.5)

    for ax, name, sf, in zip(axs, ['(a) Measured','(b) Simulated'], [SF_meas, SF_simu]):
        ax.set_xlabel('$r$ (mm)', fontsize=fontsize)
        ax.text(140,350, '$SF=$%0.1f%%' % (sf*100), fontsize=fontsize)
        ax.grid(alpha=0.5)
        ax.set_xlim(r[0], r[-1])
        ax.set_title(name, fontsize=fontsize)
    axs[0].set_ylabel('Counts', fontsize=fontsize)
    fig.legend([a,b,c,d], 
               ['Data','Fitted Primaries','Fitted Scatter','Fitted Profile'], 
               framealpha=1, ncol=2, loc='upper center', fontsize=fontsize)
    fig.subplots_adjust(top=0.74, wspace=0.05)

    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
def calc_act(scan, resolution=np.array([2.027, 4.07283, 4.07283])):
    vox_size = np.prod(resolution) # mm^3/voxel
    vox_size *= 0.001 # mL/voxel

    return np.sum(scan)*vox_size/1e6 # kBq

def calc_unif(img, spacing=np.array([2.027, 4.07283, 4.07283]), central_frac=0.6):

    # Indices of pixels
    xlocs, ylocs = np.meshgrid(np.arange(img.shape[1]), np.arange(img.shape[2]))

    # Radius of ROI in pixels
    roi_r = int(15 / spacing[1])

    # Centre locations of ROIs
    xy0 = np.repeat(np.array([img.shape])[:,1:]/2, repeats=5, axis=0)
    xy0[1,1] -= int(75 / spacing[1])
    xy0[2,0] += int(75 / spacing[1])
    xy0[3,1] += int(75 / spacing[1])
    xy0[4,0] -= int(75 / spacing[1])

    # Loop through central fraction of slices
    start_plane = int(np.rint((1-central_frac)/2*img.shape[0]))
    end_plane = int(np.rint((1-(1-central_frac)/2)*img.shape[0]))

    ROI_i = []
    for cur_slice in img[start_plane:end_plane]:

        ROI_j = []
        for i, xy in enumerate(xy0):

            # Find pixels inside circle
            inside = roi_r**2 - (xy[0]-xlocs)**2 - (xy[1]-ylocs)**2
            # Collect pixels inside circle
            ROI_j.append(np.mean(cur_slice[(inside>=0)]))
        ROI_i.append(np.array(ROI_j))

    ROI_i = np.array(ROI_i)

    # IU within slice
    IU_i = (np.max(ROI_i,axis=1) - np.min(ROI_i,axis=1)) / (np.max(ROI_i,axis=1) + np.min(ROI_i,axis=1))

    # IU for a single transverse location
    IU_j = (np.max(ROI_i,axis=0) - np.min(ROI_i,axis=0)) / (np.max(ROI_i,axis=0) + np.min(ROI_i,axis=0))
    
    return IU_i, IU_j

def plot_uniformity(pet_simu, pet_scan, 
                    spacing=np.array([2.027, 4.07283, 4.07283]), 
                    central_frac=0.6,
                    savename=None):

    # Calculate the IU values for both
    IU_i_simu, IU_j_simu = calc_unif(pet_simu, spacing, central_frac)
    IU_i_scan, IU_j_scan = calc_unif(pet_scan, spacing, central_frac)
    
    print('Simulation within slice mean and std: %0.4f, %0.4f' % (np.mean(IU_i_simu), np.std(IU_i_simu)))
    print('Measured within slice mean and std: %0.4f, %0.4f' % (np.mean(IU_i_scan), np.std(IU_i_scan)))
    
    print('Simulation within ROI mean and std: %0.4f, %0.4f' % (np.mean(IU_j_simu), np.std(IU_j_simu)))
    print('Measured within ROI mean and std: %0.4f, %0.4f' % (np.mean(IU_j_scan), np.std(IU_j_scan)))
    
    fig, axs = plt.subplots(1,2, figsize=(10,4), sharey=True)
    
    #fig.suptitle('Image Projections', x=0.43, y=0.93, fontsize=15)

    axs[0].set_title('Within Slice')
    axs[0].set_xlabel('Slice Number', fontsize=12)
    axs[0].set_ylabel('$IU$', fontsize=12)
    scan_line, = axs[0].plot(IU_i_scan, lw=1, c='k')
    simu_line, = axs[0].plot(IU_i_simu, lw=1, c='r')
    
    axs[1].set_title('Within Transverse Location')
    axs[1].set_xlabel('ROI Number', fontsize=12)
    #axs[1].set_ylabel('$IU_{\mathrm{Axial},M}$', fontsize=12)
    axs[1].plot(IU_j_scan, lw=1, c='k')
    axs[1].plot(IU_j_simu, lw=1, c='r')
    
    for ax in axs:
        ax.tick_params(labelsize=10)
        ax.grid(alpha=0.7)
    
    plt.subplots_adjust(top=0.8, wspace=0.1)
    fig.legend([scan_line, simu_line], ['Measurement', 'Simulation'], 
               loc='upper center', fontsize=15, ncol=2)
    
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
def plot_img(scan, vmin=0, vmax=75, slice_indx=[31,100,100],
              spacing=np.array([0.3125, 0.07, 0.07]), draw_lines=True,
              savename=None):
    
    mm_extent = scan.shape*spacing

    aspects = [scan.shape[2]/scan.shape[1],
               mm_extent[0]/mm_extent[1] * scan.shape[1]/scan.shape[0], 
               mm_extent[0]/mm_extent[2] * scan.shape[2]/scan.shape[0]]
    
    iz, iy, ix = slice_indx
    ix = [ix]
    iy = [iy]
    iz = [iz]
    zmax, ymax, xmax = (np.array(scan.shape)-1)

    fig = plt.figure(figsize=(8, 8))#, dpi=300)

    gs1 = gridspec.GridSpec(4, 6)
    ax1 = plt.subplot(gs1[:2,:])
    ax2 = plt.subplot(gs1[2,:3])
    ax3 = plt.subplot(gs1[2,3:])
    
    ax4 = plt.subplot(gs1[3,:2])
    ax5 = plt.subplot(gs1[3,2:4], sharey=ax4)
    ax6 = plt.subplot(gs1[3,4:], sharey=ax4)
    
    for ax in [ax1, ax2, ax3]:
        ax.tick_params(axis='both', which='both',
                       bottom=False, top=False, labelbottom=False,
                       left=False, right=False, labelleft=False)
    for ax in [ax5, ax6]:
        ax.tick_params(axis='y', labelleft=False)
        #ax.grid(True, c='w', alpha=0.2)
    
    # Plot scan
    def plot_slices(x_indx, y_indx, z_indx):
        ax1.imshow(scan[z_indx],# origin='lower',# cmap=plt.cm.bone,
                   aspect=aspects[0], vmin=vmin, vmax=vmax)
        ax2.imshow(scan[:,y_indx], origin='lower',# cmap=plt.cm.bone,
                   aspect=aspects[1], vmin=vmin, vmax=vmax)
        ax3.imshow(np.flip(scan[:,:,x_indx], 1), origin='lower',# cmap=plt.cm.bone,
                   aspect=aspects[2], vmin=vmin, vmax=vmax)
        
    def plot_profiles(x_indx, y_indx, z_indx):
        z_prof = scan[:,y_indx,x_indx]
        y_prof = scan[z_indx,:,x_indx]
        x_prof = scan[z_indx,y_indx,:]
        
        ax4.plot(np.arange(0, mm_extent[2], spacing[2]), x_prof,
                 c='r')
        ax5.plot(np.arange(0, mm_extent[1], spacing[1]), y_prof,
                 c='g')
        ax6.plot(np.arange(0, mm_extent[0], spacing[0]), z_prof,
                c='dodgerblue')
        ax4.set_ylabel('Intesnity')
        ax4.set_xlabel('x (mm)')
        ax5.set_xlabel('y (mm)')
        ax6.set_xlabel('z (mm)')
    
    def plot_lines(x_indx, y_indx, z_indx):
        ax1.plot([0,xmax],[y_indx, y_indx], lw=1, c='r')
        ax1.plot([x_indx, x_indx],[0,ymax], lw=1, c='g')
        ax1.set_xlim(0,xmax)
        ax1.set_ylim(ymax,0)
        
        ax2.plot([0,xmax],[z_indx, z_indx], lw=1, c='r')
        ax2.plot([x_indx, x_indx],[0,zmax], lw=1, c='dodgerblue')
        ax2.set_xlim(0,xmax)
        ax2.set_ylim(0,zmax)
        
        ax3.plot([0,ymax],[z_indx, z_indx], lw=1, c='g')
        ax3.plot([ymax-y_indx, ymax-y_indx],[0,zmax], lw=1, c='dodgerblue')
        ax3.set_xlim(0,ymax)
        ax3.set_ylim(0,zmax)
    
    def onclick(event):
        ax1.clear()
        ax2.clear()
        ax3.clear()
        ax4.clear()
        ax5.clear()
        ax6.clear()
        if event.inaxes == ax1:
            ix.append(event.xdata)
            iy.append(event.ydata)
        elif event.inaxes == ax2:
            ix.append(event.xdata)
            iz.append(event.ydata)
        elif event.inaxes == ax3:
            iy.append(ymax-event.xdata)
            iz.append(event.ydata)
        plot_slices(int(ix[-1]), int(iy[-1]), int(iz[-1]))
        if draw_lines:
            plot_lines(ix[-1], iy[-1], iz[-1])
        plot_profiles(int(ix[-1]), int(iy[-1]), int(iz[-1]))

    connection_id = fig.canvas.mpl_connect('button_press_event', onclick)

    plot_slices(int(ix[-1]), int(iy[-1]), int(iz[-1]))
    if draw_lines:
        plot_lines(ix[-1], iy[-1], iz[-1])
    plot_profiles(int(ix[-1]), int(iy[-1]), int(iz[-1]))
    
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    
    plt.show()
    
## Line Source functions

def determine_max(profile):
    # Select peak value and neighbouring values
    max_indx = np.argmax(profile)
    peak_vals = profile[max_indx-1:max_indx+2]
    # Fit 2nd order polynomial to the points
    coeffs = np.polyfit(np.arange(max_indx-1,max_indx+2), peak_vals, 2)
    # Parabolic function
    parab_fun = np.poly1d(coeffs)
    # Calculate points along this curve
    x_pts = np.linspace(max_indx-1,max_indx+1,1000)
    points = parab_fun(x_pts)
    # Return max value
    return np.max(points), max_indx, x_pts[np.argmax(points)]

def determine_fwhm(profile, max_val, max_indx, xps):
    # Two values on either side of half of the maximum 
    # (on the left side of the maximum)
    left_half = profile[:max_indx+1]
    left_pts = np.array([left_half[(left_half<max_val/2)][-1],
                         left_half[(left_half>=max_val/2)][0]])
    left_indices = np.array([xps[np.where(profile==left_pts[0])[0][0]],
                             xps[np.where(profile==left_pts[1])[0][0]]])

    # Linearly interpolate to find location of half-max on the left side of maximum
    left_loc = np.interp(max_val/2, left_pts, left_indices)

    # Do the same for the right hand side
    right_half = profile[max_indx:]
    right_pts = np.array([right_half[(right_half<max_val/2)][0],
                         right_half[(right_half>=max_val/2)][-1]])
    right_indices = np.array([xps[np.where(profile==right_pts[0])[0][0]],
                             xps[np.where(profile==right_pts[1])[0][0]]])
    right_loc = np.interp(max_val/2, right_pts, right_indices)
    # Return the FWHM in units of pixels
    return right_loc - left_loc

def calc_gauss(x, A, mu, sigma):
    return A*np.exp(-1/2*((x-mu)/sigma)**2)

def gauss_resid(params, x, prof):
    A = params[0]
    mu = params[1]
    sigma = params[2]
    
    return calc_gauss(x, A, mu, sigma) - prof

def fit_gauss(x, prof):
    
    
    # Perform fitting
    res = least_squares(gauss_resid, [np.max(prof), x[np.argmax(prof)], 1], 
                            args=(x, prof), method='trf')
    # Fitted params
    A = res.x[0]
    mu = res.x[1]
    sigma = res.x[2]
    
    return A, mu, sigma
    
def fwhm_supersample(scan, slice_axis, prof_axis, central_frac, threshold=1e3, pixel_size=4.07283,
                    plot=True):
   
        # Remaining axis
    if 0 not in [slice_axis, prof_axis]:
        other_axis = 0
    elif 1 not in [slice_axis, prof_axis]:
        other_axis = 1
    elif 2 not in [slice_axis, prof_axis]:
        other_axis = 2

    # Determine the slices where the line source is present
    line_slices = np.unique(np.where(scan>threshold)[slice_axis])

    # Select middle central_frac of these
    line_slices = line_slices[int(len(line_slices)/2 - len(line_slices)*central_frac/2):
                              int(len(line_slices)/2 + len(line_slices)*central_frac/2)]

    #print('Considering %i slices.' % len(line_slices))
    # Location of where to take profiles
    prof_indx = np.argmax(np.sum(np.take(scan, line_slices, axis=slice_axis), (slice_axis, other_axis)))

    # Find coordinates of the line source and the profiles
    xp = []
    yp = []
    profs = []    
    # Loop through slices and create profile
    for i, slice_indx in enumerate(line_slices):

        cur_slice = np.take(scan, slice_indx, axis=slice_axis)

        # Location of where to take profile
        indx1, indx2 = np.unravel_index(cur_slice.argmax(), cur_slice.shape)

        # Index correctly into slice
        if prof_axis==2:
            prof_indx = indx1
            prof_ax = 0
        elif prof_axis==0:
            prof_indx = indx2
            prof_ax = 1
        else:
            if slice_axis==0:
                prof_indx = indx2
                prof_ax = 1
            else:
                prof_indx = indx1
                prof_ax = 0

        # Select profile
        profs.append(np.take(cur_slice, prof_indx, axis=prof_ax))

        # Determine max of profile
        max_val, max_indx, max_indx_val = determine_max(profs[-1])

        xp.append(i)
        yp.append(max_indx_val)

    # Fit linear curve to the points
    coeffs = np.polyfit(xp, yp, 1)
    # Linear function
    line_fun = np.poly1d(coeffs)
    # Calculate points along this curve
    yp2 = line_fun(xp)

    # Shift each profile and supersample into a single profile
    xps = []
    prof_ss = []
    for y_ref, prof in zip(yp2, profs):
        xps.append(np.arange(len(prof))-y_ref)
        prof_ss.append(prof)
    xps = np.array(xps).flatten()
    prof_ss = np.array(prof_ss).flatten()

    # Take average of overlapping values
    xps_new = np.unique(np.round(xps,3))
    if len(xps_new)!=len(xps):
        prof_ss_new = []
        for x in xps:
            prof_ss_new.append(np.mean(prof_ss[xps==x]))
        prof_ss_new = np.array(prof_ss_new)
    else:
        xps_new = xps
        prof_ss_new = prof_ss

    # Sort
    order = np.argsort(xps_new)
    xps_new = xps_new[order]
    prof_ss_new = prof_ss_new[order]

    # Determine max of profile
    max_val, max_indx, _ = determine_max(prof_ss_new)

    A, mu, sigma = fit_gauss(xps_new, prof_ss_new)

    # Find the FWHM in units of pixels
    #fwhm = determine_fwhm(prof_ss_new, max_val, max_indx, xps_new)
    
    #print('FWHM=',fwhm*pixel_size)
    print('FWHM = %0.3f' % (2.35*sigma*pixel_size))
    if plot:
        plt.figure()
        plt.scatter(xps_new, prof_ss_new, c='r')
        x_gauss = np.linspace(np.min(xps_new), np.max(xps_new),1000)
        plt.plot(x_gauss, calc_gauss(x_gauss, A, mu, sigma), c='k')
        
    return 2.35*sigma*pixel_size

## NEMA phantom analyss

def plot_RCs(RC_max_scan, RC_mean_scan, RC_maxs_simu, RC_means_simu, results_scan, fontsize=15, savename=None):
    # Diameters of individual spheres based on phantom
    diams = np.array([37, 28, 22, 17, 13, 10]) # mm

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    ax = axes.ravel()
    ax[0].set_ylabel(r'$RC_{max}$', fontsize=fontsize)
    ax[1].set_ylabel(r'$RC_{mean}$', fontsize=fontsize)
    ax[0].plot(diams, RC_max_scan, 'o', label='Measurement', c='k')
    #ax[0].plot(diams, RC_max_simu, 'o', label='Simulation', c='r')
    ax[0].errorbar(diams, np.mean(RC_maxs_simu,axis=0), yerr=np.std(RC_maxs_simu,axis=0),
                   fmt='o', capsize=4, label='Simulation', c='r')
    ax[0].fill_between(diams, results_scan.RC_max_lims[:,0], results_scan.RC_max_lims[:,1],
                     color='darkgreen', alpha=0.2, label='EARL range')
    aa, = ax[1].plot(diams, RC_mean_scan, 'o', c='k')
    bb = ax[1].errorbar(diams, np.mean(RC_means_simu,axis=0), yerr=np.std(RC_means_simu,axis=0), 
                     fmt='o', capsize=4, c='r')
    cc = ax[1].fill_between(diams, results_scan.RC_mean_lims[:,0], results_scan.RC_mean_lims[:,1],
                     color='darkgreen', alpha=0.2, label='EARL range')
    for a in ax:
        a.set_xlabel(r'Sphere Diameter (mm)', fontsize=fontsize)
        a.tick_params(labelsize=12)
        #a.legend(loc='lower right', fontsize=14, fancybox=True, framealpha=0.8)
        a.grid(alpha=0.7)
        
    fig.legend([aa,bb,cc], 
           ['Measured','Simulated','EARL range'], 
           framealpha=1, ncol=3, loc='upper center', fontsize=fontsize)
    plt.tight_layout()
    plt.subplots_adjust(top=0.86, wspace=0.1)
    if savename is not None:
        plt.savefig(savename, transparent=True, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
def avg_perc_diff(RC_max_scan, RC_mean_scan, RC_max_simu, RC_mean_simu):
    max_perc_diff = ( np.abs(RC_max_simu - RC_max_scan) / 
                      np.mean([RC_max_simu, RC_max_scan], axis=0)) * 100
    mean_perc_diff = ( np.abs(RC_mean_simu - RC_mean_scan) / 
                      np.mean([RC_mean_simu, RC_mean_scan], axis=0)) * 100
    
    return np.mean([max_perc_diff, mean_perc_diff])