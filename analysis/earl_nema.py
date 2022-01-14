import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import blob_dog
import pandas as pd 

def find_spheres(pet_scan, dz, dy, dx, threshold=0.1, sigma_ratio=1.45):
    '''Use the Difference of Gaussian (DoG) method to locate the spheres.

    Parameters:
    pet_scan (numpy array): 3D PET scan.
    dz, dy, dx (floats): resolution of the PET scan.
    threshold (float, optional): The absolute lower bound for scale space maxima. Local maxima 
                                    smaller than thresh are ignored. Reduce this to detect 
                                    blobs with less intensities.
    sigma_ratio (float, optional): The ratio between the standard deviation of Gaussian Kernels 
                                    used for computing the Difference of Gaussians.

    Returns:
    blobs (numpy array): A 2d array with each row representing 3 coordinate values (z, y, x) 
                            and the radii of the ellipsoid in each direction (rz, ry, rx). 
    '''
    # Volumes of individual spheres based on phantom
    volumes = np.array([26.52, 11.49, 5.57, 2.57, 1.15, 0.52]) * 1000 # mm^3
    # Radii
    radii = (3/(4*np.pi) * volumes )**(1/3) # mm
    # Convert to pixels
    radii_zyx = np.array([radii/dz, radii/dy, radii/dx])
    # Convert to standard deviation units for Gaussian kernel
    sigma_zyx = radii_zyx / np.sqrt(3)

    # Limits for the standard deviation for Gaussian kernel for the blob detection.
    min_sigma = np.min(sigma_zyx, axis=1)
    max_sigma = np.max(sigma_zyx, axis=1)
    
    # Process PET scan by clipping the low valued pixels 
    # and the pixels above 40% of the max. Then we will
    # normalize the pixel values between 0 and 255 to agree
    # with typical greyscale values. 
    # NOTE: this processed version of the scan will NOT be 
    # used to calculate the RCs and COVs.
    image = np.copy(pet_scan)
    image = np.clip(image, 1, 0.4*np.max(image))
    image = ((image - image.min()) * (1/(image.max() - image.min()) * 255)).astype(np.uint8)

    # Find "blobs" using the Difference of Gaussian (DoG) method
    blobs = blob_dog(image, min_sigma=min_sigma, max_sigma=max_sigma, threshold=threshold, 
                     sigma_ratio=sigma_ratio, overlap=1e-5)
    # Compute radius of each blob in each direction (ie. z,y,x)
    blobs[:, 3:] = blobs[:, 3:] * np.sqrt(3) # pixels

    # Only select blobs that are within the middle 50% of scan (in z direction)
    min_slice = np.rint(0.25*pet_scan.shape[0])
    max_slice = np.rint(0.75*pet_scan.shape[0])
    indices = np.where((blobs[:,0]>=min_slice)&(blobs[:,0]<=max_slice))[0]
    blobs = blobs[indices]
        
    return blobs

def create_boolean_ellipsoid(pet_scan, z0, y0, x0, rz, ry, rx):
    ''' Create mask that locates the points within the given blob.'''
    mask = np.zeros_like(pet_scan) 
    for x in range(x0-rx, x0+rx+1):
        for y in range(y0-ry, y0+ry+1):
            for z in range(z0-rz, z0+rz+1):
                check_inside = ((x-x0)/rx)**2 + ((y-y0)/ry)**2 + ((z-z0)/rz)**2
                if check_inside<1:
                    mask[z,y,x] = 1
    return mask.astype(bool)

def create_boolean_square(pet_slice, y0, x0, sidelen):
    ''' Creates mask that locates the points within the given ROI.'''
    mask = np.zeros_like(pet_slice)    
    mask[int(np.rint(y0))-int(np.rint(sidelen/2)):int(np.rint(y0))+int(np.rint(sidelen/2)),
         int(np.rint(x0))-int(np.rint(sidelen/2)):int(np.rint(x0))+int(np.rint(sidelen/2))] = 1
    
    return mask.astype(bool)

def FDG_activity(S0, t, lam=0.0063152315):
    '''Calculate the current FDG activity.

    Parameters:
    S0 (float): Original activity (in kBq/mL) at t=0.
    t (float): Time elapsed (in minutes) since S0 was measured.
    lam (float, optional): Decay rate of FDG. Default is (ln(2)/109.758) min^(-1) (Intro to Nuclear Science, J. Bryan).

    Returns:
    float: Current activity (in kBq/mL).
    '''
    return S0 * np.exp(-lam*t) # kBq/mL

def voi_segment(pet_scan, blobs, dz, dy, dx, S0, t, RC_max_lims, RC_mean_lims):
    '''Identify positions of VOIs used to calculate the recovery coefficients.

    Parameters:
    pet_scan (numpy array): 3D PET scan.
    blobs (numpy array): A 2d array with each row representing 3 coordinate values (z, y, x) 
                            and the radii of the ellipsoid in each direction (rz, ry, rx).
    dz, dy, dx (floats): resolution of the PET scan.
    S0 (float): Original activity (in kBq/mL) at t=0.
    t (float): Time elapsed (in minutes) since S0 was measured.

    Returns:
    spheres (pandas DataFrame): Summary of information on each sphere including recovery coefficients.
    '''
    # Determine the values needed for the calculation of the recovery coefficients
    max_pixels = []
    mean_pixels = []
    VOI_sizes = []
    sphere_radii = []
    
    xy_locs = []
    # Loop through each blob
    for blob in blobs:
        # Collect location and size of current blob
        z0, y0, x0, rz, ry, rx = np.rint(blob).astype(int)
        sphere_radii.append(np.mean([rz*dz, ry*dy, rx*dx]))
        # Add two pixels to radii to insure we collect the entire sphere
        rz += 2
        ry += 2
        rx += 2
        # Create mask based on location and size of blob
        mask = create_boolean_ellipsoid(pet_scan, z0, y0, x0, rz, ry, rx)
        # Investigate only the pixels within this blob
        pet_segment = pet_scan[mask] / 1000 #convert from Bq to kBq   
       
        # Calculate maximum pixel value in kBq/mL
        max_pix = np.max(pet_segment) 
        max_pixels.append(max_pix)
        
        # # Find the VOI based on the half max
        # VOI = pet_segment[pet_segment>max_pix*1000/2]

        # Segment VOI using adapted background method: VOI A50 = (SUVmax-Background)*0.5+Background
        A50 = ((max_pix-2.0)*0.5+2.0) #This line is a placeholder, not necessary due to calc_RCs function below
        VOI_A50 = pet_segment[pet_segment>A50]

        # Calculate the mean pixel value (in kBq/mL) in this VOI
        mean_pix = np.mean(VOI_A50)
        mean_pixels.append(mean_pix)
        # Save locations to place spheres in the correct order
        xy_locs.append([x0,y0])
        
    max_pixels = np.array(max_pixels)
    mean_pixels = np.array(mean_pixels)
    sphere_radii = np.array(sphere_radii)  
    xy_locs = np.array(xy_locs)
    
    # Place correct order based on sphere locations
    order = np.zeros((6,))
    xy_half = np.array(pet_scan.shape)[1:]/2
    spheres_03 = np.argsort(np.abs(xy_locs[:,1]-xy_half[1]))[:2]
    order[0] = spheres_03[np.argmax(xy_locs[spheres_03,0])]
    order[3] = spheres_03[np.argmin(xy_locs[spheres_03,0])]
    spheres_12 = np.where(xy_locs[:,1]>xy_half[1])[0]
    spheres_12 = [loc for loc in spheres_12 if loc not in spheres_03]
    order[1] = spheres_12[np.argmax(xy_locs[spheres_12,0])]
    order[2] = spheres_12[np.argmin(xy_locs[spheres_12,0])]
    spheres_45 = np.where(xy_locs[:,1]<xy_half[1])[0]
    spheres_45 = [loc for loc in spheres_45 if loc not in spheres_03]
    order[4] = spheres_45[np.argmin(xy_locs[spheres_45,0])]
    order[5] = spheres_45[np.argmax(xy_locs[spheres_45,0])]
    order = order.astype(int)
    #order = np.argsort(-VOI_sizes)
    max_pixels = max_pixels[order]
    mean_pixels = mean_pixels[order]
    sphere_radii = sphere_radii[order]
    locations = blobs[order,:3]

    # Calculate current activity
    S_true = FDG_activity(S0, t) # kBq/mL

    # Calculate recovery coefficients
    RC_max = max_pixels/S_true
    RC_mean = mean_pixels/S_true
    
    # Check if within EARL ranges
    RC_max_EARL = []
    for min_max, RC in zip(RC_max_lims, RC_max):
        RC_max_EARL.append(((RC>=min_max[0])&(RC<=min_max[1])))
    RC_max_EARL = np.array(RC_max_EARL)
    RC_mean_EARL = []
    for min_max, RC in zip(RC_mean_lims, RC_mean):
        RC_mean_EARL.append(((RC>=min_max[0])&(RC<=min_max[1])))
    RC_mean_EARL = np.array(RC_mean_EARL)

    # Save data into pandas DataFrame
    data = np.hstack((np.expand_dims(sphere_radii, 1), locations, 
                      np.expand_dims(max_pixels, 1), np.expand_dims(mean_pixels, 1),
                      np.expand_dims(RC_max, 1), np.expand_dims(RC_max_EARL, 1),
                      np.expand_dims(RC_mean, 1), np.expand_dims(RC_mean_EARL, 1)))  
    spheres = pd.DataFrame(data, columns = ['"Blob" Radius (mm)', 'Z_loc (pixels)', 'Y_loc (pixels)', 
                                            'X_loc (pixels)', 'S_max (kBq/mL)', 'S_mean (kBq/mL)', 
                                            'RC_max', 'RC_max EARL Compatible', 'RC_mean', 'RC_mean EARL Compatible'])
    spheres['RC_max EARL Compatible'] = spheres['RC_max EARL Compatible'].astype(bool)
    spheres['RC_mean EARL Compatible'] = spheres['RC_mean EARL Compatible'].astype(bool)
    
    return spheres, order

def calc_COV(pet_scan, spheres, dz, dy, dx):
    '''Calculate the coefficient of variation for the scan.

    Parameters:
    pet_scan (numpy array): 3D PET scan.
    spheres (pandas DataFram): Summary of information on each sphere including recovery coefficients.
    dz, dy, dx (floats): resolution of the PET scan.

    Returns:
    COV (float): Coefficient of variation.
    ROI_z, ROI_y, ROI_x (floats): Center location (in pixels) of each of the 9 ROIs.
    ROI_r: Radius of ROI (in pixels) in xy plane.
    '''
    # Side length of squuare ROIs to give an area of 900mm^2
    ROI_sidelen = 30 / dx # pixels

    # Location of centers of ROIs
    ROI_x = np.array([spheres['X_loc (pixels)'][1] + 60/dx, 
                      spheres['X_loc (pixels)'][2] - 60/dx, 
                      (spheres['X_loc (pixels)'][4]+spheres['X_loc (pixels)'][5])/2])
    ROI_y = np.array([spheres['Y_loc (pixels)'][1] + 10/dy, 
                      spheres['Y_loc (pixels)'][2] + 10/dy, 
                      spheres['Y_loc (pixels)'][4] - 35/dy])
    ROI_z = np.array([spheres['Z_loc (pixels)'][1]-1,
                       spheres['Z_loc (pixels)'][1],
                       spheres['Z_loc (pixels)'][1]+1]).astype(int)
    COVs = []
    bkg = [] #background activity of each ROI
    # Loop through the 3 slices
    for z in ROI_z:
        pet_slice = np.copy(pet_scan[z])
        for y0, x0 in zip(ROI_y, ROI_x):
            # Find ROI
            mask = create_boolean_square(pet_slice, y0, x0, ROI_sidelen)
            ROI = pet_slice[mask] 
            # Calculate COV for this ROI
            COVs.append(np.std(ROI)/np.mean(ROI))
            bkg.append(np.mean(ROI))
    # Calculate average COV
    COV = np.mean(COVs)
    COV_error = np.std(COVs)
    bkg_mean = np.mean(bkg)/1000 #Background activity concentration in kBq/mL
    
    return COV, ROI_z, ROI_y, ROI_x, ROI_sidelen, COV_error, bkg_mean
    
def calc_RCs(pet_scan, blobs, dz, dy, dx, S0, B0, t, order, RC_max_lims, RC_mean_lims, bkgnd_perc=0.5):
    '''Calculate the recovery coefficients.

    Parameters:
    pet_scan (numpy array): 3D PET scan.
    blobs (numpy array): A 2d array with each row representing 3 coordinate values (z, y, x) 
                            and the radii of the ellipsoid in each direction (rz, ry, rx).
    dz, dy, dx (floats): resolution of the PET scan.
    S0 (float): Original activity (in kBq/mL) at t=0.
    t (float): Time elapsed (in minutes) since S0 was measured.

    Returns:
    spheres (pandas DataFrame): Summary of information on each sphere including recovery coefficients.
    '''
    # Determine the values needed for the calculation of the recovery coefficients
    max_pixels = []
    mean_pixels = []
    VOI_sizes = []
    sphere_radii = []
    
    xy_locs = []
    # Loop through each blob
    for blob in blobs:
        # Collect location and size of current blob
        z0, y0, x0, rz, ry, rx = np.rint(blob).astype(int)
        sphere_radii.append(np.mean([rz*dz, ry*dy, rx*dx]))
        # Add two pixels to radii to insure we collect the entire sphere
        rz += 2
        ry += 2
        rx += 2
        # Create mask based on location and size of blob
        mask = create_boolean_ellipsoid(pet_scan, z0, y0, x0, rz, ry, rx)
        # Investigate only the pixels within this blob
        pet_segment = pet_scan[mask] / 1000 #convert from Bq to kBq   
        
        # Calculate maximum pixel value in kBq/mL
        max_pix = np.max(pet_segment) 
        max_pixels.append(max_pix)
        
        # # Find the VOI based on the half max
        # VOI = pet_segment[pet_segment>max_pix*1000/2]

        # Segment VOI using adapted background method: VOI A50 = (SUVmax-Background)*0.5+Background
        A50 = ((max_pix-B0)*bkgnd_perc+B0) #isocontour at 50% SUV max, adapted for background
        VOI_A50 = pet_segment[pet_segment>A50]
        
        # Calculate the mean pixel value (in kBq/mL) in this VOI
        mean_pix = np.mean(VOI_A50)
        mean_pixels.append(mean_pix)
        # Save locations to place spheres in the correct order
        xy_locs.append([x0,y0])
        
    max_pixels = np.array(max_pixels)
    mean_pixels = np.array(mean_pixels)
    sphere_radii = np.array(sphere_radii)  
    xy_locs = np.array(xy_locs)

    max_pixels = max_pixels[order]
    mean_pixels = mean_pixels[order]
    sphere_radii = sphere_radii[order]
    locations = blobs[order,:3]

    # Calculate current activity
    S_true = FDG_activity(S0, t) # kBq/mL

    # Calculate recovery coefficients
    RC_max = max_pixels/S_true
    RC_mean = mean_pixels/S_true
    
    # Check if within EARL ranges
    RC_max_EARL = []
    for min_max, RC in zip(RC_max_lims, RC_max):
        RC_max_EARL.append(((RC>=min_max[0])&(RC<=min_max[1])))
    RC_max_EARL = np.array(RC_max_EARL)
    RC_mean_EARL = []
    for min_max, RC in zip(RC_mean_lims, RC_mean):
        RC_mean_EARL.append(((RC>=min_max[0])&(RC<=min_max[1])))
    RC_mean_EARL = np.array(RC_mean_EARL)

    # Save data into pandas DataFrame
    data = np.hstack((np.expand_dims(sphere_radii, 1), locations, 
                      np.expand_dims(max_pixels, 1), np.expand_dims(mean_pixels, 1),
                      np.expand_dims(RC_max, 1), np.expand_dims(RC_max_EARL, 1),
                      np.expand_dims(RC_mean, 1), np.expand_dims(RC_mean_EARL, 1)))  
    results = pd.DataFrame(data, columns = ['"Blob" Radius (mm)', 'Z_loc (pixels)', 'Y_loc (pixels)', 
                                            'X_loc (pixels)', 'S_max (kBq/mL)', 'S_mean (kBq/mL)', 
                                            'RC_max', 'RC_max EARL Compatible', 'RC_mean', 'RC_mean EARL Compatible'])
    results['RC_max EARL Compatible'] = results['RC_max EARL Compatible'].astype(bool)
    results['RC_mean EARL Compatible'] = results['RC_mean EARL Compatible'].astype(bool)
    
    return results

class NemaRC:
    '''EANM analysis of a PET scan.

    Parameters:
    dicom_dir (str): Path to dicom files.
    S0 (float): Original activity concentration in spheres (in kBq/mL) at t=0.
    t (float): Time elapsed (in minutes) since S0 was measured.
    threshold (float, optional): The absolute lower bound for scale space maxima. Local maxima 
                                    smaller than thresh are ignored. Reduce this to detect 
                                    blobs with less intensities. Don't reduce below 0.085.
    sigma_ratio (float, optional): The ratio between the standard deviation of Gaussian Kernels 
                                    used for computing the Difference of Gaussians.
    '''
    def __init__(self, dicom_dir, S0, t, threshold=0.085, sigma_ratio=1.45, bkgnd_perc=0.5,
                 pet_scan=None, resolution=None):
        #self.S0 = dicom_dir
        self.S0 = S0
        self.t = t

        if pet_scan is None:
            # Load PET data
            self.pet_scan, self.dz, self.dy, self.dx = load_PET(dicom_dir)
        else:
            self.pet_scan = pet_scan
            self.dz, self.dy, self.dx = resolution

        # Locate spheres
        print('Locating spheres using the Difference of Gaussian (DoG) method...')        
        self.blobs = find_spheres(self.pet_scan, self.dz, self.dy, self.dx, threshold, sigma_ratio)

        # Segment VOIs for spheres (will also be used for positioning COV ROIs)
        print('Segmenting VOI for each sphere...')
        # EARL 1 RC range
        self.RC_max_lims = np.array([[0.95,1.16],[0.91,1.13],[0.83,1.09],
                                [0.73,1.01],[0.59,0.85],[0.31,0.49]])
        self.RC_mean_lims = np.array([[0.76,0.89],[0.72,0.85],[0.63,0.78],
                                [0.57,0.73],[0.44,0.60],[0.27,0.38]])
        self.spheres, self.order = voi_segment(self.pet_scan, self.blobs, self.dz, self.dy, self.dx, 
                                             self.S0, self.t, self.RC_max_lims, self.RC_mean_lims)
        

        # Calculate coefficient of variation
        print('Calculating coefficient of variation...')
        self.COV, self.ROI_z, self.ROI_y, self.ROI_x, self.ROI_sidelen, self.COV_error, self.bkg = calc_COV(self.pet_scan, self.spheres, 
                                                                            self.dz, self.dy, self.dx)
        
        # Calculate recovery coefficients
        print('Calculating recovery coefficients for each sphere...')
        self.results = calc_RCs(self.pet_scan, self.blobs, self.dz, self.dy, self.dx, 
                        self.S0, self.bkg, self.t, self.order, self.RC_max_lims, self.RC_mean_lims,
                               bkgnd_perc=bkgnd_perc)

        self.RC_mean = self.results['RC_mean']
        self.RC_max = self.results['RC_max']
        # Calculate mean contrast recovery (MCR) for SUV mean and SUV max
        self.MCR_mean = np.mean(self.results['RC_mean'])
        self.MCR_max = np.mean(self.results['RC_max'])
        
        # Calculate curvature of RC curves (RC Mean and RC Max)
        i = 1
        mean_rmse = []
        max_rmse = []
        while i<6:
            mean_error = (self.RC_mean[0]-self.RC_mean[i])**2
            max_error = (self.RC_max[0]-self.RC_max[i])**2
            mean_rmse.append(mean_error)
            max_rmse.append(max_error)
            i+=1
        self.mean_curv = np.sqrt(np.mean(mean_rmse))
        self.max_curv = np.sqrt(np.mean(max_rmse))

        print('Procedure complete.')


    def plot_spheres(self):
        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        ax = axes.ravel()
        ax[0].set_title('PET', fontsize=20)
        ax[0].imshow(self.pet_scan[int(self.blobs[0,0])], vmax=0.5*np.max(self.pet_scan), cmap='bone_r')
        ax[1].set_title('Difference of \nGaussian Spheres', fontsize=20)
        ax[1].imshow(self.pet_scan[int(self.blobs[0,0])], vmax=0.5*np.max(self.pet_scan), cmap='bone_r')
        # Plot cirlces
        for blob in self.blobs:
            y, x = blob[1:3]
            r = blob[4]
            c = plt.Circle((x, y), r, color='red', linewidth=2, fill=False)
            ax[1].add_patch(c)
        ax[0].set_axis_off()
        ax[1].set_axis_off()
        plt.tight_layout()
        plt.show()
        
    def plot_RCs(self):
        # Diameters of individual spheres based on phantom
        diams = np.array([37, 28, 22, 17, 13, 10]) # mm

        fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
        ax = axes.ravel()
        fig.text(0.45, 0.9, r'RC Curves', fontsize=20)
        ax[0].set_ylabel(r'$RC_{max}$', fontsize=15)
        ax[1].set_ylabel(r'$RC_{mean}$', fontsize=15)
        ax[0].plot(diams, self.spheres['RC_max'].values, 'o', label='Measurements', c='navy')
        ax[0].fill_between(diams, self.RC_max_lims[:,0], self.RC_max_lims[:,1],
                         color='darkgreen', alpha=0.2, label='EARL range')
        ax[1].plot(diams, self.spheres['RC_mean'].values, 'o', label='Measurements', c='maroon')
        ax[1].fill_between(diams, self.RC_mean_lims[:,0], self.RC_mean_lims[:,1],
                         color='darkgreen', alpha=0.2, label='EARL range')
        for a in ax:
            a.set_xlabel(r'Sphere Diameter (mm)', fontsize=15)
            a.tick_params(labelsize=12)
            a.legend(loc='lower right', fontsize=14, fancybox=True, framealpha=0.8)
            a.grid(alpha=0.7)
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.3, top=0.88)
        plt.show()
        
    def plot_ROIs(self):
        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(8, 4))
        ax = axes.ravel()
        ax[0].set_title('PET', fontsize=20)
        ax[0].imshow(self.pet_scan[int(self.blobs[0,0])], vmax=0.5*np.max(self.pet_scan), cmap='bone_r')
        ax[1].set_title('Background ROIs', fontsize=20)
        ax[1].imshow(self.pet_scan[int(self.blobs[0,0])], vmax=0.5*np.max(self.pet_scan), cmap='bone_r')
        # Plot cube
        for x0, y0 in zip(self.ROI_x, self.ROI_y):

            # Plot each side of cube
            ax[1].plot([int(np.rint(x0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(x0))-int(np.rint(self.ROI_sidelen/2))],
                     [int(np.rint(y0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(y0))+int(np.rint(self.ROI_sidelen/2))], 
                     c='r', linewidth=2)
            ax[1].plot([int(np.rint(x0))+int(np.rint(self.ROI_sidelen/2)), int(np.rint(x0))+int(np.rint(self.ROI_sidelen/2))],
                     [int(np.rint(y0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(y0))+int(np.rint(self.ROI_sidelen/2))], 
                     c='r', linewidth=2)
            ax[1].plot([int(np.rint(x0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(x0))+int(np.rint(self.ROI_sidelen/2))],
                     [int(np.rint(y0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(y0))-int(np.rint(self.ROI_sidelen/2))], 
                     c='r', linewidth=2)
            ax[1].plot([int(np.rint(x0))-int(np.rint(self.ROI_sidelen/2)), int(np.rint(x0))+int(np.rint(self.ROI_sidelen/2))],
                     [int(np.rint(y0))+int(np.rint(self.ROI_sidelen/2)), int(np.rint(y0))+int(np.rint(self.ROI_sidelen/2))], 
                     c='r', linewidth=2)

        ax[0].set_axis_off()
        ax[1].set_axis_off()

        plt.tight_layout()
        plt.show()
    
    def check_COV(self):
        if self.COV>0.15:
            print('The COV for this image is %0.1f%%, which is not below the maximum of 15%%.'%(self.COV*100)+
                  ' This does not meet the EANM guidelines.'%(self.COV*100))
        else:
            print('The COV for this image is %0.1f%%, which is below the maximum of 15%%.'%(self.COV*100)+
                  ' This meets the EANM guidelines.')
