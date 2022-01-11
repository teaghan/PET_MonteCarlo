import numpy as np
import os

'''Ziegler 2015: According to the NEMA NU 2–2007 standard [3], image quality parameters of PET scanners are obtained by measuring a specific International Electrotechnical Commission (IEC) 61675–1 emission phantom [13] (NEMA image quality phantom) (Fig. 1a, PTW, Freiburg, Germany). This image quality phantom mimics the shape of an upper human body and is built of acrylic glass material. It comprises 6 hollow glass spheres (inner diameters 37, 28, 22, 17, 13, and 10 mm) which can be inserted into the large phantom compartment. Additionally a cylindrical insert containing styrofoam with an average density of 0.3 ± 0.1 g ml−1 (simulates patient lung tissue [3] (μlung-insert ~ 0.026 cm−1) and is positioned in the center of the phantom. The inner volume may vary between NEMA IQ phantoms. The volume of the tested phantom (PTW, Freiburg, Germany, Fig. 1a) was measured to be 9.5 L ± 1 % when the spheres and lung cylinder are inserted. The phantom housing has a thickness of approximately 3 mm along the phantom body, and 10 mm (in few parts 20 mm) at the lids at both ends of the phantom. The glass material (μ ~ 0.118 cm−1) of the spheres has a thickness of around 1 mm.'''

def create_phantom(phantom_size, dz, dy, dx, ph_type, spheres=True, sphere_z_off=0.8, sphere_y_off=1.4, sphere_x_off=-0.5):
    '''
    Create digital form of NEMA phantom.
    
    inputs:
        phantom_size (tuple): (z_length, y_length, x_length) in mm.
        dz, dy, dx (floats): resolution in mm/pixel for z, y, and x dimensions.
        
    returns:
        phantom (np.array): 3D array containing NEMA phantom.
    '''
    
    print('Creating NEMA phantom in a %s array, ' % (str(phantom_size)) + 
          'with a resolution of (%0.4f, %0.4f, %0.4f)mm/pixel' % (dz, dy, dx))
    
    phantom = np.zeros(phantom_size)

    # Cylindrical insert diameter
    d = 51 # mm
    # Diameters of spheres
    diams = np.array([37,28,22,17,13,10]) # mm
    # Distance from centre of spheres to centre of phantom (in plane)
    b = 114.2/2
    if ph_type=='act':
        # Length of phantom
        L = 194 # mm
        # Width of phantom
        W = 298 # mm
        # Distance from end of phantom to the sphere planes
        plane_dist = 70 # mm
    elif ph_type=='attn':
        ## Add more water to the edges of the phantom to resemble the glass housing
        # Length of phantom
        L = 214 # mm
        # Width of phantom
        W = 302 # mm
        # Distance from end of phantom to the sphere planes
        plane_dist = 80 # mm
    elif ph_type=='mu':
        ## Add more water to the edges of the phantom to resemble the glass housing
        # Length of phantom
        L = 214 # 194 # mm
        # Width of phantom
        W = 302 # 300 # mm
        # Distance from end of phantom to the sphere planes
        plane_dist = 80 # mm
        
    # Convert to pixels
    L = L/dz
    W = W/dx
    d = d/dx
    b = b/dx
    plane_dist = plane_dist/dz
    sphere_x_off = sphere_x_off/dx
    sphere_y_off = sphere_y_off/dy
    sphere_z_off = sphere_z_off/dz
    
    ## Create a 2D replicate of the D shaped cylinder
    
    plane = np.zeros(phantom.shape[1:])
    xlocs, ylocs = np.meshgrid(np.arange(plane.shape[0]), np.arange(plane.shape[0]))

    # Draw upper semi-circle

    # Centre locations of semi-circle
    x0 = np.round(plane.shape[0]/2)
    y0 = np.round(plane.shape[1]/2 + 48/dy)
    # Radius
    r = W/2
    # Find pixels inside semi-circle
    inside = r**2 - (x0-xlocs)**2 - (y0-ylocs)**2
    # Set pixels inside semi-circle to 1
    plane[(inside>=0)&(ylocs<=y0)] = 1
    
    # Draw lower two circles

    # Radius
    r = 0.52*W/2
    # Centre locations of first circle
    x0 = np.round(plane.shape[0]/2-W/2+r)
    y0 = np.round(plane.shape[1]/2 + 48/dy)
    # Find pixels inside circle
    inside = r**2 - (x0-xlocs)**2 - (y0-ylocs)**2
    # Set pixels inside semi-circle to 1
    plane[(inside>=0)] = 1
    # Centre locations of second circle
    x0 = np.round(plane.shape[0]/2+W/2-r)
    y0 = np.round(plane.shape[1]/2 + 48/dy)
    # Find pixels inside circle
    inside = r**2 - (x0-xlocs)**2 - (y0-ylocs)**2
    # Set pixels inside semi-circle to 1
    plane[(inside>=0)] = 1
    
    # Draw lower rectangle

    # Top left corner
    x0 = np.round(plane.shape[0]/2-W/2+r)
    y0 = np.round(plane.shape[1]/2 + 48/dy)
    # Width and height
    w = np.ceil(140/dx)
    h = r
    # Set pixels inside rectangle to 1
    plane[(xlocs>=x0)&(xlocs<=x0+w)&(ylocs>=y0)&(ylocs<=y0+h)] = 1
    
    # Remove cylindrical insert
    
    # Centre locations of circle
    x0 = np.round(plane.shape[0]/2)
    y0 = np.round(plane.shape[1]/2 + 13/dy)
    # Radius
    r = d/2
    # Find pixels inside circle
    inside = r**2 - (x0-xlocs)**2 - (y0-ylocs)**2
    # Set pixels inside circle to 0
    plane[(inside>=0)] = 0
    
    ## Copy this plane into the correct slices of the phantom
    
    # Place this plane throughout the phantom
    if ph_type=='mu':
        phantom[int(phantom.shape[0]/2-L/2)+3:int(phantom.shape[0]/2+L/2)+3] = plane
    else:
        phantom[int(phantom.shape[0]/2-L/2):int(phantom.shape[0]/2+L/2)] = plane
    print('The volume of the D-shaped cylinder is %0.2f L' % (len(np.where(phantom==1)[0])*dx*dy*dz*1e-6))
    
    ## Insert spheres
    if spheres:
        # Centre locations of spheres
        x0s = [np.round(plane.shape[0]/2 - b + sphere_x_off),
               np.round(plane.shape[0]/2 - b * np.cos(60 * np.pi/180) + sphere_x_off),
               np.round(plane.shape[0]/2 + b * np.cos(60 * np.pi/180) + sphere_x_off),
               np.round(plane.shape[0]/2 + b + sphere_x_off),
               np.round(plane.shape[0]/2 + b * np.cos(60 * np.pi/180) + sphere_x_off),
               np.round(plane.shape[0]/2 - b * np.cos(60 * np.pi/180) + sphere_x_off)]
        y0s = [np.round(plane.shape[1]/2 + 13/dy + sphere_y_off),
               np.round(plane.shape[1]/2 + b * np.sin(60 * np.pi/180) + 13/dy + sphere_y_off),
               np.round(plane.shape[1]/2 + b * np.sin(60 * np.pi/180) + 13/dy + sphere_y_off),
               np.round(plane.shape[1]/2 + 13/dy + sphere_y_off),
               np.round(plane.shape[1]/2 - b * np.sin(60 * np.pi/180) + 13/dy + sphere_y_off),
               np.round(plane.shape[1]/2 - b * np.sin(60 * np.pi/180) + 13/dy + sphere_y_off)]
        z0 = np.max(np.where(phantom==1)[0]) - plane_dist + sphere_z_off

        # Grid of x, y, z positions
        ylocs, zlocs, xlocs = np.meshgrid(np.arange(plane.shape[0]), np.arange(phantom.shape[0]), np.arange(plane.shape[0]))
        # Loop through each sphere
        for r,x0,y0 in zip(diams/2, x0s, y0s):    
            # Find pixels inside sphere
            inside = r**2 - ((x0-xlocs)*dx)**2 - ((y0-ylocs)*dy)**2 - ((z0-zlocs)*dz)**2
            # Set pixels inside circle to 10
            phantom[(inside>=0)] = 10
        
    return phantom

def save_phantom(phantom, data_dir, ph_type):
    ''' 
    Save phantom array as ASCII.
    
    inputs:
        phantom (np.array): 3D array containing NEMA phantom.
        data_dir (str): path to data directory where the array will be saved.
    '''
    phantom_fn = os.path.join(data_dir, 'nema_%s.bin'%ph_type)
    # Swap axes to make the simulated sinogram match the one from the scanner
    phantom = np.flip(phantom)
    phantom = np.swapaxes(phantom, 1, 2)
    print('Saving phantom array in %s' % phantom_fn)
    with open(phantom_fn, 'wb') as f:
        f.write(np.round(phantom).flatten().astype('uint16'))

def save_atten_map(data_dir, recon_dir):
    # Load attenuation phantom that consists of material indices for the Gate simulation
    attn_ph = np.fromfile(os.path.join(data_dir, 'nema_attn.bin'), dtype=np.int16).astype(np.float32)

    # Obtain image size and resolution
    hdr_fn = os.path.join(data_dir,'nema_attn.h33')
    img_size = np.array([0,0,0])
    img_res = np.array([0.0,0.0,0.0])
    with open(hdr_fn, 'r') as hdr_file:
        for line in hdr_file.readlines():
            if '!matrix size [1]' in line:
                img_size[2] = int(line.split()[-1])
            if '!matrix size [2]' in line:
                img_size[1] = int(line.split()[-1])
            if '!number of slices' in line:
                img_size[0] = int(line.split()[-1])
            if 'scaling factor (mm/pixel) [1]' in line:
                img_res[2] = float(line.split()[-1])  
            if 'scaling factor (mm/pixel) [2]' in line:
                img_res[1] = float(line.split()[-1]) 
            if 'slice thickness' in line:
                img_res[0] = float(line.split()[-1]) 

    attn_ph = attn_ph.reshape(img_size)
    # Order: z, y, x
    attn_ph = np.swapaxes(attn_ph, 1, 2)
    attn_ph = np.flip(attn_ph, (1,2))
    # MAY NEED TO FLIP IN Z AS WELL

    # Set attenuation of water and air
    attn_ph[attn_ph>=1] = 0.095996585 # cm^-1
    attn_ph[attn_ph==0] = 0.00010409389 # cm^-1
    
    # Save mu-map data
    mu_map_fn = os.path.join(recon_dir, 'umap.v')
    with open(mu_map_fn, 'wb') as f:
        f.write(attn_ph.flatten())
        
    # Save main header
    mu_map_mhdr_fn = os.path.join(recon_dir, 'umap.mhdr')
    with open(mu_map_mhdr_fn, 'w') as f:
        f.write('!INTERFILE:=\n')
        f.write('%comment:=created with code from Oct 26 2015 17:33:53\n')
        f.write('!originating system:=1104\n')
        f.write('%SMS-MI header name space:=image main header\n')
        f.write('%SMS-MI version number:=3.5\n\n')
        
        f.write('%DATA MATRIX DESCRIPTION:=\n')
        f.write('number of time frames:=1\n')
        f.write('%number of horizontal bed offsets:=1\n')
        f.write('number of time windows:=1\n\n')
        
        f.write('%DATA SET DESCRIPTION:=\n')
        f.write('!total number of data sets:=1\n')
        f.write('%data set[1]:={1000000000,umap.v.hdr,umap.v}')
        
    # Save subheader
    mu_map_hdr_fn = os.path.join(recon_dir, 'umap.v.hdr')
    with open(mu_map_hdr_fn, 'w') as f:
        f.write('!INTERFILE:=\n')
        f.write('%comment:=created with code from Oct 26 2015 17:33:53\n')
        f.write('!originating system:=1104\n')
        f.write('%SMS-MI header name space:=image subheader\n')
        f.write('%SMS-MI version number:=3.5\n\n')
        
        f.write('!GENERAL DATA:=\n')
        f.write('!name of data file:=umap.v\n\n')

        f.write('!GENERAL IMAGE DATA:=\n')
        f.write('image data byte order:=LITTLEENDIAN\n')
        f.write('%patient orientation:=HFS\n')
        f.write('%image orientation:={1,0,0,0,1,0}\n')
        f.write('!PET data type:=image\n')
        f.write('number format:=float\n')
        f.write('!number of bytes per pixel:=4\n')
        f.write('number of dimensions:=3\n')
        f.write('matrix axis label[1]:=x\n')
        f.write('matrix axis label[2]:=y\n')
        f.write('matrix axis label[3]:=z\n')
        f.write('matrix size[1]:=%i\n' % img_size[2])
        f.write('matrix size[2]:=%i\n' % img_size[1])
        f.write('matrix size[3]:=%i\n' % img_size[0])
        f.write('scale factor (mm/pixel) [1]:=%0.5f\n' % img_res[2])
        f.write('scale factor (mm/pixel) [2]:=%0.5f\n' % img_res[1])
        f.write('scale factor (mm/pixel) [3]:=%0.5f\n' % img_res[0])
        f.write('horizontal bed translation:=continuous\n')
        f.write('start horizontal bed position (mm):=%0.3f\n' % (-1170 + img_size[0]*img_res[0]/2))
        f.write('end horizontal bed position (mm):=%0.3f\n' % (-1170 - img_size[0]*img_res[0]/2))
        f.write('start vertical bed position (mm):=174\n')
        f.write('number of energy windows:=1\n')
        f.write('%energy window lower level (keV) [1]:=0\n')
        f.write('%energy window upper level (keV) [1]:=0\n\n')
        
        f.write('!IMAGE DATA DESCRIPTION:=\n')
        f.write('!total number of data sets:=1\n')
        f.write('!image duration (sec):=0\n')
        f.write('!image relative start time (sec):=0.0\n')
        f.write('%total uncorrected singles rate:=0\n')
        f.write('%image slope:=1\n')
        f.write('%image intercept:=0.0\n\n')
        
        f.write('%SUPPLEMENTARY ATTRIBUTES:=\n')
        f.write('quantification units:=1/cm\n')
        f.write('slice orientation:=transverse\n')
        f.write('%axial compression:=11')


def save_header(data_dir, phantom_size, dz, dy, dx, ph_type):
    '''Write header file.'''
    hdr_fn = os.path.join(data_dir, 'nema_%s.h33'%ph_type)
    print('Saving phantom header file in %s' % hdr_fn)

    phantom_fn = os.path.join(data_dir, 'nema_%s.bin'%ph_type)
    with open(hdr_fn,'w') as file:
        file.write('!name of data file := %s\n' % phantom_fn)
        file.write('imagedata byte order := LITTLEENDIAN\n')
        file.write('!number format := unsigned integer\n')
        file.write('!matrix size [1] := %i\n' % phantom_size[1])
        file.write('!matrix size [2] := %i\n' % phantom_size[2])
        file.write('scaling factor (mm/pixel) [1] := +%0.3f\n' % dy)
        file.write('scaling factor (mm/pixel) [2] := +%0.3f\n' % dx)
        file.write('!number of slices := %i\n' % phantom_size[0])
        file.write('slice thickness (pixels) := +%0.3f\n' % dz)

def save_ranges(data_dir, dz, dy, dx, sph_act, bgnd_act, ph_type):
    '''
    Save the activity and attenuation ranges for indexing in Gate.
    
    inputs:
        data_dir (str): path to data directory where the array will be saved.
        dz, dy, dx (floats): resolution in mm/pixel for z, y, and x dimensions.
        sph_act (float): activity concentration of spherical inserts in Bq/mL.
        bgnd_act (float): background activity concentration in phantom in Bq/mL.
    '''
    if ph_type=='attn':
        # Create file which defines the pairing between pixel value ranges and materials
        atn_rng_fn = os.path.join(data_dir, 'AttenuationRange.dat')
        print('Saving attenuation range file in %s' % atn_rng_fn)
        with open(atn_rng_fn,'w') as file:
            # First line is the number of materials
            file.write('2\n')
            # Now the index range followed by the material corresponding to that range
            # Air
            file.write('0 0 Air\n')
            # Water
            file.write('1 11 Water\n')

    if ph_type=='act':
        # Create file which defines the activity-index pairs
        # where the activities should be in Bq
        act_rng_fn = os.path.join(data_dir, 'ActivityRange.dat')

        # Volume of each voxel
        vox_vol = dx*dy*dz*1e-3 ## mL

        # Determine activity in each voxel
        sph_act = sph_act*vox_vol # Bq
        bgnd_act = bgnd_act*vox_vol # Bq

        print('Saving activity range file in %s' % act_rng_fn)
        with open(act_rng_fn,'w') as file:
            # First line is the number of indices
            file.write('3 \n')
            # Now the index range followed by the activity corresponding to that range
            # Air
            file.write('0 0 0.0\n')
            # Background
            file.write('1 1 %0.2f\n' % bgnd_act)
            # Spheres
            file.write('10 10 %0.2f' % sph_act)
