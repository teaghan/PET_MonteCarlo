import numpy as np
import os
import sys
import configparser

def run_setup(patient_name, cur_dir, config):
    # Directories
    patient_dir = os.path.join(cur_dir, '../patients/', patient_name)
    data_dir = os.path.join(patient_dir, 'phantom')
    results_dir = os.path.join(patient_dir, 'results')
    scripts_dir = os.path.join(patient_dir, 'scripts')
    done_dir = os.path.join(patient_dir, 'done')
    stderr_dir = os.path.join(patient_dir, 'stderr')
    stdout_dir = os.path.join(patient_dir, 'stdout')

    # Create directories
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    if not os.path.exists(scripts_dir):
        os.mkdir(scripts_dir)
    if not os.path.exists(done_dir):
        os.mkdir(done_dir)
    if not os.path.exists(stderr_dir):
        os.mkdir(stderr_dir)
    if not os.path.exists(stdout_dir):
        os.mkdir(stdout_dir)

    ph_type = config['Phantom']['type'].lower()
    scan_time = float(config['Simulation']['scan_time']) # s
    num_jobs = int(config['Simulation']['num_jobs'])
    bed_pos = float(config['Setup']['bed_pos']) # mm
    x_off = float(config['Setup']['x_off']) # mm
    y_off = float(config['Setup']['y_off']) # mm
    quant_eff = float(config['Simulation']['quantum_eff'])
    ener_res = float(config['Simulation']['ener_res'])
    cryst_res = eval(config['Simulation']['cryst_res'])
    multi_policy = config['Simulation']['multi_policy']
    coinc_window = float(config['Simulation']['coinc_window'])
    time_blur = float(config['Simulation']['time_blur'])

    if (ph_type=='nema_vox') or (ph_type=='nema_vox2'):
        source_mac = 'source_voxelized.mac'
        attn_mac = 'attenuation_voxelized.mac'
        skip_eq_mat = int(config['NEMA']['skip_eq_mat'])
        bed_pos_source = bed_pos + 1070 # mm
        bed_pos_atten = bed_pos + 1170 # mm
        x_off_src = x_off - 200
        y_off_src = y_off - 200
    if ph_type=='nema_vox2':
        # Axial location of sphere plane
        z_sphere_loc = bed_pos_source + 74 # mm
        # Distance from centre of spheres to centre of phantom (in plane)
        b = 114.2/2 # mm
        # Locations of each sphere in plane
        x_sphere = np.array([b, b*np.cos(60*np.pi/180), -b*np.cos(60*np.pi/180),
                             -b, -b*np.cos(60*np.pi/180), b*np.cos(60*np.pi/180)]) # mm
        y_sphere = np.array([0, -b*np.sin(60*np.pi/180), -b*np.sin(60*np.pi/180),
                             0, b*np.sin(60*np.pi/180), b*np.sin(60*np.pi/180)]) - 13 # mm
        
        # Activity of each sphere
        act_conc_sphere = float(config['NEMA']['sph_act_conc']) # Bq/mL
        r_sphere = np.array([37,28,22,17,13,10])/2 # mm
        vol_sphere = 4/3 * np.pi * r_sphere**3 * 1e-3 ## mL
        act_sphere = act_conc_sphere * vol_sphere # Bq
        
        # Outer radius of housing
        r_out_sphere = r_sphere + 2 # mm

    elif (ph_type=='xcat'):
        source_mac = 'source_voxelized.mac'
        attn_mac = 'attenuation_voxelized.mac'
        skip_eq_mat = int(config['XCAT']['skip_eq_mat'])
        num_slices = int(config['XCAT']['endslice']) - int(config['XCAT']['startslice']) + 1    
        bed_pos_source = bed_pos + 1170 - num_slices*10*float(config['XCAT']['slice_width'])/2 # mm
        bed_pos_atten = bed_pos + 1170 # mm
        x_off_src = x_off - int(config['XCAT']['array_size'])*10*float(config['XCAT']['pixel_width'])/2
        y_off_src = y_off - int(config['XCAT']['array_size'])*10*float(config['XCAT']['pixel_width'])/2

        # Data Filenames
        atn_rng_fn = os.path.join(data_dir, 'AttenuationRange.dat')
        act_rng_fn = os.path.join(data_dir, 'ActivityRange.dat')
        
        # Load breathing data
        breath_file = os.path.join(data_dir, 'patient_breathing_motion.npz')
        breath_data = np.load(breath_file, allow_pickle=True)
        breath_time = breath_data['time_stamp']
        breath_phase = breath_data['phase']
        
        # Phase associated with the different XCAT frames created
        frame_phases = np.linspace(0, 0.4, int(config['XCAT']['out_frames']))

    elif (ph_type=='uniform_vox'):
        source_mac = 'source_voxelized.mac'
        attn_mac = 'attenuation_voxelized.mac'
        skip_eq_mat = int(config['Uniform']['skip_eq_mat'])
        bed_pos_source = bed_pos + 1059 # mm
        bed_pos_atten = bed_pos + 1170 # mm
        x_off_src = x_off - 201
        y_off_src = y_off - 201
    elif ph_type=='uniform_shape':
        source_mac = 'source_cylinder.mac'
        attn_mac = 'attenuation_cylinder.mac'
        cyl_inner_rad = float(config['Cylinder']['cyl_inner_rad'])
        cyl_inner_len = float(config['Cylinder']['cyl_inner_len'])
        cyl_outer_rad = float(config['Cylinder']['cyl_outer_rad'])
        cyl_outer_len = float(config['Cylinder']['cyl_outer_len'])
        housing_material = config['Cylinder']['housing_material']
        cyl_half_len = cyl_inner_len/2
        # Volume of solution
        vol = np.pi * cyl_inner_rad**2 * cyl_inner_len * 0.001 # mL
        # Totaly activity of solution
        cyl_act = float(config['Cylinder']['cyl_act_conc']) * vol # Bq
        # Central position
        bed_pos = bed_pos + 1170 # mm
        # Central position of ends of housing
        cyl_end_len = (cyl_outer_len - cyl_inner_len)/2
        cyl_top_pos = bed_pos - cyl_half_len - cyl_end_len/2
        cyl_btm_pos = bed_pos + cyl_half_len + cyl_end_len/2
        # Rotation around x-axis
        x_rot = float(config['Setup']['x_rot']) # deg
        y_rot = float(config['Setup']['y_rot']) # deg
        # New x-axis orientation
        new_y = '%0.3f %0.3f %0.3f'%(np.cos(x_rot*np.pi/180),
                                     np.sin(y_rot*np.pi/180)*np.sin(x_rot*np.pi/180),
                                     np.cos(y_rot*np.pi/180)*np.sin(x_rot*np.pi/180))
        new_x = '0. %0.3f %0.3f'%(np.cos(y_rot*np.pi/180), -1*np.sin(y_rot*np.pi/180))
        
        # The rotation of the phantom is defined differently
        R = np.array([[np.cos(x_rot*np.pi/180), 0, -np.sin(x_rot*np.pi/180)],
             [np.sin(y_rot*np.pi/180)*np.sin(x_rot*np.pi/180), 
              np.cos(y_rot*np.pi/180),
              np.sin(y_rot*np.pi/180)*np.cos(x_rot*np.pi/180)],
             [np.cos(y_rot*np.pi/180)*np.sin(x_rot*np.pi/180),
              -1*np.sin(y_rot*np.pi/180),
             np.cos(y_rot*np.pi/180)*np.cos(x_rot*np.pi/180)]])
        rot_angle = np.arccos((np.trace(R) - 1)/2)
        rot_axis = 1 / (2*np.sin(rot_angle)) * np.array([R[2,1] - R[1,2],
                                                         R[0,2] - R[2,0],
                                                         R[1,0] - R[0,1]])
        rot_angle *= 180/np.pi # deg
        rot_axis = '%0.3f %0.3f %0.3f' % (rot_axis[0], rot_axis[1], rot_axis[2])

    if multi_policy=='all':
        multi_policy = 'takeAllGoods'
    elif multi_policy=='winner':
        multi_policy = 'takeWinnerOfGoods'

    if 'vox' in ph_type:
        # Data Filenames
        atn_rng_fn = os.path.join(data_dir, 'AttenuationRange.dat')
        atn_fn = os.path.join(data_dir, '%s_attn.h33'% ph_type.split('_')[0])
        act_rng_fn = os.path.join(data_dir, 'ActivityRange.dat')
        act_fn = os.path.join(data_dir, '%s_act.h33'%ph_type.split('_')[0])

    # Create a bash script for each simulation

    # Template script
    with open('mc_job.sh', 'rt') as old_script:
        time_slice = scan_time/num_jobs
        for i in range(1, num_jobs+1):
            # Start and end time for this simulation
            start_time = time_slice * (i-1)
            end_time = start_time + time_slice

            # Output save name
            root_fn = os.path.join(results_dir, str(i))

            # Name of script
            script_fn = os.path.join(scripts_dir, patient_name+'_'+str(i)+'.sh')  

            if ph_type=='uniform_shape':
                # Gate call for this script
                cmd = ('Gate -a "[root_fn,%s]' % (root_fn)+
                       ' [source_mac,%s] [attn_mac,%s]' % (source_mac, attn_mac)+
                       ' [TimeSlice,%0.5f] [TimeStart,%0.5f] [TimeStop,%0.5f]' % (time_slice, start_time, end_time)+
                       ' [bed_pos,%0.3f] [cyl_half_len,%0.3f] [cyl_end_len,%0.3f]' % (bed_pos, cyl_half_len, cyl_end_len)+
                       ' [x_off,%0.3f] [y_off,%0.3f]' % (x_off, y_off)+
                       ' [new_x,%0s] [new_y,%0s]' % (new_x, new_y)+
                       ' [rot_axis,%0s] [rot_angle,%0.3f]' % (rot_axis, rot_angle)+
                       ' [cyl_inner_rad,%0.3f] [cyl_inner_len,%0.3f]' % (cyl_inner_rad, cyl_inner_len)+
                       ' [cyl_outer_rad,%0.3f] [cyl_outer_len,%0.3f]' % (cyl_outer_rad, cyl_outer_len)+
                       ' [cyl_top_pos,%0.3f] [cyl_btm_pos,%0.3f] [cyl_act,%0.3f]' % (cyl_top_pos, cyl_btm_pos, cyl_act)+
                       ' [housing_material,%s] [coinc_window,%0.4f]' % (housing_material, coinc_window)+
                       ' [time_blur,%0.6f]' % (time_blur)+
                       ' [ener_res,%0.3f] [cryst_res_min,%0.3f] [cryst_res_max,%0.3f]' % (ener_res, cryst_res[0], cryst_res[1])+
                       ' [multi_policy,%s] [quant_eff,%0.3f]" main.mac\n' % (multi_policy, quant_eff))
            elif 'vox2' in ph_type:
                # Gate call for this script
                cmd = ('Gate -a "[root_fn,%s] [atn_fn,%s]' % (root_fn, atn_fn)+
                       ' [atn_rng_fn,%s] [act_fn,%s] [act_rng_fn,%s]' % (atn_rng_fn, act_fn, act_rng_fn)+
                       ' [source_mac,%s] [attn_mac,%s] [skip_eq_mat,%i]' % (source_mac, attn_mac, skip_eq_mat)+
                       ' [TimeSlice,%0.5f] [TimeStart,%0.5f] [TimeStop,%0.5f]' % (time_slice, start_time, end_time)+
                       ' [bed_pos_source,%0.3f] [bed_pos_atten,%0.3f]' % (bed_pos_source, bed_pos_atten)+
                       ' [x_off,%0.3f] [y_off,%0.3f]' % (x_off, y_off)+
                       ' [x_off_src,%0.3f] [y_off_src,%0.3f]' % (x_off_src, y_off_src)+
                       ' [z_sphere_loc,%0.3f]' % (z_sphere_loc)+
                       ' [x_sphere1,%0.3f] [y_sphere1,%0.3f]' % (x_sphere[0], y_sphere[0])+
                       ' [x_sphere2,%0.3f] [y_sphere2,%0.3f]' % (x_sphere[1], y_sphere[1])+
                       ' [x_sphere3,%0.3f] [y_sphere3,%0.3f]' % (x_sphere[2], y_sphere[2])+
                       ' [x_sphere4,%0.3f] [y_sphere4,%0.3f]' % (x_sphere[3], y_sphere[3])+
                       ' [x_sphere5,%0.3f] [y_sphere5,%0.3f]' % (x_sphere[4], y_sphere[4])+
                       ' [x_sphere6,%0.3f] [y_sphere6,%0.3f]' % (x_sphere[5], y_sphere[5])+
                       ' [act_sphere1,%0.3f] [act_sphere2,%0.3f]' % (act_sphere[0], act_sphere[1])+
                       ' [act_sphere3,%0.3f] [act_sphere4,%0.3f]' % (act_sphere[2], act_sphere[3])+
                       ' [act_sphere5,%0.3f] [act_sphere6,%0.3f]' % (act_sphere[4], act_sphere[5])+
                       ' [r_sphere1,%0.3f] [r_sphere2,%0.3f]' % (r_sphere[0], r_sphere[1])+
                       ' [r_sphere3,%0.3f] [r_sphere4,%0.3f]' % (r_sphere[2], r_sphere[3])+
                       ' [r_sphere5,%0.3f] [r_sphere6,%0.3f]' % (r_sphere[4], r_sphere[5])+
                       ' [r_out_sphere1,%0.3f] [r_out_sphere2,%0.3f]' % (r_out_sphere[0], r_out_sphere[1])+
                       ' [r_out_sphere3,%0.3f] [r_out_sphere4,%0.3f]' % (r_out_sphere[2], r_out_sphere[3])+
                       ' [r_out_sphere5,%0.3f] [r_out_sphere6,%0.3f]' % (r_out_sphere[4], r_out_sphere[5])+
                       ' [coinc_window,%0.4f]' % (coinc_window)+
                       ' [time_blur,%0.6f]' % (time_blur)+
                       ' [ener_res,%0.3f] [cryst_res_min,%0.3f] [cryst_res_max,%0.3f]' % (ener_res, cryst_res[0], cryst_res[1])+
                       ' [multi_policy,%s] [quant_eff,%0.3f]" main2.mac\n' % (multi_policy, quant_eff))
            elif 'vox' in ph_type:
                # Gate call for this script
                cmd = ('Gate -a "[root_fn,%s] [atn_fn,%s]' % (root_fn, atn_fn)+
                       ' [atn_rng_fn,%s] [act_fn,%s] [act_rng_fn,%s]' % (atn_rng_fn, act_fn, act_rng_fn)+
                       ' [source_mac,%s] [attn_mac,%s] [skip_eq_mat,%i]' % (source_mac, attn_mac, skip_eq_mat)+
                       ' [TimeSlice,%0.5f] [TimeStart,%0.5f] [TimeStop,%0.5f]' % (time_slice, start_time, end_time)+
                       ' [bed_pos_source,%0.3f] [bed_pos_atten,%0.3f]' % (bed_pos_source, bed_pos_atten)+
                       ' [x_off,%0.3f] [y_off,%0.3f]' % (x_off, y_off)+
                       ' [x_off_src,%0.3f] [y_off_src,%0.3f]' % (x_off_src, y_off_src)+
                       ' [coinc_window,%0.4f]' % (coinc_window)+
                       ' [time_blur,%0.6f]' % (time_blur)+
                       ' [ener_res,%0.3f] [cryst_res_min,%0.3f] [cryst_res_max,%0.3f]' % (ener_res, cryst_res[0], cryst_res[1])+
                       ' [multi_policy,%s] [quant_eff,%0.3f]" main.mac\n' % (multi_policy, quant_eff))
            elif ph_type=='xcat':
                
                ## USE START AND END TIME TO INDEX INTO BREATHING PATTERN
                
                # Get average time of current simulation
                simu_time = (start_time + end_time) / 2
                
                # Find nearest time in breathing data and collect phase
                breath_indx = np.argmin(np.abs(breath_time - simu_time))
                cur_phase = breath_phase[breath_indx]
                
                # Determine index of frame based on phase
                ph_indx = np.argmin(np.abs(frame_phases - cur_phase)) + 1
                
                # Use this xcat phase in the current simulation
                atn_fn = os.path.join(data_dir, '%s_attn_%i.h33'% (ph_type.split('_')[0], ph_indx))
                act_fn = os.path.join(data_dir, '%s_act_%i.h33' % (ph_type.split('_')[0], ph_indx))
                
                # Gate call for this script
                cmd = ('Gate -a "[root_fn,%s] [atn_fn,%s]' % (root_fn, atn_fn)+
                       ' [atn_rng_fn,%s] [act_fn,%s] [act_rng_fn,%s]' % (atn_rng_fn, act_fn, act_rng_fn)+
                       ' [source_mac,%s] [attn_mac,%s] [skip_eq_mat,%i]' % (source_mac, attn_mac, skip_eq_mat)+
                       ' [TimeSlice,%0.5f] [TimeStart,%0.5f] [TimeStop,%0.5f]' % (time_slice, start_time, end_time)+
                       ' [bed_pos_source,%0.3f] [bed_pos_atten,%0.3f]' % (bed_pos_source, bed_pos_atten)+
                       ' [x_off,%0.3f] [y_off,%0.3f]' % (x_off, y_off)+
                       ' [x_off_src,%0.3f] [y_off_src,%0.3f]' % (x_off_src, y_off_src)+
                       ' [time_blur,%0.6f]' % (time_blur)+
                       ' [coinc_window,%0.4f]' % (coinc_window)+
                       ' [ener_res,%0.3f] [cryst_res_min,%0.3f] [cryst_res_max,%0.3f]' % (ener_res, cryst_res[0], cryst_res[1])+
                       ' [multi_policy,%s] [quant_eff,%0.3f]" main.mac\n' % (multi_policy, quant_eff))

            # Create Compute-Canada script by adjusting template script
            with open(script_fn, 'wt') as new_script:
                for line in old_script:
                    if 'Gate' in line:
                        new_script.write(cmd)
                    elif 'cd /path/to/gate/' in line:
                        new_script.write('cd %s\n' % (os.path.join(cur_dir, '..', 'gate')))
                    #elif './run_sorting.sh path/to/output/root "delay"' in line:
                        #new_script.write('./run_sorting.sh %s "delay"\n' % (root_fn))
                    #elif 'mv path/to/output/delay path/to/output/delay_new' in line:
                        #new_script.write('mv %s %s\n' % (root_fn+'_prompt.lor', root_fn+'_delay.lor'))
                    #elif './run_sorting.sh path/to/output/root "Coincidences"' in line:
                    #    new_script.write('./run_sorting.sh %s "Coincidences"\n' % (root_fn))
                    #elif 'cd /path/to/analysis' in line:
                    #    new_script.write('cd %s\n' % (os.path.join(cur_dir, 'sorting')))
                    #elif 'cd /path/to/results/' in line:
                    #    new_script.write('cd %s\n' % (results_dir))
                    #elif '/analysis_dir/mct4r_scan_sort -L 1_promt.lor -S 1_prompt.s -d -r 42.8,0' in line:
                    #    new_script.write('%s -L %i_prompt.lor -S %i_prompt.s -r 42.8,0.7 -h' % 
                    #                     (os.path.join(cur_dir, 'sorting', 'mct4r_scan_sort'),i,i))
                    else:
                        new_script.write(line)
            old_script.seek(0)

    root_fn = os.path.join(results_dir, 'combine_root.sh')
    with open('combine_root.sh', 'rt') as old_root:
        with open(root_fn, 'w') as root_script:
                for line in old_root:
                    if 'hadd ' in line:
                        root_script.write('hadd -f combined.root ')
                        for i in range(1, num_jobs+1):
                            root_script.write('%i.root ' % i)
                    else:
                        root_script.write(line)

if __name__ == "__main__":
    
    # Command line arguments
    patient_name = sys.argv[1]
    
    cur_dir = os.path.dirname(os.path.realpath(__file__))
    config_dir = os.path.join(cur_dir, '../patients/configs/')
    
    # Patient configuration
    config = configparser.ConfigParser()
    config.read(os.path.join(config_dir, patient_name+'.ini'))
    '''
    if config['Phantom']['type'].lower() == 'xcat':
        out_frames = int(config['XCAT']['out_frames'])
        # Setup patient for each phase
        for resp_start_ph_index in np.round(np.arange(0,1,1/out_frames),3):
            print('Preparing gate scripts for %s_%i...' % (patient_name, resp_start_ph_index*1000))
            run_setup(patient_name+'_%i'%(resp_start_ph_index*1000), cur_dir, config)
    else:
    '''
    print('Preparing gate scripts for %s...' % (patient_name))
    run_setup(patient_name, cur_dir, config)
