import numpy as np
import os
import sys
import configparser

def run_setup(patient_name, cur_dir, config):
    # Directories for patient
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

    # Load the configuration
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

    # Phantom type dependent parametes
    if (ph_type=='nema_vox'):
        # NEMA voxelized phantom
        source_mac = 'source_voxelized.mac'
        attn_mac = 'attenuation_voxelized.mac'
        skip_eq_mat = int(config['NEMA']['skip_eq_mat'])
        bed_pos_source = bed_pos + 1070 # mm
        bed_pos_atten = bed_pos + 1170 # mm
        x_off_src = x_off - 200
        y_off_src = y_off - 200
    elif (ph_type=='uniform_vox'):
        # Uniform cylindrical voxelized phantom
        source_mac = 'source_voxelized.mac'
        attn_mac = 'attenuation_voxelized.mac'
        skip_eq_mat = int(config['Uniform']['skip_eq_mat'])
        bed_pos_source = bed_pos + 1059 # mm
        bed_pos_atten = bed_pos + 1170 # mm
        x_off_src = x_off - 201
        y_off_src = y_off - 201
    elif ph_type=='uniform_shape':
        # Uniform cylindrical phantom defined as a shape within GATE
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
        # Filenames for the phantoms
        atn_rng_fn = os.path.join(data_dir, 'AttenuationRange.dat')
        atn_fn = os.path.join(data_dir, '%s_attn.h33'% ph_type.split('_')[0])
        act_rng_fn = os.path.join(data_dir, 'ActivityRange.dat')
        act_fn = os.path.join(data_dir, '%s_act.h33'%ph_type.split('_')[0])

    # Split the simulation into many jobs
    # and create a bash script for each job

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

            # Create Compute-Canada script by adjusting template script
            with open(script_fn, 'wt') as new_script:
                for line in old_script:
                    if 'Gate' in line:
                        new_script.write(cmd)
                    elif 'cd /path/to/gate/' in line:
                        new_script.write('cd %s\n' % (os.path.join(cur_dir, '..', 'gate')))
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

    print('Preparing gate scripts for %s...' % (patient_name))
    run_setup(patient_name, cur_dir, config)
