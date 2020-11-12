#!/usr/bin/env python
import matplotlib, os
if not os.getenv('DISPLAY'):
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, subprocess, glob
import numpy as np
import iris
import cartopy.crs as ccrs
import shapely.geometry as sgeom

'''
Required processing
latitude/longitude names to lat/lon
regrid orog and mslp to u grid
may need u,v in same file

modules
module load libraries/gcc/6.1.0
module load mpi/mpich/3.3/gnu/6.1.0
module load hdf5/1.10.4/gnu/6.1.0
module load netcdf/4.6.2/gnu/6.1.0

'''

def write_parameter_file(outfile, cmd_detect, cmd_stitch):
    nl = '\n'
    with open(outfile, 'w') as fh:
        fh.write('detect cmd '+cmd_detect+nl)
        fh.write(nl)
        fh.write('stitch cmd '+cmd_stitch+nl)

def fix_lat_lon(fname):
    ''' 
    change latitude to lat, longitude to lon
    '''
    for old, new in zip(['latitude','longitude'], ['lat', 'lon']):
        cmd = 'ncrename -d '+old+','+new+' -v '+old+','+new+' '+fname
        subprocess.call(cmd, shell=True)

def calculate_zg_difference(fname):
    c = iris.load_cube(fname)
    cdiff = c[:,0,:,:] - c[:,1,:,:]
    print (cdiff.var_name)
    iris.save(cdiff, fname[:-3]+'_diff.nc')
    return fname[:-3]+'_diff.nc'

def MO_suite(suite, resol, zg = False):
    m_dict = {}
    m_dict['model'] = suite
    m_dict['resol'] = resol
    m_dict['model_name'] = m_dict['model']
    m_dict['psl_var'] = 'air_pressure_at_sea_level'
    m_dict['u_var'] = 'x_wind'
    m_dict['v_var'] = 'y_wind'
    if not zg:
        m_dict['t_var'] = 'ta'
    else:
        m_dict['t_var'] = 'unknown'
    m_dict['lat_var'] = 'lat'
    m_dict['lon_var'] = 'lon'
    m_dict['topo_var'] = 'surface_altitude'
    return m_dict


def getTrajectories(filename, nVars, nVars_stitch, headerDelimStr, cube):

    print ("Running getTrajectories on %s." % filename)
    print ("nVars set to %d and headerDelimStr set to '%s'" % (nVars, headerDelimStr))

    nVars_offset = nVars_stitch - 2 + 1
    # Using the newer with construct to close the file automatically.
    with open(filename) as f:
        data = f.readlines()

    # Find total number of trajectories and maximum length of trajectories
    numtraj=0
    numtraj_nh=0; numtraj_sh=0
    numPts=[]
    storms = []
    coords_index = [2, 3, 3+nVars_offset, 3+nVars_offset+1, 3+nVars_offset+2, 3+nVars_offset+3]
    coords = {'lon':coords_index[0],'lat':coords_index[1],'year':coords_index[2],'month':coords_index[3],'day':coords_index[4],'hour':coords_index[5]}
  
    for line in data:
        if headerDelimStr in line:
            headArr = line.split()
            numtraj += 1
            numPts.append(int(headArr[1]))
    maxNumPts = max(numPts) # Maximum length of ANY trajectory

    print ("Found %d number of trajectories" % numtraj)

    # Initialize storm and line counter
    stormID=-1
    lineOfTraj=-1

    prodata = np.empty((nVars,numtraj,maxNumPts))
    prodata[:] = np.NAN

    storm = {}
    for i, line in enumerate(data):
        if headerDelimStr in line:  # check if header string is satisfied
            stormID += 1      # increment storm
            lineOfTraj = 0    # reset trajectory line to zero
            headArr = line.split()
            tracklen = int(headArr[1])
            storm[stormID] = {}
            storm[stormID]['length'] = tracklen
            for coord in coords:
                storm[stormID][coord] = []
            storm[stormID]['step'] = []
        else:
            ptArr = line.split()
            if lineOfTraj <= tracklen:
                lon = ptArr[coords['lon']]
                lat = ptArr[coords['lat']]
                year = ptArr[coords['year']]
                month = ptArr[coords['month']]
                day = ptArr[coords['day']]
                hour = ptArr[coords['hour']]
                step = convert_date_to_step(cube, int(year), int(month), int(day), int(hour))
                # now check if there is a gap in the traj, if so fill it in
                if lineOfTraj > 0:
                    #print ('gap in traj ',tracklen, year, month, day, hour, step, storm[stormID]['step'][-1])
                    #print ('cube ',cube)
                    if (step - storm[stormID]['step'][-1]) > 1:
                        # add extra points before the next one in the TempExt trajectory
                        fill_traj_gap(storm, stormID, step, lon, lat, year, month, day, hour)
                for coord in coords:
                    storm[stormID][coord].append(ptArr[coords[coord]])
                storm[stormID]['step'].append(step)
            lineOfTraj += 1   # increment line
            if lineOfTraj == tracklen:
                storms.append(storm[stormID])


    for ii, zz in enumerate(range(numtraj)):
        lat = float(storms[ii]['lat'][0])
        if lat < 0.0:
            numtraj_sh += 1
        else:
            numtraj_nh += 1

    print ('numtraj tota, nh, sh ', numtraj, numtraj_nh, numtraj_sh)

    return numtraj, prodata, numtraj_nh, numtraj_sh, storms

def fill_traj_gap(storm, stormID, step, lon, lat, year, month, day, hour):
    '''
    Fill the gap by linearly interpolating the lon, lat and adding steps
    Be careful around the date line, averaging <360 and >0
    '''
    gap_length = step - storm[stormID]['step'][-1]
    if (float(storm[stormID]['lon'][-1]) > 300. and float(lon) < 20.):
        dlon = (float(storm[stormID]['lon'][-1]) - (float(lon))+360.) / float(gap_length)
        print ('dlon, ',float(storm[stormID]['lon'][-1]), float(lon), dlon)
    elif (float(storm[stormID]['lon'][-1]) < 20. and float(lon) > 300.):
        dlon = ((float(storm[stormID]['lon'][-1])+360.) - (float(lon))) / float(gap_length)
        print ('dlon, ',float(storm[stormID]['lon'][-1]), float(lon), dlon)
    else:
        dlon = (float(storm[stormID]['lon'][-1]) - float(lon)) / float(gap_length)
    dlat = (float(storm[stormID]['lat'][-1]) - float(lat)) / float(gap_length)
    for ig in range(1, gap_length):
        lon1 = float(lon) + float(ig) * dlon
        if lon1 > 360.:
            lon1 -= 360.
        if lon1 < 0.:
            lon1 += 360.
        lat1 = float(lat) + float(ig) * dlat
        storm[stormID]['lon'].append(str(lon1))
        storm[stormID]['lat'].append(str(lat1))
        storm[stormID]['step'].append(storm[stormID]['step'][-1]+1)
        storm[stormID]['year'].append(year)
        storm[stormID]['month'].append(month)
        storm[stormID]['day'].append(day)
        storm[stormID]['hour'].append(hour)

def convert_date_to_step(cube, year, month, day, hour):
    '''
    calculcate the next step, i.e. step+1, so normalise the time to get the next intege
    '''
    #from cftime import utime, datetime
    from netCDF4 import netcdftime,num2date,datetime
    import cf_units

    calendar = cube.coord('time').units.calendar
    tunit = cube.coord('time').units
    tunit1 = netcdftime.utime(str(cube.coord('time').units), calendar=calendar)
    dt_point = tunit1.date2num(datetime(year, month, day, hour))
    delta = dt_point - cube.coord('time').points[0]

    if 'hours' in str(tunit):
        return int(delta/6) + 1
    elif 'days' in str(tunit): 
        return int(delta*4) + 1
    elif 'minutes' in str(tunit):
        return int(delta/(60*6)) + 1
  
def plot_trajectories_cartopy(numtraj, prodata, trajfile, period, model, resol, storms, feature_type = 'TC', title=''):
    fig = plt.figure(figsize=(9,6), dpi=100)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-160))
    ax.set_global()
    for ii, zz in enumerate(range(numtraj)):

        # extract lat/lon from array, ditch nans, convert lon for basemap
        lon = [float(i) for i in storms[ii]['lon']]
        lat = [float(i) for i in storms[ii]['lat']]
        if lat[0] > 40:
            print ('ii ',ii, zz)
            print ('lat ', lat)
            print ('lon ', lon)
        # Shift lon since basemap needs -180 to 180
        #lon = np.where(lon > 180., lon-360., lon)

        ax.plot(lon, lat,
                     linewidth=1.2, transform=ccrs.Geodetic())

    fig.gca().coastlines() 
    if title != '':
        title_full = model+' '+resol+', '+str(period)+', '+str(numtraj)+' '+title+' '
    else:
        title_full = model+' '+resol+', '+str(period)+', '+str(numtraj)+' TempestExtremes TCs'
    print ('ii ',ii, zz)
    plt.title(title_full)
    plt.savefig(trajfile[:-4]+'.png')
    plt.show()

def read_and_plot_tracks(trackedfile, nc_file, tracking_period, suiteid, resol, title = '', feature_type = 'TC', nVars_stitch = 4):

    nVars= nVars_stitch + 6
    headerStr='start'

    cube = iris.load_cube(nc_file)
    # Load trajectories into prodata (nVar,numtraj,maxLines)
    numtraj, prodata, numtraj_nh, numtraj_sh, storms = getTrajectories(trackedfile, nVars, nVars_stitch, headerStr, cube)

    #plot_trajectories_basemap(numtraj, prodata)
    plot_trajectories_cartopy(numtraj, prodata, trackedfile, tracking_period, suiteid, resol, storms, title = title, feature_type = feature_type)

def process_input_files(data_in, data_out, runid, period):
    '''
    Set up and associate file names with variables
    
    Things to think about:
        when using levels of a variable (e.g. zg_available), how to make sure that the
        (0), (1) etc indices in the calculations refer to the right levels?

    '''
    fnames_out = []
    files = ['pslfile', 'zgfile', 'ufile', 'vfile', 'ws10mfile', 'rvfile', 'rvT63file', 'u10mfile', 'u10mfile', 'v10mfile', 'viwvefile', 'viwvnfile', 'tafile']
    filenames_remote = {}; filenames = {}

    data_orog = '/data/d04/jseddon/tempest_extremes_datafiles'
    filenames['topofile'] = os.path.join(data_orog, 'orog_HadGEM3-GC31-LM.nc')

    var_dict = {}
    var_dict['pslfile'] = {'fname':'slp', 'varname':'air_pressure_at_sea_level'}
    var_dict['zgfile'] = {'fname':'zg_available', 'varname':'unknown', 'varname_new':'zg_available'}
    var_dict['ufile'] = {'fname':'ua', 'varname':'x_wind'}
    var_dict['vfile'] = {'fname':'va', 'varname':'y_wind'}
    var_dict['ws10mfile'] = {'fname':'wsas', 'varname':'wind_speed'}
    var_dict['rvfile'] = {'fname':'rv', 'varname':'unknown', 'varname_new':'rv'}
    var_dict['rvT63file'] = {'fname':'rvT63', 'varname':'unknown', 'varname_new':'rvT63'}
    var_dict['u10mfile'] = {'fname':'u10m', 'varname':'x_wind', 'varname_new':'uas'}
    var_dict['v10mfile'] = {'fname':'v10m', 'varname':'y_wind', 'varname_new':'vas'}
    var_dict['viwvefile'] = {'fname':'viwve', 'varname':'unknown', 'varname_new':'viwve'}
    var_dict['viwvnfile'] = {'fname':'viwvn', 'varname':'unknown', 'varname_new':'viwvn'}
    var_dict['tafile'] = {'fname':'ta', 'varname':'unknown', 'varname_new':'ta'}

    fname_format = os.path.join(data_in, 'atmos_bx315a_6h_{}_pt-{}.nc')
    for fname in files:
        filenames_remote[fname] = fname_format.format(period, var_dict[fname]['fname'])

        filenames[fname] = os.path.join(data_out, os.path.basename(filenames_remote[fname]))
        if var_dict[fname]['fname'] == 'ua' or var_dict[fname]['fname'] == 'va' or var_dict[fname]['fname'] == 'wsas' or var_dict[fname]['fname'] == 'rv' or var_dict[fname]['fname'] == 'rvT63' or var_dict[fname]['fname'] == 'u10m' or var_dict[fname]['fname'] == 'v10m':

            # regrid u,v, to t grid and rename variable if necessary
            c = iris.load_cube(filenames_remote[fname])
            cref = iris.load_cube(fname_format.format(period, 'slp'))
            c_regrid = c.regrid(cref, iris.analysis.Linear())
            if 'varname_new' in var_dict[fname]:
                c_regrid.var_name = var_dict[fname]['varname_new']
            iris.save(c_regrid, filenames[fname])

        # rename variables only - could do with ncrename instead
        if var_dict[fname]['fname'] == 'ta' or var_dict[fname]['fname'] == 'zg_available':
            c = iris.load_cube(filenames_remote[fname])
            if 'varname_new' in var_dict[fname]:
                c.var_name = var_dict[fname]['varname_new']
            iris.save(c, filenames[fname])
            
        if not os.path.exists(filenames[fname]):
            cmd = 'cp '+filenames_remote[fname]+' '+filenames[fname]
            #print (cmd)
            os.system(cmd)

        # newer tempest can specify names of latitude, longitude
        #fix_lat_lon(fname_local)
 
    return filenames

def data_files(timestamp, runid, resol, data_in, data_out, m_dict):
    timestamp_day = timestamp[:8]
    print ('time_stamp_day ',timestamp_day)
    file_search = os.path.join(data_in, 'atmos_'+runid+'*'+timestamp_day+'??-*-slp.nc')
    print ('file_search ',file_search)
    files = sorted(glob.glob(file_search))
    print ('files ',files)
    #stop
    for f in files:
        period = os.path.basename(f).split('_')[3]

        filenames = process_input_files(data_in, data_out, runid, period)
   
        psl_var = 'air_pressure_at_sea_level'

        output_uv = False

    return filenames, period, output_uv, psl_var

def construct_command(nl_file, tracking_type):
    dict = {}
    with open(nl_file, 'r') as infile:
        for line in infile.readlines():
            if '&'+tracking_type in line:
                #section = line.strip().split(':')[-1][:-1]
                section = line.strip().split('_')[-1]
                dict[section] = {}
            elif '=' in line:
                line_val = line.strip().split('=',1)
                parameter = line_val[0]
                # strip the comma from the end of the namelist line
                value = line_val[1][:-1]
                dict[section][parameter] = value

    cmd_step = {}; cmds = {}
    for step in dict:
        print ('step ',step)
        cmd_step[step] = []
        for elt in dict[step]:
            if dict[step][elt] != '':
                cmd_step[step].append('--'+elt+' '+dict[step][elt])
        cmds[step] = ' '.join(cmd_step[step])

    print ('cmds ',cmds)
    return cmds

def tracking(timestamp, runid, resol, data_in, data_out, cylc_task_id, hist_file, m_dict, tc_detect_script, tc_stitch_script, zg=False, tracking_type = 'tc_slp'):

    print ('cwd ',os.getcwd())
    #nl_file = './namelist_'+tracking_type.lower()+'.nl'
    nl_file = os.path.join(os.getcwd(), tracking_type+'.nl')
    cmds = construct_command(nl_file, tracking_type)

    filenames, period, output_uv, psl_var = data_files(timestamp, runid, resol, data_in, data_out, m_dict)

    candidatefile = os.path.join(data_out, 'candidate_file_'+timestamp+'_'+tracking_type+'.txt')

    files_extend = False
    if not files_extend:
    # identify candidates
        cmd_io = '{} --in_data "{};{};{};{};{}" --out {} '.format(
               tc_detect_script, filenames['pslfile'], filenames['zgfile'], 
               filenames['ufile'], filenames['vfile'], filenames['topofile'], 
               candidatefile)
    else:
        
        cmd_io = '{} --in_data "{};{};{};{};{};{};{};{};{};{};{};{};{}" --out {} '.format(
          tc_detect_script, filenames['pslfile'], filenames['zgfile'], 
          filenames['ufile'], filenames['vfile'], filenames['rvfile'], 
          filenames['u10mfile'], filenames['v10mfile'], filenames['ws10mfile'], 
          filenames['viwvefile'], filenames['viwvnfile'], filenames['tafile'], 
          filenames['rvT63file'], filenames['topofile'], candidatefile)

    cmd_detect = cmd_io + cmds['detect']
    print ('tc command ',cmd_detect)
    # python version 2 or 3
    if int(sys.version[0]) == 2:
        sts = subprocess.check_output(cmd_detect, shell = True,  stderr=subprocess.STDOUT, universal_newlines=True)
        print ('sts ',sts)
    else:
        sts = subprocess.run(cmd_detect, shell = True,  stdout=subprocess.PIPE,  stderr=subprocess.PIPE, universal_newlines=True)
        print (sts.stdout)
        if 'EXCEPTION' in sts.stdout:
            raise Exception('Stop')

    print ('candidatefile ',candidatefile)

    #candidate_period_start = timestamps[0].split('-')[0]
    #candidate_period_end = timestamps[-1].split('-')[-1]
    tracking_period = period
    #candidatefile = './candidate_file_'+tracking_period+'_'+tracking_type+'.txt'

    #cmd = 'cat '+' '.join(candidatefiles)+' > '+candidatefile
    #print ('cat cmd ',cmd)
    #sts = subprocess.run(cmd, shell = True, stdout=subprocess.PIPE,  stderr=subprocess.PIPE, universal_newlines=True)
    #print ('sts err ',sts.stdout)

    trackedfile = os.path.join(data_out, 'track_file_'+tracking_period+'_'+tracking_type+'.txt')

    # stitch candidates together
    cmd_stitch_io = '{} --in {} --out {} '.format(tc_stitch_script, candidatefile, trackedfile)
    cmd_stitch = cmd_stitch_io + cmds['stitch']
    print ('tc command ',cmd_stitch)
    if int(sys.version[0]) == 2:
        sts = subprocess.check_output(cmd_stitch, shell = True,  stderr=subprocess.STDOUT, universal_newlines=True)
        print ('sts ',sts)
    else:
        sts = subprocess.run(cmd_stitch, shell = True, stdout=subprocess.PIPE,  stderr=subprocess.PIPE, universal_newlines=True)
        print ('sts err ',sts.stderr)

    parameter_outfile=os.path.join(data_out, m_dict['model']+'_'+resol+'_'+tracking_period+'_parameter_output_'+tracking_type+'.txt')
    write_parameter_file(parameter_outfile, cmd_detect, cmd_stitch)

    return candidatefile, trackedfile, filenames['pslfile'], tracking_period

def main():
    RESOL = 'N96'
    print ('os.environ ',os.environ)
    DATA_IN = os.environ['DATAM']
    DATA_OUT = os.path.join(DATA_IN, 'tempest_tracking')
    if not os.path.exists(DATA_OUT):
        os.makedirs(DATA_OUT)
    RUNID = os.environ['RUNID']
    # Current cycle task
    CYLC_TASK_ID = os.environ['CYLC_TASK_ID']
    CYLC_TASK_CYCLE_TIME = os.environ['CYLC_TASK_CYCLE_TIME']
    HIST_FILE = os.path.join(DATA_OUT, 'tracking_history')
    print ('CYLC_TASK_CYCLE_TIME, um_runid ',CYLC_TASK_CYCLE_TIME, RUNID)
    m_dict = MO_suite(RUNID, RESOL, True)

    nVars_stitch = 6 # lon, lat, psl, topo

    tc_detect_script = '/home/d05/hadom/tempestextremes-master/bin/DetectNodes'
    tc_stitch_script = '/home/d05/hadom/tempestextremes-master//bin/StitchNodes'
    timestamps = ['1988090106-1988100100' ,'1988100106-1988110100']

    timestamp = CYLC_TASK_CYCLE_TIME

    track_types = ['TC', 'ETC', 'EW', 'TCVORT']

    track_type = 'tc_slp'
    candidatefile, trackedfile, nc_file, tracking_period = tracking(timestamp, RUNID, RESOL, DATA_IN, DATA_OUT, CYLC_TASK_ID, HIST_FILE, m_dict, tc_detect_script, tc_stitch_script, tracking_type = track_type)

    if os.stat(candidatefile).st_size > 0:
        if os.stat(trackedfile).st_size > 0:
            title = track_type+' tracks'
            read_and_plot_tracks(trackedfile, nc_file, tracking_period, RUNID, RESOL, title=title, feature_type = track_type, nVars_stitch = nVars_stitch)
        else:
            print ('candidatefile has data but no tracks '+trackedfile)
    else:
        print ('candidatefile is empty '+candidatefile)


if __name__ == '__main__':
    main()
