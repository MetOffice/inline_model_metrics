#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

import iris, os
# Must be here to run on batch systems
import subprocess
import glob
import cf_units
from datetime import datetime
import metpy
import metpy.calc as mpcalc
from metpy.units import units
import xarray as xr
from netCDF4 import Dataset

# 10 day files, starting date/time is 01:06:00, so 00:00 is at the end of previous month. Need to make this agree with interpolating the daily data

months = {'01':'jan', '02':'feb', '03':'mar', '04':'apr', '05':'may', '06':'jun', '07':'jul', '08':'aug', '09':'sep', '10':'oct', '11':'nov', '12':'dec'}

dir_in_format = '/gws/nopw/j04/hrcm/n1280run/{suite}/6hrly'
dir_in_daily_format = '/gws/nopw/j04/hrcm/cache/malcolm/{suite}/dl3/'

stream_out = 'pt'
fname_out = '{}a.pt{}_{}.nc'
reftime = '1988-9-1,0:0:0,6hour'

stream_daily = 'pe'
variables_daily = ['sh_100', 'ta_100']

variables = ['psl', 'ua', 'va', 'zg', 'sfcWind', 'sh', 'ta', 'viwve', 'viwvn', 'sh_100', 'ta_100']
streams = ['pc', 'pc', 'pc', 'pc', 'pc', 'pc', 'pc', 'pc', 'pc', stream_daily, stream_daily]
streams_unique = list(set(streams))

stash = {}
stash['psl'] = '16222'
stash['ua'] = '30201'
stash['va'] = '30202'
stash['zg'] = '30207'
#stash['uas'] = '03209'
#stash['vas'] = '03210'
stash['sfcWind'] = '03227'
#stash['pr'] = '05216'
#stash['ts'] = '00024'
stash['sh'] = '30205'
stash['sh_100'] = '30205'
stash['viwve'] = '30407'
stash['viwvn'] = '30408'
stash['ta'] = '30204'
stash['ta_100'] = '30204'

levels = {}
levels['ua'] = [925, 850, 500, 250]
levels['va'] = [925, 850, 500, 250]
levels['zg'] = [925, 850, 700, 500, 250]
levels['ta'] = [850]
levels['sh'] = [850]

GAP = 2

def convert_to_date(date_string):
    year = date_string[0:4]
    if len(date_string) == 7:
        mon_str = date_string[4:]
        mon = datetime.strptime(mon_str.title(), '%b').month            
        day = '01'
        datestr = year+str(mon).zfill(2)+day
        date = datetime.strptime(date_string+'01', '%Y%b%d')
    else:
        mon = date_string[4:6]
        day = date_string[6:]
        datestr = year+mon+day
        date = datetime.strptime(year+mon+day, '%Y%m%d')

    return date, datestr

def latest_mass_file(suite, stream = 'p7', ens=""):
    if ens == "":
        cmd = 'moo ls -t moose:/crum/'+suite+'/a'+stream+'.pp'
    else:
        cmd = 'moo ls -t moose:/ens/'+suite+'/'+ens+'/a'+stream+'.pp'
    sts = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = sts.communicate()
    mass_list = sorted(str(out).strip().split('\\n'))
    if len(mass_list) > 0:
        print('last mass file ',mass_list[-1])
        mass_date = ((mass_list[-1].split('/')[-1]).split('.')[1][2:])
        return mass_date
    else:
        return None

def latest_local_file(suite, dirname):
    cmd1 = 'ls -l '+dirname+'/'+suite[2:]+'a.*.nc'
    sts1 = subprocess.Popen(cmd1, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = sts1.communicate()
    file_listing = (str(out).split('\\n'))
    file_list = []
    for f in file_listing:
        a = f.split(' ')
        if len(a) >= 8:
            fpath = a[-1]
            file_list.append(fpath)
    if len(file_list) > 0:
        file_date = os.path.basename(file_list[-1].split('.')[1][2:10])
        print('last local file ',file_list[-1])
        return file_date
    else:
        return None

def get_file(fin, dirname, strm, ens=""):
    if not os.path.exists(os.path.join(dirname, fin)):
        if ens == "":
            cmd = 'moo get -i moose:/crum/'+um_suiteid+'/a'+strm+'.pp/'+fin+' '+dirname
        else:
            cmd = 'moo get -i moose:/ens/'+um_suiteid+'/'+ens+'/a'+strm+'.pp/'+fin+' '+dirname            
        print(cmd)
        #os.system(cmd)
        run_cmd(cmd)

def get_stash_constraint(stashcode):
    st_pattern = 'm01s{}i{}'
    stash = st_pattern.format(stashcode[0:2], stashcode[2:])
    stash_con = iris.AttributeConstraint(STASH=stash)
    return stash_con

def run_cmd(
        cmd,
    check=True
):
    sts = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=check,
    )
    print('cmd ', sts.stdout)
    if sts.stderr:
        if 'Warning' in sts.stderr or 'warning' in sts.stderr:            
            msg = (
                "Warning found in cat output ", sts.stderr
            )
            print(msg)
        else:
            msg = (
                "Error found in cat output ", sts.stderr
            )
            raise RuntimeError(msg)
    return sts

def change_time_units(fout, fout_tmp1, fout_tmp2):
    cmd = 'cdo setreftime,'+reftime+' '+fout_tmp1+' '+fout_tmp2
    print(cmd)
    run_cmd(cmd)

    # remove the leadtime variable that cdo adds
    cmd = 'ncks -O -x -v leadtime '+ fout_tmp2+' '+fout
    print(cmd)
    run_cmd(cmd)

def get_environment_variables():
    """
    Get the required environment variables from the suite. A list and
    explanation of the required environment variables is included in the
    documentation.
    """
    global um_runid, um_suiteid, cylc_task_cycle_time, time_cycle, previous_cycle, tm2_cycle, next_cycle, startdate, enddate, current_year, current_month, current_day, period, cycleperiod, dir_out, data_freq, data_freq_str, ens, lastcycle, prevcycle, nextcycle

    #try:
    #    um_runid = os.environ["RUNID_OVERRIDE"]
    #except:
    #    um_runid = os.environ["RUNID"]
    try:
        um_suiteid = os.environ["SUITEID_OVERRIDE"]
    except:
        um_suiteid = os.environ["CYLC_SUITE_NAME"]
    try:
        ens = os.environ["ENS"]
    except:
        ens = ""

    um_runid = um_suiteid.split('-')[1]

    cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
    time_cycle = os.environ["TIME_CYCLE"]
    startdate = os.environ["STARTDATE"]
    enddate = os.environ["ENDDATE"]
    cycleperiod = os.environ["CYCLEPERIOD"]
    dir_out = os.environ["DIR_OUT"]
    data_freq_str = os.environ["DATA_FREQ_STRING"]
    data_freq = int(data_freq_str.strip("h"))
    lastcycle = os.environ["LASTCYCLE"]
    prevcycle = os.environ["PREVIOUS_CYCLE"]
    nextcycle = os.environ["NEXT_CYCLE"]

    current_year = time_cycle[0:4]
    current_month = time_cycle[4:6]
    current_day = time_cycle[6:8]
    period = str(current_year)+str(current_month)+str(current_day)

def extract_levels(pp, fout, fout_tmp1, fout_tmp2, var):
    '''
    extract levels from file and save each level separately
    '''
    fnames = []
    for level in levels[var]:
        fout_level = fout[:-3]+'_'+str(level)+'.nc'
        if not os.path.exists(fout_level):
            constraint = iris.Constraint(coord_values = {'pressure': lambda l : l == level})
            c = pp.extract(constraint)
            var_name = var+'_'+str(level)
            c.var_name = var_name
            iris.save(c, fout_tmp1, unlimited_dimensions = ['time'])
            change_time_units(fout_level, fout_tmp1, fout_tmp2)
            if os.path.exists(fout[:-3]+'_'+str(levels[var][0])+'.nc'):
                os.remove(fout_tmp1)
                os.remove(fout_tmp2)
            fnames.append(fout_level)
    return fnames

def limit_latitude(pp):
    lat_min = pp.coord('latitude').points.min()
    lat_max = pp.coord('latitude').points.max()

    if lat_min < -90.0 or lat_max > 90.0:
        pp.coord('latitude').bounds = None
        if pp.coord('latitude').points[0] > 90.0:
            pp.coord('latitude').points[0] = 90.0
        elif pp.coord('latitude').points[0] < -90.0:
            pp.coord('latitude').points[0] = -90.0

        if pp.coord('latitude').points[-1] > 90.0:
            pp.coord('latitude').points[-1] = 90.0
        elif pp.coord('latitude').points[-1] < -90.0:
            pp.coord('latitude').points[-1] = -90.0

        pp.coord('latitude').guess_bounds()

def modify_time_coord(pp, fout, fout_tmp1, fout_tmp2, st):
    '''
    Need to adjust mean time coord to have same times as instantaneous
    Note that the mean is currently wrong - offset by 3 hours
    '''
    if os.path.exists(fout):
        return
    times = pp.coord('time').points[1]
    if (times % data_freq) != 0:
        # make a new time coordinate
        calendar = pp.coord('time').units.calendar
        values = pp.coord('time').points.copy()
        values += data_freq/2
        new_time = iris.coords.DimCoord(values, standard_name = 'time', units = cf_units.Unit(pp.coord('time').units, calendar = calendar))
        print ('new_time ',new_time.points)
        pp.remove_coord('time')
        pp.add_dim_coord(new_time, 0)
        pp.var_name = st
        iris.save(pp, fout_tmp1, unlimited_dimensions = ['time'])
        change_time_units(fout, fout_tmp1, fout_tmp2)
    return [fout]

def test_file_dates(suite, dirname, period, ens=""):
    '''
    Find what data is available in MASS
    Find what data has already been tracked
    Make sure at least a 1 year gap between the most recent MASS data
    If current date fits, then go ahead and get data
    If this data already tracked then need to figure out what to do

    if cylc date is before most recent file date:
      done this already, so return?
    if cylc date after most recent file date:
      if cylc date at least 1 year before final MASS date:
        if cylc date is next date in line:
          go ahead and get data
    '''
    recent_mass_file_date = latest_mass_file(suite, ens=ens)
    if not recent_mass_file_date == None:
        datem, datestrm = convert_to_date(recent_mass_file_date)
    else:
        datem, datestrm = None, None

    recent_local_file_date = latest_local_file(suite, dirname)
    if not recent_local_file_date == None:
        datef, datestrf = convert_to_date(recent_local_file_date)
    else:
        datef, datestrf = None, None

    print('most recent file is ',datef, datestrf)
    print('most recent mass file is ',datem, datestrm)

    print('compare curr time ',period,' with file ',datestrf, 'and mass ',datestrm)
    if datestrf is None:
        if datestrm is not None:
            condition = 'GetData'
        else:
            condition = 'FailAndRetry'
    else:
        if int(period) <= int(datestrf):
            # this tracking has already been done
            condition = 'AlreadyComplete'
        elif int(period) > int(datestrf):
            # do this data if at least one year behind mass
            if int(datestrm[0:6]) - int(period[0:6]) > 12:
                condition = 'GetData'
            else:
                condition = 'FailAndRetry'
        else:
            condition = 'FailAndRetry'

    return condition

def convert_sh_to_rh(t_file, s_file, rh_file, level, t_var, s_var, rh_var):

    dst = xr.open_dataset(t_file)
    dssh = xr.open_dataset(s_file)

    dst = dst.metpy.parse_cf()
    dssh = dssh.metpy.parse_cf()

    dst = dst[t_var].metpy.quantify()
    dssh = dssh[s_var].metpy.quantify()

    rh = mpcalc.relative_humidity_from_specific_humidity(int(level)*units.hPa, dst, dssh)

    rh = rh.metpy.convert_units('percent')
    rh = rh.metpy.dequantify().drop_vars('metpy_crs')
    rh_new = rh.rename(rh_var)
    print(rh_new)
    rh_new.to_netcdf(rh_file)

    #dsin = Dataset(t_file)
    #dsout = Dataset(rh_file, 'r+')
    
    #for v_name, varin in dsin.variables.items():
    #    print('v_name ',v_name, varin)
    #    if 'latitude_longitude' in varin or 'latitude_longitude' in v_name:
    #        outVar = dsout.createVariable(v_name, varin.datatype, varin.dimensions)
    #        outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
    #        outVar[:] = varin[:]
    #dsout.close()

    #cmd = 'ncatted -O -a grid_mapping,'+rh_var+',c,c,latitude_longitude '+rh_file

def main():
    get_environment_variables()
    if ens == "":
        dir_write = os.path.join(dir_out, um_suiteid, data_freq_str)
    else:
        dir_write = os.path.join(dir_out, um_suiteid, data_freq_str, ens)

    if not os.path.exists(dir_write):
        os.makedirs(dir_write)

    dir_in = dir_in_format.format(suite=um_suiteid)

    files_to_delete = []

    files = {}
    if ens == "":
        fname_start = um_runid
    else:
        fname_start = um_runid+'-'+ens

    for strm in streams_unique:
        if strm == stream_daily:
            continue
        file_in = fname_start+'a.'+strm+''+str(current_year)+str(current_month)+'*.pp'

        print('file search ',dir_in, file_in)
        files_exist = glob.glob(os.path.join(dir_in, file_in))
        if len(files_exist) < 1:
            print('no files found ',files_exist, os.path.join(dir_in, file_in))

        files[strm] = sorted(glob.glob(os.path.join(dir_in, file_in)))

    for iv, var in enumerate(variables):
        strm = streams[iv]
        if strm == stream_daily:
            continue

        if len(files[strm]) > 0:
            stash_con = get_stash_constraint(stash[var])
            print ('process ',var, stash[var], stash_con, strm, files[strm][0])
            fname_base = os.path.basename(files[strm][0]).split('.')
            fname_pt = '.'.join([fname_base[0], stream_out+fname_base[1][2:], fname_base[2]])

            ftmp = fname_pt.replace(months[str(current_month)], str(current_month).zfill(2)+'01')
            fout = os.path.join(dir_write, ftmp[:-3]+'_'+var+'.nc')
            print('fout ',fout)
            fout_tmp1 = fout[:-3]+'_1.nc'
            fout_tmp2 = fout[:-3]+'_2.nc'
        
            if not os.path.exists(fout):
                if strm == stream_daily:
                    # need to cat files together, convert to nc, use cdo to interp, and then extract just this month
                    pass
                pp = iris.load_cube(files[strm], stash_con)
                time_coord = pp.coord('time')
                time_origin = time_coord.units.num2date(time_coord.points)[0]

                print('time_origin ',time_origin)
                if var in levels.keys():
                    files_1 = extract_levels(pp, fout, fout_tmp1, fout_tmp2, var)
                elif 'pr' in var:
                    files_1 = modify_time_coord(pp, fout, fout_tmp1, fout_tmp2, var)
                else:
                    pp.var_name = var
                    iris.save(pp, fout_tmp1, unlimited_dimensions = ['time'])
                    if var == 'sfcWind':
                        modify_time_coord(pp, fout, fout_tmp1, fout_tmp2, var)
                    else:
                        change_time_units(fout, fout_tmp1, fout_tmp2)
                    files_1 = [fout]

                print('var ',var, files_1)

                for f in files_1:
                    date = str(time_origin).split(' ')[0]
                    time = str(time_origin).split(' ')[1]
                    cmd = 'cdo selmon,'+str(int(current_month))+' '+f+' '+f[:-3]+'_3.nc'
                    print(cmd)
                    run_cmd(cmd)

                    if os.path.exists(f[:-3]+'_3.nc'):
                        cmd = 'mv '+f[:-3]+'_3.nc '+f
                        print(cmd)
                        run_cmd(cmd)
                        files_to_delete.append(f[:-3]+'_3.nc')

#                if os.path.exists(fout):
#                    os.remove(fout_tmp1)
#                    os.remove(fout_tmp2)
#            files_to_delete.append(f)

    
    strm = stream_daily
    for variable in variables_daily:
        #variable = variable_daily
        stash_code = [stash[variable]]
        stash_con = get_stash_constraint(stash[variable])
        level = int(variable.split('_')[-1])
        level_constraint = iris.Constraint(coord_values = {'pressure': lambda l : l == level})
        # daily data, need to get extra for interpolation
        # need to get month-1, month, month+1
        year = []; month = []
        year.append(prevcycle[:4])
        month.append(prevcycle[4:6])
        year.append(current_year)
        month.append(current_month)
        year.append(nextcycle[:4])
        month.append(nextcycle[4:6])
        lev = variable.split('_')[-1]
        files_daily = []; files_daily_base = []
        for i in range(len(year)):
            dir_in_daily = dir_in_daily_format.format(suite=um_suiteid)
            yearmon = year[i]+month[i]+'01'
            file_in1 = fname_start+'a.'+strm+'*.pp'
            print('file search ',dir_in_daily, yearmon, file_in1)
            files_exist = glob.glob(os.path.join(dir_in_daily, yearmon, file_in1))

            if len(files_exist) > 0:
                for f in files_exist:
                    if os.path.basename(f) not in files_daily_base:
                        files_daily.append(f)
                        files_daily_base.append(os.path.basename(f))
        print('files_daily ',files_daily)
        print('files_daily_base ',files_daily_base)

        if len(files_daily) > 1:
            for f in files_daily:
                if os.stat(f).st_size == 0:
                    raise Exception('File size is zero '+f)

            file_daily_3month = os.path.join(dir_write, fname_start+'a.'+strm+''+str(current_year)+months[str(current_month)]+'.pp')
            fname_base = os.path.basename(file_daily_3month).split('.')
            fname_pt = '.'.join([fname_base[0], stream_out+fname_base[1][2:], fname_base[2]])
            file_daily_out = fname_pt.replace(months[str(current_month)], str(current_month).zfill(2)+'01')
            fout = os.path.join(dir_write, file_daily_out[:-3]+'_'+variable+'.nc')
            if not os.path.exists(fout):
                c = iris.load_cube(files_daily, stash_con&level_constraint)
                print('c ',c)

                c.var_name = variable
                iris.save(c, fout[:-3]+'_tmp.nc', unlimited_dimensions = ['time'])
        
                files_to_delete.append(files_daily[0])
                files_to_delete.append(fout[:-3]+'.pp')
    
                date = str(time_origin).split(' ')[0]
                time = str(time_origin).split(' ')[1]
                cmd = 'cdo -inttime,'+date+','+time+','+str(data_freq)+'hour '+fout[:-3]+'_tmp.nc '+ fout[:-3]+'_tmp1.nc'
                print(cmd)
                run_cmd(cmd)
                files_to_delete.append(fout[:-3]+'_tmp.nc')
                files_to_delete.append(fout[:-3]+'_tmp1.nc')

                cmd = 'cdo selmon,'+str(int(current_month))+' '+fout[:-3]+'_tmp1.nc '+fout[:-3]+'_tmp2.nc'
                print(cmd)
                run_cmd(cmd)
                change_time_units(fout, fout[:-3]+'_tmp2.nc', fout[:-3]+'_tmp3.nc')

                files_to_delete.append(fout[:-3]+'_tmp2.nc')
                files_to_delete.append(fout[:-3]+'_tmp3.nc')

    for level in ['100', '850']:
        t_var = 'ta_'+level
        s_var = 'sh_'+level
        rh_var = 'rh_'+level
        t_file = os.path.join(dir_write, fname_start+'a.'+stream_out+str(current_year)+str(current_month)+'01_'+t_var+'.nc')
        s_file = os.path.join(dir_write, fname_start+'a.'+stream_out+str(current_year)+str(current_month)+'01_'+s_var+'.nc')
        rh_file = os.path.join(dir_write, fname_start+'a.'+stream_out+str(current_year)+str(current_month)+'01_'+rh_var+'.nc')
        if not os.path.exists(rh_file):
            convert_sh_to_rh(t_file, s_file, rh_file, level, t_var, s_var, rh_var)
            if level == '850':
                cmd = 'cdo setmisstoc,0 '+rh_file+' '+rh_file[:-3]+'_tmp.nc'
                print(cmd)
                run_cmd(cmd)
                cmd = 'mv '+rh_file[:-3]+'_tmp.nc '+rh_file
                print(cmd)
                run_cmd(cmd)
                #files_to_delete.append(rh_file[:-3]+'_tmp.nc')
            cmd = 'ncatted -O -a coordinates,'+rh_var+',c,c,"latitude longitude" '+rh_file
            print(cmd)
            run_cmd(cmd)
        

    for f in files_to_delete:
        if os.path.exists(f):
            print('delete ',f)

    #for st in stash:
    #    strm = stream[st]
    #    for f in files[strm]:
    #        print('get data, remove source ',f)
    #        if os.path.exists(f):
    #            os.remove(f)

if __name__ == '__main__':
    main()
