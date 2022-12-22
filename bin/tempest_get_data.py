#!/usr/bin/env python3
# (C) British Crown Copyright 2022, Met Office.
# Please see LICENSE for license details.

from datetime import datetime
import glob
import os
import subprocess

import cf_units
import iris

months = {
    "01": "jan",
    "02": "feb",
    "03": "mar",
    "04": "apr",
    "05": "may",
    "06": "jun",
    "07": "jul",
    "08": "aug",
    "09": "sep",
    "10": "oct",
    "11": "nov",
    "12": "dec",
}

# fname_out = 'atmos_{}a_6h_{}_pt-{}.nc'
stream_out = "pt"
fname_out = "{}a.pt{}_{}.nc"
reftime = "1988-9-1,0:0:0,6hour"

stash = {}
stash["psl"] = "16222"
stash["ua"] = "30201"
stash["va"] = "30202"
stash["zg"] = "30297"
stash["uas"] = "03225"
stash["vas"] = "03226"
stash["pr"] = "05216"
stash["ts"] = "00024"
# stash['rvT63'] = '30458'
# stash['rvT42'] = '30456'
stash["rh_850"] = "30296"
stash["viwve"] = "30462"
stash["viwvn"] = "30463"
stash["ta_850"] = "30294"

levels = {}
levels["ua"] = [925, 850]
levels["va"] = [925, 850]
levels["zg"] = [500, 250]
# levels['ta'] = [925, 850]

stream = {}
for var in stash.keys():
    stream[var] = "p7"
# stream['ta'] = 'pc'
streams = list(set(val for val in stream.values()))


def convert_to_date(date_string):
    year = date_string[0:4]
    if len(date_string) == 7:
        mon_str = date_string[4:]
        mon = datetime.strptime(mon_str.title(), "%b").month
        day = "01"
        datestr = year + str(mon).zfill(2) + day
        date = datetime.strptime(date_string + "01", "%Y%b%d")
    else:
        mon = date_string[4:6]
        day = date_string[6:]
        datestr = year + mon + day
        date = datetime.strptime(year + mon + day, "%Y%m%d")

    return date, datestr


def latest_mass_file(suite, stream="p7"):
    cmd = "moo ls -t moose:/crum/" + suite + "/a" + stream + ".pp"
    sts = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = sts.communicate()
    mass_list = sorted(str(out).strip().split("\\n"))
    if len(mass_list) > 0:
        print("last mass file ", mass_list[-1])
        mass_date = (mass_list[-1].split("/")[-1]).split(".")[1][2:]
        return mass_date
    else:
        return None


def latest_local_file(suite, dirname):
    cmd1 = "ls -l " + dirname + "/" + suite[2:] + "a.*.nc"
    sts1 = subprocess.Popen(
        cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    out, err = sts1.communicate()
    file_listing = str(out).split("\\n")
    file_list = []
    for f in file_listing:
        a = f.split(" ")
        if len(a) >= 8:
            fpath = a[-1]
            file_list.append(fpath)
    if len(file_list) > 0:
        file_date = os.path.basename(file_list[-1].split(".")[1][2:10])
        print("last local file ", file_list[-1])
        return file_date
    else:
        return None


def get_file(fin, dirname, strm):
    if not os.path.exists(os.path.join(dirname, fin)):
        cmd = (
            "moo get -i moose:/crum/"
            + um_suiteid
            + "/a"
            + strm
            + ".pp/"
            + fin
            + " "
            + dirname
        )
        print(cmd)
        os.system(cmd)


def get_stash_constraint(stashcode):
    st_pattern = "m01s{}i{}"
    stash = st_pattern.format(stashcode[0:2], stashcode[2:])
    stash_con = iris.AttributeConstraint(STASH=stash)
    return stash_con


def run_cmd(cmd, check=True):
    sts = subprocess.run(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=check,
    )
    print("cmd ", sts.stdout)
    if sts.stderr:
        if "Warning" in sts.stderr:
            msg = ("Warning found in cat output ", sts.stderr)
            print(msg)
        else:
            msg = ("Error found in cat output ", sts.stderr)
            raise RuntimeError(msg)
    return sts


def change_time_units(fout, fout_tmp1, fout_tmp2):
    cmd = "cdo setreftime," + reftime + " " + fout_tmp1 + " " + fout_tmp2
    print(cmd)
    run_cmd(cmd)
    # remove the leadtime variable that cdo adds
    cmd = "ncks -O -x -v leadtime " + fout_tmp2 + " " + fout
    print(cmd)
    run_cmd(cmd)


def get_environment_variables():
    """
    Get the required environment variables from the suite. A list and
    explanation of the required environment variables is included in the
    documentation.
    """
    global um_runid, um_suiteid, cylc_task_cycle_time, time_cycle, previous_cycle, tm2_cycle, next_cycle, startdate, enddate, current_year, current_month, current_day, period, cycleperiod, dir_out, dir_in, data_freq

    # try:
    #    um_runid = os.environ["RUNID_OVERRIDE"]
    # except:
    #    um_runid = os.environ["RUNID"]
    try:
        um_suiteid = os.environ["SUITEID_OVERRIDE"]
    except:
        um_suiteid = os.environ["CYLC_SUITE_NAME"]

    um_runid = um_suiteid.split("-")[1]

    cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
    time_cycle = os.environ["TIME_CYCLE"]
    startdate = os.environ["STARTDATE"]
    enddate = os.environ["ENDDATE"]
    cycleperiod = os.environ["CYCLEPERIOD"]
    dir_in = os.environ["DIR_IN"]
    dir_out = os.environ["DIR_OUT"]
    print("dir_in ", dir_in)
    data_freq = int(os.environ["DATA_FREQ"])

    current_year = time_cycle[0:4]
    current_month = time_cycle[4:6]
    current_day = time_cycle[6:8]
    period = str(current_year) + str(current_month) + str(current_day)


def extract_levels(pp, fout, fout_tmp1, fout_tmp2, var):
    """
    extract levels from file and save each level separately
    """
    for level in levels[var]:
        fout_level = fout[:-3] + "_" + str(level) + ".nc"
        if not os.path.exists(fout_level):
            constraint = iris.Constraint(
                coord_values={"pressure": lambda l: l == level}
            )
            c = pp.extract(constraint)
            var_name = var + "_" + str(level)
            c.var_name = var_name
            iris.save(c, fout_tmp1, unlimited_dimensions=["time"])
            change_time_units(fout_level, fout_tmp1, fout_tmp2)
            if os.path.exists(fout[:-3] + "_" + str(levels[var][0]) + ".nc"):
                os.remove(fout_tmp1)
                os.remove(fout_tmp2)


def limit_latitude(pp):
    lat_min = pp.coord("latitude").points.min()
    lat_max = pp.coord("latitude").points.max()

    if lat_min < -90.0 or lat_max > 90.0:
        pp.coord("latitude").bounds = None
        if pp.coord("latitude").points[0] > 90.0:
            pp.coord("latitude").points[0] = 90.0
        elif pp.coord("latitude").points[0] < -90.0:
            pp.coord("latitude").points[0] = -90.0

        if pp.coord("latitude").points[-1] > 90.0:
            pp.coord("latitude").points[-1] = 90.0
        elif pp.coord("latitude").points[-1] < -90.0:
            pp.coord("latitude").points[-1] = -90.0

        pp.coord("latitude").guess_bounds()


def modify_time_coord(pp, fout, fout_tmp1, fout_tmp2, st):
    """
    Need to adjust mean time coord to have same times as instantaneous
    Note that the mean is currently wrong - offset by 3 hours
    """
    if os.path.exists(fout):
        return
    times = pp.coord("time").points[1]
    if (times % data_freq) != 0:
        # make a new time coordinate
        calendar = pp.coord("time").units.calendar
        values = pp.coord("time").points.copy()
        values += data_freq / 2
        new_time = iris.coords.DimCoord(
            values,
            standard_name="time",
            units=cf_units.Unit(pp.coord("time").units, calendar=calendar),
        )
        print("new_time ", new_time.points)
        pp.remove_coord("time")
        pp.add_dim_coord(new_time, 0)
        pp.var_name = st
        iris.save(pp, fout_tmp1, unlimited_dimensions=["time"])
        change_time_units(fout, fout_tmp1, fout_tmp2)


def test_file_dates(suite, dirname, period):
    """
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
    """
    recent_mass_file_date = latest_mass_file(suite)
    if not recent_mass_file_date == None:
        datem, datestrm = convert_to_date(recent_mass_file_date)
    else:
        datem, datestrm = None, None

    recent_local_file_date = latest_local_file(suite, dirname)
    if not recent_local_file_date == None:
        datef, datestrf = convert_to_date(recent_local_file_date)
    else:
        datef, datestrf = None, None

    print("most recent file is ", datef, datestrf)
    print("most recent mass file is ", datem, datestrm)

    print("compare curr time ", period, " with file ", datestrf, "and mass ", datestrm)
    if datestrf is None:
        if datestrm is not None:
            condition = "GetData"
        else:
            condition = "FailAndRetry"
    else:
        if int(period) <= int(datestrf):
            # this tracking has already been done
            condition = "AlreadyComplete"
        elif int(period) > int(datestrf):
            # do this data if at least one year behind mass
            if int(datestrm[0:6]) - int(period[0:6]) > 12:
                condition = "GetData"
            else:
                condition = "FailAndRetry"
        else:
            condition = "FailAndRetry"

    return condition


def main():
    get_environment_variables()
    dir_write = os.path.join(dir_out, um_suiteid)
    if not os.path.exists(dir_write):
        os.makedirs(dir_write)

    # test what data is available in MASS and what data has been done so far
    condition = test_file_dates(um_suiteid, dir_write, period)
    print("condition ", condition)

    if condition == "GetData":
        pass
    elif condition == "AlreadyComplete":
        return
    else:
        raise Exception("Need to wait for data")

    files = {}
    for strm in streams:
        if "M" in cycleperiod:
            file_in = (
                um_runid
                + "a."
                + strm
                + ""
                + str(current_year)
                + months[str(current_month)]
                + ".pp"
            )
        elif "Y" in cycleperiod:
            file_in = um_runid + "a." + strm + "" + str(current_year) + "*.pp"
            file_out = file_in
        else:
            file_in = um_runid + "a." + strm + "" + period + ".pp"

        files_exist = glob.glob(os.path.join(dir_write, file_in))
        if len(files_exist) < 4:
            get_file(file_in, dir_write, strm)

        files[strm] = sorted(glob.glob(os.path.join(dir_write, file_in)))

    for st in stash:
        strm = stream[st]
        for f in files[strm]:
            stash_con = get_stash_constraint(stash[st])
            print(st, stash[st], stash_con)
            fname_base = os.path.basename(f).split(".")
            fname_pt = ".".join(
                [fname_base[0], stream_out + fname_base[1][2:], fname_base[2]]
            )

            if not "M" in cycleperiod:
                fout = os.path.join(dir_write, fname_pt[:-3] + "_" + st + ".nc")
            else:
                ftmp = fname_pt.replace(
                    months[str(current_month)], str(current_month).zfill(2) + "01"
                )
                fout = os.path.join(dir_write, ftmp[:-3] + "_" + st + ".nc")
            fout_tmp1 = fout[:-3] + "_1.nc"
            fout_tmp2 = fout[:-3] + "_2.nc"
            if not os.path.exists(fout):
                pp = iris.load_cube(f, stash_con)
                if st in levels.keys():
                    extract_levels(pp, fout, fout_tmp1, fout_tmp2, st)
                elif "pr" in st:
                    modify_time_coord(pp, fout, fout_tmp1, fout_tmp2, st)
                else:
                    pp.var_name = st
                    iris.save(pp, fout_tmp1, unlimited_dimensions=["time"])
                    change_time_units(fout, fout_tmp1, fout_tmp2)
                if os.path.exists(fout):
                    os.remove(fout_tmp1)
                    os.remove(fout_tmp2)

    for st in stash:
        strm = stream[st]
        for f in files[strm]:
            if os.path.exists(f):
                os.remove(f)


if __name__ == "__main__":
    main()
