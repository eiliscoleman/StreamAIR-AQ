#!/usr/bin/env python3
import datetime as dt
from datetime import timedelta
import glob as g
import json
import os


def get_dates(config):
    idx_date = {
        'm': slice(-9, -7),
        'd': slice(-7, -5),
        'H': slice(-5, -3)
    }
    f_paths = sorted(g.glob(os.path.join(config['input-wps'], '*')))
    f_start, f_end = (os.path.basename(f_paths[0]),
                      os.path.basename(f_paths[-1]))
    print(f_end)
    date = dt.datetime.now()
    #commented out part will create namelist that'll run simulation for whatever available input data exists
    #date_start, date_end = (date.replace(month=int(f_start[idx_date['m']]),
                                         #day=int(f_start[idx_date['d']]),
                                         #hour=int(f_start[idx_date['H']]),
                                         #minute=0, second=0, microsecond=0),
                            #date.replace(month=int(f_end[idx_date['m']]),
                                         #day=int(f_end[idx_date['d']]),
                                         #hour=int(f_end[idx_date['H']]),
                                         #minute=0, second=0, microsecond=0))
    date_start = date.replace(month=int(f_start[idx_date['m']]),
                                         day=int(f_start[idx_date['d']]),
                                         hour=int(f_start[idx_date['H']]),
                                         minute=0, second=0, microsecond=0)
    #let end date of the simulation be 4 days after start
    date_end=date_start+timedelta(days=4)
    print(date_end)
    print(date_start+timedelta(days=4))
    return date_start, date_end


def update_namelist_val(line, value, times=1):
    return '{} = {}\n'.format(line.split('=')[0].rstrip(),
                              ', '.join([value] * times))


def update_namelist_wps(config):
    max_dom = 0
    date_start, date_end = get_dates(config)
    with open(config['namelist']['wps']) as f:
        namelist = f.readlines()

    for line in namelist:
        if 'max_dom' in line:
            max_dom = int(line.split('=')[-1])
            break

    for i, line in enumerate(namelist):
        if 'start_date' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_start.strftime("'%Y-%m-%d_%H:%M:%S'"),
                times=max_dom
            )
        elif 'end_date' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_end.strftime("'%Y-%m-%d_%H:%M:%S'"),
                times=max_dom
            )

    with open(config['namelist']['wps'], 'w') as f:
        f.writelines(namelist)


def update_namelist_wrf(config):
    max_dom = 0
    date_start, date_end = get_dates(config)
    with open(config['namelist']['wrf']) as f:
        namelist = f.readlines()

    for line in namelist:
        if 'max_dom' in line:
            print(line)
            ll=line.rstrip(",\n")
            print(ll)
            max_dom = int(ll.split('=')[-1])
            break

    for i, line in enumerate(namelist):
        if 'start_year' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_start.strftime('%Y'),
                times=max_dom
            )
        elif 'start_month' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_start.strftime('%m'),
                times=max_dom
            )
        elif 'start_day' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_start.strftime('%d'),
                times=max_dom
            )
        elif 'start_hour' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_start.strftime('%H'),
                times=max_dom
            )
        elif 'end_year' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_end.strftime('%Y'),
                times=max_dom
            )
        elif 'end_month' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_end.strftime('%m'),
                times=max_dom
            )
        elif 'end_day' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_end.strftime('%d'),
                times=max_dom
            )
        elif 'end_hour' in line:
            namelist[i] = update_namelist_val(
                line=line,
                value=date_end.strftime('%H'),
                times=max_dom
            )

    with open(config['namelist']['wrf'], 'w') as f:
        f.writelines(namelist)


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    #update_namelist_wps(config)
    update_namelist_wrf(config)


if __name__ == '__main__':
    main('./config.json')
