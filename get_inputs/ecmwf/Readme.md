Getting the ECMWF data dissemnation is done in the crontab, where the IKD files are pulled from the ECMWF server at 7:35 and 19:35 daily.


# Get the ECMWF dissemination data every 12 hours
30 07 * * * rm  /mnt/raid/wrf-chem/ECMWF-op/IKD*
35 07 * * * rsync -az -e 'ssh -p 8222' ecmwf@140.203.204.134:'$(ls -t /home/ecmwf/data/IKD* | head -49)' /mnt/raid/wrf-chem/ECMWF-op/
40 07 * * * cp /mnt/raid/wrf-chem/ECMWF-op/* /mnt/raid/wrf-chem/ECMWF-op-VOLCEX/
