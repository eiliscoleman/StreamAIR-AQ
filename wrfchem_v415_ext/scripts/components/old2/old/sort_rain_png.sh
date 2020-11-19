#!/bin/bash
dir_out=/mnt/raid/rong-ming/web/temp
dir_finaldestination=/mnt/raid/rong-ming/web
hourinc=0
for fin_fil in ${dir_out}"/2017"*"rain_psl.png"; do
yyyy=${fin_fil#*temp/}
echo ${yyyy:0:13}
if [ $hourinc -lt 19 ]; then
       echo scrap this data ${yyyy:0:13} $hourinc
       mv -f ${dir_out}/${yyyy:0:13}* ${dir_out}/scrap
else
       echo keep this data ${yyyy:0:13} $hourinc
       mv -f ${dir_out}/${yyyy:0:13}* $dir_finaldestination
fi
hourinc=$((hourinc+1))
done
year=${yyyy:0:4}
echo $year
rsync -zavu -e "ssh -p 8222"  --include ${year}"*.jpg" --include ${year}"*.json"  --include ${year}"*.png" --exclude="*" /mnt/raid/rong-ming/web/* update@140.203.204.132:/home/www/html/rt/weather/
