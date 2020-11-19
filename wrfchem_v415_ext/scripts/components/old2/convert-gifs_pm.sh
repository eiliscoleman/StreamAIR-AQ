#!/usr/bin/env sh
ls ../imgs/*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
d=78
convert -delay $d -loop 0 ../imgs/pm25-ir_*.png ../imgs/pm25-ir.gif
convert -delay $d -loop 0 ../imgs/pm25-eu_*.png ../imgs/pm25-eu.gif
convert -delay $d -loop 0 ../imgs/rain-ir_*.png ../imgs/rain-ir.gif
convert -delay $d -loop 0 ../imgs/rain-e_*.png ../imgs/rain-e.gif
convert -delay $d -loop 0 ../imgs/t2-p-ir_*.png ../imgs/t2-p-ir.gif
convert -delay $d -loop 0 ../imgs/t2-p-e_*.png ../imgs/t2-p-e.gif
convert -delay $d -loop 0 ../imgs/wind-ir_*.png ../imgs/wind-ir.gif
convert -delay $d -loop 0 ../imgs/wind-e_*.png ../imgs/wind-e.gif
convert -delay $d -loop 0 ../imgs/rh-e_*.png ../imgs/rh-e.gif
convert -delay $d -loop 0 ../imgs/rh-ir_*.png ../imgs/rh-ir.gif
curl --upload-file "../imgs/{pm25-ir,pm25-eu,rain-ir,rain-e,t2-p-ir,t2-p-e,wind-ir,wind-e,rh-ir,rh-e}.gif" --user f112266:head+mace+FPT ftp://www.macehead.org/webspace/siteapps/5379/htdocs/images/rtdata/rtdata/

mv ../imgs/* ../img_pm/
