#!/usr/bin/env sh
ls ../imgs/pm25-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/pm25-eu_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/pm10-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/pm10-eu_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/so2-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/so2-eu_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/o3-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/o3-eu_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/nox-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/nox-eu_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/rain-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/rain-e_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/t2-p-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/t2-p-e_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/wind-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/wind-e_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/rh-ir_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh
ls ../imgs/rh-e_*.png | head -n 19 | sed -r 's|(.*)| mv \1 ../img_pm|' | sh

d=78
convert -delay $d -loop 0 ../imgs/pm25-ir_*.png ../imgs/pm25-ir.gif
convert -delay $d -loop 0 ../imgs/pm25-eu_*.png ../imgs/pm25-eu.gif
convert -delay $d -loop 0 ../imgs/pm10-ir_*.png ../imgs/pm10-ir.gif
convert -delay $d -loop 0 ../imgs/pm10-eu_*.png ../imgs/pm10-eu.gif
convert -delay $d -loop 0 ../imgs/so2-ir_*.png ../imgs/so2-ir.gif
convert -delay $d -loop 0 ../imgs/so2-eu_*.png ../imgs/so2-eu.gif
convert -delay $d -loop 0 ../imgs/o3-ir_*.png ../imgs/o3-ir.gif
convert -delay $d -loop 0 ../imgs/o3-eu_*.png ../imgs/o3-eu.gif
convert -delay $d -loop 0 ../imgs/nox-ir_*.png ../imgs/nox-ir.gif
convert -delay $d -loop 0 ../imgs/nox-eu_*.png ../imgs/nox-eu.gif
convert -delay $d -loop 0 ../imgs/rain-ir_*.png ../imgs/rain-ir.gif
convert -delay $d -loop 0 ../imgs/rain-e_*.png ../imgs/rain-e.gif
convert -delay $d -loop 0 ../imgs/t2-p-ir_*.png ../imgs/t2-p-ir.gif
convert -delay $d -loop 0 ../imgs/t2-p-e_*.png ../imgs/t2-p-e.gif
convert -delay $d -loop 0 ../imgs/wind-ir_*.png ../imgs/wind-ir.gif
convert -delay $d -loop 0 ../imgs/wind-e_*.png ../imgs/wind-e.gif
convert -delay $d -loop 0 ../imgs/rh-e_*.png ../imgs/rh-e.gif
convert -delay $d -loop 0 ../imgs/rh-ir_*.png ../imgs/rh-ir.gif
curl --upload-file "../imgs/{pm25-ir,pm25-eu,pm10-ir,pm10-eu,so2-ir,so2-eu,o3-ir,o3-eu,nox-ir,nox-eu,rain-ir,rain-e,t2-p-ir,t2-p-e,wind-ir,wind-e,rh-ir,rh-e}.gif" --user f112266:head+mace+FPT ftp://www.macehead.org/webspace/siteapps/5379/htdocs/images/rtdata/rtdata/

mv ../imgs/* ../img_pm/

