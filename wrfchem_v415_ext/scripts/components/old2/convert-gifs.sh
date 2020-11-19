#!/usr/bin/env sh
d=45
convert -delay $d -loop 0 ../imgs/aml-e_*.png ../imgs/aml-e.gif
convert -delay $d -loop 0 ../imgs/aml-ic_*.png ../imgs/aml-ic.gif
convert -delay $d -loop 0 ../imgs/rain-ir_*.png ../imgs/rain-ir.gif
convert -delay $d -loop 0 ../imgs/rain-e_*.png ../imgs/rain-e.gif
convert -delay $d -loop 0 ../imgs/t2-p-ir_*.png ../imgs/t2-p-ir.gif
convert -delay $d -loop 0 ../imgs/t2-p-e_*.png ../imgs/t2-p-e.gif
convert -delay $d -loop 0 ../imgs/wind-ir_*.png ../imgs/wind-ir.gif
convert -delay $d -loop 0 ../imgs/wind-e_*.png ../imgs/wind-e.gif
convert -delay $d -loop 0 ../imgs/rh-e_*.png ../imgs/rh-e.gif
convert -delay $d -loop 0 ../imgs/rh-ir_*.png ../imgs/rh-ir.gif
convert -delay $d -loop 0 ../imgs/vert_*.png ../imgs/vert.gif
curl --upload-file "../imgs/{aml-e,aml-ic,rain-ir,rain-e,t2-p-ir,t2-p-e,wind-ir,wind-e,rh-ir,rh-e,vert}.gif" --user f112266:head+mace+FPT ftp://www.macehead.org/webspace/siteapps/5379/htdocs/images/stories/

