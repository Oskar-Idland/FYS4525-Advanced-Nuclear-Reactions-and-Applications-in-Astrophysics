#!/usr/bin/sh

ln -sf ~/talys/structure/density/ground/goriely/Sn_modified.tab ~/talys/structure/density/ground/goriely/Sn.tab
talys < talys.inp > talys.out 
ln -sf ~/talys/structure/density/ground/goriely/Sn_default.tab ~/talys/structure/density/ground/goriely/Sn.tab