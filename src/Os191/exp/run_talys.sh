#!/usr/bin/sh

ln -sf ~/talys/structure/density/ground/goriely/Os_modified.tab ~/talys/structure/density/ground/goriely/Os.tab

talys < talys.inp > talys.out 

ln -sf ~/talys/structure/density/ground/goriely/Os_default.tab ~/talys/structure/density/ground/goriely/Os.tab