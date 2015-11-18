#!/bin/bash
echo "#### now in `pwd` #####################";echo
echo '#### git remote -v update #############';echo;git remote -v update;echo
echo ""
echo ""
echo '#### git status -u ####################'
git status -u
echo 
echo 
echo 
echo 
echo 
echo "--"
echo '#### git commit -a -m "`date`"#########'
echo "--"
echo 
echo 
echo 
git commit -a -m "`date`"
echo "--"
echo "--"
echo '#### git pull ######################### this makes merge'
echo "--"
echo "--"
echo "--"
git pull



