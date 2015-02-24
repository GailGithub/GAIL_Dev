#!/bin/sh

cd /home/gail/GAIL_tests/repo/gail-development/
# GIT
# Pulling the latest repository before testing
# /usr/local/bin/git --git-dir /home/gail/GAIL_tests/repo/gail-development/.git checkout .
/usr/local/bin/git --git-dir /home/gail/GAIL_tests/repo/gail-development/.git pull
# /usr/local/bin/git --git-dir /home/gail/GAIL_tests/repo/gail-development/.git checkout .

# MATLAB
# Generate files necessary for creating HTML documentation
cd GAIL_Matlab/Documentation/
g++ -std=c++11 autodoc.cpp -o autodoc
./autodoc
git add /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Documentation
git add /home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Documentation/html
git commit -m 'Update automatically generated Matlab files for HTML documentation'
rm -f autodoc
git push origin develop

# matlab < ./Developers_only/GAIL_Publish.m 

 