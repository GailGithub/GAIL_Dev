# This script can be used to remove files from the commit history. Note that the files should already be removed from the current repository by using the "git rm" command.
# To run this script in a Linux system, navigate to the folder that contains this file and run 'python PurgeFile.py'.
# Manually force updating the remote repository by running 'git push --force origin master'. 

import os

FileName = raw_input('Enter the name of the file to be purged (including the path): ')
os.system('git gc')
commits = os.popen('git log --pretty=oneline --branches -- ' + FileName).readlines()
try:
    assert len(commits), 'No such file in the commit history. The program will quit.'
except AssertionError, args:
    print '%s: %s' % (args.__class__.__name__, args)
    exit()
os.system('git filter-branch --index-filter \'git rm --cached --ignore-unmatch ' + FileName + '\'')
os.system('rm -Rf .git/refs/original')
os.system('rm -Rf .git/logs/')
os.system('git reflog expire --expire=now --all')
os.system('git gc --prune=now')
