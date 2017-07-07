import os
import sys

currDir = os.getcwd()
dir = currDir.split('/test')[0]
bindir = currDir.split('/test')[0] + '/bin'
exportText1 = 'export PYTHONPATH='+ dir +':${PYTHONPATH}'
exportText2 = 'export PYTHONPATH='+ bindir +':${PYTHONPATH}'

File = open('setup.sh', 'w')
File.writelines(exportText1 + '\n')
File.writelines(exportText2)
File.close()
os.system('chmod 777 setup.sh')
print '''sky_tiling is configured.\n 
Run the following in your terminal or put in your .bashrc'''
print exportText1
print exportText2


