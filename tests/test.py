import os
import sys

currDir = os.getcwd()
dir = currDir.split('/test')[0]
exportText = 'export PYTHONPATH='+ dir +':${PYTHONPATH}'

File = open('setup.sh', 'w')
File.writelines(exportText)
File.close()
os.system('chmod 777 setup.sh')
print '''sky_tiling is configured.\n 
Run the following in your terminal or put in your .bashrc'''
print exportText

