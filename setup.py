#! /usr/bin/env python
#
import os
import glob

DESCRIPTION = "Tools for generating ranked tiles of a square FOV telescopes"
LONG_DESCRIPTION = """\
Given a GW sky-localziation maps obtained from BAYESTAR or LALInference, this 
tool finds the ranked-tiles for observation using FOV of the telescope.
More information: A&A, 592 (2016) A82, arXiv:1511.02673
"""

DISTNAME = 'sky_tiling'
AUTHOR = 'Shaon Ghosh'
MAINTAINER = 'Shaon Ghosh' 
MAINTAINER_EMAIL = 'ghosh4@uwm.edu'
URL = 'https://github.com/shaonghosh/sky_tiling/'
LICENSE = 'GPLv3'
DOWNLOAD_URL = 'https://github.com/shaonghosh/sky_tiling/'
VERSION = '0.1.0'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
   install_requires = []

   # Make sure dependencies exist.
   try:
       import astropy
   except ImportError:
       install_requires.append('astropy')
   try:
       import healpy
   except ImportError:
       install_requires.append('healpy')
   try:
       import pickle
   except ImportError:
       install_requires.append('pickle')

   return install_requires


if __name__ == "__main__":

    install_requires = check_dependencies()
    scripts = glob.glob('bin/*') + glob.glob('utilities/*') 

    if _has_setuptools:
        packages = find_packages()
        print packages
    else:
        # This should be updated if new submodules are added
        packages = [
            'sky_tiling', 
            'sky_tiling.tile_pixel_maps', 
            'sky_tiling.tile_center_files',
            'sky_tiling.utilities']
    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
	      scripts=scripts,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          install_requires=install_requires,
          packages=packages,
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 2.7',
              'License :: OSI Approved :: GPLv3',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
      )
dir = os.getcwd()
exportText = 'export PYTHONPATH='+ dir +':${PYTHONPATH}'
print '''\n***** sky_tiling is configured *****.
Run the following in your terminal or put in your .bashrc'''
print exportText

