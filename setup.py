#!/usr/bin/env python

#
from distutils.core import setup
from setuptools import find_packages

#
setup(name='gw-phenom',
      version='1.0',
      description='Phenomenological Waveform Models',
      author='Sebastian Khan',
      author_email="KhanS22@Cardiff.ac.uk",
      packages=find_packages(),
      url='https://github.com/Cyberface/phenom',
      download_url='https://github.com/Cyberface/phenom/archive/master.zip',
      scripts = [
          'bin/phenom-ringdown-frequency',
          'bin/phenom-chirp-time'
         ]
     )
