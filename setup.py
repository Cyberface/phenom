#!/usr/bin/env python

#
from distutils.core import setup
from setuptools import find_packages

#
setup(name='phenom',
      version='1.0',
      description='Phenomenological Waveform Models',
      author='Sebastian Khan',
      author_email='sebastian,khan@LIGO.org',
      packages=find_packages(),
      url='https://github.com/Cyberface/phenom',
      download_url='https://github.com/Cyberface/phenom/archive/master.zip',
     )
