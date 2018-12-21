#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

from setuptools import find_packages
from setuptools import setup

setup(
  name = 'pyproteinsExt',
  version = '0.5',
  license='BSD',
  description = 'Extending pyproteins for bioinformatics tools&services',
  author = 'Guillaume Launay, Juliette Martin et al',
  author_email = 'pitooon@gmail.com, juliette.martin@ibcp.fr',
  url = 'https://github.com/glaunay/pyproteinsExt', # use the URL to the github repo
  packages=find_packages('src'),
  package_dir={'': 'src'},
  include_package_data=True,
  zip_safe=False,
  py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
  download_url = 'https://github.com/glaunay/pyproteinsExt/tarball/0.5', # I'll explain this in a second
  keywords = ['protein', 'sequence'], # arbitrary keywords
  classifiers = [],
  install_requires=[
          'pyproteins', 'networkx', 'bs4', 'biopython'
      ],
   package_data = {
   'pyproteinsExt': ['external/submitPsipred.sh','static/psicquicRegistryDefault.xml']
   },
  #data_files=[
  #          ('external', ['external/pathos.tar.bz']),
  #          ('bin', ['bin/module1.py']),
  #          ('conf',['conf/confModule1.json'])
  #    ]
  #  dependency_links = [
  #      "http://dev.danse.us/trac/pathos"
  #  ]
)
