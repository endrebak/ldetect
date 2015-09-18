import sys
from distutils.core import setup

if sys.version_info<(3,):
	raise Exception('Sorry, only Python v3 and higher are supported')

setup(
  name = 'ldetect',
  packages = ['ldetect', 'ldetect.baselib', 'ldetect.examples', 'ldetect.pipeline_elements'],
  version = '0.2.5',
  description = 'Package for detecting regions of linkage disequilibrium in the human genome',
  install_requires = ['numpy', 'scipy', 'matplotlib', 'commanderline'],
  author = 'Tomaz Berisa',
  author_email = 'tomaz.berisa@gmail.com',
  url = 'https://bitbucket.org/nygcresearch/ldetect', 
  download_url = 'https://bitbucket.org/nygcresearch/ldetect/get/0.2.5.tar.gz', 
  keywords = ['LD', 'linkage disequilibrium', 'genome', 'bioinformatics', 'ldetect'], 
  classifiers = [],
)
