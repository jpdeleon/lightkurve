#!/usr/bin/env python
import os
import sys
from setuptools import setup

# Prepare and send a new release to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/lightkurve*")
    sys.exit()


# Load the __version__ variable without importing the package already
exec(open('lightkurve/version.py').read())

setup(name='lightkurve2',
      version=__version__,
      description="scripts written based on forked lightkurve",
      long_description=open('README_orig.rst').read(),
      author='KeplerGO and Jerome de Leon',
      author_email = 'jpdeleon@astron.s.u-tokyo.ac.jp',
      url = 'https://github.com/jpdeleon/exofop',
      license='MIT',
      #package_dir={"lightkurve": "lightkurve"},
      packages=['lightkurve','K2tools'],
      package_data={"data": "Yu18candidates"},
      scripts=['scripts/make_mask','scripts/tpf2lc'],
      include_package_data=True,
      install_requires=['numpy>=1.11', 'astropy>=1.3', 'scipy>=0.19.0',
                        'matplotlib>=1.5.3', 'tqdm', 'oktopus', 'bs4',
                        'requests', 'astroquery>=0.3.7', 'pandas', 'scipy',
                        'sklearn', 'photutils', 'skimage'],
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'pytest-cov', 'pytest-remotedata'],
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
    )
