#!/usr/bin/env python
# encoding: UTF-8

from setuptools import setup, find_packages

setup(name='pyGaussDCA', version='0.1',
      description='',
      url='',
      author='David Menéndez Hurtado',
      author_email='davidmenhur@gmail.com',
      license='BSD',
      packages=find_packages(),
      package_data={'gaussdca.data': ['gaussdca/test/data/*.a3m']},
      include_package_data=True,
      requires=['numpy', 'scipy', 'pythran'],
      classifiers=['Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Intended Audience :: Science/Research'],
      zip_safe=False)
