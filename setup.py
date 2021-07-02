#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: setup.py
Author: syliu @ liusongyu@cau.edu.cn
Created on: 2021-07-1 10:33:30
'''

from setuptools import setup

setup(name='MRBIGR',
    python_requires='>3.7',
    version='1.0',
    description='MRBIGR: Mendelian Randomization-Based Inference of Genetic Regulation',
    url='https://github.com/liusy-jz/MRBIGR.git',
    author='Songyu Liu',
    author_email='',
    license='GPLv3',
    packages=['mrbigr'],
    scripts=['MRBIGR.py'],
    install_requires=[
        'cython',
        'pandas_plink',
        'numpy',
        'pandas',
        'scipy',
        'scikit-learn',
    ]
)
