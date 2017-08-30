#!/usr/bin/env python3

from setuptools import setup, find_packages

desc = '''Removes specific viral sequences from a 
          host genome and generates a SNAP index'''

setup(
        name='host_snap',
        version='0.1',
        packages=find_packages(),
        install_requires=['biopython', 'pandas', 'numpy'],
        author='Dillon Barker',
        author_email='dillon.barker@canada.ca',
        description=desc,
        license='GPL 3+',
        keywords='bioinformatics science biology',
        url='https://github.com/dorbarker/host_snap',

        entry_points={
            'console_scripts': [
                'host-snap = host_snap.host_snap:main'
                ]}
)
