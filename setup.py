from setuptools import setup

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='RCMAP',
    version='1.2.0',
    packages=['RCMAP'],
    package_dir={'': 'src'},
    url='https://github.com/fplewniak/RCMAP',
    license='CeCILL FREE SOFTWARE LICENSE AGREEMENT Version 2.1 dated 2013-06-21',
    description='Residue Conservation in Multiple Alignment of Proteins ',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Development Status :: 4 - Beta'
        ],
    python_requires='>=3.6', install_requires=['biopython', 'scipy'],
    entry_points={
        "console_scripts": [
            "evaluate_seq = RCMAP.evaluation_seq:main"
            ]
        }
    )
