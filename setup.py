from setuptools import setup

VERSION = '0.0.1' 
DESCRIPTION = 'A spectrum viewing tool'
LONG_DESCRIPTION = 'A gui based tool to view and bin spectra of planetary surfaces by selecting regions on a map.'

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='specview',
    url='add later',
    author='M. Ryleigh Davis',
    author_email='rdavis@caltech.edu',
    # Needed to actually package something
    packages=['specview'],
    # Needed for dependencies
    install_requires=['numpy', 'scipy', 'matplotlib', 'pillow', 'cartopy'],
    #package_data={'spectools': ['linear_models/*.hdf5']}
    # *strongly* suggested for sharing
    version=VERSION,
    # The license can be anything you like
    license='TBD',
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
        ]
    )
