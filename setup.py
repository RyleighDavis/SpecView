from setuptools import setup

VERSION = '0.0.1' 
DESCRIPTION = 'A spectrum viewing tool'
LONG_DESCRIPTION = 'A gui based tool to view and bin spectra of planetary surfaces by selecting regions on a map.'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="specview", 
        version=VERSION,
        url="https://github.com/RyleighDavis/specview",
        author="M. Ryleigh Davis",
        author_email="<rdavis@caltech.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        install_requires=[numpy, scipy, matplotlib, pillow, cartopy], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'specview', 'spectral viewing tool'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
        ]
)
