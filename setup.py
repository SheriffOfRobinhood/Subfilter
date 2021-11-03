from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='Subfilter',
    url='https://github.com/ReadingClouds/Subfilter',
    author='Peter Clark',
    author_email='p.clark@reading.ac.uk',
    contributors='Todd Jones',
    # Needed to actually package something
    packages=['subfilter', 
              'subfilter/cli',
              'subfilter/io',
              'subfilter/thermodynamics',
              'subfilter/utils',
              ],
    # Needed for dependencies
    install_requires=['numpy', 'scipy', 'dask', 'xarray'],
    # *strongly* suggested for sharing
    version='0.5.2',
    # The license can be anything you like
    license='MIT',
    description='python code to compute sub-filter quantities from MONC output.',
    # We will also need a readme eventually (there will be a warning)
    long_description=open('README.md').read(),
)