from setuptools import setup, find_packages

setup(
    name='uniqprimer',
    version='0.5.3',
    description='A Python tool for finding primers unique to a genome',
    author='John Herndon, Diego Gutierrez',
    author_email='johnlherndon@gmail.com, diego.ernesto.gutierrez@colostate.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'biopython>=1.78',
    ],
    entry_points={
        'console_scripts': [
            'uniqprimer=uniqprimer.main:main',
        ],
    },
    python_requires='>=3.6',
)
