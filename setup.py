# setup.py

from setuptools import setup, find_packages

setup(
    name='social_network_analysis',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'networkx',
        'matplotlib',
        'numpy',
        'scipy',
        'openpyxl',  # Per leggere Excel
        'requests',
    ]
)