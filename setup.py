from setuptools import setup

package = 'dnds'
version = '2.1'

setup(
    name=package,
    version=version,
    description="Calculate dN/dS ratio precisely (Ka/Ks) using a codon-by-codon counting method.",
    long_description=open("README.rst").read(),
    author="Adel Qalieh",
    author_email="adelq@med.umich.edu",
    url="https://github.com/adelq/dnds",
    license="MIT",
    py_modules=['dnds', 'codons'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ])
