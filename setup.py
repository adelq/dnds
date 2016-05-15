from setuptools import setup

package = 'dnds'
version = '1.0.1'

setup(
    name=package,
    version=version,
    description="Calculate dN/dS ratio precisely (Ka/Ks) using a codon-by-codon counting method.",
    author="Adel Qalieh",
    author_email="adelq@sas.upenn.edu",
    url="https://github.com/adelq/dnds",
    license="MIT",
    py_modules=['dnds', 'codons'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ])
