dN/dS Calculator
================

.. image:: https://travis-ci.org/adelq/dnds.svg?branch=master
   :target: https://travis-ci.org/adelq/dnds

.. image:: https://img.shields.io/pypi/v/dnds.svg?maxAge=2592000?style=plastic
    :target: https://pypi.python.org/pypi/dnds

Calculate dN/dS ratio precisely (Ka/Ks) using a codon-by-codon counting
method.

Usage
-----

.. code:: python

    >>> sequence_1 = "ATGCTTTTGAAATCG"
    >>> sequence_2 = "ATGCGTTCGAAGTCG"
    >>> dnds(sequence_1, sequence2)
    Fraction(38, 71)
    >>> round(float(dnds(sequence_1, sequence2)), 3)
    0.535
