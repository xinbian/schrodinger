===============================
schrodinger
===============================


.. image:: https://img.shields.io/pypi/v/schrodinger.svg
        :target: https://pypi.python.org/pypi/schrodinger

.. image:: https://img.shields.io/travis/xinbian/schrodinger.svg
        :target: https://travis-ci.org/xinbian/schrodinger

.. image:: https://readthedocs.org/projects/schrodinger/badge/?version=latest
        :target: https://schrodinger.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/xinbian/schrodinger/shield.svg
     :target: https://pyup.io/repos/github/xinbian/schrodinger/
     :alt: Updates

.. image:: https://coveralls.io/repos/github/xinbian/schrodinger/badge.png?branch=master
      :target: https://coveralls.io/github/xinbian/schrodinger?branch=master

This is a project of Schrodinger equation

* Free software: MIT license
* Documentation: https://schrodinger.readthedocs.io.


HOW TO USE
---------------
* Main code: run the /schrodinger/schrodinger.py use python
* Unit test: use the following script in folder /tests/ : nosetests --with-cov --cov-config .coveragerc --cover-html
* Some important variables are period, period of the input function; time, independable variable domain; y, input wave function; resol, basis function set length; vx, potential energy; c, constant; basis, basis type.

Features
--------

* The code mainly solves two problems. First, given a function, it can output the function operated by Hamilton operator. Second, the code uses variational method to find the ground state and plot the ground state (corresponding to smallest eigenvalue)
* Fourier series and Legendre polynomials are used as basis functions
* For simplicity, potential energy function can only be constant
* If potential function is a constant, the ground state should be 1 (after normalization), the default output has the same result. However, if we change the resol (resolution), the ground state of Legendre will change. Sometimes it can output the right results, sometimes not. I think the right results can be extracted when an infinite resolution is used. 

Methodology
-----------------
* Variational method is used. The problem can be simplified to solve the eigenvalue equation 


   **H** • **C** =λ • **S** • **C**

, where H_nm=<φ_i | H |φ_j >, S_nm=< φ_i | φ_j >. 

More details can be found in `variational method <http://www.physics.metu.edu.tr/~hande/teaching/741-lectures/lecture-01.pdf>`_.


TODO
--------
* Modify potential energy to a function

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
