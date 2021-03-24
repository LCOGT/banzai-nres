Development
===========
BANZAI-NRES is developed on GitHub and follows a standard open-source project contribution workflow. We follow a very
similar workflow to Astropy which is documented `here <https://www.astropy.org/contribute.html>`_.

Testing
-------
To run the unit tests locally run::

  tox -e test

To run the code style checks run::

  tox -e codestyle

To test building the documentation run::

  tox -e build_docs

The package infrastructure is based on the `Astropy Package Template <https://github.com/astropy/package-template>`_
so the tox testing infrastructure to test different versions of dependencies is provided by the template and
documented there.

Testing on Real Data
--------------------
Generally, the simplest way to reduce data with BANZAI-NRES is to follow the example reduction on the next page.
For more complex setups that use RabbitMQ and celery, see the helm chart of deployment details.
