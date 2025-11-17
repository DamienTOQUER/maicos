Basic usage
===========

Loading library
-----------------

.. include:: res-tabs/splitted-usage-python-basic.rst
  :start-after: # .. start_basic_initial
  :end-before: # .. end_basic_initial

.. tabs::

  .. group-tab:: Python

    .. include:: res-tabs/splitted-usage-python-basic.rst
      :start-after: # .. start_basic_import_py
      :end-before: # .. end_basic_import_py

  .. group-tab:: Bash

    .. include:: res-tabs/splitted-usage-bash.rst
      :start-after: .. start_basic_import_sh
      :end-before: # .. end_basic_import_sh

Opening trajectory
------------------

For this tutorial we use a system consisting of a 2D slab with 1176 water molecules
confined in a 2D slit made of NaCl atoms, where the two water/solid interfaces are
normal to the axis :math:`z` as shown in the snapshot below:

.. image:: ../../static/slit-flow-dark.png
  :alt: Snapshot Slit Flow System
  :class: only-dark
.. image:: ../../static/slit-flow-light.png
  :alt: Snapshot Slit Flow System
  :class: only-light

An acceleration :math:`a = 0.05\,\text{nm}\,\text{ps}^{-2}` was applied to the water
molecules in the :math:`\boldsymbol{e}_x` direction parallel to the NaCl wall, and the
atoms of the wall were maintained frozen along :math:`\boldsymbol{e}_x`.

You can download the required :download:`topology <../examples/basics/slit_flow.tpr>`
and the :download:`trajectory <../examples/basics/slit_flow.trr>` from our website.

.. tabs::

  .. group-tab:: Python

    .. include:: res-tabs/splitted-usage-python-basic.rst
      :start-after: .. start_basic_traj_py
      :end-before: .. end_basic_traj_py

  .. group-tab:: Bash

    .. include:: res-tabs/splitted-usage-bash.rst
      :start-after: .. start_basic_traj_sh
      :end-before: # .. end_basic_traj_sh

Now, we define an atom group containing the oxygen and the hydrogen atoms (of the water molecules).

.. tabs::

  .. group-tab:: Python

    .. include:: res-tabs/splitted-usage-python-basic.rst
      :start-after: .. start_basic_group_py
      :end-before: .. end_basic_group_py

  .. group-tab:: Bash

    .. include:: res-tabs/splitted-usage-bash.rst
      :start-after: # .. start_basic_group_sh
      :end-before: # .. end_basic_group_sh

Density Profiles
----------------

Let us use the :class:`maicos.DensityPlanar` class to extract the density profile of
the ``group_H2O`` along the (default) :math:`z` axis by running the analysis.

.. tabs::

  .. group-tab:: Python

    .. include:: res-tabs/splitted-usage-python-basic.rst
      :start-after: .. start_basic_run_py
      :end-before: .. end_basic_run_py

  .. group-tab:: Bash

    .. include:: res-tabs/splitted-usage-bash.rst
      :start-after: # .. start_basic_run_sh
      :end-before: # .. end_basic_run_sh

.. include:: res-tabs/splitted-usage-python-basic.rst
  :start-after: .. start_basic_plot
  :end-before: .. end_basic_plot

Help function
-------------

.. tabs::

  .. group-tab:: Python

    .. include:: res-tabs/splitted-usage-python-basic.rst
      :start-after: .. start_basic_help_py
      :end-before: .. end_basic_help_py

  .. group-tab:: Bash

    .. include:: res-tabs/splitted-usage-bash.rst
      :start-after: # .. start_basic_help_sh
      :end-before: # .. end_basic_help_sh
