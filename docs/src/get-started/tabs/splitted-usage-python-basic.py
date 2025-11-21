#!/usr/bin/env python
#
# Copyright (c) 2025 Authors and contributors
# (see the AUTHORS.rst file for the full list of names)
#
# Released under the GNU Public Licence, v3 or any higher version
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Basic usage - Python interpreter - Splitted version
===================================================
# .. start_basic_initial
To follow this tutorial, it is assumed that MAICoS has been :ref:`installed
<label_installation>` on your computer.

MAICoS heavily depends on the `MDAnalysis`_ infrastructure for trajectory loading and
atom selection. Here we will only cover a small aspects of the capabilities of
`MDAnalysis`_. If you want to learn more about the library, take a look at their
`documentation <https://docs.mdanalysis.org/stable/index.html>`_.

.. _`MDAnalysis`: https://www.mdanalysis.org

Only one MAICoS analysis modules is used in this tutorial
:class:`maicos.DensityPlanar` but all modules follow the same structure:

1. load your simulation data into an :class:`MDAnalysis.core.universe.Universe`
2. define analysis parameters like bin width or the direction of the analysis
3. after the analysis was succesful, access all results in a
   :class:`MDAnalysis.analysis.base.Results` of the analysis object.

Note that some of the calculations may contain pitfall, such as dielectric profiles
calculation. Potential pitfalls and best practices are listed in the
:ref:`userdoc-how-to` section.
# .. end_basic_initial
# .. start_basic_import_py
To start, let us first import Matplotlib, MDAnalysis and MAICoS
"""  # noqa: D415

# %%

import matplotlib.pyplot as plt
import MDAnalysis as mda

import maicos

# .. end_basic_import_py
# %%
# Load Simulation Data
# --------------------
#
# For this tutorial we use a system consisting of a 2D slab with 1176 water molecules
# confined in a 2D slit made of NaCl atoms, where the two water/solid interfaces are
# normal to the axis :math:`z` as shown in the snapshot below:
#
# .. image:: ../../../static/slit-flow-dark.png
#   :alt: Snapshot Slit Flow System
#   :class: only-dark
#
# .. image:: ../../../static/slit-flow-light.png
#   :alt: Snapshot Slit Flow System
#   :class: only-light
#
# An acceleration :math:`a = 0.05\,\text{nm}\,\text{ps}^{-2}` was applied to the water
# molecules in the :math:`\boldsymbol{e}_x` direction parallel to the NaCl wall, and the
# atoms of the wall were maintained frozen along :math:`\boldsymbol{e}_x`.
# .. start_basic_traj_py
# We first create an :class:`MDAnalysis.core.universe.Universe` by loading a topology
# and trajectory from disk.

u = mda.Universe("slit_flow.tpr", "slit_flow.trr")
# %%
# Let us print a few information about the trajectory:

print(f"Number of frames in the trajectory is {u.trajectory.n_frames}.")
timestep = round(u.trajectory.dt, 2)
print(f"Time interval between two frames is {timestep} ps.")
total_time = round(u.trajectory.totaltime, 2)
print(f"Total simulation time is {total_time} ps.")
# %%
# .. end_basic_traj_py
#
# Now, we define an atom group containing the oxygen and the hydrogen atoms:
# .. start_basic_group_py
group_H2O = u.select_atoms("type OW HW")
# %%
# Let us print a few information about the groups

print(f"Number of water molecules is {group_H2O.n_atoms // 3}.")

# %%
# .. end_basic_group_py
#
# Density Profiles
# ----------------
#
# Let us use the :class:`maicos.DensityPlanar` class to extract the density profile of
# the ``group_H2O`` along the (default) :math:`z` axis by running the analysis:
# .. start_basic_run_py
dplan = maicos.DensityPlanar(group_H2O).run()

# %%
# The warning starting with *Unwrapping* is perfectly normal and can be ignored for now.
# Let us extract the bin coordinates :math:`z`, the averaged density profile and its
# uncertainty estimated by MAICoS from the ``results`` attribute:

zcoor = dplan.results.bin_pos
dens = dplan.results.profile
uncertainity = dplan.results.dprofile

# %%
# The density profile is given as a 1D array, let us look at the 10 first lines:

print(dens[:10])

# %%
# By default the ``bin_width`` is 1 Å, and the unit is atomic mass per :math:`Å^3`
# (:math:`\text{u}/\text{Å}^3`).
# .. end_basic_run_py
# .. start_basic_plot_py
# Using Matplotlib:

fig, ax = plt.subplots()

ax.errorbar(zcoor, dens, 5 * uncertainity)

ax.set_xlabel(r"z coordinate ($\rm Å$)")
ax.set_ylabel(r"density H2O ($\rm u \cdot Å^{-3}$)")

fig.show()

# For this example we scale the error by 5 to be visible in the plot. More detail on the
# uncertainty estimation can be found in the advanced usages section.

# %%
# .. end_basic_plot_py
# .. start_basic_help_py
#
# The general help of MAICoS can be accessed using

help(maicos)
# %%
# Package-specific page can also be accessed

help(maicos.DensityPlanar)

# %%
# .. end_basic_help_py
