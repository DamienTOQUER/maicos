# Basic usage - Command line
# ==========================
# .. start_basic_import_sh
# 
# MAICoS can be used directly from the command line (cli). Using cli instead of a Jupyter
# notebook can sometimes be more comfortable, particularly for lengthy analysis. The cli
# in particular is handy because it allows for updating the analysis results during the
# run. You can specify the number of frames after the output is updated with the
# ``-concfreq`` flag.
#
# Note that in this documentation, we almost exclusively describe the use of MAICoS from
# the python interpreter, but all operations can be equivalently performed from the cli.
#
# Any maicos module can be accessed simply by using the module name, for example in this tutorial we will use :class:`maicos.DensityPlanar`:

maicos densityplanar
# .. end_basic_import_sh
# .. start_basic_traj_sh
# %%
# We can specify the topology and trajectory using the '-s' and '-f' flags.

maicos densityplanar -s slit_flow.tpr\
                     -f slit_flow.trr
# .. end_basic_traj_sh
# .. start_basic_group_sh
# %%
# Using the '-atomgroup' flag.

maicos densityplanar -s slit_flow.tpr\
                     -f slit_flow.trr\
                     -atomgroup 'type OW HW'

# .. end_basic_group_sh
# .. start_basic_run_sh
# %%
#
# The density profile has been written in a file named ``density.dat`` in the current
# directory. The written file starts with the following lines

head -n 20 density.dat

# %%
#
# For lengthy analysis, use the ``concfreq`` option to update the result during the run

maicos densityplanar -s slit_flow.tpr \
                     -f slit_flow.trr \
                     -atomgroup 'type OW HW' \
                     -concfreq '10'
# .. end_basic_run_sh
# .. start_basic_plot_sh
# %%
# Using gnuplot:

echo " \
    set xlabel 'z coordinate, (Å)'; \
    set ylabel 'density H2O (u.Å⁻³)'; \
    plot 'density.dat' using (column(1)):(column(2)):(5*column(3)) with yerrorlines title '' \
" | gnuplot -persist || true
# .. end_basic_plot_sh
# .. start_basic_help_sh
# %%
#
# The general help of MAICoS can be accessed using

maicos -h

# %%
#
# Package-specific page can also be accessed from the cli

maicos densityplanar -h
# .. end_basic_help_sh