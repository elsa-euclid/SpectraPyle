Quick Start
===========

Installation
------------

Install the latest stable release from PyPI:

.. code-block:: bash

   pip install spectraPyle

This installs SpectraPyle together with the dependencies required by its
interactive GUI.

If you also need a local JupyterLab environment and Voilà, install the
``notebook`` extras:

.. code-block:: bash

   pip install "spectraPyle[notebook]"

For development or an editable installation from the source repository, see
the development documentation.

Try SpectraPyle immediately
---------------------------

SpectraPyle includes a ready-to-use sample of 100 Euclid NISP Q1 spectra,
together with a pre-configured ``default.gui``. This quick-start sample is
intended for users who want to see the complete stacking workflow and its
output immediately, before configuring their own data.

In a Jupyter notebook, launch the GUI with:

.. code-block:: python

   from spectraPyle.gui import start

   start()

When SpectraPyle is launched for the first time in a working directory, it
creates:

.. code-block:: text

   configs/GUI/default.gui

Load this configuration from the GUI and run SpectraPyle. The bundled catalogue
and example spectra are located automatically, so no input paths need to be
configured.

The example performs a complete stack of the supplied Euclid spectra and
generates the standard SpectraPyle outputs, including the interactive
stacked-spectrum visualization.

The example data are installed as read-only package resources and are not
copied into your project directory. Output files are instead written to the
local ``./outputs`` directory.

Once you have verified the workflow, you can modify the GUI configuration for
your own catalogue and spectra and save project-specific GUI, JSON, or YAML
configurations under ``./configs/``.

Jupyter GUI
-----------

The GUI can be launched directly from any Jupyter environment:

.. code-block:: python

   from spectraPyle.gui import start

   start()

The interface is rendered directly inside the notebook using ``ipywidgets``;
no local server or separate browser window is required.

ESA DataLabs users can use the same interface directly from a notebook. For
Euclid internal data and the DataLabs-specific workflow, see :doc:`datalabs`.

Local GUI via Voilà
-------------------

For a standalone local Voilà interface, the source repository also provides:

.. code-block:: bash

   python notebooks/run_gui.py

This requires SpectraPyle to be installed from the source repository with the
``notebook`` extras.

The launcher starts a Voilà server and opens the GUI in your default browser.
If the browser does not open automatically, use the local URL printed in the
terminal. Press ``Ctrl+C`` in the terminal to stop the server.

Running via CLI
---------------

After installation, the standard command-line entry point is:

.. code-block:: bash

   spectrapyle --config path/to/config.yaml

When working from the source repository, a helper script is also available:

.. code-block:: bash

   python notebooks/run_cli.py --config path/to/config.yaml [--log-level INFO]

The helper provides automatic timestamped logging and supports YAML/JSON
configuration files. Available log levels include ``DEBUG``, ``INFO``, and
``WARNING``.

Override individual configuration keys at runtime, for example:

.. code-block:: bash

   spectrapyle --config config.yaml --instrument.grisms '["red","blue"]'

Post-stacking Analysis
----------------------

**Stacked spectrum** — SpectraPyle automatically generates an interactive
Plotly visualization of the final stack after a successful run.

The plotting function can also be called directly:

.. code-block:: python

   from spectraPyle.plot.plot import plotting

   plotting("result_STACKING.fits", width=1200, height=700)

**2D spectral distribution** — The intermediate HDF5 file stores the full
matrix of individual resampled spectra before combination. It can be explored
with ``plot_h5_heatmap()``:

.. code-block:: python

   from spectraPyle.plot.plot import plot_h5_heatmap

   fig = plot_h5_heatmap(
       h5_path="result_array.h5",
       fits_path="result_STACKING.fits",
       metric="specMedian",
       mode="heatmap",
   )

   fig.show()

When working from the source repository, the supplied
``notebooks/plot_helper.ipynb`` provides an interactive interface to the same
visualization tools.

**Spectral Line Manager** — When working from the source repository,
``notebooks/line_manager.ipynb`` can be used to enable or disable emission
lines and absorption features displayed by the plotting routines.

Configuration Pipeline
----------------------

All inputs pass through a strict validation pipeline::

   Widgets / JSON / YAML / CLI
       ↓
   normalize_raw_config()        runtime/runtime_adapter.py
       ↓
   StackingConfig (Pydantic v2)  schema/schema.py
       ↓
   StackingConfigResolver        schema/schema.py
       ↓
   flatten_schema_model()        runtime/runtime_adapter.py
       ↓
   Stacking(flat_dict).run()     stacking/stacking.py

Optional continuum normalization is configured inside the ``norm`` block and
works with every instrument driver, including ``generic``. See
:doc:`normalization` for the polynomial fitting, feature masking, and EW
caveats.
