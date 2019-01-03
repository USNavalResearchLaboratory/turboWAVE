Tools Install for RHEL 7.5
==========================

.. caution::

	We assume core TurboWAVE has already been installed.

.. note::

	CentOS and Scienfific Linux are clones of RHEL, and so should behave the same way.

Python 3 Conda Environment
--------------------------

#. In this instruction we will use conda virtual environments.  You can also use the native venv system.
#. Download Miniconda3 installer from internet
#. Navigate to downloaded file
#. :samp:`bash {downloaded_file}`
#. Respond with defaults to prompts.  Open a new terminal window when finished.
#. :samp:`conda update conda`
#. Choose a name for the environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME}`
#. :samp:`conda activate {NAME}`

	* This is the new way to activate the environment.  If prompted to modify your login files follow the instructions and repeat.

#. You are now in an isolated conda environment.  The environment must be activated each time you open a new terminal window.
#. :samp:`conda install scipy matplotlib jupyter`
#. Jupyter may default to the Konqueror browser.  It is recommended to switch the default to Firefox.

	* :samp:`jupyter notebook --generate-config`
	* Edit :samp:`~/.jupyter/jupyter_notebook_config.py`, uncomment the :samp:`c.NotebookApp.browser` definition, and replace the empty string with :samp:`'firefox'`.

#. If there are problems with Jupyter notebooks any or all of the following may be tried:

	* Try adding ``-c conda-forge`` at any install step
	* :samp:`conda install widgetsnbextension={n}`, where :samp:`{n}` is some preferred version.
	* :samp:`conda install ipywidgets`
	* :samp:`jupyter nbextension install --py --sys-prefix widgetsnbextension`
	* :samp:`jupyter nbextension enable --py --sys-prefix widgetsnbextension`

TurboWAVE Python Packages
-------------------------

#. If this is a new terminal session, activate the virtual environment (see above)
#. Navigate to the :samp:`{twroot}/tools/twutils` directory
#. Do **not** descend into the second :samp:`twutils` directory within.
#. :samp:`pip install --upgrade pip`
#. :samp:`pip install .`
#. Your python programs should now have access to twutils and sub-packages.

Python DataViewer
-----------------

#. The Python DataViewer has to be run in a Jupyter Notebook.
#. Copy :samp:`{twroot}/tools/DataViewer.ipynb` to :samp:`~/bin`
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{twroot}/tools/config-files/custom.css` to the new directory.
