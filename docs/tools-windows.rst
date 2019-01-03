Tools Install for Windows 10
============================

.. caution::

	We assume core TurboWAVE has already been installed.

.. caution::

	Anaconda warns that spaces in the install path should be avoided. If your account name has spaces, it may be worth creating a new account without spaces for this work.

Python 3 via Anaconda
---------------------

#. If you already have Anaconda3 installed, skip the next step.  If you have enough packages in your conda environment, it is possible you can skip this entire section: but to be safe create a new environment as detailed below.
#. Run Miniconda3 installer from internet
#. Open the Anaconda Prompt from the Start menu
#. :samp:`conda update conda`
#. Choose a name for your environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME}`
#. :samp:`conda activate {NAME}`
#. You are now in an isolated conda environment.  The environment must be activated each time you open a new terminal window.
#. :samp:`conda install scipy matplotlib jupyter`
#. If there are problems with Jupyter notebooks any or all of the following may be tried:

	* Try adding ``-c conda-forge`` at any install step
	* :samp:`conda install widgetsnbextension={n}`, where :samp:`{n}` is some preferred version.
	* :samp:`conda install ipywidgets`
	* :samp:`jupyter nbextension install --py --sys-prefix widgetsnbextension`
	* :samp:`jupyter nbextension enable --py --sys-prefix widgetsnbextension`


TurboWAVE Python Packages
-------------------------

#. If this is a new terminal session, activate the conda environment with :samp:`activate {NAME}`
#. Remember to use the Anaconda prompt, not the usual command prompt
#. Navigate to the :samp:`{twroot}\\tools\\twutils` directory
#. Do **not** descend into the second :samp:`twutils` directory within.
#. :samp:`pip install --upgrade pip`
#. :samp:`pip install .`
#. Your python programs should now have access to twutils and sub-packages.

Native DataViewer
-----------------

#. For Windows there is a native DataViewer application
#. You should be able to immediately run the :samp:`DataViewer.exe` application in :samp:`{twroot}\\tools`.
#. This was written for Windows XP and we have lost the source, but it mostly still works.

Python DataViewer
-----------------

#. The Python DataViewer may also be useful since you can modify the source
#. Copy :samp:`{twroot}\\tools\\DataViewer.ipynb` to some convenient place, such as :samp:`{Run}`.
#. Create a directory :samp:`{C}:\\Users\\{account_name}\\.jupyter\\custom\\` and copy :samp:`{twroot}\\core\\documentation\\config-files\\custom.css` to the new directory.
