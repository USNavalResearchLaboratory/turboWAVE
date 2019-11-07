Tools Install for MacOS 10.15
=============================

.. caution::

	We assume core TurboWAVE has already been installed.

X Windows Display Manager
-------------------------

The graphical output from some python packages is displayed using X Windows.  On MacOS this is typically emulated using XQuartz.  Find XQuartz via internet search and install.

Python 3 via Anaconda
---------------------

#. If you already have Anaconda3 installed, skip the next 4 steps.  If you have enough packages in your conda environment, it is possible you can skip this entire section: but to be safe create a new environment as detailed below.
#. Download Miniconda3 installer from internet
#. Navigate to downloaded file
#. :samp:`bash {filename}`, where :samp:`{filename}` is the file that you just downloaded
#. Respond with defaults to prompts.  Open a new terminal window when finished.
#. :samp:`conda update conda`
#. :samp:`conda init`
#. Choose a name for your environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME} scipy matplotlib pillow jupyter`
#. :samp:`conda activate {NAME}`
#. You are now in an isolated conda environment.  The environment must be activated each time you open a new terminal window.
#. If there are problems with Jupyter notebooks any or all of the following may be tried:

	* Try adding ``-c conda-forge`` at any install step
	* :samp:`conda install widgetsnbextension={n}`, where :samp:`{n}` is some preferred version.
	* :samp:`conda install ipywidgets`
	* :samp:`jupyter nbextension install --py --sys-prefix widgetsnbextension`
	* :samp:`jupyter nbextension enable --py --sys-prefix widgetsnbextension`

TurboWAVE Python Packages
-------------------------

#. If this is a new terminal session, activate the conda environment as above.
#. Navigate to the :samp:`{twroot}/tools/twutils` directory
#. Do **not** descend into the second :samp:`twutils` directory within.
#. :samp:`pip install --upgrade pip`
#. :samp:`pip install .`
#. Your python programs should now have access to twutils and sub-packages.

Native DataViewer
-----------------

#. For MacOS there is a native DataViewer application
#. Double-click on :samp:`{twroot}/tools/DataViewer.dmg`
#. Open the disk image and copy the DataViewer application to :samp:`Applications` or wherever you like.

Python DataViewer
-----------------

#. The Python DataViewer may also be useful since you can modify the source
#. Copy :samp:`{twroot}/tools/DataViewer.ipynb` to :samp:`~/bin`
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{twroot}/tools/config-files/custom.css` to the new directory.

Input File Syntax Highlights
----------------------------

You can add syntax highlights for ``vim`` and ``Atom`` editors.  Syntax highlights assign different colors to different input file elements, such as comments, macros, keywords, etc..  This often makes the file easier to read and helps identify errors.

#. To enable turboWAVE input file syntax highlights with the :samp:`vim` editor

	* Copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim/`
	* Copy :samp:`{twroot}/tools/config-files/turbowave.vim` to :samp:`~/.vim/syntax/`
	* Files with extension ``.tw`` or the name ``stdin`` will be highlighted

#. To enable turboWAVE input file syntax highlights with the :samp:`Atom` editor

	* Create directory :samp:`~/.atom/packages/language-turbowave/grammars/` and copy :samp:`{twroot}/tools/config-files/turbowave.cson` to the new directory
	* Create directory :samp:`~/.atom/packages/language-turbowave/` and copy :samp:`{twroot}/tools/config-files/package.json` to the new directory
	* Restart Atom
	* Files with the extension ``.tw`` will be highlighted.  Using the ``Select Grammar`` menu item and choosing ``turbowave`` allows any file to be highlighted.
