Tools Install for Windows 10
============================

.. caution::

	We assume core TurboWAVE has already been installed.

.. caution::

	Anaconda warns that spaces in the install path should be avoided. If your account name has spaces, it may be worth creating a new account without spaces for this work.

Python 3 via Anaconda
---------------------

#. If you already a sufficiently recent Anaconda3 installed, you may be able to skip the next step.  The danger is that older versions have trouble with PowerShell.
#. Run Miniconda3 installer from internet, accept defaults.
#. Open a new PowerShell window
#. If you get an error you likely need to update the Execution Policy

	* :samp:`Set-ExecutionPolicy Bypass -Scope CurrentUser`
	* Respond affirmatively to the prompt, close and reopen the PowerShell

#. :samp:`conda update conda`
#. Choose a name for your environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME} scipy matplotlib pillow jupyter ipympl`
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

#. If this is a new PowerShell session, activate the conda environment with :samp:`conda activate {NAME}`
#. :samp:`cd {twroot}/tools/twutils`
#. Do **not** descend into the second :samp:`twutils` directory within.
#. :samp:`pip install --upgrade pip`
#. :samp:`pip install .`
#. Your python programs should now have access to twutils and sub-packages.

Native DataViewer
-----------------

#. For Windows there is a native DataViewer application
#. You should be able to immediately run the :samp:`DataViewer.exe` application in :samp:`{twroot}/tools`.
#. This was written for Windows XP and we have lost the source, but it mostly still works.

Python DataViewer
-----------------

#. The Python DataViewer may also be useful since you can modify the source
#. Copy :samp:`{twroot}/tools/DataViewer.ipynb` to :samp:`~/bin`.
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{twroot}/tools/config-files/custom.css` to the new directory.

Input File Syntax Highlights
----------------------------

You can add syntax highlights for ``vim`` and ``Atom`` editors.  Syntax highlights assign different colors to different input file elements, such as comments, macros, keywords, etc..  This often makes the file easier to read and helps identify errors.

#. To enable turboWAVE input file syntax highlights with the :samp:`vim` editor

	* Create :samp:`~/vimfiles/` and :samp:`~/vimfiles/syntax/` if they do not already exist
	* Copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/vimfiles/`
	* Copy :samp:`{twroot}/tools/config-files/turbowave.vim` to :samp:`~/vimfiles/syntax/`
	* Files with extension ``.tw`` or the name ``stdin`` will be highlighted

#. To enable turboWAVE input file syntax highlights with the :samp:`Atom` editor, go to the package installation screen and search for the :samp:`language-turbowave` package.  Press the button to install the package.
