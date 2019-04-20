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
#. :samp:`conda init powershell`
#. Put away the terminal and open a user PowerShell
#. If you get an error you likely need to update the Execution Policy

	* :samp:`Set-ExecutionPolicy Bypass -Scope CurrentUser`
	* Respond affirmatively to the prompt, close and reopen the PowerShell

#. Choose a name for your environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME} scipy matplotlib jupyter`
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
