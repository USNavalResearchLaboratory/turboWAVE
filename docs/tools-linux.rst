Tools Install for Linux
=======================

Install Anaconda
----------------

#. If you already have Anaconda3 installed, skip to the next step.
#. Download Miniconda3 installer from internet
#. Navigate to downloaded file
#. :samp:`bash {filename}`, where :samp:`{filename}` is the file that you just downloaded
#. Respond with defaults to prompts.  Open a new terminal window when finished.

TurboWAVE Python Environment
----------------------------

#. :samp:`conda update conda`
#. :samp:`conda init`
#. Open a new terminal window.
#. Choose a name for your environment, denoted :samp:`{NAME}`
#. :samp:`conda create -n {NAME} -c dfxgordon twutils`
#. :samp:`conda activate {NAME}`
#. You are now in an isolated conda environment.  The environment must be activated each time you open a new terminal window.

Python DataViewer
-----------------

#. The Python DataViewer has to be run in a Jupyter Notebook.
#. Copy :samp:`{twroot}/tools/DataViewer.ipynb` to :samp:`~/bin`
#. Create a directory :samp:`~/.jupyter/custom/` and copy :samp:`{twroot}/tools/config-files/custom.css` to the new directory.

Input File Syntax Highlights
----------------------------

You can add syntax highlights for ``vim`` and ``Atom`` editors.  Syntax highlights assign different colors to different input file elements, such as comments, macros, keywords, etc..  This often makes the file easier to read and helps identify errors.

#. To enable turboWAVE input file syntax highlights with the :samp:`vim` editor

	* Copy :samp:`{twroot}/tools/config-files/filetype.vim` to :samp:`~/.vim/`
	* Copy :samp:`{twroot}/tools/config-files/turbowave.vim` to :samp:`~/.vim/syntax/`
	* Files with extension ``.tw`` or the name ``stdin`` will be highlighted

#. To enable turboWAVE input file syntax highlights with the :samp:`Atom` editor, go to the package installation screen and search for the :samp:`language-turbowave` package.  Press the button to install the package.
