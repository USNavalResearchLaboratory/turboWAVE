Easy Install
============

Install
-------

The following should work in any operating system, assuming compilers are installed in standard locations.  The shell must be aware of Anaconda.

#. Open a terminal
#. Choose a name for the turboWAVE environment, here denoted :samp:`{NAME}`
#. Type :samp:`conda create -n {NAME} -c dfxgordon twutils`
#. Type :samp:`conda activate {NAME}`

	* This activates the environment. Each time a new terminal session is started, the environment needs to be activated.

#. Run the installer

	* Type ``twinstall`` for installations on a local machine.
	* Type ``twinstall --terminal`` for installation on a remote machine.

	.. note::

		As of this writing, the installer's ``--terminal`` option is not compatible with Windows PowerShell.  As a result, if you want to install on a remote non-Windows server from a Windows terminal, use the old-style command prompt, or a third party terminal emulator.

#. Use the installer to complete the sequence of steps in the ``Tasks`` area.

	* You can usually accept default responses.
	* On Cray, the platform and vector type popups should be set explicitly.
	* For now use OpenMP for the accelerator.
	* The installer can configure for GPGPU, but you may need to fulfill some prerequisites as root for the compiler to succeed.

Upgrade
-------

#. Open a terminal
#. :samp:`conda activate {NAME}`
#. :samp:`conda update -c dfxgordon twutils`
#. Run the installer

	* If you delete the old local repository first, or select a new location for ``Get Components``, the process is identical to a new installation.
	* Alternatively you can point ``Get Components`` to the old local repository and let the installer pull the latest from upstream.  Note the local repository has to be clean in this case.

Uninstall
---------

#. Open a terminal
#. To remove the whole turboWAVE environment type :samp:`conda remove -n {NAME} --all`
#. Delete the ``turboWAVE`` local repository (the location was chosen by you)
#. Delete ``DataViewer.ipynb`` (the location was chosen by you)
#. Remove ``~/.vim/filetype.vim`` and ``~/.vim/syntax/turbowave.vim`` (if you installed them)
#. If you installed ``language-turbowave`` in Atom remove it using Atom's package manager
