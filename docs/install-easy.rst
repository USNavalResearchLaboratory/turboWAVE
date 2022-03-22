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

	* Type ``twinstall`` for installations on a local machine, or a remote machine with tolerable X forwarding.
	* Type ``twinstall --terminal`` for installation on a remote machine where X forwarding fails or is too slow.

	.. note::

		As of this writing, the ``--terminal`` option is not compatible with Windows PowerShell.  You can use a third party terminal emulator (e.g. PuTTY) instead.  Correct rendering of lines can depend on terminal settings.

#. Use the installer to complete the sequence of steps in the ``Tasks`` area.

	* You can usually accept default responses.
	* On Cray, the platform and vector type popups should be set manually.
	* For now use OpenMP for the accelerator.
	* The installer can configure for GPGPU, but you may need to fulfill some prerequisites as root for the compiler to succeed.

.. note::

	ANSI escape sequences provide color and formatting.  This works out of the box with most shells.  As of this writing, one way to get the correct formatting in PowerShell (mostly) is by issuing ``Set-ItemProperty HKCU:\Console VirtualTerminalLevel -Type DWORD 1``.

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
#. If you installed syntax highlights and you want to remove them:

	* Delete ``~/.vim/syntax/turbowave.vim``, ``~/.turbowave.nanorc``, and ``~/.config/micro/syntax/turbowave.micro.yaml``, or their Windows counterparts
	* Delete ``~/.vim/filetype.vim`` and ``~/.nanorc``, or their Windows counterparts, unless you need them for other purposes.  In the latter case, search for lines containing ``turbowave`` and delete them.
	* If you installed ``language-turbowave`` in Atom remove it using Atom's package manager
