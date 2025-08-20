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

#. Run the installer (for v5.x this may be dicey until further notice)

	* Type ``twinstall`` for installations on a local machine, or a remote machine with tolerable X forwarding.
	* Type ``twinstall --terminal`` for installation on a remote machine where X forwarding fails or is too slow (not available in PowerShell).

#. Use the installer to complete the sequence of steps in the ``Tasks`` area.

.. note:: Pretty PowerShell

	To get color working try ``Set-ItemProperty HKCU:\Console VirtualTerminalLevel -Type DWORD 1``.  To get Unicode characters working try ``$OutputEncoding = [console]::InputEncoding = [console]::OutputEncoding = New-Object System.Text.UTF8Encoding`` (put this in the ``$PROFILE`` file to make it persist). If unicode is still not displayed correctly, start cycling fonts (right-click title bar and choose ``Properties``).

Runtime
-------

When ``twinstall`` is used, the executable installs to the conda environment.  As a result, when outside the environment, the executable will not be found.  However, one should understand that the executable itself does not depend on Anaconda or Python in any way.  It can be manually installed anywhere.

Upgrade
-------

#. Open a terminal
#. :samp:`conda activate {NAME}`
#. :samp:`conda update -c dfxgordon twutils`
#. Run the installer (for v5.x this may be dicey until further notice)

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
