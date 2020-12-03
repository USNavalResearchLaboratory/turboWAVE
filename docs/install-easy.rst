Easy Install
============

The following should work in any operating system, assuming compilers are installed in standard locations.  The shell must be aware of Anaconda.

#. Open a terminal
#. Choose a name for the turboWAVE environment, here denoted :samp:`{NAME}`
#. Type :samp:`conda create -n {NAME} -c dfxgordon twutils`
#. Type :samp:`conda activate {NAME}`

	* This activates the environment. Each time a new terminal session is started, the environment needs to be activated.

#. Run the installer

	* Type ``twinstall`` for installations on a local machine.
	* Type ``twinstall --terminal`` for installation on a remote terminal.

#. Use the installer to complete the sequence of steps in the ``Tasks`` area.

	* You can usually accept default responses.
	* For now stick with OpenMP for the accelerator.
	* The installer can configure for GPGPU, but you may need to fulfill some prerequisites as root for the compiler to succeed.
