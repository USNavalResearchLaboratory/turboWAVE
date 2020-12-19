Parameter Scans
===============

There is a Python script ``twscan`` which carries out parameter scans, taking advantage of parallelism.  Scans can be run in two modes:

.. glossary::

  parallel-serial scan
    Each simulation employs parallelism, but the set of simulations are run serially

  parallel-parallel scan
    Each simulation may employ parallelism, while the whole set also runs in parallel

Generally, parallel-serial scans are used on the desktop, while parallel-parallel scans are used on HPC systems or large local clusters.

Creating the Scan File
----------------------

The basic premise of the ``twscan`` tool is that a parameter scan can be carried out by letting some set of ``#define`` constants vary from one simulation to another.  The particular constants and how they vary are defined by a scan file (a serialized Python dictionary), which the user must create.  A simple example follows::

  import json
  import numpy as np

  # set up a 2D parameter space
  # we could have any number of parameters varying, but will use 2 in this example
  dims = (5,5) # parameter space is a 5x5 grid
  param1_key = '$a0' # some amplitude referenced in the input file
  param2_key = '$t0' # some pulse duration referenced in the input file
  param1_vals = np.linspace(0,1,dims[0])
  param2_vals = np.linspace(3,10,dims[1])

  # create the dictionary
  # N.b. numpy arrays are converted to lists so they can be serialized
  scan = {}
  scan[param1_key] = np.outer(param1_vals,np.ones(dims[1])).tolist()
  scan[param2_key] = np.outer(np.ones(dims[0]),param2_vals).tolist()

  # serialize the dictionary and write to disk
  # you may replace the filename with one of your own choice
  with open('tw_scan_file.json','w') as f:
  	json.dump(scan,f)

Running Parallel-Serial Scans
-----------------------------

Navigate to the directory where the scan file is located.  Copy the input file to the same directory.  The input file should have the parameter keys as ``#define`` constants.  These constants can be used as many times in the input file as desired.  Then type

:samp:`twscan {scan_file} -i {input_file} -n {mpi_tasks} -c {threads} --dry-run`

The ``--dry-run`` argument causes the tool to create the set of run directories without actually running the jobs.  Check inside a sampling of the directories to see if the input file ``#define`` constants have been set up as expected.  If everything appears correct, type

:samp:`twscan {scan_file} -i {input_file} -n {mpi_tasks} -c {threads}`

The jobs should run in sequence, and diagnostic outputs for a given job will appear in its directory.  Post-processing of the data set is up to the user.  However, most post-processing scripts will want to load the scan file::

  import json
  import numpy as np

  # Keys we expect to find in this example
  param1_key = '$a0' # some amplitude referenced in the input file
  param2_key = '$t0' # some pulse duration referenced in the input file

  # read the data into a dictionary
  # you may replace the filename with one of your own choice
  with open('tw_scan_file.json','w') as f:
  	scan = json.load(f)

  param1_vals = np.array(scan[param1_key])
  param2_vals = np.array(scan[param2_key])

Running Parallel-Parallel Scans
-------------------------------

The process for parallel-parallel scans is similar to that of parallel-serial, except that command line arguments are different, and the working directory also needs a job script.  As of this writing, the tool understands PBS and SLURM job scripts.  In order to test the scan, type

:samp:`twscan {scan_file} -i {input_file} -l {launch} --submit {sub} --script {script} --dry-run`

Here, :samp:`{launch}` is an MPI launcher such as ``aprun``, :samp:`{sub}` is a job submission command such as ``qsub``, and :samp:`{script}` is the name of the job script.  Note that the parallel parameters are not given as command line options, instead they are expected to be in the job script.  Inspect the job directories, and if everything looks correct, type

:samp:`twscan {scan_file} -i {input_file} -l {launch} --submit {sub} --script {script}`

The jobs will be put into the queue in order.
