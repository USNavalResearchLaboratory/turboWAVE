Core Install for Windows 10
===========================

We will use the free version of Visual Studio from Microsoft to compile turboWAVE.
Most turboWAVE text files, such as input file examples, have UNIX line feeds.  This is no problem for WordPad (set word wrap to no wrap), but Notepad will not display them properly.  Installing a developer-oriented text editor, such as Atom, might be useful.

#. Put the turboWAVE components (:samp:`core` and :samp:`tools`) into some directory, denoted :samp:`{twroot}`.
#. Edit :samp:`{twroot}\\core\\source\\definitions.h`
#. In the definitions file, you must comment/uncomment lines to select platform and acceleration options.  In a C++ file, comments are preceded by :samp:`//`, and :samp:`#` is **not** a comment.  For this installation, only :samp:`#define USE_DESKTOP` and :samp:`#define USE_OPENMP` should be uncommented.
#. Install Visual Studio Community Edition (find installer via internet search)

	* Install the :samp:`Desktop Development with C++` and :samp:`Python Development` components.

#. When installation of Visual Studio is complete, open the :samp:`Developer Command Prompt for VS {YYYY}`.  You can find this in the start menu under :samp:`Visual Studio {YYYY}`.

	* The usual Windows command prompt may not be used.

#. :samp:`cd` :samp:`{twroot}`:samp:`\\core\\source`
#. :samp:`nmake /F winmakefile`
#. Create a directory called :samp:`{Run}` wherever you prefer in the file system.
#. Move :samp:`tw3d.exe` to :samp:`{Run}`.
