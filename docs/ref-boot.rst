
.. toctree::
	:maxdepth: 2

.. highlight:: none

Understanding Multiboot Workstations
------------------------------------

Overview of File system
,,,,,,,,,,,,,,,,,,,,,,,

Let us refer to a physical storage device generically as a drive.
A drive can be split into *partitions*.  There are various partitioning schemes
that can be utilized by Linux.
The partitions behave as if they were separate drives for most purposes.
From within a Linux file system, drives and partitions are somewhat hidden.
Instead, the user perceives *mount points*, which appear simply as elements of the
directory structure.

.. tip::
	If you are accustomed to DOS or Windows drive designations such as ``C:\``, ``D:\``, etc.,
	you may be surprised that file system paths in Linux have no analog of this.  Linux views files as being in an
	abstract directory space that is detached from physical devices.

Although drives and partitions are not explicitly referenced in most Linux file operations,
there is a conventional notation that designates them in cases where this is needed, such as
in creating mount points.
For a typical drive, the notation is ``sd`` (presumably stands for SATA drive) followed by an
alphabetic sequence.  For example, the first drive listed would be ``sda``, the
next ``sdb``, etc..  To identify a partition a digit is appended, in sequence, starting
with one.  For example, the first partition of the second drive is ``sdb1``. The drives
and partitions that are recognized by the OS will appear as "files" in the ``/dev`` directory of
the Linux file system.  These "files" are *not* a substitute for the mount points.

.. note::
	The naming scheme given above is typical, but it can change depending on the type of
	drive, the partitioning scheme, and the volume management layer used within the partition.
	Furthermore, lower level software such as a bootloader does not have to follow the same scheme.


UEFI / EFI Firmware
,,,,,,,,,,,,,,,,,,,,,,,,,,

When the computer starts there is a program built into the Motherboard that runs first.
The modern variety of program is called UEFI, or sometimes EFI, as opposed to the legacy BIOS.
We will focus on the modern variety, UEFI.  The UEFI program uses information from a disk partition called the EFI System Partition (ESP).
There can be multiple ESP's on a given disk.  Moreover, there can be multiple disks with any number of ESP's per disk.  The role of the ESP is to contain at least one :dfn:`bootloader`.
The bootloader is a piece of software designed to load one or more operating systems.
When the UEFI program first loads, it typically runs a default bootloader, unless interrupted by a keystroke.
If there is a keystroke interruption, UEFI may offer the user a menu to choose from different bootloaders, something like this::

	ubuntu (M.2_1 : TS512GMTS600)
	centos (M.2_1 : TS512GMTS600)
	windows (M.2_1 : TS512GMTS600)

The menu may also list physical drives.  Here we have a scenario where there is an ESP on an M.2 solid state drive,
containing bootloaders associated with 3 operating systems.  When one of the menu items is selected,
the associated bootloader is executed.

The connections and associations allowed by UEFI are very flexible, but this has a downside,
namely, that things may not associate in the way you expect.  In particular,

	1.	Selecting the bootloader associated with a given OS does not necessarily boot that OS,
		because the bootloader often has its own multiboot capability.
	2.	The OS does not necessarily mount the ESP that booted it.  This can cause major
		confusion if you are trying to make changes to the bootloader from within the OS.

If the workstation has a permanent drive configuration, one can avoid confusion by using a single ESP for all the bootloaders.
If drives are to be swapped in and out, one must either designate a permanent drive to contain the ESP, or put an ESP on each drive.

.. sidebar:: Legacy BIOS

	The legacy BIOS firmware is usually supported by UEFI.  Legacy BIOS firmware
	searches for a Master Boot Record (MBR).  The MBR is like the ESP, but less flexible.

In order to discover which partitions are being mounted upon startup by a Linux OS, use:

	:samp:`cat /etc/fstab`

To identify the ESP, look for the mount point ``/boot/efi``.  In order to discover all the available partitions:

	:samp:`sudo fdisk -l`

The ESP partitions can be identified by looking at the ``Type`` column.  To decode a UUID, use:

 	:samp:`blkid -U <UUID>`

Bootloader Software
,,,,,,,,,,,,,,,,,,,

The primary bootloader used with Linux is GRUB2.  When UEFI runs the bootloader it
typically presents the user with a menu like this::

	Ubuntu 16.04
	Advanced Options for Ubuntu 16.04
	CentOS 7.4
	Advanced Options for CentOS 7.4
	Windows Boot Loader

If everything is working, the user simply selects one of the items and the corresponding OS
will start up.  This can go wrong, however, if the GRUB2 configuration file, usually called ``grub.cfg``,
contains errors, or is missing.

.. caution::
	The name of the configuration file can vary across different GRUB2 installations. Another common name
	for the configuration file is ``menu.lst``.

You are not supposed to directly edit :file:`grub.cfg`, but you still need to know where it is (in some cases).
On UEFI systems, the location is :file:`/boot/efi/EFI/{osname}/` where :file:`{osname}` depends on the OS.
You can determine :file:`{osname}` with

:samp:`sudo ls /boot/efi/EFI/`

In order to change :file:`grub.cfg` without direct editing, there is a utility that scans the
system and overwrites it, accounting for each OS it can find and understand.
On Ubuntu, this utility is invoked using:

	:samp:`sudo update-grub`

.. tip::
	Utilities like ``update-grub`` update the ``grub.cfg`` on the ESP that is mounted.
	If you are booting from a different ESP it can seem as if the update had no effect.
	Linux determines which partitions to mount by reading ``/etc/fstab``.

More generically, one has to resort to

	:samp:`sudo grub2-mkconfig -o /boot/efi/EFI/{osname}/grub.cfg`

.. tip::
	When invoking ``grub2-mkconfig`` directly, be sure to
	correctly specify the path of the ``grub.cfg`` that goes with the OS you are working in.

Modifying the GRUB2 menu
,,,,,,,,,,,,,,,,,,,,,,,,

There are situations where you may want to customize or repair the GRUB2 menu.
For example, you may need to pass customized arguments to the linux kernel.

The menu that is presented by GRUB2 is formed from entries in ``grub.cfg``.
The menu entries can be modified on the fly by highlighting a menu item with the arrow keys
and pressing ``e``.  This will display the section of ``grub.cfg`` which you can edit.  **The
changes made in this way are not permanent**.  After making the changes, you can press :samp:`Control-x` to
load and run the selected OS.

You can make permanent changes to the GRUB2 menu in two ways.  The simple way is to edit ``/etc/default/grub``.  If all you want to do is add or subtract kernel arguments, simply edit the string assigned to ``GRUB_CMDLINE_LINUX``, and run the grub update utility.

You can make more elaborate changes by editing the set of files in
``/etc/grub.d/``, followed by running the grub update utility.  If you edit ``10_linux``
you will be modifying how you load the OS you are currently
working in.  If you edit ``30_os-prober`` you will be modifying how other OS's are loaded.
The file ``40_custom`` is for manually adding your own specialized menu items.

.. note::
	Although you are potentially affecting how multiple OS's are loaded, you are doing this
	with respect to a single bootloader.  The other OS's may have their own bootloaders
	which are not affected until you overwrite their associated ``grub.cfg`` files.

In modifying the configuration files,
the primary piece of knowledge you need is that the GRUB2 command to load the Linux kernel is ``linux``.
Look for this in the GRUB2 configuration files.  The argument immediately following ``linux`` is the
name of the kernel file to load (Linux typically keeps past versions of its kernel on the disk).
The arguments after the kernel file are the arguments to pass to the kernel itself.

.. tip::
	In addition to the ``linux`` command, there is ``linuxefi`` for loading an OS in secure boot mode.
	If ``linux`` fails try turning off secure-boot in UEFI setup, or using ``linuxefi`` to
	load the kernel in the GRUB2 configuration.

Synchronizing Clocks
,,,,,,,,,,,,,,,,,,,,

You may find that the time of day is off by several hours when using a multiboot configuration.  This happens because different operating systems make different assumptions about the meaning of the time given by the Real Time Clock (RTC) built into the motherboard.  Linux usually assumes the RTC is giving universal time, while Windows assumes that it is giving the local time.

To make all operating systems treat the RTC as a universal time clock, proceed as follows.

	#. For each Linux operating system:

		* :samp:`sudo timedatectl -set-local-rtc 0`

	#. For Windows, use ``regedit`` to change the RTC to universal time.

		* Descend through the folders to find the key :samp:`HKEY_LOCAL_MACHINE\\SYSTEM\\CurrentControlSet\\Control\\TimeZoneInformation`
		* Select the ``TimeZoneInformation`` key and use the ``Edit`` menu to add a ``DWORD`` type with name ``RealTimeIsUniversal`` and value ``1``.

Summary
,,,,,,,

You may want to have a multiboot workstation for
running turboWAVE.  It is useful to understand how this works, especially when
problems arise.  We assume the modern UEFI system is used in place of the
legacy BIOS and MBR.  The basic boot sequence is:

	1.	UEFI firmware selects a bootloader, which it looks for in any of the EFI
		System Partitions (ESP) that are installed in the system.

	2.	The bootloader selects an operating system (OS), and starts it, possibly
		passing arguments to the kernel and running an initialization program.

	3.	The OS handles login and all subsequent computing.
