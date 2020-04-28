## mt2mseed
Magnetotellurics: IRIS/OSU developed software to convert from nimsread *.bin files to miniseeds

# Overview

Convert MT time series data to Mini-SEED. Currently, only NIMSread binary output files are supported.

Software developed for EarthScope MT data archiving in 2006-2008.

Authors: Chad Trabant, Anna Kelbert.

# Building/Installing 

In most environments a simple 'make' will build the program.

In the Win32 environment the Makefile.wat can be used with Open
Watcom's wmake program.

Using GCC, running 'make static' will compile a static version
if possible.

For further installation simply copy the resulting binary and man page
(in the 'doc' directory) to appropriate system directories.

# Licensing 

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License (GNU-GPL) for more details.  The GNU-GPL and
further information can be found here: http://www.gnu.org/
