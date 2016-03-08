GAPPATH = ../..
#
# makefile for the genss package                             Max Neunhoeffer
#
##  Copyright (C) 2009  Max Neunhoeffer
##  This file is free software, see license information at the end.
#
default: doc

doc:
	(echo 'Read("makedoc.g");' | $(GAPPATH)/bin/gap.sh -A -q -b)

clean:
	(cd doc ; ./clean)

.PHONY: default doc clean

##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; version 3 of the License.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
##

