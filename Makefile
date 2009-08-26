GAPPATH = ../..
#
# makefile for the genss package                             Max Neunhoeffer
#
##  Copyright (C) 2009  Max Neunhoeffer
##  This file is free software, see license information at the end.
#
default: doc

doc: doc/manual.six

doc/manual.six: doc/genss.xml doc/intro.xml doc/install.xml doc/stabchains.xml \
		doc/backtrack.xml doc/examples.xml VERSION
	($(GAPPATH)/bin/gap.sh makedoc.g)

clean:
	(cd doc ; ./clean)

archive: doc
	(cd .. ; tar czvf genss-`cat genss/VERSION`.tar.gz --exclude ".svn" --exclude .git genss)

WEBPOS=/scratch/neunhoef/mywebpage.pub/Computer/Software/Gap/genss

towww: archive
	echo '<?xml version="1.0" encoding="ISO-8859-1"?>' >${WEBPOS}.version
	echo '<mixer>' >>${WEBPOS}.version
	cat VERSION >>${WEBPOS}.version
	echo '</mixer>' >>${WEBPOS}.version
	cp PackageInfo.g ${WEBPOS}
	cp README ${WEBPOS}/README.genss
	cp doc/manual.pdf ${WEBPOS}/genss.pdf
	cp ../genss-`cat VERSION`.tar.gz ${WEBPOS}
	@echo
	@echo Now please adjust the archive link in 
	@echo "    " ${WEBPOS}.mixer !

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

