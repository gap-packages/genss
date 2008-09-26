#############################################################################
##
##  read.g                genss package
##                                                          Max Neunhoeffer
##                                                             Felix Noeske
##
##  Copyright 2006 Lehrstuhl D für Mathematik, RWTH Aachen
##
##  Reading the implementation part of the genss package.
##
#############################################################################

ReadPackage("genss","gap/genss.gi");
ReadPackage("genss","gap/setwise.gi");

if IsBound(IO_PackageIsLoaded) then
    ReadPackage("genss","gap/picklers.gi");
else
    if not(IsBound(IO_PkgThingsToRead)) then
        IO_PkgThingsToRead := [];
    fi;
    Add(IO_PkgThingsToRead,["genss","gap/picklers.gi"]);
fi;
