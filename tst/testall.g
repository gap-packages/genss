#
# This file runs package tests. It is also referenced in the package
# metadata in PackageInfo.g.
#
LoadPackage("genss");
d := DirectoriesPackageLibrary("genss", "tst");

TestDirectory(d[1], rec(exitGAP := true));

FORCE_QUIT_GAP(1);
