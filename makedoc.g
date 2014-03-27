##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##

SetPackagePath("genss", ".");
PrintTo("VERSION", PackageInfo("genss")[1].Version);

LoadPackage("GAPDoc");

MakeGAPDocDoc("doc", "genss", [], "genss");
CopyHTMLStyleFiles("doc");
GAPDocManualLab("genss");

QUIT;

