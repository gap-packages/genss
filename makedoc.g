##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  
##  Call this with GAP.
##

LoadPackage("GAPDoc");

MakeGAPDocDoc("doc", "genss", [], "genss", "../../..");

CopyHTMLStyleFiles("doc");

GAPDocManualLab("genss");

QUIT;

