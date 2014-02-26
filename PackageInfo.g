#############################################################################
##  
##  PackageInfo.g for the package `genss'
##                                                            Max Neunhoeffer
##                                                               Felix Noeske
##  (created from Frank Lübeck's PackageInfo.g template file)
##  

SetPackageInfo( rec(

##  This is case sensitive, use your preferred spelling.
PackageName := "genss",

##  This may be used by a default banner or on a Web page, should fit on
##  one line.
Subtitle := "genss - generic Schreier-Sims",

##  See '?Extending: Version Numbers' in GAP help for an explanation
##  of valid version numbers. For an automatic package distribution update
##  you must provide a new version number even after small changes.
Version := "1.5.dev",
##  Please adjust also the VERSION file in the package directory when
##  changing this.

##  Release date of the current version in dd/mm/yyyy format.
Date := "31/05/2012",

##  Information about authors and maintainers.
Persons := [
  rec( 
    LastName      := "Neunhoeffer",
    FirstNames    := "Max",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "neunhoef@mcs.st-and.ac.uk",
    WWWHome       := "http://www-groups.mcs.st-and.ac.uk/~neunhoef",
    PostalAddress := Concatenation( [
                       "School of Mathematics and Statistics\n",
                       "Mathematical Institute\n",
                       "North Haugh\n",
                       "St Andrews, Fife KY16 9SS\n",
                       "Scotland, UK" ] ),
    Place         := "St Andrews",
    Institution   := "University of St Andrews"
  ),
  rec( 
    LastName      := "Noeske",
    FirstNames    := "Felix",
    IsAuthor      := true,
    IsMaintainer  := true,
    Email         := "felix.noeske@math.rwth-aachen.de",
    WWWHome       := "http://www.math.rwth-aachen.de/~Felix.Noeske",
    PostalAddress := Concatenation( [
                       "Felix Noeske\n",
                       "Lehrstuhl D fuer Mathematik, RWTH Aachen\n",
                       "Templergraben 64\n",
                       "52056 Aachen\n",
                       "Germany" ] ),
    Place         := "Aachen",
    Institution   := "RWTH Aachen"
  ),
],

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "deposited"     for packages for which the GAP developers agreed 
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages 
##    "other"         for all other packages
##
# Status := "accepted",
Status := "deposited",

##  You must provide the next two entries if and only if the status is 
##  "accepted" because is was successfully refereed:
# format: 'name (place)'
# CommunicatedBy := "Mike Atkinson (St. Andrews)",
#CommunicatedBy := "",
# format: mm/yyyy
# AcceptDate := "08/1999",
#AcceptDate := "",

BaseURL := "http://www-groups.mcs.st-and.ac.uk/~neunhoef/Computer/Software/Gap/",

PackageWWWHome := Concatenation( ~.BaseURL, "genss.html" ),
ArchiveURL     := Concatenation( ~.BaseURL, "genss/genss-", ~.Version ),
README_URL     := Concatenation( ~.BaseURL, "genss/README.genss" ),
PackageInfoURL := Concatenation( ~.BaseURL, "genss/PackageInfo.g" ),

ArchiveFormats := ".tar.gz",

##  Here you  must provide a short abstract explaining the package content 
##  in HTML format (used on the package overview Web page) and an URL 
##  for a Webpage with more detailed information about the package
##  (not more than a few lines, less is ok):
##  Please, use '<span class="pkgname">GAP</span>' and
##  '<span class="pkgname">MyPKG</span>' for specifing package names.
##  
AbstractHTML := 
  "The <span class=\"pkgname\">genss</span> package implements the \
   randomised Schreier-Sims algorithm to compute a stabiliser chain \
   and a base and strong generating set for arbitrary finite groups.",

PackageDoc := rec(
  BookName  := "genss",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "genss - generic Schreier-Sims",
  Autoload  := false
),

Dependencies := rec(
  GAP := ">=4.5",
  NeededOtherPackages := [
    ["GAPDoc", ">= 1.5"],
    ["IO", ">= 4.2"],
    ["orb", ">= 4.5"],
  ],
  SuggestedOtherPackages := [],
  ExternalConditions := []
),

AvailabilityTest := ReturnTrue,

Autoload := false,

##  *Optional*, but recommended: path relative to package root to a file which 
##  contains as many tests of the package functionality as sensible.
#TestFile := "tst/testall.g",

##  *Optional*: Here you can list some keyword related to the topic 
##  of the package.
Keywords := ["Schreier-Sims", "Schreier", "Sims", "Stabilizer chain"]

));


