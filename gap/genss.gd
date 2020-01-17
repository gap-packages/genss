#############################################################################
##
##  genss.gd              genss package           
##                                                           Max Neunhoeffer
##                                                              Felix Noeske
##
##  Copyright 2006 Lehrstuhl D für Mathematik, RWTH Aachen
##
##  Declaration stuff for generic Schreier-Sims
##
#############################################################################

#############################################################################
# Our Info class:
#############################################################################

DeclareInfoClass( "InfoGenSS" );
SetInfoLevel(InfoGenSS,1);

#############################################################################
# The following global record contains default values for options for the 
# main function "StabilizerChain":
#############################################################################

BindGlobal( "GENSS", rec() );


#############################################################################
# Our main type:
#############################################################################

BindGlobal( "StabilizerChainFamily", NewFamily( "StabilizerChainFamily" ) );
DeclareCategory( "IsStabilizerChain", IsComponentObjectRep );
DeclareRepresentation( "IsStabilizerChainByOrb", IsStabilizerChain,
  [ "size",        # false or a size if known
    "orb",         # the orbit
    "cand",        # a record for base point candidates (shared between layers)
    "stab",        # the next layer
    "base",        # the global list of base points (shared between layers)
    "layer",       # the layer, 1 is top
    "opt",         # the options record (shared between layers)
  ] );
BindGlobal( "StabChainByOrbType", 
  NewType(StabilizerChainFamily,IsStabilizerChain and IsStabilizerChainByOrb));


#############################################################################
# A few helper functions needed elsewhere:
#############################################################################

DeclareGlobalFunction( "GENSS_CopyDefaultOptions" );
DeclareGlobalFunction( "GENSS_MapBaseImage" );
DeclareGlobalFunction( "GENSS_FindVectorsWithShortOrbit" );
DeclareGlobalFunction( "GENSS_FindShortOrbit" );
DeclareGlobalFunction( "GENSS_IsOneProjective" );
DeclareGlobalFunction( "GENSS_RandomElementFromAbove" );
DeclareGlobalFunction( "GENSS_OpFunctionMaker" );

#############################################################################
# Now to the heart of the method, the Schreier-Sims:
#############################################################################

DeclareOperation( "FindBasePointCandidates", 
  [ IsGroup, IsRecord, IsInt, IsObject ] );
DeclareGlobalFunction( "GENSS_NextBasePoint" );
DeclareGlobalFunction( "GENSS_CreateStabChainRecord" );
DeclareGlobalFunction( "GENSS_StabilizerChainInner" );
DeclareGlobalFunction( "GENSS_DeriveCandidatesFromStabChain" );
DeclareGlobalFunction( "GENSS_TrivialOp" );
DeclareOperation( "StabilizerChain", [ IsGroup ] );
DeclareOperation( "StabilizerChain", [ IsGroup, IsRecord ] );
DeclareOperation( "AddGeneratorToStabilizerChain",
                  [IsStabilizerChain,IsObject] );
DeclareOperation( "SiftGroupElement", [ IsStabilizerChain, IsObject ] );
DeclareOperation( "SiftGroupElementSLP", [ IsStabilizerChain, IsObject ] );
DeclareOperation( "StrongGenerators", [IsStabilizerChain] );
DeclareOperation( "NrStrongGenerators", [IsStabilizerChain] );
DeclareOperation( "GENSS_CreateSchreierGenerator",
  [ IsStabilizerChain, IsPosInt, IsPosInt ] );
DeclareGlobalFunction( "GENSS_FindGensStabilizer" );
DeclareGlobalFunction( "GENSS_FindShortGensStabilizerOld" );
DeclareGlobalFunction( "GENSS_FindShortGensStabilizer" );
DeclareGlobalFunction( "SLPChainStabilizerChain" );
DeclareGlobalFunction( "GENSS_Prod" );
DeclareOperation( "VerifyStabilizerChainMC", [ IsStabilizerChain, IsInt ] );
DeclareGlobalFunction( "VerifyStabilizerChainTC" );
DeclareGlobalFunction( "VerifyStabilizerChainTC2" );
DeclareGlobalFunction( "GENSS_ImageElm" );
DeclareGlobalFunction( "GENSS_PreImagesRepresentative" );
DeclareGlobalFunction( "GroupHomomorphismByImagesNCStabilizerChain" );
DeclareOperation( "AddNormalizingGenToLayer", 
                  [ IsStabilizerChain, IsObject, IsPosInt ] );

#############################################################################
# The following operations apply to stabilizer chains:
#############################################################################

DeclareAttribute( "Size", IsStabilizerChain );
DeclareOperation( "IsProved", [IsStabilizerChain] );
DeclareOperation( "StabChainOp", [IsPermGroup, IsStabilizerChain] );
DeclareOperation( "GroupIteratorByStabilizerChain", [IsStabilizerChain] );
DeclareGlobalFunction( "GENSS_GroupNextIterator" );
DeclareGlobalFunction( "GENSS_GroupIsDoneIterator" );
DeclareGlobalFunction( "GENSS_GroupShallowCopy" );
DeclareGlobalFunction( "GENSS_MakeIterRecord" );

DeclareOperation( "BaseStabilizerChain", [IsStabilizerChain] );
DeclareOperation( "SiftBaseImage", [IsStabilizerChain, IsList] );

DeclareOperation( "SetStabilizerChain", [IsGroup,IsStabilizerChain] );
DeclareAttribute( "StoredStabilizerChain", IsGroup );

#############################################################################
# The following operations are about methods to compute stabilisers:
#############################################################################

DeclareOperation( "Stab", [IsGroup, IsObject, IsFunction, IsRecord] );
DeclareOperation( "Stab", [IsList, IsObject, IsFunction, IsRecord] );
DeclareOperation( "Stab", [IsGroup, IsObject, IsFunction] );
DeclareOperation( "Stab", [IsList, IsObject, IsFunction] );
DeclareOperation( "StabMC", [IsGroup, IsObject, IsFunction, IsRecord] );
DeclareOperation( "StabMC", [IsGroup, IsObject, IsFunction] );
DeclareOperation( "StabMC", [IsList, IsObject, IsFunction, IsRecord] );
DeclareOperation( "StabMC", [IsList, IsObject, IsFunction] );

#######################################################
# The following operations are for backtrack searches: 
#######################################################

DeclareGlobalFunction( "BacktrackSearchStabilizerChainElement" );
DeclareGlobalFunction( "ComputeSuborbitsForStabilizerChain" );
DeclareGlobalFunction( "BacktrackSearchStabilizerChainSubgroup" );

