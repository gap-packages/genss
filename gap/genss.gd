# A little try to implement a generic Schreier-Sims using the orb package

DeclareInfoClass( "InfoGenSS" );

DeclareGlobalVariable( "GENSS" );
InstallValue(GENSS,rec());

BindGlobal( "StabilizerChainFamily", NewFamily( "StabilizerChainFamily" ) );
DeclareCategory( "IsStabilizerChain", IsComponentObjectRep );
DeclareRepresentation( "IsStabilizerChainByOrb", IsStabilizerChain,
  [ "gens", "gensi", "size", "orb", "cand", "stab" ] );
BindGlobal( "StabChainByOrbType", 
  NewType(StabilizerChainFamily,IsStabilizerChain and IsStabilizerChainByOrb));

DeclareGlobalFunction( "GENSS_CreateStabChainRecord" );
DeclareGlobalFunction( "GENSS_CollectSchreierGens" );
DeclareGlobalFunction( "GENSS_ConsiderNewSchreierGenerator" );
DeclareGlobalFunction( "GENSS_StartNewOrbitEnumeration" );
DeclareGlobalFunction( "GENSS_MapBaseImage" );
DeclareGlobalFunction( "GENSS_FindVectorsWithShortOrbit" );
DeclareGlobalFunction( "GENSS_FindShortOrbit" );
DeclareGlobalFunction( "GENSS_CompleteStabilizerChain" );

DeclareGlobalFunction( "GENSS_NextBasePoint" );
DeclareOperation( "StabilizerChain", [ IsGroup ] );
DeclareOperation( "StabilizerChain", [ IsGroup, IsRecord ] );
DeclareOperation( "FindBasePointCandidates", [ IsGroup, IsInt ] );
DeclareOperation( "AddGeneratorToStabilizerChain",
  [IsStabilizerChain,IsObject,IsBool] );
DeclareOperation( "SiftGroupElement", [ IsStabilizerChain, IsObject ] );
DeclareOperation( "SiftGroupElement2", [ IsStabilizerChain, IsObject ] );
DeclareAttribute( "Size", IsStabilizerChain );

DeclareGlobalFunction( "GENSS_CreateStabChainRecord2" );
DeclareGlobalFunction( "GENSS_StabilizerChain2Inner" );
DeclareOperation( "AddGeneratorToStabilizerChain2",
  [IsStabilizerChain,IsObject] );
DeclareGlobalFunction( "GENSS_CopyDefaultOptions" );
DeclareOperation( "StabilizerChain2", [ IsGroup ] );
DeclareOperation( "StabilizerChain2", [ IsGroup, IsRecord ] );

