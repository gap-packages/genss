# A little try to implement a generic Schreier-Sims using the orb package

DeclareInfoClass( "InfoGenSS" );
SetInfoLevel(InfoGenSS,1);

DeclareGlobalVariable( "GENSS" );
InstallValue(GENSS,rec());

BindGlobal( "StabilizerChainFamily", NewFamily( "StabilizerChainFamily" ) );
DeclareCategory( "IsStabilizerChain", IsComponentObjectRep );
DeclareRepresentation( "IsStabilizerChainByOrb", IsStabilizerChain,
  [  ] );
BindGlobal( "StabChainByOrbType", 
  NewType(StabilizerChainFamily,IsStabilizerChain and IsStabilizerChainByOrb));

DeclareGlobalFunction( "GENSS_MapBaseImage" );
DeclareGlobalFunction( "GENSS_FindVectorsWithShortOrbit" );
DeclareGlobalFunction( "GENSS_FindShortOrbit" );

DeclareGlobalFunction( "GENSS_NextBasePoint" );
DeclareOperation( "FindBasePointCandidates", [ IsGroup, IsInt ] );
DeclareOperation( "SiftGroupElement", [ IsStabilizerChain, IsObject ] );
DeclareAttribute( "Size", IsStabilizerChain );

DeclareGlobalFunction( "GENSS_CreateStabChainRecord" );
DeclareGlobalFunction( "GENSS_StabilizerChainInner" );
DeclareOperation( "AddGeneratorToStabilizerChain",
  [IsStabilizerChain,IsObject] );
DeclareGlobalFunction( "GENSS_CopyDefaultOptions" );
DeclareOperation( "StabilizerChain", [ IsGroup ] );
DeclareOperation( "StabilizerChain", [ IsGroup, IsRecord ] );

