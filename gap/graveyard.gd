DeclareGlobalFunction( "GENSS_CreateStabChainRecord" );
DeclareGlobalFunction( "GENSS_CollectSchreierGens" );
DeclareGlobalFunction( "GENSS_ConsiderNewSchreierGenerator" );
DeclareGlobalFunction( "GENSS_StartNewOrbitEnumeration" );
DeclareGlobalFunction( "GENSS_CompleteStabilizerChain" );
DeclareOperation( "StabilizerChain", [ IsGroup ] );
DeclareOperation( "StabilizerChain", [ IsGroup, IsRecord ] );
DeclareOperation( "AddGeneratorToStabilizerChain",
  [IsStabilizerChain,IsObject,IsBool] );
DeclareOperation( "SiftGroupElement", [ IsStabilizerChain, IsObject ] );

