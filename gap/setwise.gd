#############################################################################
##
##  setwise.gd              genss package           
##                                                           Max Neunhoeffer
##                                                              Felix Noeske
##
##  Copyright 2008 by the authors
##
##  Declaration stuff for setwise stabilizer code.
##
#############################################################################

DeclareGlobalFunction( "GENSS_FindElmMappingBaseSubsetIntoSet" );
DeclareGlobalFunction( "GENSS_ComputeFoundElm" );
DeclareOperation( "SetwiseStabilizer", [IsGroup,IsFunction,IsList] );
DeclareGlobalFunction( 
     "GENSS_FindElmMappingBaseSubsetIntoSetPartitionBacktrack" );
DeclareOperation( "SetwiseStabilizerPartitionBacktrack",
   [IsGroup,IsFunction,IsList] );

