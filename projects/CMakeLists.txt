cmake_minimum_required (VERSION 2.8)

# Standard LeMonADE
add_subdirectory(SimpleSimulatorStandard)

# Branched Catalysis
add_subdirectory(CreatorChainWalkingTertiaryBondWalking)

# Rouse Matrix evaluation of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_RouseMatrix)

# Topological properties of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_TopologicalProperties)

# Creational properties of BraCat structures
add_subdirectory(CreatorChainWalkingTertiaryBondWalking_CreationProperties)

# BGRW_CreationProperties
add_subdirectory(BGRW_CreationProperties)

# Rouse Matrix evaluation of BGRW
add_subdirectory(BGRW_RouseMatrix)

# Rouse Matrix evaluation of Hyperstar
add_subdirectory(Hyperstar_RouseMatrix)

# Rouse Matrix evaluation of Hyperdendrimer
add_subdirectory(Hyperdendrimer_RouseMatrix)

# Creator of BarabasiAlbert for m=1 and initial m0=2
# without excluded volume nor restricted functionality
add_subdirectory(BarabasiAlbert_Creator)

# Rouse Matrix evaluation of BarabasiAlbert for m=1 and initial m0=2
add_subdirectory(BarabasiAlbert_RouseMatrix)

# Rouse Matrix evaluation of BarabasiAlbert for m=1 and initial m0=2
# under excluded volume and restricted functionality f<=3
add_subdirectory(BarabasiAlbertRestrictedF3_RouseMatrix)

# Creator of BarabasiAlbert for m=1 and initial m0=2
# under excluded volume and restricted functionality f<=3
add_subdirectory(BarabasiAlbertRestrictedF3_Creation)

# Radius of gyration analyzer for single tree-like structure
add_subdirectory(ChainWalking_Analyzer_RG2)

# Rouse Matrix evaluation of SlowGrowth
add_subdirectory(SlowGrowth_RouseMatrix)

# Perfect Dendrimer generation algorithm
add_subdirectory(CreateDendrimerGeneral_FGS)

# Analyzer for extracting the topological properties of the tree (not generic graphs)
add_subdirectory(Analyzer_TopologicalPropertiesTree)

