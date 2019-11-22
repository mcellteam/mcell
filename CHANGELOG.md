CHANGELOG
===============================================================================

Changes in MCell 3.5
-------------------------------------------------------------------------------

* Major changes:
  * NFSim reaction engine gets the same seed as mcell, no longer it generates its seed from current time 

* Bug Fixes:
  * Fixed several issues in compartments support for BNGL style reactions
  * Fixed case when memory usage grew indefinitely when logging BNGL style reactions
  
* Optimizations:
  * Optimized search in nfSim reaction cache (up to 2x speedup for some MCellR models)

