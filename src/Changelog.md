MCell Changelog
===============

**MCell v3.4.0 (release date 09-29-2016)**

  * Major changes

    * Support for periodic boundary conditions was added and comes in two
      different forms: traditional and mirrored geometry.

      With traditional PBCs, molecules simply wrap around when they hit the
      edge of the boundary. For instance, if a molecule hits the the positive X
      boundary of the periodic box, then it will wrap around to the negative X
      boundary and continue diffusing from there.

      With the mirrored geometry form, the geometry is mirrored in adjacent
      boxes and molecules do not wrap around. When a molecule hits a boundary,
      it will continue diffusing in the mirrored adjacent geometry.

    * The experimental and undocumented macromolecule feature was removed.
      Note: this will eventually be replaced with an improved version in a
      future release.

    * Support for DReAMM output was completely removed, since DReAMM has been
      completely superseded by CellBlender. Any users who still need support
      for DReAMM output should use MCell 3.3

  * Minor changes

    * More consistent handling of warnings and errors.

    * Removed REACTION_GROUPs. This was never documented anywhere, nor did it
      serve much of a purpose, so its removal shouldn't have much impact.
    
    * Fixed multiple memory leaks.

    * A number of bugs were fixed. Here's a list of several notable ones:

      * Removed false positives when checking for overlapping triangles.

      * Explicitly disallow trimolecular reactions and microscopic
        reversibility. Previously this would cause MCell to crash.

      * Fixed cases where molecules could leak out of meshes with sharp
        corners.

      * Fixed counting bug when using ALL_ENCLOSED keyword.

    * MCell now supports the CMAKE build system. The autoconf build system is
      now deprecated.
