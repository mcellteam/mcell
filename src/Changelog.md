MCell Changelog
===============

**MCell v3.3.0 (release date 04-23-2015)**

  * Changed default visualization mode to `CELLBLENDER`.

    **NOTE:** This may break models for users who still use DREAMM
    combined with this default setting to produce DREAMM viz data.

  * Added exact concentration based volume releases for single
    regions.

    Previously, concentration based volume releases would happen
    by picking molecule positions in the closed volume's bounding box
    (at the proper concentration) and only molecules that would fall
    within the actual region would be placed. However, since the
    computation of the number of molecules would happen based on the
    volume of the bounding box, the number of released molecules
    within the region would fluctuate from seed to seed and within
    a given seed for different regions of identical volume. This
    behaviour seemed counterintuitive. Regions of identical
    volume should receive identical numbers of molecules (placed
    at random positions).

    The current commit enables these exact region releases based
    on region volume for single closed regions without any boolean
    operators. For regions composed via boolean operators the code
    still uses the previous bounding box based release. This is due
    to the fact that the computation of volume for regions composed
    via boolean operators requires sophisticated CSG methods
    currently missing from MCell. This may be addressed in the
    future.

  * Fixed several issues in the checkpointing code. Previously,
    the checkpointing code assumed that past (checkpointed)
    simulations were run with the same timestep as the resumed
    run and scheduled events based on this assumption. Thus, if
    the resumed similation was run with a different timestep,
    events were scheduled at the incorrect time (either too early
    or too late). This MCell version corrects this issue for
    release events (such as release patterns) and unimolecular
    reaction events by re-scheduling release events appropriately.

    **NOTE 1:** Starting with version 3.3, MCell checkpoint files
    contain an API tag and MCell 3.3 and newer should be able to
    parse any checkpoint file, including legacy ones created with
    pre 3.3 versions.

    **NOTE 2:** The lifetimes for unimolecular reactions are always
    recomputed when continuing from a checkpoint. As such, results
    generated when using checkpointing will not exactly match results
    without checkpointing. Since lifetimes are always recomputed, it
    is now possible to change unimolecular rate constants between
    checkpoints.

  * Added ability to continue checkpointed simulations. Previously,
    simulations that were checkpointed via `CHECKPOINT_ITERATIONS`
    would stop after the specified number of checkpoint iterations
    and had to be restarted manually. In contrast, simulations that
    were checkpointed via `CHECKPOINT_REALTIME` could continue
    via the `NOEXIT` keyword. The `NOEXIT` keyword is now also  
    available via `CHECKPOINT_ITERATIONS`.

        CHECKPOINT_ITERATIONS = <integer> <EXIT|NOEXIT>

    with `EXIT` being the default for backward compatibility.

    In addition, MCell now has the keyword

        KEEP_CHECKPOINT_FILES = <TRUE|FALSE>

    When enabled, this keyword triggers MCell keeping backups of
    all checkpoint files generated along the way by tagging
    files with their respective iteration number. The default
    for this option is `FALSE`, in which case subsequent
    checkpoints will overwrite previous checkpoint files
    consistent with the current behavior.

  * Adressed issue in which volume molecules created by surface
    reactions could be released on the wrong side of a closed
    object.

  * Enable counting or release inside (closed) regions which are
    part of an overall non-manifold object (e.g. a cube with
    triangles sticking out at various angles). Previously, MCell
    was not able to deal with these cases.

  * Fixed incorrect `CONCENTRATION_CLAMP` behavior when multiple
    clamps were applied to individual surface classes.

  * Fixed problem in which reactive surface boundaries would not
    be properly detected in the case of multiple overlapping
    reactive surface classes.

  * Fixed bug in which volume molecules would not be removed
    when using a release site with negative molecule count and a
    closed region as release shape.

  * Disallow negative molecule releases via `SHAPES` and `LISTS`.
    Negative molecule releases (to remove them) are only
    supported on or within regions.

  * The maximum number of iterations is now limited by the size of
    long long instead of int previously.

  * Removed the legacy DX visualization output format.

  * Added several improvements to floating point math within MCell.

  * Fixed several memory leaks and segmentation faults.

  * The previous python based MCell test-suite was removed and
    replaced by completely rewritten test-framework
    [nutmeg](https://github.com/haskelladdict/nutmeg).

  * This release includes a major refactoring across the complete
    MCell codebase. The code was significantly cleaned up and
    simplified with a focus on working toward an API.
