# //////////////////////////////////////////////////////////////////////////////
docs_main = dict(
    main = "Description for program",
    calc = "description...",
    bse = "description...",
    cluster = "description...",
    cmap = "description...",
    pca = "description...",
    pyinteraph = "description...",
    rama = "description...",
    rgyr = "description...",
    rmsd = "description...",
    rmsf = "description...",
    sasa = "description...",
)

# //////////////////////////////////////////////////////////////////////////////
docs_wald = dict(
    main = """
  +===============================================+
  | Water Average Local Densities (WALD) pipeline |
  +===============================================+
    * Go through step 0 by selecting the reference frame of interest and the allowed RMSD treshold of "similarity".
    * Execute steps 1 and 2.
    * BEFORE step 3:
      + Open the VMD file generated in step 1.
      + Reset the view and toggle off the axes.
      + Select the reference frame as indicated in the console.
      + File -> Render... -> "Render the current scene using:" Wavefront (OBJ and MTL).
      + Rename the file with the same run prefix as the VMD file, but with an OBJ extension (e.g. mt2_rep0.obj).
      + Start Rendering.
    * Execute step 3.
    * Visualize the result with step 4 as many times as needed.
""",
    rmsd = "description...",
    calc = "description...",
    vmd = "description...",
    meshes = "description...",
    plot = "description...",
)

# //////////////////////////////////////////////////////////////////////////////
