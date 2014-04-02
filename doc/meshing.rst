.. _sumb_meshing:

Meshing
============

*1. Load or open geometry file and block file*

We generally use ANSYS ICEM-CFD to generate multiblock structured mesh for ``SUmb``.
If you already have the .tin file ( geometry file ) and .blk file ( gird block file), use ``File->Geometry->Open Geometry`` and ``File->Blocking->Open Blocking`` to load them. 
If you have to generate a grid from a geometry file (CATIA .model file or .iges file ), use ``File->Import Geometry`` and select the type of the geometry file. After that, a .tin file would be automatically generated in the current working folder. 

Using ``Blocking->Create Block (type:3D Bounding Box)`` to generate a bounding box around the farfield geometry. The farfield geometry is usually a box or hemisphere, which could be easily generated using ``Geometry->Create/Modify Surface->Standard Shapes``.

*2. Create the topology and specify the grid nodes on each edge*

Use ``Blocking->Split Blocks`` to generate the topology that would well define the geometry shape. Use ``Blocking->Ogrid Block``
to generate a bounding O-grid to simulate the boundary layer. We sometimes use nested O-grid inside another O-grid to improve the local mesh quality.
Use ``Meshing Parameters->Edge Params`` to specify the spacing and distribution law of each edge.

*3. Pre-mesh and check grid quality, check grid distribution*

In the hierarchy tree window on the left,  click on ``Blocking->Pre-Mesh``, then use ``Blocking->Pre-Mesh Quality Histoframs`` to check the mesh
quality. If negative volume exists somewhere, you have to fix them (that's another very long story..)

In the hierarchy tree window, use ``Blocking->Pre-Mesh->Scan Planes`` to scan the mesh on the slices of the grid. If the grid distribution
is not smooth, you have to change the spacing and distribution law of the related edges to fix that. 

*4. Ouput the grid to CGNS format for SUmb*

1. In the hierarchy tree window, right click ``Blocking->Pre-Mesh->Convert to Multi-Block Mesh`` , then select ``Volume`` and ``Yes``.

2. If the unit is different from what you really want, use ``Edit Mesh->Struct Domain Transformation`` to scale the mesh to the right
   unit. Then save the mesh using ``File->Mesh->Save Mesh as``.

3. Use ``Output->Select Solver`` to specify the ``Output Solver`` type as ``CGNS``.
   Use ``Output->Boundary Condition`` to specify the boundary condition type on each surfaces.
   ``BCWallViscous`` for the wall boundary such as 'wing'
   ``BCFarField`` for the farfield boundary such as 'farfiled'
   ``BCSymmetryPlane`` for the symmetry boundary such as 'symmetry'
   Then select accept.

4. Finally use ``Output->Write Input`` to output the grid as CGNS format file.
   Please use the following settings:
   ``Input Grid Type : Structured``
   ``Create BC Patches ? : No``
   ``CGNS File Output Version : 2.4``
   And you are done!
