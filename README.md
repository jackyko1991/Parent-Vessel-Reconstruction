# Parent-Vessel-Reconstruction
An objective approach to digital removal of saccular aneurysms.

This is the C++ implementation of [VMTK PARENT VESSEL RECONSTRUCTION](http://www.vmtk.org/tutorials/ParentVesselReconstruction.html). The main objective of the software is to reconstruct parent vessel from aneurysm one. Detail theory can be found from [Ford's publication](https://www.birpublications.org/doi/epub/10.1259/bjr/67593727)

**Note: Clipping section now is identified manually with given clipping center and length for robustness**

## Usage

## Workflow
1. Compute forward centerline with user defined source and target points
2. Pick the section to clip
3. Clip centerline and corresponding Voronoi diagram
4. Smooth Vornoi diagram
4. Perform centerline and Voronoi diagram interpolation
5. (Optional) End points extension
6. Reconstruct vessel surface

## References
- [Ford et al, An objective approach to digital removal of saccular aneurysm: techniques and applications. BJR, 2009, ss55-61](https://www.birpublications.org/doi/epub/10.1259/bjr/67593727)
- [VMTK PARENT VESSEL RECONSTRUCTION](http://www.vmtk.org/tutorials/ParentVesselReconstruction.html)
