## **DXtool**
###### Robert Alberstein, 2020-2022

DXtool is a Python package for manipulating DataExplorer (DX) volumetric maps. It is capable of reading and writing DX files, perfoming arithmetic manipulations, combining and splitting volumes from/to subvolumes, generating 1D and 2D projections (which can be plotted and/or saved to disk), as well as some multi-volume operations.

Heavily inspired by, a partial port of, and designed to be complementary to the volutil program in VMD (Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics 1996, 14.1, 33-38) - now dubbed "voltool" as of v1.9.4).

Originally developed for "Discrete Orientations of Interfacial waters direct crystallization of mica-binding proteins". Robert G. Alberstein, Jesse L. Prelesnik, Elias Nakouzi, Shuai Zhang, James J. De Yoreo, Jim Pfaendtner, F. Akif Tezcan, Christopher J. Mundy. ChemRxiv (2022) https://dx.doi/.....

### Example usage
COMING SOON

### **Library contents**

**Classes:**

	Volume : Primary container class for DX volume data.

**Volume functions**:
```
General functions
	update_info(self) -> None
		Updates all Volume attributes based on the underlying density.
	write_dxfile(self, dxoutname:str) -> Bool
		Writes the current Volume to disk with filename dxoutname. Returns True if successful.
	clone(self) -> Volume
		Returns a new Volume with identical attributes and contents.

Getter functions
	origin(self) -> [x,y,z]
		Getter function for the origin of the Volume.
	delta(self) -> [x,y,z]
		Getter function for the deltas (voxel edge lengths) of the Volume.
	size(self) -> [nx,ny,nz]
		Computes the size of the Volume (in voxels) along each direction
	total_voxels(self) -> int
		Computes the total number of voxels in the Volume.
	footer(self) -> str
		Getter function for the footer of the Volume.
	xlabels(self) -> list
		Getter function for list of x-axis labels in the Volume.
	ylabels(self) -> list
		Getter function for list of y-axis labels in the Volume.
	zlabels(self) -> list
		Getter function for list of z-axis labels in the Volume.
	value_at(self, xyzlist, byindex=False) -> float
		Gets density value at specific xyz label [x,y,z] list.

Setter functions
	TBD

Single-volume operations
	sadd(self, value:numeric) -> None
		Adds value to all voxels in the Volume (scalar add). Add negative numbers to subtract.
	smult(self, value:numeric) -> None
		Multiplies all voxels in the Volume by value (scalar multiply). Multiply fractions to divide.
	trim(self, [options]) -> None
		Trims the current Volume to the range specified by provided x/y/zmin + x/y/zmax variables (default to None, which doesn't trim at all). Use byindex flag to define whether using xyz label (str) or voxel indices (int). This operation trims the Volume in-place, so if you don't want to overwrite the data, make a new copy with clone() and trim that.
	trim_bydrop(self, [options]) -> None
		Experimental version of trim that relies more heavily on pandas operations. Need to test further to see if it is more memory/speed efficient or not.
```

**Library functions**:
```
	volume_from_dxfile(filename:str) -> Volume
		Read a DX file from disk and return a new Volume.

	volslice(Volume, position:str, axis:str) -> DataFrame
		2D slice through density along axis ("x", "y", or "z") at a given position (label for that dimension).
	project(Volume, axis:str) -> Series
		1D projected density along axis ("x", "y", or "z").

	plot_volslice(volslice:DataFrame, [options]) -> plot of volslice density
		Plot the given volslice data using Matplotlib. Can optionally explort plot to png at a given dpi.
	write_volslice(volslice:DataFrame, outname:str, [options]) -> None
		Write numeric text of volslice density to output tab-separated file.
	plot_projection(proj:Series, [options]) -> plot of volslice density
		Plot the given density projection using Matplotlib. Can optionally explort plot to png at a given dpi.
	write_projection(proj:Series, outname:str, [options]) -> None
		Write numeric text of projected volume density to output tab-separated file.

	join(vol1:Volume, vol2:Volume, axis:str) -> Volume
		Takes in two volumes and joins them along the provided axis. This function will determine the correct origin automatically. Volumes MUST BE adjacent in space to avoid a break in labels. Future versions of DXtool will check for this automatically.
	join2D(vol1:Volume, vol2:Volume, vol3:Volume, vol4:Volume) -> Volume
		Takes in four volumes and joins them as vol1+vol2 (along x axis), then vol3+vol4 (along x axis), then combine all along y axis. This allows one to take 4 subvolumes (all at the same range of z positions) and merge them into a single larger Volume without loss of information. Future versions of DXtool may remove this explicit ordering for more flexibility. If necessary, users can use the join() function as necessary to achieve the desired operations.
	join3D(vollist1:list, vollist2:list) -> Volume
		Takes in two lists each containing four Volumes. Runs join2D on each set of four Volumes independently, then joins those merged Volumes along the z axis. This allows one to take 8 subvolumes (that span a continuous range of xyz coordinates) and merge them into a single larger Volume without loss of information.

	voladd(vol1:Volume, vol2:Volume, [options]) -> None
		Takes in two Volumes and adds vol2 to vol1 voxel-wise, then returns vol1. This method WILL modify vol1 in-place, so make a new Volume object first if you do not want to overwrite the original vol1 data. This method DOES NOT do interpolation, so this method requires that Volumes have at least some matching coordinate indices. Can sum as either intersection (method=="inner") or the union (method="outer") of the volumes, where the latter requires a default value to fill into any empty voxels produced by incomplete overlap of the Volume coordinates. Set trim_intersect=True to trim vol1 to the exact dimensions of the intersection, else intersecting voxels will be summed while leaving the rest untouched.
```
