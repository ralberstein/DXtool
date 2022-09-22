"""
DXtool: a Python library for manipulation of Data Explorer (DX) volume files.

Classes:
	Volume

Volume functions:
	update_info(self) -> None
		Updates all Volume attributes based on the underlying density.
	write_dxfile(self, dxoutname:str) -> Bool
		Writes the current Volume to disk with filename dxoutname. Returns True if successful.
	clone(self) -> Volume
		Returns a new Volume with identical attributes and contents.

	# Getter functions
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

	# Setter functions
	TBD

	# Single-volume operations
	sadd(self, value:numeric) -> None
		Adds value to all voxels in the Volume (scalar add). Add negative numbers to subtract.
	smult(self, value:numeric) -> None
		Multiplies all voxels in the Volume by value (scalar multiply). Multiply fractions to divide.
	trim(self, [options]) -> None
		Trims the current Volume to the range specified by provided x/y/zmin + x/y/zmax variables (default to None, which doesn't trim at all). Use byindex flag to define whether using xyz label (str) or voxel indices (int). This operation trims the Volume in-place, so if you don't want to overwrite the data, make a new copy with clone() and trim that.
	trim_bydrop(self, [options]) -> None
		Experimental version of trim that relies more heavily on pandas operations. Need to test further to see if it is more memory/speed efficient or not.


Library functions:
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

TODO:
1) Convert volslice, project to Volume class functions
2) More detailed documentation for functions
3) Add more explicit error handling & unit tests
"""

#__future__ goes here if it was ever needed
__author__ = "Robert G. Alberstein"
__version__ = "0.1.0"

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
##### CONSTANTS #####
VERSION = __version__

##### CLASSES #####
class Volume:
	"""
	The Volume class holds all requisite information to define a specific DX volume, such as the origin, voxel spacing, size, footer info, and the density itself.

	Attributes:
	-----------
	origin:	starting position of the density map
	delta:	xyz voxel edge widths (typically equal, but not necessary)
	density:actual 3D pandas dataframe structure holding the density values
	footer:	output string for the bottom of the DX file, optional.

	Methods:
	--------
	TBD

	"""
	def __init__(self, origin, delta, density, footer=""):
		# Input parameters
		self._origin = origin
		self._delta = delta
		self._density = density
		self._footer = footer

		# Calculated parameters
		self._numx = len(self._density.iloc[0].columns)
		self._numy = len(self._density.iloc[0].index)
		self._numz = len(self._density)
		self._total_voxels = self._numx * self._numy * self._numz
		self._xlabels = self._density.iloc[0].columns
		self._ylabels = self._density.iloc[0].index
		self._zlabels = self._density.index

	def update_info(self):
		"""
		Recalculates the origin and size of the current Volume based on it's actual volumetric data. Called after modifications such as trimming or moving.
		"""
		# Update origin
		self._origin = [self._density.iloc[0].columns.min(), self._density.iloc[0].index.min(), self._density.index.min()]
		# Update voxel positions
		self._xlabels = self._density.iloc[0].columns
		self._ylabels = self._density.iloc[0].index
		self._zlabels = self._density.index
		# Update size
		self._numx = len(self._density.iloc[0].columns)
		self._numy = len(self._density.iloc[0].index)
		self._numz = len(self._density)
		self._total_voxels = self._numx * self._numy * self._numz


	### Volume I/O functions ###
	def write_dxfile(self, dxoutname):
		"""
		Writes the current Volume to disk.

		Args:
			dxoutname (string): filename of the DX file to write out. Full path must be included if target folder is not current directory.

		Returns:
			Bool: True if successful, False otherwise. Prints message to STDOUT.

		Raises:
			TBD

		Usage:
			>>> write_dxfile("./myoutputvolume.dx")
			Successfully wrote file ./myoutputvolume.dx to disk!
			True
		"""
		### Prepare for output ###
		# Needed to avoid truncation of array
		np.set_printoptions(threshold=np.inf)

		# Check for non-integer multiples of 3 (necessary for reshaping for output)
		self.num_extra_vals = self._total_voxels % 3
		# Calculate number of output rows
		self.num_rows = (self._total_voxels // 3)

		# Convert density dataframe to 3D np array & transpose for output
		self.np_density = np.array([xyvals.to_numpy() for xyvals in self._density]).T
		# Flatten array to 1D for output (can also overwrite np_density variable)
		self.density_flat = np.reshape(self.np_density, -1)

		# If already divisble by 3, can directly reshape to Nx3 array
		if (self.num_extra_vals == 0):
			self.outarray = np.reshape(self.density_flat, (self.num_rows, -1))
		# Else append extra values to the footer for output and reshape the remainder
		else:
			self._footer = " ".join([str(v) for v in self.density_flat[-self.num_extra_vals:]]) + "\n" + self._footer
			self.outarray = np.reshape(self.density_flat[:-self.num_extra_vals], (self.num_rows, -1))

		### Write the file using savetxt function to handle formatting
		# Define the header (separated into individual lines here for clarity only)
		self.outheader = f"# Processed by PyDXtool v{VERSION}\n"
		self.outheader += f"object 1 class gridpositions counts {self._numx} {self._numy} {self._numz}\n"
		self.outheader += f"origin {self._origin[0]} {self._origin[1]} {self._origin[2]}\n"
		self.outheader += f"delta {self._delta[0]} 0 0\n"
		self.outheader += f"delta 0 {self._delta[1]} 0\n"
		self.outheader += f"delta 0 0 {self._delta[2]}\n"
		self.outheader += f"object 2 class gridconnections counts {self._numx} {self._numy} {self._numz}\n"
		self.outheader += f"object 3 class array type double rank 0 items {self._total_voxels} data follows"

		# Write out array as string. comments="" -> prints header/footer w/o leading '#'
		if (self._footer != None):
			np.savetxt(f"{dxoutname}", self.outarray, fmt="%s", header=self.outheader, footer=self._footer, comments="")
		else:
			np.savetxt(f"{dxoutname}", self.outarray, fmt="%s", header=self.outheader, comments="")

	def clone(self):
		"""
		Returns an exact duplicate of the current Volume by initializing with same origin, delta, and footer, and using pandas .copy() to preserve data.
		"""
		return Volume(self._origin, self._delta, self._density.copy(), self._footer)


	### Getter functions ###
	def origin(self):
		"""
		Getter alias for volume origin.
		"""
		return self._origin

	def delta(self):
		"""
		Getter alias for volume delta (voxel edge lengths).
		"""
		#return self._delta
		# Return list of floats instead of strings
		return [float(val) for val in self._delta]

	def size(self):
		"""
		Getter alias for volume size in number of voxels [x,y,z].
		"""
		return [self._numx, self._numy, self._numz]

	def total_voxels(self):
		"""
		Getter alias for total number of voxels.
		"""
		return self._total_voxels

	def footer(self):
		"""
		Getter alias for footer of this Volume instance.
		"""
		return self._footer

	def xlabels(self):
		"""
		Getter for list of x-axis labels in this Volume.
		"""
		return self._xlabels

	def ylabels(self):
		"""
		Getter for list of y-axis labels in this Volume.
		"""
		return self._ylabels

	def zlabels(self):
		"""
		Getter for list of z-axis labels in this Volume.
		"""
		return self._zlabels

	def value_at(self, xyzlist, byindex=False):
		"""
		Gets data value at specific xyz label [x,y,z] list.
		"""
		x,y,z = xyzlist
		if byindex:
			x,y,z = self._xlabels[x], self._ylabels[y], self._zlabels[z]
		return self._density.at[z].at[y,x]


	### Setter functions ###


	### Single-volume operations ###
	def sadd(self, value):
		"""
		Loops over all XY slices in dataframe and adds scalar value to all positions via broadcasting.

		Args:
			value (float or int): value to add to every voxel in Volume.

		Raises:
			TBD
		"""
		for z in self._density:
			z += value

	def smult(self, value):
		"""
		Loops over all XY slices in dataframe and multiplies all positions by value via broadcasting.

		Args:
			value (float or int): value to multiply every voxel in Volume by.

		Raises:
			TBD
		"""
		for z in self._density:
			z *= value

	def trim(self, xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, byindex=False):
		"""
		Trims current volume to specified ranges (values inclusive). Defaults to label selection, but can use index with byindex=True. If you don't want to overwrite this data directly, make a new Volume with clone() and trim that.
		"""
		# Do trimming by overwriting volume with slicing
		if (byindex == True):
			self._density = pd.Series([pd.DataFrame(self._density[z].iloc[ymin:ymax,xmin:xmax]) for z in self._density[zmin:zmax].index], index=self._density[zmin:zmax].index)
		else:
			self._density = pd.Series([pd.DataFrame(self._density[z].loc[ymin:ymax,xmin:xmax]) for z in self._density[zmin:zmax].index], index=self._density[zmin:zmax].index)
		# Update volume info to reflect new size
		self.update_info()

	def trim_bydrop(self, xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0, byindex=False):
		"""
		Trims in-place using pandas operations.
		"""
		# Define lists of labels to drop for each axis
		toremove_z = list(self._density[:zmin].index[:-1]) + list(self._density[zmax:].index[1:])
		toremove_y = list(self._density[zmin].loc[:ymin].index[:-1]) + list(self._density[zmin].loc[ymax:].index[1:])
		toremove_x = list(self._density[zmin].loc[:,:xmin].columns[:-1]) + list(self._density[zmin].loc[:,xmax:].columns[1:])

		# Drop z (if inplace=False, returns a copy)
		self._density.drop(toremove_z, inplace=True)
		# Drop xy for each layer of volume
		for z in self._density:
			z.drop(index=toremove_y, columns=toremove_x, inplace=True)
		# Update volume info to reflect new size
		self.update_info()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
##### MODULE FUNCTIONS #####
def volume_from_dxfile(dxfile, savefooter=True):
	"""
	Takes the filename to a DX file and returns a Volume class.

	Args:
		dxfile (string): filename of the DX file to read in
		savefooter (bool, optional): True (default) keeps footer of input file. If False, will use generic footer information when writing to disk. Footer does not affect volume in any way during operations in-memory.

	Returns:
		Volume: new Volume class holding DX information if successful.

	Usage:
		>>> volume_from_dxfile("myfile.dx", savefooter=False)
		<PyDXtool.Volume at 0x7fd02a87abe0>
	"""
	# Open the file and store all lines in list
	fin = open(dxfile, "r")
	dxcontent = fin.readlines()
	fin.close()

	### PROCESS HEADER ###
	# Read the header, removing any initial comments
	while(dxcontent[0].startswith("#")):
		dxcontent.pop(0)

	# Get gridcounts
	currline = dxcontent.pop(0).split()
	numx, numy, numz = int(currline[5]), int(currline[6]), int(currline[7])

	# Get origin
	currline = dxcontent.pop(0).split()
	ox, oy, oz = float(currline[1]), float(currline[2]), float(currline[3])

	# Get delta values (store as string for precise output, cast to float for calculating labels)
	dx_str = dxcontent.pop(0).split()[1]
	dy_str = dxcontent.pop(0).split()[2]
	dz_str = dxcontent.pop(0).split()[3]
	dx = float(dx_str)
	dy = float(dy_str)
	dz = float(dz_str)

	# Finish out header (object values we don't need)
	dxcontent.pop(0)	#object 2 class gridconnections (same as object 1 counts)
	dxcontent.pop(0)	#object 3 class array (can calc numx*numy*numz = data follows)

	### PROCESS FOOTER ###
	# Clear out any empty lines at the end of the file (so dealing with data only)
	while (dxcontent[-1].startswith("\n")):
		dxcontent.pop(-1)

	# Process and save footer only if exists & specified
	if (savefooter == True):
		# For multi-line footers (e.g. APBS), keep all lines
		footline = -1
		# Loop backwards from end
		while (dxcontent[footline].startswith("\n") != True):
			# Checks if we are in density data by attempting cast to float
			try:
				# If we are in density, undo one line counter & store footer
				if (type(float(dxcontent[footline].split()[0])) == float):
					footline += 1
					# If there was NO footer (footline=0) use generic footer
					if (footline == 0):
						footer = "object \"processed with PyDXtool\" class field"
						break
					# Otherwise save all footer lines as a single string
					else:
						footer = "".join(dxcontent[footline:])
						break
			# If cast to float fails, we are still dealing with footer strings
			except ValueError:
				footline -= 1

	# Else we use a generic footer, just so programs like VMD don't complain
	else:
		footer = "object \"processed with PyDXtool\" class field"

	# And finally remove any final trailing newlines
	while(dxcontent[-1] == "\n"):
		dxcontent.pop(-1)

	### READ IN DENSITY ###
	# Should only have datapoints left, so can count # of lines remaining
	numlines = len(dxcontent)

	# Convert to 1D list of 3x value strings, then single space-separated string
	dens_str = " ".join([" ".join(val.split()) for val in dxcontent])
	# Convert to 1D np array
	dens = np.fromstring(dens_str, dtype=float, sep=" ")
	# Reshape to 3D array of Z,Y,X
	density = np.reshape(dens, (numx, numy, numz)).T
	#Note: possible that reshape order='F' would work too, by changing the first
	# element change fastest instead of the last, but .T achieves the same goal

	# Create dataframe axis labels as list of coordinates, rounded to 0.001 A precision
	x_coords = [ox + (dx*x) for x in range(numx)] #All x_axis labels
	y_coords = [oy + (dy*y) for y in range(numy)] #All y_axis labels
	z_coords = [oz + (dz*z) for z in range(numz)] #All z_axis labels
	x_coords = [round(x, 3) for x in x_coords]
	y_coords = [round(y, 3) for y in y_coords]
	z_coords = [round(z, 3) for z in z_coords]

	# Convert volume to dataframe
	volume_df = pd.Series([pd.DataFrame(density[z], index=y_coords, columns=x_coords) for z in range(len(z_coords))], index=z_coords)

	# Return the Volume object (footer already empty string if not saved)
	return Volume([ox,oy,oz], [dx_str,dy_str,dz_str], volume_df, footer)

def volslice(vol, pos, axis):
	"""
	Returns a DataFrame corresponding to a volume slice at the given position along the specified axis. Axes are: x (YZ plane), y (XZ plane), or z (XY plane). Slices retain dataframe convention of top-left value being the minimum for both axes.
	"""
	# Ensure axis selected properly
	if (axis.lower() not in "xyz"):
		print("Invalid axis selection. Please specify x, y, or z.")
		return False
	# Generate the slices
	if (axis.lower() == "x"):
		xslice = [zslice.loc[:,pos:pos].transpose().to_numpy()[0] for zslice in vol._density]
		return pd.DataFrame(xslice, index=vol._density.index, columns=vol._density.iloc[0].loc[:,pos:pos].index)
	elif (axis.lower() == "y"):
		yslice = [zslice.loc[pos:pos,:].to_numpy()[0] for zslice in vol._density]
		return pd.DataFrame(yslice, index=vol._density.index, columns=vol._density.iloc[0].loc[pos:pos,:].columns)
	else:
		# axis == "z"
		return vol._density.loc[pos]

def project(vol, axis):
	"""
	Calculates a 1D projection for the entire volume along the specified axis and returns it as a pandas Series.
	"""
	# Ensure axis selected properly
	if (axis.lower() not in "xyz"):
		print("Invalid axis selection. Please specify x, y, or z.")
		return False
	# Generate the projections
	if (axis.lower() == "x"):
		all_xval = vol._density.iloc[0].columns
		return pd.Series([volslice(vol, xval, "x").mean().mean() for xval in all_xval], index=all_xval)
	elif (axis.lower() == "y"):
		all_yval = vol._density.iloc[0].index
		return pd.Series([volslice(vol, yval, "y").mean().mean() for yval in all_yval], index=all_yval)
	else:
		# axis == "z"
		all_zval = vol._density.index
		return pd.Series([volslice(vol, zval, "z").mean().mean() for zval in all_zval], index=all_zval)

def plot_volslice(volslice, origin="lower", cmap="binary", vmin=0.0, vmax=0.1, outname=None, dpi=100):
	"""
	Plots a given volume slice as a fullscreen image (no axes). Since image convention starts with the lowest XY pixel at the bottom-left (not top-left as it is in dataframes), origin flag is set to lower by default ([xmin,ymin] at bottom-left), though this can be prevented by setting origin=upper ([xmin,ymin] at top-left).

	TODO:
	1) Make sure having 400 DPI thing doesn't mess up my images at all vs 100 dpi (check that 4x4 pixel values are all identical)
	1B) Change interpolation back to 'nearest' if necessary

	2) Add flag to plot +/- axes for those just wanting a plot rather than a standalone image.
	3) Add plot defaults from plot_projections if axes are shown.
	"""
	num_cols = len(volslice.columns)
	num_rows = len(volslice.index)

	# Typically dataframe will come with upper origin
	if (volslice.index[0] == volslice.index.min()):
		if (origin == "upper"):
			pass
		# For lower origin, need to flip the slice
		elif (origin == "lower"):
			volslice = volslice.iloc[::-1]
			# If it didn't work, something went wrong
			if (volslice.index[0] == volslice.index.max()):
				pass
			else:
				print("Volume slice orientation error. Could not flip slice to bottom-left origin.")
				return False
		# If it didn't pass these statements, probably origin checks failed
		else:
			print("Error plotting volume slice. Did you correctly specify origin as 'lower' or 'upper'?")

	# Division by 100.0 ensures 1:1 pixel:voxel ratio
	fig = plt.figure(frameon=False, figsize=(num_cols/100.0,num_rows/100.0))
	# Axes ensure full-screen plotting
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	ax.set_axis_off()
	fig.add_axes(ax)

	# Plot the figure
	ax.imshow(volslice, aspect="equal", interpolation='none', cmap=cmap, vmin=vmin, vmax=vmax)
	if (outname != None):
		# 400 DPI a nice 2x2 expansion of pixels
		fig.savefig(f"{outname}", dpi=dpi)
	else:
		plt.show()
	plt.close()

def write_volslice(volslice, outname, origin="lower", writelabels=False):
	"""
	Writes the provided volume slice to a text file as a tab-separated 2D matrix. Origin flag follows matplotlib imshow convention: upper = [xmin, ymin] at top-left, lower = [xmin, ymin] at bottom-left.
	"""
	num_cols = len(volslice.columns)
	num_rows = len(volslice.index)

	# Typically dataframe will come with upper origin
	if (volslice.index[0] == volslice.index.min()):
		if (origin == "upper"):
			pass
		# For lower origin, need to flip the slice
		elif (origin == "lower"):
			volslice = volslice.iloc[::-1]
			# If it didn't work, something went wrong
			if (volslice.index[0] == volslice.index.max()):
				pass
			else:
				print("Volume slice orientation error. Could not flip slice to bottom-left origin.")
				return False
		# If it didn't pass these statements, probably origin checks failed
		else:
			print("Error writing volume slice file. Did you correctly specify origin as 'lower' or 'upper'?")

	# Write slice to disk
	volslice.to_csv(f"{outname}", sep="\t", header=writelabels, index=writelabels)
	print(f"Successfully wrote {outname} to disk!")

def plot_projection(proj, dims=None, xlabel=None, ylabel=None, title=None, figsize=(8,6), outname=None):
	"""
	Wrapper for making a basic plot of a calculated 1D projection.
	"""
	# Initialize & configure the plot
	fig, ax = plt.subplots(1,1, figsize=figsize)
	if (title != None):
		ax.set_title(title)
	if (dims != None):
		ax.axis(dims)
	xtick_spacing = int(len(proj.index) / 5)

	# Make the plot
	proj.plot(color="Blue", linewidth=2.0, xlabel=xlabel, ylabel=ylabel, ax=ax)
	plt.tight_layout()

	# Save or show the plot, then close it
	if (outname != None):
		fig.savefig(f"{outname}", dpi=300)
	else:
		plt.show()
	plt.close()

def write_projection(proj, outname, xlabel="distance", ylabel="projection"):
	"""
	Writes the provided 1D projection to a text file as a tab-separated pair of columns.
	"""
	proj.to_csv(f"{outname}", sep="\t", index_label=xlabel, header=[ylabel])

##### MULTI-VOLUME OPERATIONS ###
def join(vol1, vol2, axis):
	"""
	Joins two volumes along a specified axis. Automatically ensures that volumes are merged left->right along x, and bottom->top along y/z to preserve volume origin. Provided volumes MUST match along merging dimension (index for x/z, columns for y) AND have identical voxel sizes to be merged. At the moment no check is performed to ensure that the volumes are adjacent, so the user must ensure this to avoid breaks in the voxel labels.

	TODO: remove print functions?
	"""
	# Check compatibility of volumes
	if (vol1._delta != vol2._delta):
		print("These volumes do not have matching voxel sizes!! Join operation failed.")
		return False
	# If footers match, keep for consistency, else use first volume's as default
	if (vol1._footer != vol2._footer):
		print("The volume footers do not match! Defaulting to first volume's footer.")

	# Store indices and columns for convenience
	vol1_zindex = vol1._density.index
	vol2_zindex = vol2._density.index
	vol1_indices = vol1._density.iloc[0].index
	vol2_indices = vol2._density.iloc[0].index
	vol1_columns = vol1._density.iloc[0].columns
	vol2_columns = vol2._density.iloc[0].columns

	# Ensure compatibility along join axis and define origin
	if (axis.lower() == "x"):
		if (vol1_indices.equals(vol2_indices) != True):
			print("Volume indices do not line up! Cannot merge along this axis.")
			return False
		# Pick the left-most volume to start with
		if (vol1_columns.min() < vol2_columns.min()):
			start_vol = vol1
			next_vol = vol2
		else:
			start_vol = vol2
			next_vol = vol1
		# Do the merging
		joined_df = pd.Series([pd.concat([start_vol._density[z], next_vol._density[z]], axis=1, join="outer") for z in start_vol._density.index], index=start_vol._density.index)
	elif (axis.lower() == "y"):
		if (vol1_columns.equals(vol2_columns) != True):
			print("Volume columns do not line up! Cannot merge along this axis.")
			return False
		# Pick the bottom-most volume to start with
		if (vol1_indices.min() < vol2_indices.min()):
			start_vol = vol1
			next_vol = vol2
		else:
			start_vol = vol2
			next_vol = vol1
		# Do the merging
		joined_df = pd.Series([pd.concat([start_vol._density[z], next_vol._density[z]], axis=0, join="outer") for z in start_vol._density.index], index=start_vol._density.index)
	# Special treatment for z because it's a series instead
	elif (axis.lower() == "z"):
		# Pick the left-most volume to start with
		if (vol1_zindex.min() < vol2_zindex.min()):
			start_vol = vol1
			next_vol = vol2
		else:
			start_vol = vol2
			next_vol = vol1
		# Do the merging
		joined_df = pd.concat([start_vol._density, next_vol._density], axis=0, join="outer")
	else:
		print("Invalid axis selection. Please specify x, y, or z.")
		return False

	# Return the joined Volume object
	#print(f"Volumes successfully joined along {axis.lower()}-axis.")
	return Volume(start_vol._origin, start_vol._delta, joined_df, vol1._footer)

def join2D(vol1, vol2, vol3, vol4):
	"""
	Joins 4 Volumes along xy to form a single larger 2x2x1 layer. The first two and last two volumes should be joinable along the x-axis, however the order is not important. The join function will identify the correct origin before concatenation.
	"""
	join1 = join(vol1, vol2, "x")
	join2 = join(vol3, vol4, "x")
	return join(join2,join1, "y")

def join3D(vollist1, vollist2):
	"""
	Joins 8 Volumes into a single larger 2x2x2 cube. Takes in two lists of four Volumes. The first two and last two Volumes of each list should be joinable along the x-axis. The order of the rest is not important, as the functions will identify the correct origin during each concatenation operation.

	TODO: technically still needs to be tested from lists -> 3D, but all component functions work.
	"""
	join1 = join2D([vollist1[0], vollist1[1], vollist1[2], vollist1[3]])
	join2 = join2D([vollist2[0], vollist2[1], vollist2[2], vollist2[3]])
	return join(join1, join2, "z")

def voladd(vol1, vol2, method="inner", fill_value=0.0, trim_intersect=False):
	"""
	Takes two Volumes and adds vol2 to vol1 voxel-wise, then returns vol1. This method WILL modify vol1 in-place, so make a new Volume object first if you do not want to overwrite the original vol1 data. This method DOES NOT do interpolation, so this method requires that Volumes have at least some matching coordinate indices. Can sum as either intersection (method=="inner") or the union (method="outer") of the volumes, where the latter requires a default value to fill into any empty voxels produced by incomplete overlap of the Volume coordinates. Set trim_intersect=True to trim vol1 to the exact dimensions of the intersection, else intersecting voxels will be summed while leaving the rest untouched.

	NOTE: df.add unions indices by default, so must trim first for intersection
	NOTE: when later add volsub, volmult, etc should write a separate function that does the intersection calculations to call from each. After that the only difference is the pandas syntax of pd.add, pd.sub, pd.mul, pd.div, pd.pow ...
	"""
	# For intersection, need to define overlapping regions
	if (method == "inner"):
		# Get True/False list for which labels are matching
		match_zlabels = (vol1.zlabels() == vol2.zlabels())
		match_ylabels = (vol1.ylabels() == vol2.ylabels())
		match_xlabels = (vol1.xlabels() == vol2.xlabels())

		# Then check that there is at least one set of matching indices in volumes
		if ((True in match_zlabels) and (True in match_ylabels) and (True in match_xlabels)):
			# Define range of matching labels
			inner_zlabels = vol1.zlabels()[match_zlabels]
			inner_ylabels = vol1.ylabels()[match_ylabels]
			inner_xlabels = vol1.xlabels()[match_xlabels]

			# Do the addition only for intersecting labels (add 0.0 to non-matching)
			for z in inner_zlabels:
				vol1._density.at[z] = vol1._density.at[z].add(vol2._density.at[z].loc[inner_ylabels,inner_xlabels], fill_value=0.0)

			# Return the volume (trimmed or not)
			if (trim_intersect):
				vol1.trim(inner_xlabels.min(), inner_xlabels.max(),
						  inner_ylabels.min(), inner_ylabels.max(),
						  inner_zlabels.min(), inner_zlabels.max())
			print("Successfully added vol1 to vol2 by intersection.")
			return True
		# Else there were no overlapping voxel labels
		else:
			print("No matches between any voxel indices in the two volumes! Cannot add these volumes via intersection (method=\"inner\").")
			return False

	# For addition by union, just add datframes over the series
	elif (method == "outer"):
		print("Successfully added vol1 to vol2 by union.")
		vol1._density = vol1._density.add(vol2._density, axis=0, fill_value=fill_value)
		return True
	else:
		print("Unexpected method argument! Please use \"inner\" for intersection or \"outer\" for union.")
		return False
