# 3dThesis

Heat transfer code utilizing a nondimensionalized semi-analytic solution to moving heat sources with a 3D Gaussian power density

A detailed explanation of the mathematics can be found in [Stump and Plotkowski](CITATION.bib).

## Citing

If you use 3dThesis in your work, please cite the [Stump and Plotkowski](CITATION.bib).
The original release is available on [DOE Code](https://doi.org/10.11578/dc.20200909.3). Subsequent releases are updated on [Zenodo](https://doi.org/10.5281/zenodo.13686743).

## Acknowledgements
3dThesis has been authored by UT-Battelle, LLC under Contract No. DE-AC05-00OR22725 with the U.S. Department of Energy.

3dThesis was originally co-sponsored by the U.S. Department of Energy, Office of Energy Efficiency and Renewable Energy, Advanced Manufacturing Office and the Office of Electricity Delivery and Energy Reliability (OE) - Transformer Resilience and Advanced Components (TRAC) Program.

## License

3dThesis is distributed under an [open source 3-clause BSD license](LICENSE).

## Build
3dThesis uses `make` to build. Options inside the `makefile` can be changed as needed for the local hardware. Run `make` from the commandline in order to build 3dThesis.

## Run
By default, 3dThesis will be built within `./build/application`. The input files described below can be modified as needed and then 3dThesis can be run on the commandline: `./build/application/3dThesis`

## Inputs
This section shows how to create a custom simulation with all tunable parameters.
If a necessary keyword is not included, often the program will default to values found in the test case. Note when editing files, it is important to begin each group with a brace `{`, as well as close each group `}` as in the example files.

When run from the command line, an input file is optional (e.g. `./3DThesis ParamInput.txt`), defaulting to `TestInputs/ParamInput.txt`. If `Name` is set to `TestSim`, then the data for the simulation can be found at `TestSim/Data`. All the files under `Simulation` are necessary, the rest are optional.

### Simulation Files
This set of files dictate everything having to do with the physics of the simulation, such as the material, the heat source, and the path of the heat source. All the inputs in this section are necessary.

- Simulation Name
- Mode file
- Material file
- Beam file
- Path file

#### Mode File
This file contains all information about how the simulation is run. There are two types of modes: `Solidification` and `Snapshots`. Only one mode can be run at a time.

Solidification
- `Tracking`
  - `None` always calculates all points. 
  - `Volume` always calculates all molten points.
  - `Surface` always calculates the surface of the molten pool. This tracking is fastest if only simulation solidification conditions.
- `Timestep `: interval between when the temperatures are evaluated. Setting this value too low results in unnecessary computational cost; too high and solidification will be missed.
- `OutputFrequency`: How many timesteps between subsequent outputs. The final result is always output so if intermediate results are not desired, set this value to be a large number.
- `Secondary `: Calculates the secondary solidification characteristics. This is experimental and should never be used but has been left in for research purposes.

Snapshots
- `Times` Times at which to calculate the temperature snapshots (can either choose times OR scanFracs but not both) 
- `ScanFracs`: Times, in terms of percentage of the length of a scan path, at which to calculate the temperature snapshots
- `Tracking`
  - `None` always calculates all points. 
  - `Volume` always calculates all molten points.
  - `Surface` always calculates the surface of the molten pool. For temperature snapshots, Volume and Surface behave identically.

#### Material File
This file contains all the material constants to be used by the simulation. The CET parameters are not necessary and are used to calculate the equiaxed grain fraction according to the model by Gaumann et al., Acta Materialia, 2001.
- `T_0`: Initial/background temperature (K)
- `T_L`: Liquidus temperature (K)
- `k`: Thermal conductivity (W/(m*K))
- `c`: Specific Heat (J/(kg*K))
- `p`: Density (kg/m^3 )
- CET
  - `N0`
  - `n`
  - `a`

#### Beam File
This file contains information on the energy source, which is a volumetric gaussian. It should be noted that everything input here is equivalent to `√6 σ` as defined by the `D4σ` beam diameter. For example, if the standard deviation of an actual beam is `20 µm` and radially symmetric, both `Width_X` and `Width_Y` should be set to 48.9898e-6.

Shape
- `Width_X`: Width of the beam in the X direction (m)
- `Width_Y`: Width of the beam in the Y direction (m)
- `Depth_Z`: Penetration depth of the beam (m)

Intensity
- `Power`: Power of energy source (W)
- `Efficiency`: Absorption efficiency of beam; refer to literature for accurate values. Typically it is 0.35 for laser powder bed fusion (PBF) and 0.85 for electron beam PBF.

#### Path File
This file dictate where and how the heat source travels. This file and all its components are necessary. An incomplete path file will cause failures. The format is different than other files and is in the form:
`Mode	X(mm)	Y(mm)	Z(mm)	Pmod	Vel(m/s)/Time(s)`

The Mode dictates how the heat source moves.
- `Mode` = 0 is a line melt. X, Y, Z controls where the beam travels TO and the last parameter controls the constant speed of this melt.
- `Mode` = 1 is a spot melt. X, Y, Z control the location of the spot melt and the last parameter controls the duration of this melt.
	
- `Pmod` is a power multiplier to the heat source. This is typically 0 or 1 to control the beam turning on and off. More advanced scan strategies have variable powers. In these cases, it’s generally best to let `Pmod` control the power and set the `Power` in the beam file to 1.

Note: Make sure that before a line melt, the starting point is correct! A raster pattern is best represented by alternating mode 1 and 0 where mode 1 is only there to set the start point and would have a Pmod of 0 and a short value for Time(s).

### Option Files
This set of files dictate everything having to do with the numerics of the simulation, such as the domain considered as well as various settings. These files are optional but it is highly recommended to understand how to change these for specific uses as they have a large impact on the speed and accuracy of the simulation. Default behavior is specified.

#### Domain File
This file contains information about the domain over which to calculate the temperature solution. It contains ways to control the resolution and bounds of the simulation domain but highly recommended to use. Boundary conditions can be set to provide somewhat nearly insulative boundaries using the method of images (just 1 iteration). For very thins walls (<1mm), this may not be enough to provide a completely insulative effect. Keep in mind that each point is calculated independently, so increasing the resolution by a factor of 2 in each direction will slow down the simulation by a factor of 8.

- `X`
  - `Min`: Minimum X value
  - `Max`: Maximum X value
  - `Res`: Resolution in the X-direction of the domain (either Res OR Num can be used)
  - `Num`: Number of unique X-values in the grid. If this is set to 1, only the Max value is used. (either Res OR Num can be used)
- `Y`
	`Min`
	`Max`
	`Res`
	`Num`
- `Z`
	`Min`
	`Max`: Typically set to 0 (where the top surface is)
	`Res`
	`Num`
- `BoundaryConditions`
  - `X_min`
  - `X_max`
  - `Y_min`
  - `Y_max`
- `Custom`
  - `File`: File which has a specific set (x,y,z) coordinates to be used for the domain Note: Using a point file negates any tracking modes (`Volume` and `Surface`), thus all points will always calculated.

Note: Not specifying Min or Max value results in the unspecified value being calculed via the path file(s) with a domain file 500µm buffer on each side and a 250µm simulated depth. Not specifying a Res or Num results in a default resolution of 50 µm.

#### Output File

This file contains all variables which can be output. A value of 0 indicated to not output the variable whereas a value of 1 indicates that variable should be output. Most variables default to a value of 0. The memory required to run a simulation increases when more outputs are selected; therefore, it is good practice to output only the necessary or desired information.

- Grid
  - `x`: x-coordinate
  - `y`: y-coordinate
  - `z`: z-coordinate
- Temperature
  - `T`: Temperature
  - `T_hist`: Temperature history. Note: this will create a separate file for each point. DO NOT USE if many points are in the domain
- Solidification
  - `tSol`: Solidification time
  - `G`: Magnitude of the thermal gradient at solidification
  - `Gx`: x-component of normalized gradient
  - `Gy`: y-component of normalized gradient
  - `Gz`: z-component of normalized gradient
  - `V`: Velocity of solidification front
  - `dTdt`: Cooling rate
  - `eqFrac`: Equiaxed fraction
  - `RDF`: Export results in “Reduced Data Format” compatible with ExaCA
  - `numMelt`: Number of times a point melted and solidified
- Solidification+
  - `H`: Magnitude of the orthogonal differential change in the solidification gradient in the direction of the solidification gradient
  - `Hx`: x-component of normalized H
  - `Hy`: y-component of normalized H
  - `Hz`: z-component of normalized H

#### Settings File
This file contains all the tunable parameters for how the simulation is run. Generally, only `MaxThreads` (which speeds up the simulation if the hardware being used has at least that many available threads) needs to be changed.

- Temperature
  - `Sol_Tol`: Controls how close the temperature must be to the liquidus temperature (as a fraction) before the search algorithm stops. As an example, if `T_L = 10K`, `timestep=0.1s`, `T(t=1.7s)=12K`, `T(t=1.8s)=9K`, and `Sol_Tol=1e-2` THEN the search algorithm won’t stop trying to find the EXACT time the point solidified until it finds a time where `9.9K < T < 10.1K`. Defaults to `1e-3`
  - `Sol_Iter`: Controls the maximum number of iterations for the algorithm used above. Defaults to 10.
  - `Cutoff_Peak`: Helps control how far back in time the integration is done. Defaults to 1e-9. For example, if a simulation runs for 30s, maybe the first 5s have so little temperature contribution that they aren’t worth integrating. This setting rarely gets used since the introduction of the path compression algorithm.
  - `Cutoff_T0TL`: Another setting to help control how far back in time the integration is done. Defaults to 1e-9.
- Path
  - `Buffer`: Control how far outside the domain the scan path is considered (in meters). For example, if a scan file contains information for 6 cubes, but you only want to analyze one of them. Set the domain to include the 1 cube and set the buffer to 0 or 0.001 (1mm) or something small. This will result in a large speedup.
  - `Compression`: Binary toggle (0-off, 1-on) for using a path compression algorithm. This is primary useful for LARGE scan paths (1 Mb or larger) with smooth-ish movement of the heat source (ex: raster or point raster) rather than discontinuous movement (ex: random point fill). It does so by compressing multiple nearby diffuse heat sources into a single source. There is a small loss in accuracy but a large speedup for certain simulations.
- Compute
  - `MaxThreads`: Number of threads to use for the simulation.

## 4. OutputFiles
Analysis of 3dThesis results can be done with a variety of methods (e.g. Python), but quick visualization is possible with Paraview.

First ensure that the run finished successfully for `TestInputs/ParamInput.txt` the data should be found in `TestSim/Data/TestSim.Final.csv`. To view this file in Paraview, open it through the folder icon, in the top left, and then click the apply button. A table of data should appear on the right side of the screen. To convert this data to something visual, click on the file name on the right side of the screen to highlight it, then apply the `Table to Points` filter under `Filters->Alphabetical`.

On the left side, three dropdown menus should appear titled “X Column”, “Y Column”, and “Z Column”. Simply change these to be `x`, `y`, and `z`, click apply, and exit out of the tabular view of the data. Now at the top of the screen, locate the dropdown menu that says, “Solid Color.” Change this to `G` (the thermal solidification gradient). The point size and color scale can be changed using other options in the left menu (under the “Coloring” and “Styling” groups respectively).
