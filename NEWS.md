# _slendr_ 0.1 (upcoming first public release)

## New functionality

- Geographic landscapes and coordinate systems are no longer tied to a geographic location on Earth and can be fully abstract (including completely arbitrary shapes of the "continents", geographic barriers etc.).

- New set operations for manipulating spatial regions: `join()`, `overlap()`, `subtract()`. Users can use these do define more complex larger population boundaries and geographic regions from smaller ones.

- New [shiny](http://shiny.rstudio.com)-based interactive exploration "app": `explore()`. The user can provide a compiled model object (i.e. formal configuration of spatial population dynamics and gene flow events) and the package will spawn a browser app with a time slider. This can be used to inspect the model spatial dynamics interactively before they are sent over to the SLiM side.

- A new function `distance()` for computing distances between geographic objects in the package (either between their borders or centroids) has been implemented. A small utility function `dimension()` which calculates the dimension of the world map in the "real" projected units (i.e. meters) has been also added.

- It is now possible to define models in forward or backwards direction, depending on what is more convenient for the user and the scenario that is being modeled. The direction of time is automatically detected and translated to SLiM's units of generations in the forward direction. The implementation of this features might be a bit wonky in some aspects and more testing is required.

- Simulations of marine species are now possible (this required adding a single argument `marine` to the `population()` call which instructs the population range rendering procedure to "invert' intersection with landscape features).

- Any simulation can now have multiple "ancestral" populations (i.e. populations without an immediate ancestor created by `sim.addSubpop()`). Until now, all populations had to trace their ancestry to a single ancestor (a leftover of a hard-coded requirement from the very first version of my code).

- A function `resize()` for scheduling changes in population size (in a single step, or as an exponential growth/shrinking) is implemented.

- A new `shrink()` function (analogous to `expand()`) is implemented.

## Changes to the R interface

- Renamed `admixture()` to `geneflow()`. No ugly historical baggage and a more general name for this concept.

- The number of intermediate spatial snapshots used in `move()` and `expand()` is now determined iteratively instead of forcing the user to specify this by themselves (which is something they would have to do manually anyway and was very annoying to deal with). If the user specifies `snapshots`, the search is not performed and the specified value is used to generate intermediate spatial maps.

- Spatial interaction distances (translated to `maxDistance` on the SLiM side) as well as offspring distances from parents (i.e. the standard deviations of their normal distributions) can now be specified for each population individually. The users can still provide default values for all populations (or just those which did not have their own dedicated parameter values) in the main `compile()` call.

## Changes to the SLiM backend

- Significant overhaul of the entire [SLiM layer](https://github.com/bodkan/slendr/blob/main/inst/extdata/backend.slim). Most importantly, complex hacks approximating dataframe behavior using SLiM matrices and other were completely replaced by Dictionaries. The intention of the SLiM backend is to be as easy to understand as possible and as transparent as possible for the users to be able to easily build on top of it. If something could be written more elegantly or in a simpler way, please file a GitHub issue as I considered this a bug.

## Other changes

- More than a hundred new unit tests checking the consistency of the R interface (misspecified event times, etc.) as well as the correctness of the SLiM simulation runs programmed by _slendr_.
