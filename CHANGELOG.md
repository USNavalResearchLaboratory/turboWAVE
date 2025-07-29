# Changelog

## [Unreleased]

This is the largest rewrite turboWAVE has ever undergone.
The two biggest changes are C++20 modules, and handling of field topology.

### Physics Content

* Strong field QED processes can be included
    - Quantum synchrotron
    - Nonlinear Breit-Wheeler

### User Interface

* The input file parser is far more sophisticated
    - More strict syntax, but...
    - ...also better error messages
    - Faster: no more reliance on the comically slow std::regex

### Physics Structure

* Topology and metric classes have been further broken down
    - `StaticSpace` represents the full topology and can be counted on to stay constant during a simulation
    - `DynSpace` includes elements related to an evolving window
    - `MetricSpace` includes the metric data
* The `Field` class derives from `StaticSpace` so there is no danger of its underlying data becoming invalidated.
* The `Field` class no longer changes the interpretation of inherited data
* The `Field` class and supporting classes are more relativistically covariant
    - Every `StaticSpace` has 5 dimensions with the universally exposed ordering (t,x,y,z,c), where c is an internal dimension.  The addition of the time dimension allows us to more easily handle scenarios where a few prior time levels are of interest to retain, or where time needs to be treated on an equal footing.
* Integral transforms can be applied more flexibly
* The particle handling layer is easier to follow thanks to the new C++20 module structure
* The concept of an evolving simulation window ("moving window") is generalized to 4 dimensions
    - Ordinarily the time axis has 1 cell that "moves" by the time-step each time-level
    - More generally, the time-level can be viewed as parameterizing arbitrary spacetime structures
    - This is reflected in 4-dimensional tuples that appear in the input file

### Code Structure

* C++20 modules are used throughout
    - This pushes the envelope of where compilers and build tools stand as of this writing
    - Header files are only used for unavoidable macros or external dependencies
    - The benefits are many, not the least of which is the ability of the compiler to deter circular dependencies
* New iterator layer to handle processing cells and strips in a 5-dimensional world
* Tree-sitter parsing library is embedded as source
* `Slice` objects are updated to be 5-dimensional
* Eliminated `AutoField` template
    - `ScalarField` and `ComplexField` derive directly from `Field` and live in the `fields:aggregates` partition module
    - Aggregated field components just work (no more caveats about storage patterns or stride units)

### About Inclusive and Exclusive Bounds

* Meaning of bounds has changed in some places
* Lower bounds are always inclusive, for upper bounds it depends
    - The notation `beg` and `end` means `[beg,end)` (upper bound excluded)
    - Everything else is inclusive, e.g. `[start,stop]`, `[first,last]`, `[LFG,UFG]`, `[lb,ub]`

### Modernization

* std::array is used for index nodes, e.g. `typedef std::array<tw::Int,4> node4`
* smart pointer usage