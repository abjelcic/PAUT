# PAUT_AxialWedge

* <code>data</code> folder contains file <code>data.in</code> with wedge/specimen/inspection parameter written in format given by <code>format.in</code> file.

* Upon execution, <code>out</code> folder will contain .law file.

* How to compile: see <code>Makefile</code>

# Remarks

* Wedge geometry parameters follow the usual convention. Nevertheless, <code>doc/WedgeGeometry.pdf</code> defines the wedge geometry uniquely. The only new
parameter introduced is the "Wedge radius".

* Notice that parameter: "Height at the middle of the first element" might be defined in many different ways. The one selected is depicted in <code>doc/WedgeGeometry.pdf</code>.

* Only dual linear-array probes are supported for now, i.e. matrix probes are not supported although one can easily generalize the code to matrix probes.

* For a given sector scan inspection parameters ("Refracted angle start", "Refracted angle stop", "Refracted angle resolution"), the code for every 
refracted angle finds the optimal focal depth. That corresponds to "Auto-focusing" option in UltraVision .cal file. Therefore, the user should not try
using the flat wedge (Wedge radius = +Inf) with zero roof angle, otherwise the natural focal depth will be infinite.


abjelcic@phy.hr
