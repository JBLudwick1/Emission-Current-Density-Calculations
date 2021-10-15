HOW TO GENERATE A TABLE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

1)	Begin with an existing numerical function to calculate a desired quantity (such as local emission current density) and store the value into a lookup table.
	This function must be able to allow each input parameter in 1D array form and have an N-dimensional matrix output where N is the number of input parameters (input parameters kept constant do not count towards N)
	a.	Determine the number of input parameters for this function.
	b.	Determine the range of each input parameter for this function.
	c.	Determine a maximum acceptable error when using a lookup table to approximate this function.

2)	In the matlab file “Main.m”:
	a.	Input the lower and upper bounds for the ranges of each input parameter in the LBs and UBs arrays, respectively:
		i.	Ensure the parameters are listed in the same order for both the LBs and UBs arrays and follow the input parameter order for the numerical function.
		ii.	If a parameter is to be held constant, enter the same constant value for the lower and upper bound.
			These constants will not be stored as lookup table input parameters, and will reduce the dimensionality of the lookup table (expediting build time).

		EXAMPLE
		LBs = [0.5,300,20,2,5]; %Lower Bounds
		UBs = [10,3000,5000,6.5,5]; %Upper Bounds

		This above example shows input parameter range bounds for 5 parameters in determining local emission current density. 
		In order from left to right, these arrays contain quantities for electric field (V/nm), emitter temperature (K), radius of curvature (nm),
		emitter work function (eV), and Fermi energy (eV). Fermi energy is held constant at 5 eV in this example.

	b.	Determine which input parameters are to be stored in special funciton space
		i.	These cell arrays must be the same size and parameter order as UBs/LBs
		ii.	a zero value uses parameter in standard linear space, an anonymous function input uses function in special parameter space
		iii.	If a special anonymous function is used in the variable 'special_funcs', the corresponding inverse of that function must be supplied in 'special_funcs_inv'
		iv.	These anonymous functions must be able to handle 1D array inputs.
		v.	It is useful for certain parameters to be stored in special funciton space, such as emitter radius of curvature when calculating emission current.
			This parameter varies greatly at small values and varies extremely small amounts at large values, thus it is better conducted in log10 space.

		EXAMPLE
		special_funcs = {0,0,@(x)log10(x),0,0};
		special_funcs_inv = {0,0,@(x)10.^x,0,0};
		
		This above example conducts parameters 3 in log10 space, and all other parameters in standard linear space. Note that the inverse of log10(x) (10^x) is stored in the inverse function variable

	c.	Set the maximum acceptable error (in percent error) relative to the numerical function when approximating a function with the lookup table in the 'error_threshold' variable

	d.	Determine if table should be built in natural logspace. Functions where linear interpolation is easier in natural logspace (such as local emission current density functions)
		should be marked true here, as this will greatly expedite lookup table build time. Note, the resulting finished table will store quantities in natural logspace if value is marked true.
		When this is the case, when using the table, ensure linear interpolation is conducted first, then the interpolated lookup table value must be converted out of natural logspace
		i.e. true value = exp(interpolated quantity).

	e.	Set the output full path to store the finished lookup table in the variable 'output_full_path'. This output will be saved each lookup table build iteration in case of runtime failure.

	f.	Set the numerical function to have a lookup table approximate following the specifications in step 1.
		i.	The function name must be inserted as string in the variable 'num_fun_str'
		ii.	It is recommended this function is in the same folder as these other MATLAB functions

3)	Run “Main.m” and wait for lookup table to finish building. This process may take multiple hours or more depending on the function being calculated, the quantity and range of input parameters,
	and the maximum acceptable error of the lookup table. Too stringent of specifications may induce memory errors.
	a.	As the lookup table is building, the result of each iteration will print out in the MATLAB command window:
		Iteration: 1
		Approximated Maximum Error: 200.00%
		Parameter: 1, size: 5, max error: 200.00%
		Parameter: 2, size: 3, max error: 200.00%
		Parameter: 3, size: 3, max error: 171.85%
		Parameter: 4, size: 3, max error: 174.57%
		This Iteration Run Time (min): 0.12
		Total Run Time (min): 0.12
		
		Iteration: 2
		Approximated Maximum Error: 200.00%
		Parameter: 1, size: 9, max error: 200.00%
		Parameter: 2, size: 3, max error: 179.82%
		Parameter: 3, size: 3, max error: 109.11%
		Parameter: 4, size: 3, max error: 81.39%
		This Iteration Run Time (min): 0.08
		Total Run Time (min): 0.20
		
		Iteration: 3
		Approximated Maximum Error: 200.00%
		Parameter: 1, size: 17, max error: 200.00%
		Parameter: 2, size: 3, max error: 198.16%
		Parameter: 3, size: 3, max error: 191.72%
		Parameter: 4, size: 3, max error: 128.45%
		This Iteration Run Time (min): 0.13
		Total Run Time (min): 0.33

	b.	In this above example, the approximated errors from each parameter, along with the number of quantities for each parameter is shown, as well as a total lookup table approximated error.
		Due to the relative error equation, some errros above 200% are falsely listed as a 200% error. In each iteration, the input parameter with the highest associated error has additional stored calculations 
		added to it in areas of high error, ultimately lowering the net lookup table error each iteration. In tables with many input parameters, such as the above example,
		it may take a few iterations before the total lookup table error drops below 200%. A lookup table .mat file is saved after each iteration.

HOW TO USE A GENERATED TABLE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Once a table is built, a .mat file is save at a specified location. That .mat file contains the variables:
	•	“error_threshold”
		o	Maximum acceptable error in lookup table approximation
	•	“max_errors”
		o	Approximate max error caused by each input parameter
	•	“dvars_linspace”
		o	A cell array, where each cell contains a 1D array for each input parameter. These arrays align with the lookup table quantities.
		o	If an input parameter used a special function, it is converted back to standard linear space when stored here.
	•	“table”
		o	N-Dimensional (where N is the number of input parameters) table of stored numerical quantities from desired numerical function.
	•	“total_time”
		o	Total time used to generate lookup table. 

To use the lookup table to approximate a queried quantity, it is recommended to use the “linear_interp.m” MATLAB function. See additional description of function in “Brief Description of MATLAB files”.
After loading the .mat file, an example of interpolating a quantity is shown below:

J = exp(linear_interp.interpn(dvars_linspace{:},table,F_queries,T_queries,R_queries,phi_queries));

In this example, J is local emission current density stored in a 4-dimensional lookup table. The input query parameters used to approximate J are F, T, R, and phi representing electric field magnitudes,
temperatures, radii of curvature, and work functions, respectively. “table” was stored in natural logspace, so “exp” is needed to convert the table to A/cm2.

It is recommended when using this interpolation call, to only load the table once if possible. The interpolation command is very fast, but loading the table many times can create unnecessary overhead.

BRIEF DESCRIPTION OF MATLAB FILES----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
•	Main.m
	o	Main file to generate a lookup table. All user inputs are contained in this file.
•	LUT_gen.m
	o	Function to build lookup table
•	NUMERICAL_FUNC.m
	o	Example numerical function used to calculate local emission current density values with the quantum mechanical wave impedance approach. Output units are in A/cm^2. For more details, see journal paper at https://doi.org/10.1063/5.0065612
•	linear_interp.m
	o	Function used to rapidly linearly interpolate values from a table. This code was modified from: Jeffrey Wu (2021). Faster linear Interpolation (https://www.mathworks.com/matlabcentral/fileexchange/28376-faster-linear-interpolation), MATLAB Central File Exchange.
	o	Input 1-dimensional arrays X1, X2, ... , Xn, a n-dimensional array V, and x1, x2, ..., xn. Query values assume Xi's are in increasing order. Can query multiple points, but x1, x2, ..., xn must be the same size
