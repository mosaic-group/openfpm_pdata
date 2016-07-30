/*
 * ### WIKI 1 ###
 *
 * ## Simple example for graph plotting with google charts
 *
 * In this example we show how to plot several different graphs type
 *
 * ### WIKI END ###
 *
 */

#include "VCluster.hpp"
#include "Plot/GoogleChart.hpp"

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library and we define Ghost size
	// and non-periodic boundary conditions
	//
	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	// This test is not made to run in parallel
	if (v_cl.getProcessingUnits() > 1)
	{
		std::cerr << "Error: only one processor is allowed" << "\n";
		return 1;
	}

	//
	// ### WIKI 3 ###
	//
	// Here we define vectors that contain information about the x value and y value
	//

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<size_t>> y;
	openfpm::vector<std::string> yn;

	//
	// ### WIKI 4 ###
	//
	// The x values can contain numbers but also strings
	//

	x.add("colum1");
	x.add("colum2");
	x.add("colum3");
	x.add("colum4");
	x.add("colum5");
	x.add("colum6");

	//
	// ### WIKI 5 ###
	//
	// Each x value can have different values (or dataset, lines) in this case 4.
	// This mean that for each x value we have to define 4 y values
	//
	//

	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Each x has 4 values
	y.add({2,3,5,6});
	y.add({5,6,1,6});
	y.add({2,1,6,9});
	y.add({1,6,3,2});
	y.add({3,3,0,6});
	y.add({2,1,4,6});

	//
	// ### WIKI 6 ###
	//
	// We can specify several options for the graphs.
	//
	// * Title of the graph
	// * Title of the y axis
	// * Title of the x axis
	// * stype is the type of graph. Possible values are "bars" and "lines"
	// * stypeext is the extension property, in this case we say that the third value
	//            in y must be reppresented as a line instead of a bar
	//
	// * more In this section can be used non default options from Google Charts
	//
	// For more options refer to doxygen and Google Charts
	//

	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.stype = std::string("bars");
	options.lineWidth = 5;

	// it say that the colum4 must me represented with a line
	options.stypeext = std::string("{3: {type: 'line'}}");

	//
	// ### WIKI 7 ###
	//
	// We create the object to create plots with Google Charts
	//
	// A writer can produce several graphs interleaved with HTML code
	//
	// Hist Graph is a graph with histograms for each x values
	// Lines Graph is the typical graph with lines
	//

	GoogleChart cg;
	//
	cg.addHTML("<h2>First graph</h2>");
	cg.AddHistGraph(x,y,yn,options);
	cg.addHTML("<h2>Second graph</h2>");
	cg.AddLinesGraph(x,y,yn,options);
	cg.write("gc_out.html");

	//
	// ### WIKI 8 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
