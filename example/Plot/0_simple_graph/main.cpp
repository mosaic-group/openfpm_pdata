 /*! \page Plot Plot
 *
 * \subpage Plot_0_cg
 *
 */

/*!
 * \page Plot_0_cg Plot 0 Google Chart
 *
 * # Simple example for plotting 2D graph with google charts # {#e0_pl}
 *
 * In this example we show how to plot several different 2D graphs type with Google
 *
 * ## Inclusion ## {#e0_pl_inc}
 *
 * To use the Plot features we have to include GoogleChart
 *
 * \snippet Plot/0_simple_graph/main.cpp include
 *
 */

//! \cond [include] \endcond

#include "Plot/GoogleChart.hpp"
#include "VCluster.hpp"

//! \cond [include] \endcond

int main(int argc, char* argv[])
{
	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## Initialization ##
	 *
	 * Here we Initialize the library, and we check the we are on a single processor. GoogleChart
	 * cannot do parallel IO or write big files. So or we collect all data on one processor, or each
	 * processor write a distinct file. In this particular example we simply stop if the program  start
	 * on more than one processor
	 *
	 * \snippet Plot/0_simple_graph/main.cpp initialize
	 *
	 */

	//! \cond [initialize] \endcond

	openfpm_init(&argc,&argv);
	auto & v_cl = create_vcluster();

	// Google chart is only single processor
	if (v_cl.getProcessingUnits() > 1)
	{
		std::cerr << "Error: only one processor is allowed" << "\n";
		return 1;
	}

	//! \cond [initialize] \endcond

	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## Graph data ##
	 *
	 * Here we have the vectors that will contain the information about the graph.
	 *
	 * \snippet Plot/0_simple_graph/main.cpp datas vector
	 *
	 */

	//! \cond [datas vector] \endcond

	openfpm::vector<std::string> x;
	openfpm::vector<openfpm::vector<double>> y;
	openfpm::vector<std::string> yn;

	//! \cond [datas vector] \endcond

	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * We will try now to produce the following situation. Six values on **x** each of them having 4 values on **y**
	 *
	 * This mean that for each x value we have to define 4 y values. Having multiple values on on x can be used for
	 * several purpose.
	 *
	 * * Define multiple lines. For example if we connect all the points # we obtain one line. If we connect
	 *   all the @ points we obtain another line, an so on ... (figure below)
	 *
	 * * Define error bands
	 *
	 * * Visualize different observables/parameters for the same value x
	 *
	 *
	 * \verbatim

		  y  ^                          $ dataset1
		     |                          * dataset2
		 0.9 |                          # dataset3
		     |       @                  @ dataset4
		     |   #
		 0.6 |       *   *   @       *
		     |   $   #   @   #   #
		     |   @       $   $   @   @
		 0.3 |           #   *   $   #
		     |       $           *
		     |   *                   $
		  0  |_________________________________
				 o   t   t   f   f   s          x
				 n   w   h   o   i   i
				 e   o   r   u   v   x
						 e   r   e
						 e



	  \endverbatim
	 *
	 * We start from the first case (Define multiple lines)
	 *
	 * \snippet Plot/0_simple_graph/main.cpp data fill
	 *
	 */

	//! \cond [data fill] \endcond

	// Fill the x values
	x.add("one");
	x.add("two");
	x.add("three");
	x.add("four");
	x.add("five");
	x.add("six");

	// we have 4 dataset  or lines
	yn.add("dataset1");
	yn.add("dataset2");
	yn.add("dataset3");
	yn.add("dataset4");

	// Because we have 6 points on x each containing 4 lines or dataset, we have to provides
	// 6 point with 4 values at each x point
	y.add({2,3,5,6});
	y.add({5,6,1,6});
	y.add({2,1,6,9});
	y.add({1,6,3,2});
	y.add({3,3,0,6});
	y.add({2,1,4,6});

	//! \cond [data fill] \endcond

	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## Graph options ##
	 *
	 * We can specify several options for the graphs.
	 *
	 * * Title of the graph
	 * * Title of the y axis
	 * * Title of the x axis
	 *
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart
	 *
	 */

	//! \cond [google chart] \endcond

	GCoptions options;

	options.title = std::string("Example");
	options.yAxis = std::string("Y Axis");
	options.xAxis = std::string("X Axis");
	options.lineWidth = 5;

	//! \cond [google chart] \endcond

	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## Graph write ##
	 *
	 * We create the object to create plots with Google Charts
	 *
	 * A writer can produce several graphs optionally interleaved with HTML code.
	 * Here we write in HTML a description of the graph, than we output the graph
	 *
	 * AddLinesGraph create a typical graph with lines
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart write1
	 *
	 * \htmlonly
		<div id="chart_div0" style="width: 900px; height: 500px;"></div>
	   \endhtmlonly
	 *
	 */

	//! \cond [google chart write1] \endcond

	GoogleChart cg;
	//
	cg.addHTML("<h2>First graph</h2>");
	cg.AddLinesGraph(x,y,yn,options);

	//! \cond [google chart write1] \endcond


	/*!
	 *
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## Hist graph ##
	 *
	 * Hist graph is instead a more flexible Graph writer. In particular we can specify
	 * how to draw each dataset. With the option
	 *
	 * * **stype** specify how to draw each dataset
	 * * **stypeext** we can override the default stype option. In this case we say that the third dataset
	 *            in must be reppresented as a line instead of a bars
	 *
	 * To note that we can reuse the same Google chart writer to write multiple
	 * Graph on the same page, interleaved with HTML code
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart write2
	 *
	 * \htmlonly
		<div id="chart_div1" style="width: 900px; height: 500px;"></div>
	   \endhtmlonly
	 *
	 *
	 */

	//! \cond [google chart write2] \endcond

	options.stype = std::string("bars");

	// it say that the dataset4 must me represented with a line
	options.stypeext = std::string("{3: {type: 'line'}}");

	cg.addHTML("<h2>Second graph</h2>");
	cg.AddHistGraph(x,y,yn,options);

	//! \cond [google chart write2] \endcond

	/*!
	 *
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## %Error bars ##
	 *
	 * Here we show how to draw error bars. %Error bars are drawn specifying intervals with a min and a max.
	 * Intervals in general does not have to encapsulate any curve. First we construct the vector y with 3
	 *  values the first value contain the curve points, the second and third contain the min,max interval.
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart write3
	 *
	 * \htmlonly
		<div id="chart_div2" style="width: 900px; height: 500px;"></div>
	   \endhtmlonly
	 *
	 *
	 */

	//! \cond [google chart write3] \endcond

	cg.addHTML("<h2>Third graph</h2>");

	// The first colum are the values of a line while the other 2 values
	// are the min and max of an interval, as we can see interval does not
	// have to encapsulate any curve
	y.clear();
	y.add({0.10,0.20,0.19});
	y.add({0.11,0.21,0.18});
	y.add({0.12,0.22,0.21});
	y.add({0.15,0.25,0.20});
	y.add({0.09,0.29,0.25});
	y.add({0.08,0.28,0.27});

	// Here we mark that the the colum 2 and 3 are intervals
	yn.clear();
	yn.add("line1");
	yn.add("interval");
	yn.add("interval");

	cg.AddLinesGraph(x,y,yn,options);

	//! \cond [google chart write3] \endcond

	/*!
	 *
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * The style of each interval can be controlled, and the definition of intervals can be interleaved with definition of
	 * other lines. In this example we show how to define 3 lines and 3 intervals, controlling the style of the last interval
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart write4
	 *
	 * \htmlonly
		<div id="chart_div3" style="width: 900px; height: 500px;"></div>
	  \endhtmlonly
	 *
	 *
	 */

	//! \cond [google chart write4] \endcond

	cg.addHTML("<h2>Four graph</h2>");

	// again 6 point but 9 values
	y.clear();
	y.add({0.10,0.20,0.19,0.22,0.195,0.215,0.35,0.34,0.36});
	y.add({0.11,0.21,0.18,0.22,0.19,0.215,0.36,0.35,0.37});
	y.add({0.12,0.22,0.21,0.23,0.215,0.225,0.35,0.34,0.36});
	y.add({0.15,0.25,0.20,0.26,0.22,0.255,0.36,0.35,0.37});
	y.add({0.09,0.29,0.25,0.30,0.26,0.295,0.35,0.34,0.36});
	y.add({0.08,0.28,0.27,0.29,0.275,0.285,0.36,0.35,0.37});

	// colum  0 and 1 are lines
	// colums 2-3 and 4-5 are intervals
	// colum 6 is a line
	// colum 7-8 is an interval
	yn.add("line1");
	yn.add("line2");
	yn.add("interval");
	yn.add("interval");
	yn.add("interval");
	yn.add("interval");
	yn.add("line3");
	yn.add("interval");
	yn.add("interval");

	// Intervals are enumerated with iX, for example in this case with 3 intervals we have i0,i1,i2
	// with this line we control the style of the intervals. In particular we change from the default
	// values
	options.intervalext = std::string("{'i2': { 'color': '#4374E0', 'style':'bars', 'lineWidth':4, 'fillOpacity':1 } }");

	cg.AddLinesGraph(x,y,yn,options);

	//! \cond [google chart write4] \endcond

	/*!
	 *
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 * ## More options ##
	 *
	 * In this last example we also show how to:
	 *
	 *
	 * * Make the graph bigger, setting **width** and **height** options
	 * * Give the possibility to to zoom-in and zoom-out with **GC_EXPLORER**
	 * * Use lines instead a smooth function to connect points
	 * * Use logaritmic scale
	 *
	 * \note For more options refer to doxygen and Google Charts
	 *
	 * \snippet Plot/0_simple_graph/main.cpp google chart write5
	 *
	 *
	 * \htmlonly
		<div id="chart_div4" style="width: 1280px; height: 700px;"></div>
	  \endhtmlonly
	 *
	 */

	//! \cond [google chart write5] \endcond

	openfpm::vector<double> xn;

	xn.add(1.0);
	xn.add(2.0);
	xn.add(3.0);
	xn.add(4.0);
	xn.add(5.0);
	xn.add(6.0);

	options.intervalext = "";
	options.width = 1280;
	options.heigh = 720;
	options.curveType = "line";
	options.more = GC_ZOOM + "," + GC_X_LOG + "," + GC_Y_LOG;

	cg.AddLinesGraph(xn,y,yn,options);

	cg.write("gc_out.html");

	//! \cond [google chart write5] \endcond


	/*!
	 * \page Plot_0_cg Plot 0 Google Chart
	 *
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Plot/0_simple_graph/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();


	//! \cond [finalize] \endcond
}

////////////// FOR DOXYGEN DOCUMENTATION ///////////////////////////

/*!
 * \page Plot_0_cg Plot 0 Google Chart
 *
 * \htmlonly
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {'packages':['corechart']});
      google.charts.setOnLoadCallback(drawVisualization);


      function drawVisualization() {
var data0 = new google.visualization.DataTable();
data0.addColumn('string','X Axis');
data0.addColumn('number','dataset1');
data0.addColumn('number','dataset2');
data0.addColumn('number','dataset3');
data0.addColumn('number','dataset4');
data0.addRows([
['one',2,3,5,6],
['two',5,6,1,6],
['three',2,1,6,9],
['four',1,6,3,2],
['five',3,3,0,6],
['six',2,1,4,6],
]);
var data1 = new google.visualization.DataTable();
data1.addColumn('string','X Axis');
data1.addColumn('number','dataset1');
data1.addColumn('number','dataset2');
data1.addColumn('number','dataset3');
data1.addColumn('number','dataset4');
data1.addRows([
['one',2,3,5,6],
['two',5,6,1,6],
['three',2,1,6,9],
['four',1,6,3,2],
['five',3,3,0,6],
['six',2,1,4,6],
]);
var data2 = new google.visualization.DataTable();
data2.addColumn('string','X Axis');
data2.addColumn('number','line1');
data2.addColumn({id:'i0', type:'number', role:'interval'});
data2.addColumn({id:'i0', type:'number', role:'interval'});
data2.addRows([
['one',0.1,0.2,0.19],
['two',0.11,0.21,0.18],
['three',0.12,0.22,0.21],
['four',0.15,0.25,0.2],
['five',0.09,0.29,0.25],
['six',0.08,0.28,0.27],
]);
var data3 = new google.visualization.DataTable();
data3.addColumn('string','X Axis');
data3.addColumn('number','line1');
data3.addColumn({id:'i0', type:'number', role:'interval'});
data3.addColumn({id:'i0', type:'number', role:'interval'});
data3.addColumn('number','line1');
data3.addColumn('number','line2');
data3.addColumn({id:'i1', type:'number', role:'interval'});
data3.addColumn({id:'i1', type:'number', role:'interval'});
data3.addColumn({id:'i2', type:'number', role:'interval'});
data3.addColumn({id:'i2', type:'number', role:'interval'});
data3.addRows([
['one',0.1,0.2,0.19,0.22,0.195,0.215,0.35,0.34,0.36],
['two',0.11,0.21,0.18,0.22,0.19,0.215,0.36,0.35,0.37],
['three',0.12,0.22,0.21,0.23,0.215,0.225,0.35,0.34,0.36],
['four',0.15,0.25,0.2,0.26,0.22,0.255,0.36,0.35,0.37],
['five',0.09,0.29,0.25,0.3,0.26,0.295,0.35,0.34,0.36],
['six',0.08,0.28,0.27,0.29,0.275,0.285,0.36,0.35,0.37],
]);
var data4 = new google.visualization.DataTable();
data4.addColumn('number','X Axis');
data4.addColumn('number','line1');
data4.addColumn({id:'i0', type:'number', role:'interval'});
data4.addColumn({id:'i0', type:'number', role:'interval'});
data4.addColumn('number','line1');
data4.addColumn('number','line2');
data4.addColumn({id:'i1', type:'number', role:'interval'});
data4.addColumn({id:'i1', type:'number', role:'interval'});
data4.addColumn({id:'i2', type:'number', role:'interval'});
data4.addColumn({id:'i2', type:'number', role:'interval'});
data4.addRows([
[1,0.1,0.2,0.19,0.22,0.195,0.215,0.35,0.34,0.36],
[2,0.11,0.21,0.18,0.22,0.19,0.215,0.36,0.35,0.37],
[3,0.12,0.22,0.21,0.23,0.215,0.225,0.35,0.34,0.36],
[4,0.15,0.25,0.2,0.26,0.22,0.255,0.36,0.35,0.37],
[5,0.09,0.29,0.25,0.3,0.26,0.295,0.35,0.34,0.36],
[6,0.08,0.28,0.27,0.29,0.275,0.285,0.36,0.35,0.37],
]);
var options0= {
title : 'Example',
vAxis: {title: 'Y Axis'},
hAxis: {title: 'X Axis'},
curveType: 'function',
lineWidth: 5,
intervals: { 'style':'area' },
};
var options1= {
title : 'Example',
vAxis: {title: 'Y Axis'},
hAxis: {title: 'X Axis'},
seriesType: 'bars',
series: {3: {type: 'line'}},
};
var options2= {
title : 'Example',
vAxis: {title: 'Y Axis'},
hAxis: {title: 'X Axis'},
curveType: 'function',
lineWidth: 5,
intervals: { 'style':'area' },
};
var options3= {
title : 'Example',
vAxis: {title: 'Y Axis'},
hAxis: {title: 'X Axis'},
curveType: 'function',
lineWidth: 5,
intervals: { 'style':'area' },
interval: {'i2': { 'color': '#4374E0', 'style':'bars', 'lineWidth':4, 'fillOpacity':1 } }
,
};
var options4= {
title : 'Example',
vAxis: {title: 'Y Axis'},
hAxis: {title: 'X Axis'},
curveType: 'line',
lineWidth: 5,
intervals: { 'style':'area' },
explorer: {actions: ['dragToZoom', 'rightClickToReset'],axis: 'horizontal,vertical',keepInBounds: true, maxZoomIn: 128.0},hAxis: { logscale: true },vAxis: { logscale: true }};
var chart = new google.visualization.ComboChart(document.getElementById('chart_div0'));chart.draw(data0, options0);
var chart = new google.visualization.ComboChart(document.getElementById('chart_div1'));chart.draw(data1, options1);
var chart = new google.visualization.ComboChart(document.getElementById('chart_div2'));chart.draw(data2, options2);
var chart = new google.visualization.ComboChart(document.getElementById('chart_div3'));chart.draw(data3, options3);
var chart = new google.visualization.ComboChart(document.getElementById('chart_div4'));chart.draw(data4, options4);
}</script>
\endhtmlonly
 *
 */
