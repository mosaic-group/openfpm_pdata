#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_before_map0csv = CSVReader(FileName=['Vector/vector_before_map0.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1020, 495]

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_before_map0csv)
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.ColorArrayName = [None, '']

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# create a new 'CSV Reader'
vector_before_map1csv = CSVReader(FileName=['Vector/vector_before_map1.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=vector_before_map1csv)
tableToPoints2.XColumn = 'x[0]'
tableToPoints2.YColumn = 'x[1]'
tableToPoints2.a2DPoints = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)
# trace defaults for the display properties.
tableToPoints2Display.ColorArrayName = [None, '']

# change representation type
tableToPoints2Display.SetRepresentationType('Points')

# change solid color
tableToPoints2Display.AmbientColor = [1.0, 1.0, 1.0]

# create a new 'CSV Reader'
vector_before_map2csv = CSVReader(FileName=['Vector/vector_before_map2.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=vector_before_map2csv)
tableToPoints3.XColumn = 'x[0]'
tableToPoints3.YColumn = 'x[1]'
tableToPoints3.a2DPoints = 1

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)
# trace defaults for the display properties.
tableToPoints3Display.ColorArrayName = [None, '']

# change representation type
tableToPoints3Display.SetRepresentationType('Points')

# change solid color
tableToPoints3Display.AmbientColor = [1.0, 1.0, 1.0]

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 10000.0]
renderView1.CameraFocalPoint = [0.50000391295, 0.499549132735, 0.0]
renderView1.CameraParallelScale = 0.7067808453124975

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...). 

# destroy tableToPoints1
Delete(tableToPoints1)
del tableToPoints1
Delete(tableToPoints2)
del tableToPoints2
Delete(tableToPoints3)
del tableToPoints3

# destroy vector_before_map0csv
Delete(vector_before_map0csv)
del vector_before_map0csv
Delete(vector_before_map1csv)
del vector_before_map1csv
Delete(vector_before_map2csv)
del vector_before_map2csv

# create a new 'CSV Reader'
vector_before_map0csv = CSVReader(FileName=['Vector/vector_after_map0.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1020, 495]

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_before_map0csv)
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)
# trace defaults for the display properties.
tableToPoints1Display.ColorArrayName = [None, '']

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# create a new 'CSV Reader'
vector_before_map1csv = CSVReader(FileName=['Vector/vector_after_map1.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=vector_before_map1csv)
tableToPoints2.XColumn = 'x[0]'
tableToPoints2.YColumn = 'x[1]'
tableToPoints2.a2DPoints = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)
# trace defaults for the display properties.
tableToPoints2Display.ColorArrayName = [None, '']

# change representation type
tableToPoints2Display.SetRepresentationType('Points')

# change solid color
tableToPoints2Display.AmbientColor = [1.0, 1.0, 1.0]

# create a new 'CSV Reader'
vector_before_map2csv = CSVReader(FileName=['Vector/vector_after_map2.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=vector_before_map2csv)
tableToPoints3.XColumn = 'x[0]'
tableToPoints3.YColumn = 'x[1]'
tableToPoints3.a2DPoints = 1

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)
# trace defaults for the display properties.
tableToPoints3Display.ColorArrayName = [None, '']

# change representation type
tableToPoints3Display.SetRepresentationType('Points')

# change solid color
tableToPoints3Display.AmbientColor = [1.0, 1.0, 1.0]

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 10000.0]
renderView1.CameraFocalPoint = [0.50000391295, 0.499549132735, 0.0]
renderView1.CameraParallelScale = 0.7067808453124975

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...). 

WriteImage("generated/vector_particles.jpg")
