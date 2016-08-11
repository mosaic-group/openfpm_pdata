#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_before_map0csv = CSVReader(FileName=['Vector/vector_before_map_0.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [1000, 500]

# get layout
viewLayout1 = GetLayout()

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_before_map0csv)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.ZColumn = 'column_2_[0][1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 10000.0]
renderView1.CameraFocalPoint = [0.50000391295, 0.499549132735, 0.0]

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)

# Properties modified on calculator1
calculator1.Function = 'iHat*column_1_[0]+jHat*column_1_[1]'

# show data in view
calculator1Display = Show(calculator1, renderView1)

# hide data in view
Hide(tableToPoints1, renderView1)

# set active source
SetActiveSource(tableToPoints1)

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# set active source
SetActiveSource(calculator1)

# create a new 'Glyph'
glyph1 = Glyph(Input=calculator1,
    GlyphType='Arrow')

# Properties modified on glyph1
glyph1.ScaleFactor = 0.03
glyph1.MaximumNumberOfSamplePoints = 7000

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'column0'
column0LUT = GetColorTransferFunction('column0')

# get opacity transfer function/opacity map for 'column0'
column0PWF = GetOpacityTransferFunction('column0')

# Properties modified on glyph1
glyph1.Scalars = ['POINTS', 'column_1_[0]']

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'column_1_[0]'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# get color transfer function/color map for 'column10'
column10LUT = GetColorTransferFunction('column10')

# get opacity transfer function/opacity map for 'column10'
column10PWF = GetOpacityTransferFunction('column10')

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 9999.999999999998]
renderView1.CameraFocalPoint = [0.50000391295, 0.499549132735, 0.0]
renderView1.CameraViewUp = [0.024239745600625473, 0.9997061741998082, 0.0]
renderView1.CameraParallelScale = 0.18611748651565127

# save screenshot
#SaveScreenshot('Vector/vector_vectors.jpg', magnification=1, quality=100, view=renderView1)
WriteImage("generated/vector_vectors.jpg");

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
tableToPoints1 = FindSource('TableToPoints1')

# set active source
SetActiveSource(tableToPoints1)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [982, 495]

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# find source
glyph1 = FindSource('Glyph1')

# set active source
SetActiveSource(glyph1)

# get color transfer function/color map for 'column10'
column10LUT = GetColorTransferFunction('column10')

# get opacity transfer function/opacity map for 'column10'
column10PWF = GetOpacityTransferFunction('column10')

# find source
calculator1 = FindSource('Calculator1')

# set active source
SetActiveSource(calculator1)

# hide data in view
Hide(glyph1, renderView1)

# show data in view
calculator1Display = Show(calculator1, renderView1)

# destroy glyph1
Delete(glyph1)
del glyph1

# set active source
SetActiveSource(tableToPoints1)

# hide data in view
Hide(calculator1, renderView1)

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# destroy calculator1
Delete(calculator1)
del calculator1

# find source
cSVReader1 = FindSource('CSVReader1')

# set active source
SetActiveSource(cSVReader1)

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(tableToPoints1)

# set active source
SetActiveSource(cSVReader1)

# hide data in view
Hide(tableToPoints1, renderView1)

# destroy tableToPoints1
Delete(tableToPoints1)
del tableToPoints1

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=cSVReader1)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 10000.0]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]

# set scalar coloring
ColorBy(tableToPoints1Display, ('POINTS', 'column_1_[0]'))

# rescale color and/or opacity maps used to include current data range
tableToPoints1Display.RescaleTransferFunctionToDataRange(True)

# show color bar/color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, True)

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 5.0

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.50000391295, 0.499549132735, 10000.0]
renderView1.CameraFocalPoint = [0.50000391295, 0.499549132735, 0.0]
renderView1.CameraParallelScale = 0.15381610455838945

# hide color bar/color legend
tableToPoints1Display.SetScalarBarVisibility(renderView1, False)

# save screenshot
#SaveScreenshot('Vector/vector_vectors.jpg', magnification=1, quality=100, view=renderView1)
WriteImage("generated/vector_scalar.jpg");

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
