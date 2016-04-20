#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_move0_1csv = CSVReader(FileName=['Vector/vector_move0_1.csv'])

# create a new 'CSV Reader'
vector_move1_1csv = CSVReader(FileName=['Vector/vector_move1_1.csv'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [982, 495]

# get layout
viewLayout1 = GetLayout()

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_move0_1csv)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=vector_move1_1csv)

# Properties modified on tableToPoints2
tableToPoints2.XColumn = 'x[0]'
tableToPoints2.YColumn = 'x[1]'
tableToPoints2.a2DPoints = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)

# change solid color
tableToPoints2Display.DiffuseColor = [1.0, 0.0, 0.0]

# change representation type
tableToPoints2Display.SetRepresentationType('Points')

# change solid color
tableToPoints2Display.AmbientColor = [1.0, 0.0, 0.0]

# Properties modified on tableToPoints2Display
tableToPoints2Display.PointSize = 10.0

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 10.0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.8363000195118148, 0.3097736377565421, 10000.0]
renderView1.CameraFocalPoint = [0.8363000195118148, 0.3097736377565421, 0.0]
renderView1.CameraParallelScale = 0.036959899933429505

WriteImage("generated/particle_map1.jpg")

# create a new 'Calculator'
calculator1 = Calculator(Input=tableToPoints1)

# Properties modified on calculator1
calculator1.Function = 'iHat*column_1_[0]+jHat*column_1_[1]'

# show data in view
calculator1Display = Show(calculator1, renderView1)

# create a new 'Glyph'
glyph1 = Glyph(Input=calculator1,
    GlyphType='Arrow')

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'column10'
column10LUT = GetColorTransferFunction('column10')

# get opacity transfer function/opacity map for 'column10'
column10PWF = GetOpacityTransferFunction('column10')

# Properties modified on glyph1
glyph1.ScaleFactor = 0.00800218604000001

# create a new 'Calculator'
calculator2 = Calculator(Input=tableToPoints2)

# Properties modified on calculator2
calculator2.Function = 'iHat*column_1_[0]+jHat*column_1_[1]'

# show data in view
calculator2Display = Show(calculator2, renderView1)

# create a new 'Glyph'
glyph2 = Glyph(Input=calculator2,
    GlyphType='Arrow')

# Properties modified on glyph2
glyph2.ScaleFactor = 0.008

# show data in view
glyph2Display = Show(glyph2, renderView1)

# show color bar/color legend
glyph2Display.SetScalarBarVisibility(renderView1, True)

# hide data in view
Hide(glyph1, renderView1)

# hide data in view
Hide(calculator2, renderView1)

# hide data in view
Hide(calculator1, renderView1)

# hide data in view
Hide(glyph2, renderView1)

# show data in view
glyph2Display = Show(glyph2, renderView1)

# show color bar/color legend
glyph2Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# create a new 'CSV Reader'
vector_move_before_map0_2csv = CSVReader(FileName=['Vector/vector_move_before_map0_2.csv'])

# create a new 'CSV Reader'
vector_move_before_map1_2csv = CSVReader(FileName=['Vector/vector_move_before_map1_2.csv'])

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=vector_move_before_map0_2csv)

# Properties modified on tableToPoints3
tableToPoints3.XColumn = 'x[0]'
tableToPoints3.YColumn = 'x[1]'
tableToPoints3.a2DPoints = 1

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)

# change representation type
tableToPoints3Display.SetRepresentationType('Points')

# Properties modified on tableToPoints3Display
tableToPoints3Display.PointSize = 10.0

# set active source
SetActiveSource(glyph1)

# Properties modified on glyph1
glyph1.GlyphMode = 'All Points'

# hide data in view
Hide(tableToPoints3, renderView1)

# create a new 'Table To Points'
tableToPoints4 = TableToPoints(Input=vector_move_before_map1_2csv)

# Properties modified on tableToPoints4
tableToPoints4.XColumn = 'x[0]'
tableToPoints4.YColumn = 'x[1]'
tableToPoints4.a2DPoints = 1

# show data in view
tableToPoints4Display = Show(tableToPoints4, renderView1)

# change representation type
tableToPoints4Display.SetRepresentationType('Points')

# Properties modified on tableToPoints4Display
tableToPoints4Display.PointSize = 10.0

# change solid color
tableToPoints4Display.AmbientColor = [1.0, 0.0, 0.0]

# hide data in view
Hide(tableToPoints4, renderView1)

# set active source
SetActiveSource(glyph2)

# Properties modified on glyph2
glyph2.GlyphMode = 'All Points'

# create a new 'CSV Reader'
vector_move1_2csv = CSVReader(FileName=['Vector/vector_move1_2.csv'])

# create a new 'CSV Reader'
vector_move0_2csv = CSVReader(FileName=['Vector/vector_move0_2.csv'])

# create a new 'Table To Points'
tableToPoints5 = TableToPoints(Input=vector_move1_2csv)

# Properties modified on tableToPoints5
tableToPoints5.XColumn = 'x[0]'
tableToPoints5.YColumn = 'x[1]'
tableToPoints5.a2DPoints = 1

# show data in view
tableToPoints5Display = Show(tableToPoints5, renderView1)

# change representation type
tableToPoints5Display.SetRepresentationType('Points')

# change solid color
tableToPoints5Display.AmbientColor = [1.0, 0.0, 0.0]

# Properties modified on tableToPoints5Display
tableToPoints5Display.PointSize = 10.0

# hide data in view
Hide(glyph1, renderView1)

# hide data in view
Hide(glyph2, renderView1)

# hide data in view
Hide(tableToPoints5, renderView1)

# hide data in view
Hide(tableToPoints4, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.9223207232553731, 0.3210044613442961, 10000.0]
renderView1.CameraFocalPoint = [0.9223207232553731, 0.3210044613442961, 0.0]
renderView1.CameraParallelScale = 0.0305453718458095

WriteImage("generated/particle_map1.jpg")

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(glyph2)

# show data in view
glyph2Display = Show(glyph2, renderView1)

# show color bar/color legend
glyph2Display.SetScalarBarVisibility(renderView1, True)

# save screenshot
WriteImage("generated/particle_map2.jpg")

# set active source
SetActiveSource(tableToPoints3)

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)

# set active source
SetActiveSource(tableToPoints4)

# show data in view
tableToPoints4Display = Show(tableToPoints4, renderView1)

# hide data in view
Hide(glyph1, renderView1)

# hide data in view
Hide(glyph2, renderView1)

# hide data in view
Hide(tableToPoints1, renderView1)

# hide data in view
Hide(tableToPoints2, renderView1)

# save screenshot
WriteImage("generated/particle_map3.jpg")

# set active source
SetActiveSource(vector_move0_2csv)

# create a new 'Table To Points'
tableToPoints6 = TableToPoints(Input=vector_move0_2csv)

# Properties modified on tableToPoints6
tableToPoints6.XColumn = 'x[0]'
tableToPoints6.YColumn = 'x[1]'
tableToPoints6.a2DPoints = 1

# show data in view
tableToPoints6Display = Show(tableToPoints6, renderView1)

# change representation type
tableToPoints6Display.SetRepresentationType('Points')

# Properties modified on tableToPoints6Display
tableToPoints6Display.PointSize = 10.0

# hide data in view
Hide(tableToPoints6, renderView1)

# hide data in view
Hide(tableToPoints3, renderView1)

# hide data in view
Hide(tableToPoints4, renderView1)

# set active source
SetActiveSource(tableToPoints5)

# show data in view
tableToPoints5Display = Show(tableToPoints5, renderView1)

# set active source
SetActiveSource(tableToPoints6)

# show data in view
tableToPoints6Display = Show(tableToPoints6, renderView1)

# hide data in view
Hide(tableToPoints5, renderView1)

# hide data in view
Hide(tableToPoints6, renderView1)

# show data in view
tableToPoints5Display = Show(tableToPoints5, renderView1)

# show data in view
tableToPoints6Display = Show(tableToPoints6, renderView1)

# hide data in view
Hide(tableToPoints1, renderView1)

# hide data in view
Hide(tableToPoints2, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.9223207232553731, 0.3210044613442961, 10000.0]
renderView1.CameraFocalPoint = [0.9223207232553731, 0.3210044613442961, 0.0]
renderView1.CameraParallelScale = 0.0305453718458095

# save screenshot
WriteImage("generated/particle_map4.jpg")

# hide data in view
Hide(tableToPoints5, renderView1)

# hide data in view
Hide(tableToPoints6, renderView1)

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4989760698, 0.18753136299999998, 2.059347991872556]
renderView1.CameraFocalPoint = [0.4989760698, 0.18753136299999998, 0.0]
renderView1.CameraParallelScale = 0.5329984807902486

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...). 
