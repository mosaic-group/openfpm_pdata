#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_ghost_fill = CSVReader(FileName=['Vector/vector_ghost_fill_0.csv'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
renderView1.ViewSize = [400, 400]

# get layout
viewLayout1 = GetLayout()

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
vector_ghost_fillDisplay = Show(vector_ghost_fill, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_ghost_fill)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
#tableToPoints1Display = Show(tableToPoints1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.499922575, 0.18749537000000002, 10000.0]
renderView1.CameraFocalPoint = [0.499922575, 0.18749537000000002, 0.0]

# create a new 'Legacy VTK Reader'
vect_decompositionexternal_ghost_0vtk = LegacyVTKReader(FileNames=['Vector/vect_decompositionexternal_ghost_0.vtk'])

# destroy vect_decompositionexternal_ghost_0vtk
Delete(vect_decompositionexternal_ghost_0vtk)
del vect_decompositionexternal_ghost_0vtk

# create a new 'CSV Reader'
vector_after_map0csv = CSVReader(FileName=['Vector/vector_after_map_0.csv'])

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.BlockSize = 1024L
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

# place view in the layout
viewLayout1.AssignView(2, spreadSheetView1)

# show data in view
vector_after_map0csvDisplay = Show(vector_after_map0csv, spreadSheetView1)

# destroy spreadSheetView1
Delete(spreadSheetView1)
del spreadSheetView1

# close an empty frame
viewLayout1.Collapse(2)

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=vector_after_map0csv)

# Properties modified on tableToPoints2
tableToPoints2.XColumn = 'x[0]'
tableToPoints2.YColumn = 'x[1]'
tableToPoints2.a2DPoints = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)

tableToPoints2Display.PointSize = 4.0

# create a new 'Legacy VTK Reader'
vect_decompositionsubdomains_0vtk = LegacyVTKReader(FileNames=['Vector/vect_decompositionsubdomains_0.vtk'])

# show data in view
vect_decompositionsubdomains_0vtkDisplay = Show(vect_decompositionsubdomains_0vtk, renderView1)

# show color bar/color legend
vect_decompositionsubdomains_0vtkDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'data'
dataLUT = GetColorTransferFunction('data')

# get opacity transfer function/opacity map for 'data'
dataPWF = GetOpacityTransferFunction('data')

# hide data in view
Hide(tableToPoints1, renderView1)

# change representation type
vect_decompositionsubdomains_0vtkDisplay.SetRepresentationType('Points')

# change representation type
vect_decompositionsubdomains_0vtkDisplay.SetRepresentationType('Wireframe')

# turn off scalar coloring
ColorBy(vect_decompositionsubdomains_0vtkDisplay, None)

# set active source
SetActiveSource(tableToPoints2)

# change representation type
tableToPoints2Display.SetRepresentationType('Points')

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.45889311126344623, 0.2866467748466169, 10000.0]
renderView1.CameraFocalPoint = [0.45889311126344623, 0.2866467748466169, 0.0]
renderView1.CameraParallelScale = 0.1742200032569875

# create a new 'Legacy VTK Reader'
vect_decompositionexternal_ghost_0vtk = LegacyVTKReader(FileNames=['Vector/vect_decompositionexternal_ghost_0.vtk'])

# show data in view
vect_decompositionexternal_ghost_0vtkDisplay = Show(vect_decompositionexternal_ghost_0vtk, renderView1)

# show color bar/color legend
vect_decompositionexternal_ghost_0vtkDisplay.SetScalarBarVisibility(renderView1, False)

# hide data in view
Hide(vect_decompositionexternal_ghost_0vtk, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4322940474199577, 0.3362326839869475, 10000.0]
renderView1.CameraFocalPoint = [0.4322940474199577, 0.3362326839869475, 0.0]
renderView1.CameraParallelScale = 0.08127491729954842

# save screenshot
WriteImage('generated/vector_one_p_zoom.jpg')

# set active source
SetActiveSource(vect_decompositionexternal_ghost_0vtk)

# show data in view
vect_decompositionexternal_ghost_0vtkDisplay = Show(vect_decompositionexternal_ghost_0vtk, renderView1)

# show color bar/color legend
vect_decompositionexternal_ghost_0vtkDisplay.SetScalarBarVisibility(renderView1, True)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4322940474199577, 0.3362326839869475, 10000.0]
renderView1.CameraFocalPoint = [0.4322940474199577, 0.3362326839869475, 0.0]
renderView1.CameraParallelScale = 0.08127491729954842

# save screenshot
WriteImage('generated/vector_one_p_zoom_ghost.jpg')

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# change representation type
#tableToPoints1Display.SetRepresentationType('Points')

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 4.0

# change solid color
tableToPoints1Display.AmbientColor = [1.0, 0.0, 0.0]

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.4322940474199577, 0.3362326839869475, 10000.0]
renderView1.CameraFocalPoint = [0.4322940474199577, 0.3362326839869475, 0.0]
renderView1.CameraParallelScale = 0.08127491729954842

#### saving camera placements for all active views

# save screenshot
WriteImage('generated/vector_one_p_zoom_ghost_part.jpg')


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...). 
