#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_move0_ = CSVReader(FileName=['Vector/vector_move0_0.csv', 'Vector/vector_move0_1.csv', 'Vector/vector_move0_2.csv', 'Vector/vector_move0_3.csv', 'Vector/vector_move0_4.csv', 'Vector/vector_move0_5.csv', 'Vector/vector_move0_6.csv', 'Vector/vector_move0_7.csv', 'Vector/vector_move0_8.csv', 'Vector/vector_move0_9.csv', 'Vector/vector_move0_10.csv', 'Vector/vector_move0_11.csv', 'Vector/vector_move0_12.csv', 'Vector/vector_move0_13.csv', 'Vector/vector_move0_14.csv', 'Vector/vector_move0_15.csv', 'Vector/vector_move0_16.csv', 'Vector/vector_move0_17.csv', 'Vector/vector_move0_18.csv', 'Vector/vector_move0_19.csv', 'Vector/vector_move0_20.csv', 'Vector/vector_move0_21.csv', 'Vector/vector_move0_22.csv', 'Vector/vector_move0_23.csv', 'Vector/vector_move0_24.csv', 'Vector/vector_move0_25.csv', 'Vector/vector_move0_26.csv', 'Vector/vector_move0_27.csv', 'Vector/vector_move0_28.csv', 'Vector/vector_move0_29.csv', 'Vector/vector_move0_30.csv', 'Vector/vector_move0_31.csv', 'Vector/vector_move0_32.csv', 'Vector/vector_move0_33.csv', 'Vector/vector_move0_34.csv', 'Vector/vector_move0_35.csv', 'Vector/vector_move0_36.csv', 'Vector/vector_move0_37.csv', 'Vector/vector_move0_38.csv', 'Vector/vector_move0_39.csv', 'Vector/vector_move0_40.csv', 'Vector/vector_move0_41.csv', 'Vector/vector_move0_42.csv', 'Vector/vector_move0_43.csv', 'Vector/vector_move0_44.csv', 'Vector/vector_move0_45.csv', 'Vector/vector_move0_46.csv', 'Vector/vector_move0_47.csv', 'Vector/vector_move0_48.csv', 'Vector/vector_move0_49.csv', 'Vector/vector_move0_50.csv', 'Vector/vector_move0_51.csv', 'Vector/vector_move0_52.csv', 'Vector/vector_move0_53.csv', 'Vector/vector_move0_54.csv', 'Vector/vector_move0_55.csv', 'Vector/vector_move0_56.csv', 'Vector/vector_move0_57.csv', 'Vector/vector_move0_58.csv', 'Vector/vector_move0_59.csv', 'Vector/vector_move0_60.csv', 'Vector/vector_move0_61.csv', 'Vector/vector_move0_62.csv', 'Vector/vector_move0_63.csv', 'Vector/vector_move0_64.csv', 'Vector/vector_move0_65.csv', 'Vector/vector_move0_66.csv', 'Vector/vector_move0_67.csv', 'Vector/vector_move0_68.csv', 'Vector/vector_move0_69.csv', 'Vector/vector_move0_70.csv', 'Vector/vector_move0_71.csv', 'Vector/vector_move0_72.csv', 'Vector/vector_move0_73.csv', 'Vector/vector_move0_74.csv', 'Vector/vector_move0_75.csv', 'Vector/vector_move0_76.csv', 'Vector/vector_move0_77.csv', 'Vector/vector_move0_78.csv', 'Vector/vector_move0_79.csv', 'Vector/vector_move0_80.csv', 'Vector/vector_move0_81.csv', 'Vector/vector_move0_82.csv', 'Vector/vector_move0_83.csv', 'Vector/vector_move0_84.csv', 'Vector/vector_move0_85.csv', 'Vector/vector_move0_86.csv', 'Vector/vector_move0_87.csv', 'Vector/vector_move0_88.csv', 'Vector/vector_move0_89.csv', 'Vector/vector_move0_90.csv', 'Vector/vector_move0_91.csv', 'Vector/vector_move0_92.csv', 'Vector/vector_move0_93.csv', 'Vector/vector_move0_94.csv', 'Vector/vector_move0_95.csv', 'Vector/vector_move0_96.csv', 'Vector/vector_move0_97.csv', 'Vector/vector_move0_98.csv', 'Vector/vector_move0_99.csv'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1020, 495]

# get layout
viewLayout1 = GetLayout()

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=vector_move0_)

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x[0]'
tableToPoints1.YColumn = 'x[1]'
tableToPoints1.a2DPoints = 1

# show data in view
tableToPoints1Display = Show(tableToPoints1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.500257785, 0.1870064525, 10000.0]
renderView1.CameraFocalPoint = [0.500257785, 0.1870064525, 0.0]

# create a new 'CSV Reader'
vector_move1_ = CSVReader(FileName=['Vector/vector_move1_0.csv', 'Vector/vector_move1_1.csv', 'Vector/vector_move1_2.csv', 'Vector/vector_move1_3.csv', 'Vector/vector_move1_4.csv', 'Vector/vector_move1_5.csv', 'Vector/vector_move1_6.csv', 'Vector/vector_move1_7.csv', 'Vector/vector_move1_8.csv', 'Vector/vector_move1_9.csv', 'Vector/vector_move1_10.csv', 'Vector/vector_move1_11.csv', 'Vector/vector_move1_12.csv', 'Vector/vector_move1_13.csv', 'Vector/vector_move1_14.csv', 'Vector/vector_move1_15.csv', 'Vector/vector_move1_16.csv', 'Vector/vector_move1_17.csv', 'Vector/vector_move1_18.csv', 'Vector/vector_move1_19.csv', 'Vector/vector_move1_20.csv', 'Vector/vector_move1_21.csv', 'Vector/vector_move1_22.csv', 'Vector/vector_move1_23.csv', 'Vector/vector_move1_24.csv', 'Vector/vector_move1_25.csv', 'Vector/vector_move1_26.csv', 'Vector/vector_move1_27.csv', 'Vector/vector_move1_28.csv', 'Vector/vector_move1_29.csv', 'Vector/vector_move1_30.csv', 'Vector/vector_move1_31.csv', 'Vector/vector_move1_32.csv', 'Vector/vector_move1_33.csv', 'Vector/vector_move1_34.csv', 'Vector/vector_move1_35.csv', 'Vector/vector_move1_36.csv', 'Vector/vector_move1_37.csv', 'Vector/vector_move1_38.csv', 'Vector/vector_move1_39.csv', 'Vector/vector_move1_40.csv', 'Vector/vector_move1_41.csv', 'Vector/vector_move1_42.csv', 'Vector/vector_move1_43.csv', 'Vector/vector_move1_44.csv', 'Vector/vector_move1_45.csv', 'Vector/vector_move1_46.csv', 'Vector/vector_move1_47.csv', 'Vector/vector_move1_48.csv', 'Vector/vector_move1_49.csv', 'Vector/vector_move1_50.csv', 'Vector/vector_move1_51.csv', 'Vector/vector_move1_52.csv', 'Vector/vector_move1_53.csv', 'Vector/vector_move1_54.csv', 'Vector/vector_move1_55.csv', 'Vector/vector_move1_56.csv', 'Vector/vector_move1_57.csv', 'Vector/vector_move1_58.csv', 'Vector/vector_move1_59.csv', 'Vector/vector_move1_60.csv', 'Vector/vector_move1_61.csv', 'Vector/vector_move1_62.csv', 'Vector/vector_move1_63.csv', 'Vector/vector_move1_64.csv', 'Vector/vector_move1_65.csv', 'Vector/vector_move1_66.csv', 'Vector/vector_move1_67.csv', 'Vector/vector_move1_68.csv', 'Vector/vector_move1_69.csv', 'Vector/vector_move1_70.csv', 'Vector/vector_move1_71.csv', 'Vector/vector_move1_72.csv', 'Vector/vector_move1_73.csv', 'Vector/vector_move1_74.csv', 'Vector/vector_move1_75.csv', 'Vector/vector_move1_76.csv', 'Vector/vector_move1_77.csv', 'Vector/vector_move1_78.csv', 'Vector/vector_move1_79.csv', 'Vector/vector_move1_80.csv', 'Vector/vector_move1_81.csv', 'Vector/vector_move1_82.csv', 'Vector/vector_move1_83.csv', 'Vector/vector_move1_84.csv', 'Vector/vector_move1_85.csv', 'Vector/vector_move1_86.csv', 'Vector/vector_move1_87.csv', 'Vector/vector_move1_88.csv', 'Vector/vector_move1_89.csv', 'Vector/vector_move1_90.csv', 'Vector/vector_move1_91.csv', 'Vector/vector_move1_92.csv', 'Vector/vector_move1_93.csv', 'Vector/vector_move1_94.csv', 'Vector/vector_move1_95.csv', 'Vector/vector_move1_96.csv', 'Vector/vector_move1_97.csv', 'Vector/vector_move1_98.csv', 'Vector/vector_move1_99.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints2 = TableToPoints(Input=vector_move1_)

# Properties modified on tableToPoints2
tableToPoints2.XColumn = 'x[0]'
tableToPoints2.YColumn = 'x[1]'
tableToPoints2.a2DPoints = 1

# show data in view
tableToPoints2Display = Show(tableToPoints2, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(tableToPoints1)

# change representation type
tableToPoints1Display.SetRepresentationType('Points')

# set active source
SetActiveSource(tableToPoints2)

# change representation type
tableToPoints2Display.SetRepresentationType('Points')

# create a new 'CSV Reader'
vector_move2_ = CSVReader(FileName=['Vector/vector_move2_0.csv', 'Vector/vector_move2_1.csv', 'Vector/vector_move2_2.csv', 'Vector/vector_move2_3.csv', 'Vector/vector_move2_4.csv', 'Vector/vector_move2_5.csv', 'Vector/vector_move2_6.csv', 'Vector/vector_move2_7.csv', 'Vector/vector_move2_8.csv', 'Vector/vector_move2_9.csv', 'Vector/vector_move2_10.csv', 'Vector/vector_move2_11.csv', 'Vector/vector_move2_12.csv', 'Vector/vector_move2_13.csv', 'Vector/vector_move2_14.csv', 'Vector/vector_move2_15.csv', 'Vector/vector_move2_16.csv', 'Vector/vector_move2_17.csv', 'Vector/vector_move2_18.csv', 'Vector/vector_move2_19.csv', 'Vector/vector_move2_20.csv', 'Vector/vector_move2_21.csv', 'Vector/vector_move2_22.csv', 'Vector/vector_move2_23.csv', 'Vector/vector_move2_24.csv', 'Vector/vector_move2_25.csv', 'Vector/vector_move2_26.csv', 'Vector/vector_move2_27.csv', 'Vector/vector_move2_28.csv', 'Vector/vector_move2_29.csv', 'Vector/vector_move2_30.csv', 'Vector/vector_move2_31.csv', 'Vector/vector_move2_32.csv', 'Vector/vector_move2_33.csv', 'Vector/vector_move2_34.csv', 'Vector/vector_move2_35.csv', 'Vector/vector_move2_36.csv', 'Vector/vector_move2_37.csv', 'Vector/vector_move2_38.csv', 'Vector/vector_move2_39.csv', 'Vector/vector_move2_40.csv', 'Vector/vector_move2_41.csv', 'Vector/vector_move2_42.csv', 'Vector/vector_move2_43.csv', 'Vector/vector_move2_44.csv', 'Vector/vector_move2_45.csv', 'Vector/vector_move2_46.csv', 'Vector/vector_move2_47.csv', 'Vector/vector_move2_48.csv', 'Vector/vector_move2_49.csv', 'Vector/vector_move2_50.csv', 'Vector/vector_move2_51.csv', 'Vector/vector_move2_52.csv', 'Vector/vector_move2_53.csv', 'Vector/vector_move2_54.csv', 'Vector/vector_move2_55.csv', 'Vector/vector_move2_56.csv', 'Vector/vector_move2_57.csv', 'Vector/vector_move2_58.csv', 'Vector/vector_move2_59.csv', 'Vector/vector_move2_60.csv', 'Vector/vector_move2_61.csv', 'Vector/vector_move2_62.csv', 'Vector/vector_move2_63.csv', 'Vector/vector_move2_64.csv', 'Vector/vector_move2_65.csv', 'Vector/vector_move2_66.csv', 'Vector/vector_move2_67.csv', 'Vector/vector_move2_68.csv', 'Vector/vector_move2_69.csv', 'Vector/vector_move2_70.csv', 'Vector/vector_move2_71.csv', 'Vector/vector_move2_72.csv', 'Vector/vector_move2_73.csv', 'Vector/vector_move2_74.csv', 'Vector/vector_move2_75.csv', 'Vector/vector_move2_76.csv', 'Vector/vector_move2_77.csv', 'Vector/vector_move2_78.csv', 'Vector/vector_move2_79.csv', 'Vector/vector_move2_80.csv', 'Vector/vector_move2_81.csv', 'Vector/vector_move2_82.csv', 'Vector/vector_move2_83.csv', 'Vector/vector_move2_84.csv', 'Vector/vector_move2_85.csv', 'Vector/vector_move2_86.csv', 'Vector/vector_move2_87.csv', 'Vector/vector_move2_88.csv', 'Vector/vector_move2_89.csv', 'Vector/vector_move2_90.csv', 'Vector/vector_move2_91.csv', 'Vector/vector_move2_92.csv', 'Vector/vector_move2_93.csv', 'Vector/vector_move2_94.csv', 'Vector/vector_move2_95.csv', 'Vector/vector_move2_96.csv', 'Vector/vector_move2_97.csv', 'Vector/vector_move2_98.csv', 'Vector/vector_move2_99.csv'])

# set active view
SetActiveView(renderView1)

# create a new 'Table To Points'
tableToPoints3 = TableToPoints(Input=vector_move2_)

# Properties modified on tableToPoints3
tableToPoints3.XColumn = 'x[0]'
tableToPoints3.YColumn = 'x[1]'
tableToPoints3.a2DPoints = 1

# show data in view
tableToPoints3Display = Show(tableToPoints3, renderView1)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.500684285, 0.49876095249999997, 2.7250155738876063]
renderView1.CameraFocalPoint = [0.500684285, 0.49876095249999997, 0.0]
renderView1.CameraParallelScale = 0.7052859287230878

# save animation images/movie
WriteAnimation('generated/particles_mooving.ogv', Magnification=1, FrameRate=25.0, Compression=True)
