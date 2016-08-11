#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
vector_move0_ = CSVReader(FileName=['Vector/vector_move_0_0.csv', 'Vector/vector_move_0_1.csv', 'Vector/vector_move_0_2.csv', 'Vector/vector_move_0_3.csv', 'Vector/vector_move_0_4.csv', 'Vector/vector_move_0_5.csv', 'Vector/vector_move_0_6.csv', 'Vector/vector_move_0_7.csv', 'Vector/vector_move_0_8.csv', 'Vector/vector_move_0_9.csv', 'Vector/vector_move_0_10.csv', 'Vector/vector_move_0_11.csv', 'Vector/vector_move_0_12.csv', 'Vector/vector_move_0_13.csv', 'Vector/vector_move_0_14.csv', 'Vector/vector_move_0_15.csv', 'Vector/vector_move_0_16.csv', 'Vector/vector_move_0_17.csv', 'Vector/vector_move_0_18.csv', 'Vector/vector_move_0_19.csv', 'Vector/vector_move_0_20.csv', 'Vector/vector_move_0_21.csv', 'Vector/vector_move_0_22.csv', 'Vector/vector_move_0_23.csv', 'Vector/vector_move_0_24.csv', 'Vector/vector_move_0_25.csv', 'Vector/vector_move_0_26.csv', 'Vector/vector_move_0_27.csv', 'Vector/vector_move_0_28.csv', 'Vector/vector_move_0_29.csv', 'Vector/vector_move_0_30.csv', 'Vector/vector_move_0_31.csv', 'Vector/vector_move_0_32.csv', 'Vector/vector_move_0_33.csv', 'Vector/vector_move_0_34.csv', 'Vector/vector_move_0_35.csv', 'Vector/vector_move_0_36.csv', 'Vector/vector_move_0_37.csv', 'Vector/vector_move_0_38.csv', 'Vector/vector_move_0_39.csv', 'Vector/vector_move_0_40.csv', 'Vector/vector_move_0_41.csv', 'Vector/vector_move_0_42.csv', 'Vector/vector_move_0_43.csv', 'Vector/vector_move_0_44.csv', 'Vector/vector_move_0_45.csv', 'Vector/vector_move_0_46.csv', 'Vector/vector_move_0_47.csv', 'Vector/vector_move_0_48.csv', 'Vector/vector_move_0_49.csv', 'Vector/vector_move_0_50.csv', 'Vector/vector_move_0_51.csv', 'Vector/vector_move_0_52.csv', 'Vector/vector_move_0_53.csv', 'Vector/vector_move_0_54.csv', 'Vector/vector_move_0_55.csv', 'Vector/vector_move_0_56.csv', 'Vector/vector_move_0_57.csv', 'Vector/vector_move_0_58.csv', 'Vector/vector_move_0_59.csv', 'Vector/vector_move_0_60.csv', 'Vector/vector_move_0_61.csv', 'Vector/vector_move_0_62.csv', 'Vector/vector_move_0_63.csv', 'Vector/vector_move_0_64.csv', 'Vector/vector_move_0_65.csv', 'Vector/vector_move_0_66.csv', 'Vector/vector_move_0_67.csv', 'Vector/vector_move_0_68.csv', 'Vector/vector_move_0_69.csv', 'Vector/vector_move_0_70.csv', 'Vector/vector_move_0_71.csv', 'Vector/vector_move_0_72.csv', 'Vector/vector_move_0_73.csv', 'Vector/vector_move_0_74.csv', 'Vector/vector_move_0_75.csv', 'Vector/vector_move_0_76.csv', 'Vector/vector_move_0_77.csv', 'Vector/vector_move_0_78.csv', 'Vector/vector_move_0_79.csv', 'Vector/vector_move_0_80.csv', 'Vector/vector_move_0_81.csv', 'Vector/vector_move_0_82.csv', 'Vector/vector_move_0_83.csv', 'Vector/vector_move_0_84.csv', 'Vector/vector_move_0_85.csv', 'Vector/vector_move_0_86.csv', 'Vector/vector_move_0_87.csv', 'Vector/vector_move_0_88.csv', 'Vector/vector_move_0_89.csv', 'Vector/vector_move_0_90.csv', 'Vector/vector_move_0_91.csv', 'Vector/vector_move_0_92.csv', 'Vector/vector_move_0_93.csv', 'Vector/vector_move_0_94.csv', 'Vector/vector_move_0_95.csv', 'Vector/vector_move_0_96.csv', 'Vector/vector_move_0_97.csv', 'Vector/vector_move_0_98.csv', 'Vector/vector_move_0_99.csv'])

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
vector_move1_ = CSVReader(FileName=['Vector/vector_move_1_0.csv', 'Vector/vector_move_1_1.csv', 'Vector/vector_move_1_2.csv', 'Vector/vector_move_1_3.csv', 'Vector/vector_move_1_4.csv', 'Vector/vector_move_1_5.csv', 'Vector/vector_move_1_6.csv', 'Vector/vector_move_1_7.csv', 'Vector/vector_move_1_8.csv', 'Vector/vector_move_1_9.csv', 'Vector/vector_move_1_10.csv', 'Vector/vector_move_1_11.csv', 'Vector/vector_move_1_12.csv', 'Vector/vector_move_1_13.csv', 'Vector/vector_move_1_14.csv', 'Vector/vector_move_1_15.csv', 'Vector/vector_move_1_16.csv', 'Vector/vector_move_1_17.csv', 'Vector/vector_move_1_18.csv', 'Vector/vector_move_1_19.csv', 'Vector/vector_move_1_20.csv', 'Vector/vector_move_1_21.csv', 'Vector/vector_move_1_22.csv', 'Vector/vector_move_1_23.csv', 'Vector/vector_move_1_24.csv', 'Vector/vector_move_1_25.csv', 'Vector/vector_move_1_26.csv', 'Vector/vector_move_1_27.csv', 'Vector/vector_move_1_28.csv', 'Vector/vector_move_1_29.csv', 'Vector/vector_move_1_30.csv', 'Vector/vector_move_1_31.csv', 'Vector/vector_move_1_32.csv', 'Vector/vector_move_1_33.csv', 'Vector/vector_move_1_34.csv', 'Vector/vector_move_1_35.csv', 'Vector/vector_move_1_36.csv', 'Vector/vector_move_1_37.csv', 'Vector/vector_move_1_38.csv', 'Vector/vector_move_1_39.csv', 'Vector/vector_move_1_40.csv', 'Vector/vector_move_1_41.csv', 'Vector/vector_move_1_42.csv', 'Vector/vector_move_1_43.csv', 'Vector/vector_move_1_44.csv', 'Vector/vector_move_1_45.csv', 'Vector/vector_move_1_46.csv', 'Vector/vector_move_1_47.csv', 'Vector/vector_move_1_48.csv', 'Vector/vector_move_1_49.csv', 'Vector/vector_move_1_50.csv', 'Vector/vector_move_1_51.csv', 'Vector/vector_move_1_52.csv', 'Vector/vector_move_1_53.csv', 'Vector/vector_move_1_54.csv', 'Vector/vector_move_1_55.csv', 'Vector/vector_move_1_56.csv', 'Vector/vector_move_1_57.csv', 'Vector/vector_move_1_58.csv', 'Vector/vector_move_1_59.csv', 'Vector/vector_move_1_60.csv', 'Vector/vector_move_1_61.csv', 'Vector/vector_move_1_62.csv', 'Vector/vector_move_1_63.csv', 'Vector/vector_move_1_64.csv', 'Vector/vector_move_1_65.csv', 'Vector/vector_move_1_66.csv', 'Vector/vector_move_1_67.csv', 'Vector/vector_move_1_68.csv', 'Vector/vector_move_1_69.csv', 'Vector/vector_move_1_70.csv', 'Vector/vector_move_1_71.csv', 'Vector/vector_move_1_72.csv', 'Vector/vector_move_1_73.csv', 'Vector/vector_move_1_74.csv', 'Vector/vector_move_1_75.csv', 'Vector/vector_move_1_76.csv', 'Vector/vector_move_1_77.csv', 'Vector/vector_move_1_78.csv', 'Vector/vector_move_1_79.csv', 'Vector/vector_move_1_80.csv', 'Vector/vector_move_1_81.csv', 'Vector/vector_move_1_82.csv', 'Vector/vector_move_1_83.csv', 'Vector/vector_move_1_84.csv', 'Vector/vector_move_1_85.csv', 'Vector/vector_move_1_86.csv', 'Vector/vector_move_1_87.csv', 'Vector/vector_move_1_88.csv', 'Vector/vector_move_1_89.csv', 'Vector/vector_move_1_90.csv', 'Vector/vector_move_1_91.csv', 'Vector/vector_move_1_92.csv', 'Vector/vector_move_1_93.csv', 'Vector/vector_move_1_94.csv', 'Vector/vector_move_1_95.csv', 'Vector/vector_move_1_96.csv', 'Vector/vector_move_1_97.csv', 'Vector/vector_move_1_98.csv', 'Vector/vector_move_1_99.csv'])

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
vector_move2_ = CSVReader(FileName=['Vector/vector_move_2_0.csv', 'Vector/vector_move_2_1.csv', 'Vector/vector_move_2_2.csv', 'Vector/vector_move_2_3.csv', 'Vector/vector_move_2_4.csv', 'Vector/vector_move_2_5.csv', 'Vector/vector_move_2_6.csv', 'Vector/vector_move_2_7.csv', 'Vector/vector_move_2_8.csv', 'Vector/vector_move_2_9.csv', 'Vector/vector_move_2_10.csv', 'Vector/vector_move_2_11.csv', 'Vector/vector_move_2_12.csv', 'Vector/vector_move_2_13.csv', 'Vector/vector_move_2_14.csv', 'Vector/vector_move_2_15.csv', 'Vector/vector_move_2_16.csv', 'Vector/vector_move_2_17.csv', 'Vector/vector_move_2_18.csv', 'Vector/vector_move_2_19.csv', 'Vector/vector_move_2_20.csv', 'Vector/vector_move_2_21.csv', 'Vector/vector_move_2_22.csv', 'Vector/vector_move_2_23.csv', 'Vector/vector_move_2_24.csv', 'Vector/vector_move_2_25.csv', 'Vector/vector_move_2_26.csv', 'Vector/vector_move_2_27.csv', 'Vector/vector_move_2_28.csv', 'Vector/vector_move_2_29.csv', 'Vector/vector_move_2_30.csv', 'Vector/vector_move_2_31.csv', 'Vector/vector_move_2_32.csv', 'Vector/vector_move_2_33.csv', 'Vector/vector_move_2_34.csv', 'Vector/vector_move_2_35.csv', 'Vector/vector_move_2_36.csv', 'Vector/vector_move_2_37.csv', 'Vector/vector_move_2_38.csv', 'Vector/vector_move_2_39.csv', 'Vector/vector_move_2_40.csv', 'Vector/vector_move_2_41.csv', 'Vector/vector_move_2_42.csv', 'Vector/vector_move_2_43.csv', 'Vector/vector_move_2_44.csv', 'Vector/vector_move_2_45.csv', 'Vector/vector_move_2_46.csv', 'Vector/vector_move_2_47.csv', 'Vector/vector_move_2_48.csv', 'Vector/vector_move_2_49.csv', 'Vector/vector_move_2_50.csv', 'Vector/vector_move_2_51.csv', 'Vector/vector_move_2_52.csv', 'Vector/vector_move_2_53.csv', 'Vector/vector_move_2_54.csv', 'Vector/vector_move_2_55.csv', 'Vector/vector_move_2_56.csv', 'Vector/vector_move_2_57.csv', 'Vector/vector_move_2_58.csv', 'Vector/vector_move_2_59.csv', 'Vector/vector_move_2_60.csv', 'Vector/vector_move_2_61.csv', 'Vector/vector_move_2_62.csv', 'Vector/vector_move_2_63.csv', 'Vector/vector_move_2_64.csv', 'Vector/vector_move_2_65.csv', 'Vector/vector_move_2_66.csv', 'Vector/vector_move_2_67.csv', 'Vector/vector_move_2_68.csv', 'Vector/vector_move_2_69.csv', 'Vector/vector_move_2_70.csv', 'Vector/vector_move_2_71.csv', 'Vector/vector_move_2_72.csv', 'Vector/vector_move_2_73.csv', 'Vector/vector_move_2_74.csv', 'Vector/vector_move_2_75.csv', 'Vector/vector_move_2_76.csv', 'Vector/vector_move_2_77.csv', 'Vector/vector_move_2_78.csv', 'Vector/vector_move_2_79.csv', 'Vector/vector_move_2_80.csv', 'Vector/vector_move_2_81.csv', 'Vector/vector_move_2_82.csv', 'Vector/vector_move_2_83.csv', 'Vector/vector_move_2_84.csv', 'Vector/vector_move_2_85.csv', 'Vector/vector_move_2_86.csv', 'Vector/vector_move_2_87.csv', 'Vector/vector_move_2_88.csv', 'Vector/vector_move_2_89.csv', 'Vector/vector_move_2_90.csv', 'Vector/vector_move_2_91.csv', 'Vector/vector_move_2_92.csv', 'Vector/vector_move_2_93.csv', 'Vector/vector_move_2_94.csv', 'Vector/vector_move_2_95.csv', 'Vector/vector_move_2_96.csv', 'Vector/vector_move_2_97.csv', 'Vector/vector_move_2_98.csv', 'Vector/vector_move_2_99.csv'])

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
