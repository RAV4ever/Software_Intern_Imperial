import logging
import copy
import ctk
import numpy as np
import math as math
import qt
import slicer
import vtk
import os.path
import stl
from stl import mesh
from pyquaternion import Quaternion



from slicer.ScriptedLoadableModule import *
from Vector_Plotting import Vector_Plotting_Module as vpm
from numpy import linalg as LA

distance_map = []
start_target_fid = []
path = []
points = []
graph = []
solutions = []
safe_paths = []
final_paths = []
#
# EDEN_tract_planner
#

class EDEN_tract_planner(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "EDEN_tract_planner" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["Alberto Favaro (Politecnico di Milano)"]
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
...
""" # replace with organization, grant and thanks.

#
# EDEN_tract_plannerWidget
#

class EDEN_tract_plannerWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    self.orthogonal_visibility = False
    self.logic = EDEN_tract_plannerLogic()

    # Instantiate and connect widgets ...

    #
    # Data Area widget
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Tractography area selection"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)



    #
    # LOAD DATA button
    #
    self.loadButton = qt.QPushButton("Load data")
    self.loadButton.enabled = True
    parametersFormLayout.addRow(self.loadButton)


    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLModelNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Tractography: ", self.inputSelector)


    #
    # Skull volume selector
    #
    self.inputSkullSelector = slicer.qMRMLNodeComboBox()
    self.inputSkullSelector.nodeTypes = ["vtkMRMLModelNode"]
    self.inputSkullSelector.selectNodeUponCreation = True
    self.inputSkullSelector.addEnabled = False
    self.inputSkullSelector.removeEnabled = False
    self.inputSkullSelector.noneEnabled = False
    self.inputSkullSelector.showHidden = False
    self.inputSkullSelector.showChildNodeTypes = False
    self.inputSkullSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSkullSelector.setToolTip( "Pick the skull model." )
    parametersFormLayout.addRow("Skull: ", self.inputSkullSelector)




    #
    # combo box for FIDUCIALS
    #
    self.inputFiducialsNodeSelector = slicer.qMRMLNodeComboBox()
    self.inputFiducialsNodeSelector.objectName = 'inputFiducialsNodeSelector'
    self.inputFiducialsNodeSelector.nodeTypes = ['vtkMRMLMarkupsFiducialNode']
    self.inputFiducialsNodeSelector.toolTip = "Select the marker"
    self.inputFiducialsNodeSelector.selectNodeUponCreation = True
    self.inputFiducialsNodeSelector.noneEnabled = False
    self.inputFiducialsNodeSelector.addEnabled = False
    self.inputFiducialsNodeSelector.removeEnabled = False
    self.inputFiducialsNodeSelector.setMRMLScene(slicer.mrmlScene)
    parametersFormLayout.addRow("Fiducial:", self.inputFiducialsNodeSelector)

    self.OrthogonalBox = qt.QCheckBox()
    self.OrthogonalBox.checked = True
    parametersFormLayout.addRow("Orthogonal entries", self.OrthogonalBox)

    #
    # Sphere radius value
    #
    self.sphereRadiusSliderWidget = ctk.ctkSliderWidget()
    self.sphereRadiusSliderWidget.singleStep = 1
    self.sphereRadiusSliderWidget.minimum = 5
    self.sphereRadiusSliderWidget.maximum = 20
    self.sphereRadiusSliderWidget.value = 6
    self.sphereRadiusSliderWidget.setToolTip("Set the sphere radius.")
    parametersFormLayout.addRow("Sphere radius", self.sphereRadiusSliderWidget)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Direction")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)

    #
    # Path Button
    #
    self.PathButton = qt.QPushButton("Path")
    self.PathButton.toolTip = "Compute the path."
    self.PathButton.enabled = False
    parametersFormLayout.addRow(self.PathButton)

    #
    # Reset Start Button + Opposite area Button
    #

    self.oppositeStartButton = qt.QPushButton("Opposite Entry Point")
    self.oppositeStartButton.toolTip = "Opposite Entry Area"
    self.oppositeStartButton.enabled = False
    parametersFormLayout.addRow(self.oppositeStartButton)


    #
    # Reset Button
    #
    self.resetButton = qt.QPushButton("Reset")
    self.resetButton.toolTip = "Reset Everything"
    self.resetButton.enabled = False
    parametersFormLayout.addRow(self.resetButton)

    #
    # Trajectory Area widget
    #
    trajectoryCollapsibleButton = ctk.ctkCollapsibleButton()
    trajectoryCollapsibleButton.text = "Vectors"
    self.layout.addWidget(trajectoryCollapsibleButton)

    # Layout within the dummy collapsible button
    trajectoryCollapsibleLayout = qt.QFormLayout(trajectoryCollapsibleButton)

    self.DirectionBox=qt.QCheckBox()
    self.DirectionBox.checked = False
    trajectoryCollapsibleLayout.addRow("Fiber Direction", self.DirectionBox)

    self.TractBox=qt.QCheckBox()
    self.TractBox.checked = True
    trajectoryCollapsibleLayout.addRow("Tractography", self.TractBox)

    self.UsedTractBox = qt.QCheckBox()
    self.UsedTractBox.checked = False
    trajectoryCollapsibleLayout.addRow("Tractography Used", self.UsedTractBox)

    self.skullBox=qt.QCheckBox()
    self.skullBox.checked = True
    trajectoryCollapsibleLayout.addRow("skull", self.skullBox)

    self.SphereBox = qt.QCheckBox()
    self.SphereBox.checked = False
    trajectoryCollapsibleLayout.addRow("Sphere", self.SphereBox)

    self.EntryBox = qt.QCheckBox()
    self.EntryBox.checked = False
    trajectoryCollapsibleLayout.addRow("Entry Area", self.EntryBox)

    self.PathBox = qt.QCheckBox()
    self.PathBox.checked = False
    trajectoryCollapsibleLayout.addRow("Path", self.PathBox)

    #self.OrthogonalBox = qt.QCheckBox()
    #self.OrthogonalBox.checked = False
    #trajectoryCollapsibleLayout.addRow("Orthogonal entries", self.OrthogonalBox)


    # CONNECTIONS
    self.loadButton.connect('clicked(bool)', self.onLoadButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.PathButton.connect('clicked(bool)', self.onPathButton)
    self.resetButton.connect('clicked(bool)', self.onResetButton)
    self.oppositeStartButton.connect('clicked(bool)', self.onOppositeStartButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputSkullSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.sphereRadiusSliderWidget.connect('valueChanged(double)', self.onRadiusVariation)

    # Input fiducials connection
    self.inputFiducialsNodeSelector.connect('currentNodeChanged(bool)', self.onSelect)
    self.parent.connect('mrmlSceneChanged(vtkMRMLScene*)', self.inputFiducialsNodeSelector,'setMRMLScene(vtkMRMLScene*)')
    self.DirectionBox.connect('clicked(bool)', self.onRedVectorBoxClicked)
    self.TractBox.connect('clicked(bool)', self.onTractographyBoxClicked)
    self.UsedTractBox.connect('clicked(bool)', self.onUsedTractographyBoxClicked)
    self.skullBox.connect('clicked(bool)', self.onSkullBoxClicked)
    self.EntryBox.connect('clicked(bool)', self.onEntryBoxClicked)
    self.SphereBox.connect('clicked(bool)', self.onSphereBoxClicked)
    self.PathBox.connect('clicked(bool)', self.onPathBoxClicked)
    self.OrthogonalBox.connect('clicked(bool)', self.onOrthogonalBoxClicked)


    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()


  def onSelect(self):
    if self.inputFiducialsNodeSelector.currentNode():
      self.applyButton.enabled = True
      self.resetButton.enabled = True
      self.logic.plotSphere(self.inputFiducialsNodeSelector.currentNode(), self.sphereRadiusSliderWidget.value)

  def onRadiusVariation(self):
    self.logic.plotSphere(self.inputFiducialsNodeSelector.currentNode(), self.sphereRadiusSliderWidget.value)
    fid=False
    self.logic.resetAll(fid, self.inputFiducialsNodeSelector.currentNode())
    self.oppositeStartButton.enabled=False
    self.PathButton.enabled = False
    self.applyButton.enabled = True
    self.skullBox.checked = True
    #self.onSkullBoxClicked()
    #mod = slicer.util.getNode('skull')
    #disp = mod.GetDisplayNode()
    #disp.SetOpacity(0.4)

  def onLoadButton(self):
    global distance_map
    slicer.util.loadVolume(os.path.dirname(__file__) + '/Resources/3DT1s601a1006.nii', returnNode=True)
    slicer.modules.models.logic().AddModel(os.path.dirname(__file__) + '/Resources/sheep_CST_to_T13D.vtk')
    slicer.modules.models.logic().AddModel(os.path.dirname(__file__) + '/Resources/tubes.vtk')
    slicer.modules.models.logic().AddModel(os.path.dirname(__file__) + '/Resources/skull.vtk')
    slicer.modules.models.logic().AddModel(os.path.dirname(__file__) + '/Resources/Skull_rotated_hardened_top_decimated_transf.stl')
    #slicer.modules.models.logic().AddModel(os.path.dirname(__file__) + '/Resources/brain.vtk')
    # slicer.modules.models.logic().AddModel('/Users/albertofavaro/Dropbox/projects/Slicer-projects/EDEN2020 extension/Model_17_artery.vtk')
    # center the 3D view on the screen
    mod = slicer.util.getNode('sheep_CST_to_T13D')
    disp = mod.GetDisplayNode()
    disp.SetVisibility(False)
    mod = slicer.util.getNode('tubes')
    disp = mod.GetDisplayNode()
    disp.SetColor([0.8,0.1,0.9])
    mod = slicer.util.getNode('skull')
    disp = mod.GetDisplayNode()
    disp.SetColor([1,1,0.8])
    disp.SetOpacity(0.4)
    mod = slicer.util.getNode('brain')
    disp = mod.GetDisplayNode()
    disp.SetColor([0.7,0.9,0.6])


    layoutManager = slicer.app.layoutManager()
    threeDWidget = layoutManager.threeDWidget(0)
    threeDView = threeDWidget.threeDView()
    threeDView.resetFocalPoint()
  def onRedVectorBoxClicked(self):
    Display = slicer.util.getNode('RedVector').GetDisplayNode()
    Display.SetVisibility(self.DirectionBox.isChecked())

  def onTractographyBoxClicked(self):
    Display = self.inputSelector.currentNode()
    Display = Display.GetDisplayNode()
    Display.SetVisibility(self.TractBox.isChecked())

  def onUsedTractographyBoxClicked(self):
    Display = slicer.util.getNode('Used Tract').GetDisplayNode()
    Display.SetVisibility(self.UsedTractBox.isChecked())

  def onSkullBoxClicked(self):
    Display = self.inputSkullSelector.currentNode()
    Display = Display.GetDisplayNode()
    Display.SetVisibility(self.skullBox.isChecked())

  def onEntryBoxClicked(self):
    Display = slicer.util.getNode('acceptedRegion_poly').GetDisplayNode()
    Display.SetVisibility(self.EntryBox.isChecked())

  def onSphereBoxClicked(self):
    Display = slicer.util.getNode('SphereModelNode').GetDisplayNode()
    Display.SetVisibility(self.SphereBox.isChecked())

  def onOrthogonalBoxClicked(self):
    model = slicer.util.getNode('PlaneModelNode')
    if model != None:
      Display = slicer.util.getNode('PlaneModelNode').GetDisplayNode()
      Display.SetVisibility(self.OrthogonalBox.isChecked())
    else:
      print "No orthogonal directions to show, first define a target on the tract"


  def onPathBoxClicked(self):
    Display = slicer.util.getNode('direction_line').GetDisplayNode()
    Display.SetVisibility(self.PathBox.isChecked())


  def onApplyButton(self):
    #logic = EDEN_tract_plannerLogic()
    sphereRadius = self.sphereRadiusSliderWidget.value
    self.logic.run(self.inputSelector.currentNode(), self.inputSkullSelector.currentNode(), sphereRadius, self.inputFiducialsNodeSelector.currentNode())
    #self.vector.setCurrentNodeID(self.logic.OutputModelRed.GetID())
    #mod = slicer.util.getNode('skull')
    #disp = mod.GetDisplayNode()
    #disp.SetOpacity(1)
    self.DirectionBox.checked = True
    self.UsedTractBox.checked = True
    self.SphereBox.checked = True
    self.EntryBox.checked = True
    self.skullBox.checked = True
    self.onSkullBoxClicked()
    self.oppositeStartButton.enabled = True
    self.PathButton.enabled = True
    self.applyButton.enabled = False
    self.SphereBox.checked = False
    self.onSphereBoxClicked()

  def onResetButton(self):
    fid=True
    self.logic.resetAll(fid, self.inputFiducialsNodeSelector.currentNode())
    self.oppositeStartButton.enabled = False
    self.PathButton.enabled = False
    self.applyButton.enabled = True
    self.skullBox.checked = False
    self.onskullBoxClicked()
    mod = slicer.util.getNode('skull')
    disp = mod.GetDisplayNode()
    disp.SetOpacity(0.4)
    print("Test test test reset button")


  def onOppositeStartButton(self):
    fiducials = self.inputFiducialsNodeSelector.currentNode()
    if fiducials.GetNumberOfFiducials() == 2:
      fiducials.RemoveMarkup(1)

    ok = False
    while ok == False:
      model = slicer.util.getNode('acceptedRegion_poly')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    self.logic.oppositeEntry(self.inputFiducialsNodeSelector.currentNode())




  def onPathButton(self):
    fiducials= self.inputFiducialsNodeSelector.currentNode()
    target = fiducials.GetMarkupPointVector(0,0)
    start = fiducials.GetMarkupPointVector(1,0)
    fiducials.SetNthFiducialLabel(1, "Start")
    self.logic.drawLine(start, target, [0, 0, 1], 'direction_line', 1.25, 1)
    self.PathBox.checked = True
    mod = slicer.util.getNode('skull')
    disp = mod.GetDisplayNode()
    disp.SetOpacity(1)
    self.oppositeStartButton.enabled = False

#
# EDEN_tract_plannerLogic
#

class EDEN_tract_plannerLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    self.OutputModelRed = None
    self.OutputModelBlue = None
    self.OutputModelGreen = None

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qpixMap = qt.QPixmap().grabWidget(widget)
    qimage = qpixMap.toImage()
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def plotSphere(self, fiducial, slidersphere):

    while slicer.util.getNode('SphereModelNode') is not None:
        node = slicer.util.getNode('SphereModelNode')
        slicer.mrmlScene.RemoveNode(node)

    selectedPt = [0, 0, 0]
    fiducial.GetNthFiducialPosition(0, selectedPt)

    modelNode = slicer.vtkMRMLModelNode()
    dispNode = slicer.vtkMRMLModelDisplayNode()
    transform = slicer.vtkMRMLLinearTransformNode()

    # Display node characteristics
    dispNode.SetVisibility(True)
    dispNode.SetSliceIntersectionVisibility(True)
    dispNode.SetOpacity(0.6)
    dispNode.SetColor(1, 1, 0)
    dispNode.SetScene(slicer.mrmlScene)

    # generate sphere data
    sphere = vtk.vtkSphereSource()
    sphere.SetCenter(10, 10, 10)
    sphere.SetRadius(slidersphere)
    sphere.SetCenter(selectedPt[0], selectedPt[1], selectedPt[2])
    sphere.Update()

    # adding necessary nodes to the Scene
    slicer.mrmlScene.AddNode(dispNode)
    slicer.mrmlScene.AddNode(transform)
    slicer.mrmlScene.AddNode(modelNode)

    # model node name and associations!
    modelNode.SetName("SphereModelNode")
    modelNode.SetScene(slicer.mrmlScene)
    modelNode.SetAndObserveTransformNodeID(transform.GetID())
    modelNode.SetAndObserveDisplayNodeID(dispNode.GetID())

    apd = vtk.vtkAppendPolyData()
    apd.AddInputData(sphere.GetOutput())
    apd.Update()

    # adding model node poly data! Here there are  sphere's data
    modelNode.SetAndObservePolyData(apd.GetOutput())
    # update the scene
    slicer.mrmlScene.Modified()


  def run(self, tractography_data, skull, sphereRadius, Fiducials):
    """
    Run the actual algorithm
    """
    sourceVolumeNode = tractography_data
    Fiducials.SetNthFiducialLabel(0, "Target")
    F1 = Fiducials.GetMarkupPointVector(0, 0)  # Point on the Fiber
    print(F1)
    print("line 542")
    F1 = [F1[0], F1[1], F1[2]]
    print(F1)
    RadiusSphere = sphereRadius
    RadiusNotValid = True
    PointNotGood = False
    # print("test tests test from 'run' command.")
    # your_mesh = mesh.Mesh.from_file(os.path.dirname(__file__) + r'/Resources/skull.vtk')
    # print(your_mesh)


    while RadiusNotValid:

      # Sphere definition (creation of the main direction of the fiber)
      sphere = vtk.vtkSphere()
      sphere.SetCenter(F1[0], F1[1], F1[2])
      sphere.SetRadius(RadiusSphere)

      # Clip tractography data
      clip = vtk.vtkClipPolyData()
      clip.SetInputData(sourceVolumeNode.GetPolyData())
      clip.SetClipFunction(sphere)
      clip.InsideOutOn()
      clip.Update()
      clipped_old = clip.GetOutput()

      # Algotirhm to find the main direction
      fid_pos = np.array([[F1[0]], [F1[1]], [F1[2]]])

      points = clipped_old.GetPoints()
      points = points.GetData()

      new_temporary_pointArray = vtk.vtkFloatArray()
      new_temporary_pointArray.SetNumberOfComponents(3)
      new_temporary_pointArray.DeepCopy(points)

      np_tensors = vtk.util.numpy_support.vtk_to_numpy(new_temporary_pointArray)

      if len(
              np_tensors) == 0:  # ALI sometimes it is necessary that the radius is bigger then 0.5, for example on the tumor so we increment it
        print "the second fiducual is not positioned on the fiber"
        RadiusNotValid = False
        PointNotGood = True
      else:
        mean_x = np.sum(np_tensors[:, 0]) / np_tensors.shape[0]
        mean_y = np.sum(np_tensors[:, 1]) / np_tensors.shape[0]
        mean_z = np.sum(np_tensors[:, 2]) / np_tensors.shape[0]

        np_tensors[:, 0] = np_tensors[:, 0] - mean_x
        np_tensors[:, 1] = np_tensors[:, 1] - mean_y
        np_tensors[:, 2] = np_tensors[:, 2] - mean_z

        covariance = np.cov(np_tensors.transpose())
        eigval, eigvect = np.linalg.eig(covariance)

        # ALI this part is added to be sure that the max eigValue is >> then the others, so that the result is reliable.
        maxV = max(eigval)
        minV = min(eigval)
        midV = [i for i, j in enumerate(eigval) if j != maxV and j != minV]
        a = maxV
        b = eigval[midV[0]]
        threshold = (abs(a - b) / ((abs(a) + abs(b)) / 2.0)) * 100.0

        if threshold < 60 and RadiusSphere >= 1:
          RadiusSphere = RadiusSphere - 0.5
        elif threshold < 60 and RadiusSphere < 1:
          RadiusNotValid = False
          PointNotGood = True
        else:
          RadiusNotValid = False

    if PointNotGood == True:
      print("select another point")
    else:
      # creating the used tract model
      tuber = vtk.vtkTubeFilter()
      tuber.SetInputData(clipped_old)
      tuber.SetRadius(0.5)
      tuber.SetNumberOfSides(5)
      tuber.Update()
      tubes = tuber.GetOutputDataObject(0)
      scene = slicer.mrmlScene
      modelt = slicer.vtkMRMLModelNode()
      modelt.SetScene(scene)
      modelt.SetName("Used Tract")
      modelt.SetAndObservePolyData(tubes)

      modelDisplayt = slicer.vtkMRMLModelDisplayNode()
      modelDisplayt.SetColor(0.8, 1, 1)
      modelDisplayt.SetBackfaceCulling(0)
      modelDisplayt.SetScene(scene)
      scene.AddNode(modelDisplayt)
      modelt.SetAndObserveDisplayNodeID(modelDisplayt.GetID())
      scene.AddNode(modelt)

      logging.info('Processing started')

      print("RadiusSphere", RadiusSphere)
      print("tot punti", len(np_tensors))
      print ("autovalori", eigval)
      print ("reliability of:", int(100 - (threshold / 2)), '%')
      idMaxEigValue = [i for i, j in enumerate(eigval) if j == maxV]
      idMinEigValue = [i for i, j in enumerate(eigval) if j == minV]
      # Comp1 is the main direction of the fiber
      comp1 = [eigvect[0, idMaxEigValue[0]], eigvect[1, idMaxEigValue[0]], eigvect[2, idMaxEigValue[0]]]
      comp2 = [eigvect[0, idMinEigValue[0]], eigvect[1, idMinEigValue[0]], eigvect[2, idMinEigValue[0]]]
      # to plot comp1
      vpm.plot_vector(comp1, 1, [1, 0, 0], fid_pos, 10, 'RedVector', 1)

      self.Accepted_points(Fiducials, skull, F1, comp1)
      self.drawPlane(comp1, F1)
      # to see the plane from all the views
      #comp1_flipped = (comp1[0], -comp1[1], -comp1[2]);
      #self.drawPlane(comp1_flipped, F1)

      your_mesh = mesh.Mesh.from_file(os.path.dirname(__file__) + r'\Resources\Skull_rotated_hardened_top_decimated_transf.stl')        #From skull file, import matrix of vertices and normals.
      target_direction = np.array(comp1) * -1                                                                                           #Invert vector of target point. Save as 'target direction.'
      target_point = np.array([ fid_pos[0],fid_pos[1],fid_pos[2], 1.0 ])                                                                #Create 1x4 vector of target point.
      print(target_point)
      norm_target_direction = target_direction/np.linalg.norm(target_direction)                                                         #Normalise 'target_direction' vector.
      normals = your_mesh.normals * -1                                                                                                  #Invert normals on skull.
      length = len(your_mesh.normals)
      ave_coord = np.zeros([length,3])                                                                                                  #Creates matrix of 0s to save average coordinates of triangles.
      ref_axis = np.array([1,0,0])                                                                                                      #Set axis of reference as the x-axis.
      step = 1                                                                                                                          #Set iteration step of 'for' loop = 1
      count = 0                                                                                                                         # Set count = 0
      R = 70                                                                                                                            #Set Radius of curvature = 70

      for i in range(0,length,step):
          #Calculates the average x,y,z component of the face to find central coordinates.
          ave_coord[i,0] = (your_mesh.vectors[i,0][0] + your_mesh.vectors[i,1][0] + your_mesh.vectors[i,2][0]) / 3     #Compute centre of triangle by averaging vertices.
          ave_coord[i,1] = (your_mesh.vectors[i,0][1] + your_mesh.vectors[i,1][1] + your_mesh.vectors[i,2][1]) / 3
          ave_coord[i,2] = (your_mesh.vectors[i,0][2] + your_mesh.vectors[i,1][2] + your_mesh.vectors[i,2][2]) / 3
          start_point = np.array([[ave_coord[i,0]],[ave_coord[i,1]],[ave_coord[i,2]],[1.0]])                           #Add 4th element = 1 (for matrix multiplication later on).
          start_direction = normals[i]                                                                                    #Save normal vector for face as 'start_direction'
          norm_start_direction = start_direction / np.linalg.norm(start_direction)                                     #Computes normalised start_direction vector
          v_angle1 = np.arccos(np.transpose(ref_axis).dot(norm_start_direction))                                    #Computes angle between norm_start_direction vector and reference axis.
          axis = np.cross(ref_axis,norm_start_direction) / np.linalg.norm(np.cross(ref_axis,norm_start_direction))      #Computes axis by which to apply quaternion function.
          q1 = Quaternion(axis = axis , angle = v_angle1)                                                                 #Compute quaternion given angle of rotation and new axis.
          T1 = q1.transformation_matrix                                                                                  #Turn T1 into a 4x4 transformation matrix.
          T1[:,3] = np.transpose(start_point)                                                                           #Make the last column equal to the transposed "start_point" vector.
          T1_inv = np.linalg.inv(T1)                                                                                   #Compute inverse of T1
          target_point_inverted = T1_inv.dot(target_point)                                                                #Matrix multiplication of T1_inv and target_point
          xx = np.array([target_point_inverted[0],target_point_inverted[1],target_point_inverted[2]])                  #Take 3 first elements of target_point_inverted and rename as xx.
          flag1 = False                                                                                                   #Set flag1 as False.
          if xx[0] > 0:                                                                                                   #Is greater than 0 to take the point inside the skull, not out.
              if (R - (xx[1]**2 + xx[2]**2)**0.5)**2 + xx[0]**2 > R**2:                                                   #if point is within torous.
                  flag1 = True                                                                                            #Set flag1 = 1

          v_angle2 = np.arccos(np.transpose(ref_axis).dot(norm_target_direction))                                         #Compute angle between reference axis and normalised target direction vector.
          axis = np.cross(ref_axis,norm_target_direction) / np.linalg.norm(np.cross(ref_axis,norm_target_direction))
          q2 = Quaternion(axis = axis , angle = v_angle2)                                                                 #Compute quaternion given angle of rotation and reference axis.
          T2 = q2.transformation_matrix                                                                                   #Turn T2 into a 4x4 transformation matrix.
          T2[:,3] = np.transpose(target_point)                                                                           #Set last column of T2 equal to transpose of "target_point" vector.
          T2_inv = np.linalg.inv(T2)                                                                                   #Compute inverse of T2
          start_point_inverted = T2_inv.dot(start_point)                                                                  #Matrix multiplication of T2_inv and start_point.
          xx = np.array([start_point_inverted[0],start_point_inverted[1],start_point_inverted[2]])                     #Take 3 first elements of start_point_inverted and rename as xx.
          flag2 = False                                                                                                   #Set flag2 as False.
          if xx[0] < 0:                                                                                                   #Is greater than 0 to take the point inside the skull, not out.
              if (R - (xx[1]**2 + xx[2]**2)**0.5)**2 + xx[0]**2 > R**2:                                                    #if point is within torous.
                  flag2 = True                                                                                            #Set flag2 = 1

          if flag1 and flag2 == True:                                                                                     #If both conditions are satisfied then plot vector on skull.
              ##plot_vector(vector, axisDiameter, color, translation=None, axisLength=None, Name=None, showTip=None)
              vpm.plot_vector(start_direction, 1.5, [0,0,1], start_point, 15.0, None, 1)
              count = count + 1                                                                                           #Increment 'count' to count total number of plotted vectors.

      print('count equals to:' ,count)                                                                                    #Display total number of plotted vectors.

  def drawPlane(self, comp1, F1):

      # to calculate the plane
      self.Direction = [comp1[0], comp1[1], comp1[2]]
      self.OutputModelRed = slicer.util.getNode('RedVector')

      # Square Plane definition
      self.plane_source = vtk.vtkPlaneSource()
      self.plane_source.SetCenter(F1[0], F1[1], F1[2])
      self.plane_source.SetNormal(self.Direction[0], self.Direction[1], self.Direction[2])
      self.plane_source.Update()

      dist = math.sqrt(8000)
      O = self.plane_source.GetOrigin()
      P1 = self.plane_source.GetPoint1()
      P2 = self.plane_source.GetPoint2()

      OF = [O[0] - F1[0], O[1] - F1[1], O[2] - F1[2]]
      norm = math.sqrt(OF[0] * OF[0] + OF[1] * OF[1] + OF[2] * OF[2])
      OF = [OF[0] * dist / norm, OF[1] * dist / norm, OF[2] * dist / norm]
      O = [OF[0] + F1[0], OF[1] + F1[1], OF[2] + F1[2]]

      PPF = [P2[0] - F1[0], P2[1] - F1[1], P2[2] - F1[2]]
      PPF = [PPF[0] * dist / norm, PPF[1] * dist / norm, PPF[2] * dist / norm]
      P2 = [PPF[0] + F1[0], PPF[1] + F1[1], PPF[2] + F1[2]]

      PF = [P1[0] - F1[0], P1[1] - F1[1], P1[2] - F1[2]]
      PF = [PF[0] * dist / norm, PF[1] * dist / norm, PF[2] * dist / norm]
      P1 = [PF[0] + F1[0], PF[1] + F1[1], PF[2] + F1[2]]

      self.plane_source.SetOrigin(O)
      self.plane_source.SetPoint1(P1)
      self.plane_source.SetPoint2(P2)
      self.plane_source.SetResolution(200,200);
      self.plane_source.Update()

      ####  PLANE

      modelNode = slicer.vtkMRMLModelNode()
      dispNode = slicer.vtkMRMLModelDisplayNode()
      transform = slicer.vtkMRMLLinearTransformNode()
      # Display node characteristics
      dispNode.SetVisibility(True)
      dispNode.SetSliceIntersectionVisibility(True)
      dispNode.SetOpacity(0.8)
      dispNode.SetColor(1, 1, 0)
      dispNode.SetScene(slicer.mrmlScene)

      # adding necessary nodes to the Scene
      slicer.mrmlScene.AddNode(dispNode)
      slicer.mrmlScene.AddNode(transform)
      slicer.mrmlScene.AddNode(modelNode)

      # model node name and associations!
      modelNode.SetName("PlaneModelNode")
      modelNode.SetScene(slicer.mrmlScene)
      modelNode.SetAndObserveTransformNodeID(transform.GetID())
      modelNode.SetAndObserveDisplayNodeID(dispNode.GetID())

      apd = vtk.vtkAppendPolyData()
      apd.AddInputData(self.plane_source.GetOutput())
      apd.Update()

      # adding model node poly data! Here there are  sphere's data
      modelNode.SetAndObservePolyData(apd.GetOutput())
      # update the scene
      slicer.mrmlScene.Modified()

  def Accepted_points(self, Fiducials, Accepted_skull, F1, comp1):

    Accepted_skull_poly = Accepted_skull.GetPolyData()

    cylinder = vtk.vtkCylinder()
    cylinder.SetCenter(F1[0], F1[1], F1[2])
    cylinder.SetRadius(2.5)
    cylinder.SetAxis(comp1)

    clip = vtk.vtkClipPolyData()
    clip.SetInputData(Accepted_skull_poly)
    clip.SetClipFunction(cylinder)
    clip.InsideOutOn()
    clip.Update()
    clipped_withCylinder = clip.GetOutput()

    pointsNode = clipped_withCylinder.GetPoints()
    pointsNode = pointsNode.GetData()
    temporary_pointArray_ofNode = vtk.vtkFloatArray()
    temporary_pointArray_ofNode.SetNumberOfComponents(3)
    temporary_pointArray_ofNode.DeepCopy(pointsNode)
    pointArray_ofNode1 = vtk.util.numpy_support.vtk_to_numpy(temporary_pointArray_ofNode)

    if pointArray_ofNode1.size==0:
      print "WARNING: No intersection with the skull available. Change target location"
      return

    space=200
    TempPoint = ([F1[0] + (comp1[0] *-1* space), F1[1] + (comp1[1] *-1*space), F1[2] + (comp1[2] * -1*space)])
    tempArray = []
    for p in pointArray_ofNode1:
      distPoint1 = math.sqrt((TempPoint[0] - p[0]) ** 2 + (TempPoint[1] - p[1]) ** 2 + (TempPoint[2] - p[2]) ** 2)
      tempArray.append(distPoint1)
    maxPoint1 = max(tempArray)
    maxPoint2 = min(tempArray)

    P1 = ([TempPoint[0] + (comp1[0] * maxPoint1), TempPoint[1] + (comp1[1] * maxPoint1), TempPoint[2] + (comp1[2] * maxPoint1)])

    P2 = ([TempPoint[0] + (comp1[0] * maxPoint2), TempPoint[1] + (comp1[1] * maxPoint2), TempPoint[2] + (comp1[2] * maxPoint2)])

    dist1 = math.sqrt((F1[0] - P1[0]) ** 2 + (F1[1] - P1[1]) ** 2 + (F1[2] - P1[2]) ** 2)
    dist2 = math.sqrt((F1[0] - P2[0]) ** 2 + (F1[1] - P2[1]) ** 2 + (F1[2] - P2[2]) ** 2)

    global F2, F3

    if dist1 <= dist2:
      F2 = P1
      F3 = P2
    else:
      F2 = P2
      F3 = P1

    global entryNow
    entryNow= 'F2'

    slicer.modules.markups.logic().AddFiducial(F2[0], F2[1], F2[2])
    Fiducials.SetNthFiducialLabel(1, "Start")

    # Sphere definition (creation of the entry Area)
    sphereSurface = vtk.vtkSphere()
    sphereSurface.SetCenter(F2[0], F2[1], F2[2])
    sphereSurface.SetRadius(2.5)

    # Clip tractography data
    clip = vtk.vtkClipPolyData()
    clip.SetInputData(Accepted_skull_poly)
    clip.SetClipFunction(sphereSurface)
    clip.InsideOutOn()
    clip.Update()
    clipped_Surface = clip.GetOutput()

    #green Area
    sceneAcc = slicer.mrmlScene
    modelAcc = slicer.vtkMRMLModelNode()
    modelAcc.SetScene(sceneAcc)
    modelAcc.SetName("acceptedRegion_poly")
    modelAcc.SetAndObservePolyData(clipped_Surface)
    modelDisplayAcc = slicer.vtkMRMLModelDisplayNode()
    modelDisplayAcc.SetColor(0, 1, 0)
    modelDisplayAcc.SetOpacity(0.9)
    modelDisplayAcc.SetBackfaceCulling(0)
    modelDisplayAcc.SetScene(sceneAcc)
    sceneAcc.AddNode(modelDisplayAcc)
    modelAcc.SetAndObserveDisplayNodeID(modelDisplayAcc.GetID())
    sceneAcc.AddNode(modelAcc)


  def drawLine(self, p1, p2, color, name, radius, opacity):
    line = vtk.vtkLineSource()
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    line.Update()

    tuber = vtk.vtkTubeFilter()
    tuber.SetInputData(line.GetOutput())
    tuber.SetRadius(radius)
    tuber.SetNumberOfSides(15)
    tuber.Update()
    tube = tuber.GetOutputDataObject(0)
    scene = slicer.mrmlScene
    modelt = slicer.vtkMRMLModelNode()
    modelt.SetScene(scene)
    modelt.SetName(name)
    modelt.SetAndObservePolyData(tube)

    modelDisplayt = slicer.vtkMRMLModelDisplayNode()
    modelDisplayt.SetColor(color)
    modelDisplayt.SetOpacity(opacity)
    modelDisplayt.SetBackfaceCulling(0)
    modelDisplayt.SetScene(scene)
    scene.AddNode(modelDisplayt)
    modelt.SetAndObserveDisplayNodeID(modelDisplayt.GetID())
    scene.AddNode(modelt)

  def resetAll(self, fid, fiducials):

    ok = False
    while ok == False:
      model = slicer.util.getNode('RedVector')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    ok = False
    while ok == False:
      model = slicer.util.getNode('PlaneModelNode')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    ok = False
    while ok == False:
      model = slicer.util.getNode('acceptedRegion_poly')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    ok = False
    while ok == False:
      model = slicer.util.getNode('direction_line')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    ok = False
    while ok == False:
      model = slicer.util.getNode('Used Tract')
      if model != None:
        slicer.mrmlScene.RemoveNode(model)
      else:
        ok = True

    if fiducials.GetNumberOfFiducials()>0:
      fiducials.RemoveMarkup(1)

    if fid==True:
      sphere = slicer.util.getNode('SphereModelNode')
      slicer.mrmlScene.RemoveNode(sphere)

      while slicer.util.getNode('F') is not None:
        node = slicer.util.getNode('F')
        slicer.mrmlScene.RemoveNode(node)

  def oppositeEntry(self,Fiducials):

    if entryNow=='F2':
      F=F3
      global entryNow
      entryNow = 'F3'
    else:
      F=F2
      global entryNow
      entryNow = 'F2'

    Accepted_skull = slicer.util.getNode('skull')
    Accepted_skull_poly = Accepted_skull.GetPolyData()

    slicer.modules.markups.logic().AddFiducial(F[0], F[1], F[2])
    Fiducials.SetNthFiducialLabel(1, "Start")

    # Sphere definition (creation of the entry Area)
    sphereSurface = vtk.vtkSphere()
    sphereSurface.SetCenter(F[0], F[1], F[2])
    sphereSurface.SetRadius(2.5)

    # Clip tractography data
    clip = vtk.vtkClipPolyData()
    clip.SetInputData(Accepted_skull_poly)
    clip.SetClipFunction(sphereSurface)
    clip.InsideOutOn()
    clip.Update()
    clipped_Surface = clip.GetOutput()

    # green Area
    sceneAcc = slicer.mrmlScene
    modelAcc = slicer.vtkMRMLModelNode()
    modelAcc.SetScene(sceneAcc)
    modelAcc.SetName("acceptedRegion_poly")
    modelAcc.SetAndObservePolyData(clipped_Surface)
    modelDisplayAcc = slicer.vtkMRMLModelDisplayNode()
    modelDisplayAcc.SetColor(0, 1, 0)
    modelDisplayAcc.SetOpacity(0.9)
    modelDisplayAcc.SetBackfaceCulling(0)
    modelDisplayAcc.SetScene(sceneAcc)
    sceneAcc.AddNode(modelDisplayAcc)
    modelAcc.SetAndObserveDisplayNodeID(modelDisplayAcc.GetID())
    sceneAcc.AddNode(modelAcc)



class EDEN_tract_plannerTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_EDEN_tract_planner1()

  def test_EDEN_tract_planner1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #

    slicer.util.loadVolume('/Users/albertofavaro/Desktop/clinical dataset/EDEN Dataset/PTNT01/20170323_USR_MRI_3T_Philips/Tractography/3DT1s601a1006.nii', returnNode=True)
    # slicer.modules.models.logic().AddModel('/Users/albertofavaro/Desktop/Model_1.vtk')

    # clippingModel=slicer.util.getNode('Model_1')
    # volume = slicer.util.getNode('3DT1s601a1006')


    self.delayDisplay('Test passed!')
