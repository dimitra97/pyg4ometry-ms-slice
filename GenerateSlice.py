# Load pyg4ometry
import json
import pyg4ometry
import pyg4ometry.geant4 as pyg4
import math
import numpy

# overlap protection
ovt = 0.001

x = [1, 0, 0]
y = [0, 1, 0]
z = [0, 0, 1]
idMat = numpy.matrix([x, y, z]) 
idRot = pyg4ometry.transformation.matrix2tbxyz(idMat)
nSectorsLarge=8
print(idMat)
print(idRot)

""" create mulit-wire chamber logical volume  

paranmeters are : name (of the chamber), reg (registry), mother (volume), 
                  rWire (wire radius), rTube (tube radius), hlTube (half length of tube)
                  nTubes (number of tubes), nLayers (number of mlayer in a multi layer)
                  nMultiLayers (number of multi layers per station)
                  yGap the gap betweeb tge multi layers
    
"""
def createMulitTubeChamber(name, reg, mother, rWire=1, rTube = 20, tTube = 0.1, hlTube = 350, nTubes = 40, nLayers = 4, nMulitLayers = 2, yGap = 300, structy=3) :
   
    rPack = rTube * math.sqrt(3);
    # Create the tube solids 
    wireSolid = pyg4ometry.geant4.solid.Tubs(name+'WireSolid', 0, rWire, 2*hlTube-ovt, 0., 2*math.pi, reg)
    gasSolid =  pyg4ometry.geant4.solid.Tubs(name+'GasSolid', rWire+ovt, rTube-tTube-ovt, 2*hlTube, 0., 2*math.pi, reg)
    skinSolid = pyg4ometry.geant4.solid.Tubs(name+'SkinSolid', rTube-tTube, rTube, 2*hlTube, 0., 2*math.pi, reg)
   
    # Create the multi layer solids
    #put the tubes along x-direction like ATLAS MS
    #the tubes are perpendicual to beam-direction (z-axis)
    mlZ=  (2 * nTubes + 1) * rTube
    mlY =  2 * rTube + (nLayers-1) * rPack
    print(mlY)
    mlX=2*(hlTube+ovt)
    mlSolid = pyg4ometry.geant4.solid.Box(name+'MultiLayerSolid', mlX, mlY,mlZ, reg)
    # Create the chamber solids
    chY = 2 * mlY + yGap+2*structy
    print(chY)
    chSolid = pyg4ometry.geant4.solid.Box(name+'ChamberSolid', mlX + 2 * ovt, chY + 2 * ovt, mlZ + 2 * ovt, reg)
      
    # Logical volumes for the tube
    Air = pyg4ometry.geant4.nist_material_2geant4Material('G4_AIR')    
    Cu = pyg4ometry.geant4.nist_material_2geant4Material('G4_Cu')
    Ar = pyg4ometry.geant4.nist_material_2geant4Material('G4_Ar')
    Al = pyg4ometry.geant4.nist_material_2geant4Material('G4_Al')
    CO2 = pyg4ometry.geant4.nist_material_2geant4Material('G4_CARBON_DIOXIDE')
    gas_mixture = pyg4ometry.geant4.MaterialCompound('ArCO2',1.8223e-03,2)
    gas_mixture.add_material(Ar,0.9)
    gas_mixture.add_material(CO2,0.1)
    W = pyg4ometry.geant4.nist_material_2geant4Material('G4_W')
    Re = pyg4ometry.geant4.nist_material_2geant4Material('G4_Re')
    wire_material = pyg4ometry.geant4.MaterialCompound('W/Re',19.3,2)
    wire_material.add_material(W,0.97)
    wire_material.add_material(Re,0.03)    
    wireLog = pyg4ometry.geant4.LogicalVolume(wireSolid, wire_material, name+'WireLog', reg)
    gasLog =  pyg4ometry.geant4.LogicalVolume(gasSolid, gas_mixture, name+'GasLog', reg)
    skinLog =  pyg4ometry.geant4.LogicalVolume(skinSolid, Al, name+'SkinLog', reg)
    
    # Logical volume for the multi layer
    mlLog = pyg4ometry.geant4.LogicalVolume(mlSolid, Air, name+'MultiLayerLog', reg)
    if (nLayers%2) == 0 :
      cY = -(nLayers/2 - 0.5) * rPack
    else :
      cY = -(nLayers-1)/2 * rPack  
    # place the volumes
    tX = 0
    for layer in range(nLayers) : 
      tY = cY + layer * rTube * math.sqrt(3)
      offsetZ = (layer%2) * rTube
      for tube in range(nTubes) :          
          tZ = -0.5 * mlZ + (2 * tube + 1) * rTube + offsetZ
          # Place the wire, gas, skin
          pos = [ tX, tY, tZ ]
          print("tube:")
          print(tube)
          print(pos)
          wirePhys = pyg4ometry.geant4.PhysicalVolume([0,math.pi/2,0], pos, wireLog, name+'Wire_l'+str(layer)+'_t'+str(tube), mlLog, reg)
          gasPhys = pyg4ometry.geant4.PhysicalVolume([0,math.pi/2,0], pos, gasLog, name+'Gas_l'+str(layer)+'_t'+str(tube), mlLog, reg)
          skinPhys = pyg4ometry.geant4.PhysicalVolume([0,math.pi/2,0], pos, skinLog, name+'Skin_l'+str(layer)+'_t'+str(tube), mlLog, reg)
                    
    # change log 
    chLog = pyg4ometry.geant4.LogicalVolume(chSolid, Air, name+'ChamberLog', reg)

    # place the multi layers twice
    mllPhys = pyg4ometry.geant4.PhysicalVolume(idRot, [0, cY-0.5*yGap, 0], mlLog, name+'MuliWire_l', chLog, reg)
    mlhPhys = pyg4ometry.geant4.PhysicalVolume(idRot, [0, -cY+0.5*yGap, 0], mlLog, name+'MuliWire_h', chLog, reg)    
    
    return chLog
    
        
# registry to store gdml data
reg = pyg4ometry.geant4.Registry()

# world solid and logical
world_box = pyg4ometry.geant4.solid.Box("world_box", 50000, 50000, 50000, reg)
Air = pyg4ometry.geant4.nist_material_2geant4Material('G4_AIR')
world = pyg4ometry.geant4.LogicalVolume(world_box, Air, "world", reg)


# Create 3 logical volumes
#inner chamber
chInnerLog  = createMulitTubeChamber("Inner_", reg, world, rWire=0.025, rTube = 15, tTube = 0.4, hlTube = 1350, nTubes = 36, nLayers = 4, nMulitLayers = 2, yGap = 170, structy=15)
#middle chamber
chMiddleLog = createMulitTubeChamber("Middle_" ,reg, world, rWire=0.025, rTube = 15, tTube = 0.4, hlTube = 1790, nTubes = 56, nLayers = 3, nMulitLayers = 2, yGap = 317, structy=15)
#outer chamber
chOuterLog  = createMulitTubeChamber("Outer_", reg, world, rWire=0.025, rTube = 15, tTube = 0.4, hlTube = 2495, nTubes = 72, nLayers = 3, nMulitLayers = 2, yGap = 317, structy=15)

# The radial position of the stations
innerR = 4949
middleR = 7139 
outerR = 9500
chInnerPhys = pyg4ometry.geant4.PhysicalVolume(idRot, [0, innerR, 0], chInnerLog, 'InnerChamber', world, reg)    
chMiddlePhys = pyg4ometry.geant4.PhysicalVolume(idRot, [0, middleR, 0], chMiddleLog, 'MiddleChamber', world, reg)    
chOuterPhys = pyg4ometry.geant4.PhysicalVolume(idRot, [0, outerR, 0], chOuterLog, 'OuterChamber', world, reg)    

# Visualiztion options
voAIR = pyg4ometry.visualisation.VisualisationOptions()
voAIR.visible = False

voCu = pyg4ometry.visualisation.VisualisationOptions()
voCu.colour = [0.3, 0.5, 1.0]
voCu.alpha = 1.0

v = pyg4ometry.visualisation.VtkViewerColoured(defaultColour="random")
v.addMaterialVisOption("G4_AIR", voAIR)

# Obj write out
v.addLogicalVolume(world)
v.exportOBJScene("MuonChamber")
v.addAxes(2000)
v.view()


# GDML wirte out
reg.setWorld(world.name)
w = pyg4ometry.gdml.Writer()
w.addDetector(reg)
w.write('MuonChamber.gdml')
