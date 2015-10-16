using Types

mySettings = UserSettings()

mySettings.nIterations = 100

mySettings.meshFileName = "../input/sphere2.ply"

mySettings.domainSizeX = 200.0
mySettings.domainSizeY = 200.0
mySettings.domainSizeZ = 200.0

mySettings.nCellsPerBlockX = 5
mySettings.nCellsPerBlockY = 5
mySettings.nCellsPerBlockZ = 5

mySettings.nMaxRefinementLevel = 3
mySettings.nMaxTrianglesPerCell = 1
mySettings.nNewParticlesPerIteration = 3500
