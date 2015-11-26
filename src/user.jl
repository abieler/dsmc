using Types

mySettings = UserSettings()
body = MeshBody

mySettings.nIterations = 50
mySettings.meshFileName = "../input/sphere2.ply"

mySettings.domainSizeX = 5000.0
mySettings.domainSizeY = 5000.0
mySettings.domainSizeZ = 5000.0

mySettings.nCellsPerBlockX = 5
mySettings.nCellsPerBlockY = 5
mySettings.nCellsPerBlockZ = 5

mySettings.nMaxRefinementLevel = 5
mySettings.nMaxTrianglesPerCell = 1
mySettings.nNewParticlesPerIteration = 90000
