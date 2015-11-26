using Types

mySettings = UserSettings()
body = MeshBody

mySettings.nIterations = 30
mySettings.meshFileName = "../input/sphere2.ply"

mySettings.domainSizeX = 50000.0
mySettings.domainSizeY = 50000.0
mySettings.domainSizeZ = 50000.0

mySettings.nCellsPerBlockX = 5
mySettings.nCellsPerBlockY = 5
mySettings.nCellsPerBlockZ = 5

mySettings.nMaxRefinementLevel = 11
mySettings.nMaxTrianglesPerCell = 1
mySettings.nNewParticlesPerIteration = 4000000
