double cd138Threshold = 100.0
String cd138Measurement = "Cell: CD138 mean"
String outputMeasurement = "Nearest CD138 distance (µm)"

def imageData = getCurrentImageData()
def hierarchy = imageData.getHierarchy()

// Get cells
def cells = getCellObjects()

println "Total cells: ${cells.size()}"

// Select CD138+ cells
def cd138Cells = cells.findAll { cell ->
    def value = cell.getMeasurementList().get(cd138Measurement)
    value != null && !Double.isNaN(value) && value > cd138Threshold
}

println "CD138+ cells: ${cd138Cells.size()}"

if (cd138Cells.size() < 2) {
    println "Not enough CD138+ cells"
    return
}

// Coordinates
def coords = cd138Cells.collect { cell ->
    def roi = cell.getROI()
    [
        cell: cell,
        x: roi.getCentroidX(),
        y: roi.getCentroidY()
    ]
}

// Calibration
def cal = imageData.getServer().getPixelCalibration()
double pixelSize = (cal.getPixelWidthMicrons() + cal.getPixelHeightMicrons()) / 2.0

// Nearest neighbour calculation
coords.eachWithIndex { p, i ->

    double minDist = Double.MAX_VALUE

    coords.eachWithIndex { q, j ->

        if (i != j) {

            double dx = p.x - q.x
            double dy = p.y - q.y
            double d = Math.sqrt(dx*dx + dy*dy)

            if (d < minDist)
                minDist = d
        }
    }

    p.cell.getMeasurementList()
        .put(outputMeasurement, minDist * pixelSize)
}

fireHierarchyUpdate()

println "Done"