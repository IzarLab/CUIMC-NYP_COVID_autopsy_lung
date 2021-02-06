setImageType('FLUORESCENCE');
setChannelNames(
     'DAPI',
     'CD4',
     'CD19',
     'GzmB',
     'CD103',
     'CD8',
     'CD163',
     'AF')
createSelectAllObject(true)
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "DAPI",  "requestedPixelSizeMicrons": 0.5,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 1.5,  "watershedPostProcess": true,  "cellExpansionMicrons": 4.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');
