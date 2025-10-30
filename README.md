// TerraBloom Source Code

var index_palettes = require('users/gena/packages:palettes');
var jrc = ee.Image("JRC/GSW1_4/GlobalSurfaceWater");
var waterOccurrence = jrc.select('occurrence');
var permanentWater = waterOccurrence.gt(10);
var s2Sr = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');

function getClearSentinel2Composite(aoi, startDate, endDate) {
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));
  
  function maskClouds(image) {
    var scl = image.select('SCL');
    var qa60 = image.select('QA60');
    
    var cloudBitMask = ee.Number(2).pow(10).int();
    var cirrusBitMask = ee.Number(2).pow(11).int();
    var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0)
                     .and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
    
    var sclMask = scl.neq(3)
                    .and(scl.neq(7))
                    .and(scl.neq(8))
                    .and(scl.neq(9))
                    .and(scl.neq(10))
                    .and(scl.neq(11));
    
    var water = scl.eq(6);
    var combinedMask = qaMask.and(sclMask).and(water);
    return image.updateMask(combinedMask);
  }
  
  var s2Masked = s2.map(maskClouds);
  var composite = s2Masked.median().clip(aoi);
  var filled = composite.unmask(composite.focal_mean(3, 'square', 'pixels', 1));
  
  return filled;
}

function getClearMODISComposite(aoi, startDate, endDate) {
  var modis = ee.ImageCollection('MODIS/006/MOD09GA')
    .filterBounds(aoi)
    .filterDate(startDate, endDate);
  
  function maskClouds(image) {
    var qa = image.select('state_1km');
    var cloudMask = qa.bitwiseAnd(3).eq(0);
    return image.updateMask(cloudMask);
  }
  
  var modisMasked = modis.map(maskClouds);
  var composite = modisMasked.median().clip(aoi);
  var filled = composite.unmask(composite.focal_mean(3, 'square', 'pixels', 1));
  
  return filled;
}

var leftMap = ui.Map();
var rightMap = ui.Map();
var splitPanel = ui.SplitPanel({firstPanel: leftMap, secondPanel: rightMap, orientation: 'horizontal', wipe: true});
var linkPanel = ui.Map.Linker([leftMap, rightMap]);

var drawingTools = leftMap.drawingTools();
drawingTools.setShown(false);

var aoi = null;
var binaryMap, binaryMap_right, indexMap, indexMap_right, chlaMap, chlaMap_right;
var collection, collection_right;
var x_left = ['2017-08-01', '2017-08-31'];
var x_right = ['2017-08-01', '2017-08-31'];

var visualization_S2 = {min: 0, max: 3000, bands: ['B4', 'B3', 'B2']};
var visualization_L8 = {min: 7272, max: 18182, bands: ['SR_B4', 'SR_B3', 'SR_B2']};
var visualization_L7 = {min: 7272, max: 18182, bands: ['SR_B3', 'SR_B2', 'SR_B1']};
var visualization_L5 = {min: 7272, max: 18182, bands: ['SR_B3', 'SR_B2', 'SR_B1']};
var visualization_VIIRS = {min: 0, max: 0.15, bands: ['M7', 'M4', 'M3']};
var visualization_MODIS = {min: 0, max: 3000, bands: ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03']};

function clearGeometry() {
  while (drawingTools.layers().length() > 0) drawingTools.layers().remove(drawingTools.layers().get(0));
}

function reset(side) {
  if (side === 'left') {
    binaryMap = null;
    indexMap = null;
    chlaMap = null;
  } else {
    binaryMap_right = null;
    indexMap_right = null;
    chlaMap_right = null;
  }
}

function getCollection(dataType, aoi, date) {
  var collection;
  if (dataType === 'S2') {
    collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 95));
  } else if (dataType === 'L8') {
    collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filter(ee.Filter.lt('CLOUD_COVER', 95));
  } else if (dataType === 'L7') {
    collection = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filter(ee.Filter.lt('CLOUD_COVER', 95));
  } else if (dataType === 'L5') {
    collection = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filter(ee.Filter.lt('CLOUD_COVER', 95));
  } else if (dataType === 'VIIRS') {
    collection = ee.ImageCollection('NASA/VIIRS/002/VNP09GA')
      .select(['M1','M2','M3','M4','M5','M7','M8','M10','M11','QF1']);
  } else if (dataType === 'MODIS') {
    collection = ee.ImageCollection('MODIS/006/MOD09GA');
  }
  return collection.filterBounds(aoi).filterDate(date[0], date[1]);
}

function getData(side) {
  var dateRange = (side === 'left') ? x_left : x_right;
  var panel = (side === 'left') ? left.dataDate.data : right.dataDate.data;
  var start = ee.Date(dateRange[0]), end = ee.Date(dateRange[1]);
  if (end.difference(start, 'day').getInfo() <= 0) return null;
  var studyA = (aoi === null || (aoi.type().getInfo() !== 'Polygon' && aoi.type().getInfo() !== 'Point')) ? ee.Geometry.Point([-83.5379, 41.6528]) : aoi;
  return getCollection(panel.getValue() || 'S2', studyA, dateRange);
}

function createWaterMask(image, dataType) {
  var scaled;
  
  if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') {
    scaled = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var green = scaled.select('SR_B3');
    var swir = scaled.select('SR_B6');
    var swir2 = scaled.select('SR_B7');
    var nir = scaled.select('SR_B5');
    var red = scaled.select('SR_B4');
    
    var mndwi = green.subtract(swir).divide(green.add(swir));
    var ndvi = nir.subtract(red).divide(nir.add(red));
    var ndwi = green.subtract(nir).divide(green.add(nir));
    
    var waterMask = mndwi.gt(0.2)
                    .and(ndvi.lt(0.2))
                    .and(ndwi.gt(0.0))
                    .and(swir2.lt(0.15));
    
    return waterMask;
    
  } else if (dataType === 'S2') {
    scaled = image.select(['B.*']).multiply(0.0001);
    var greenS2 = scaled.select('B3');
    var swirS2 = scaled.select('B11');
    var swir2S2 = scaled.select('B12');
    var nirS2 = scaled.select('B8');
    var redS2 = scaled.select('B4');
    
    var mndwiS2 = greenS2.subtract(swirS2).divide(greenS2.add(swirS2));
    var ndviS2 = nirS2.subtract(redS2).divide(nirS2.add(redS2));
    var ndwiS2 = greenS2.subtract(nirS2).divide(greenS2.add(nirS2));
    
    var waterMaskS2 = mndwiS2.gt(0.2)
                      .and(ndviS2.lt(0.2))
                      .and(ndwiS2.gt(0.0))
                      .and(swir2S2.lt(0.15));
    
    return waterMaskS2;
    
  } else if (dataType === 'VIIRS') {
    scaled = image.select(['M1','M2','M3','M4','M5','M7','M8','M10','M11']).multiply(0.0001);
    var greenV = scaled.select('M4');
    var swirV = scaled.select('M11');
    var nirV = scaled.select('M7');
    var redV = scaled.select('M5');
    
    var mndwiV = greenV.subtract(swirV).divide(greenV.add(swirV));
    var ndviV = nirV.subtract(redV).divide(nirV.add(redV));
    var ndwiV = greenV.subtract(nirV).divide(greenV.add(nirV));
    
    var waterMaskV = mndwiV.gt(0.2)
                     .and(ndviV.lt(0.2))
                     .and(ndwiV.gt(0.0));
    
    return waterMaskV;
    
  } else if (dataType === 'MODIS') {
    scaled = image.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07']).multiply(0.0001);
    var greenM = scaled.select('sur_refl_b04');
    var swirM = scaled.select('sur_refl_b06');
    var swir2M = scaled.select('sur_refl_b07');
    var nirM = scaled.select('sur_refl_b02');
    var redM = scaled.select('sur_refl_b01');
    
    var mndwiM = greenM.subtract(swirM).divide(greenM.add(swirM));
    var ndviM = nirM.subtract(redM).divide(nirM.add(redM));
    var ndwiM = greenM.subtract(nirM).divide(greenM.add(nirM));
    
    var waterMaskM = mndwiM.gt(0.2)
                     .and(ndviM.lt(0.2))
                     .and(ndwiM.gt(0.0))
                     .and(swir2M.lt(0.15));
    
    return waterMaskM;
  }
  
  return ee.Image.constant(1);
}

function maskLandsatClouds(image) {
  var qa = image.select('QA_PIXEL');
  
  var cloudBit = 1 << 3;
  var cloudShadowBit = 1 << 4;
  var snowBit = 1 << 5;
  var cirrusBit = 1 << 2;
  
  var cloudMask = qa.bitwiseAnd(cloudBit).eq(0);
  var cloudShadowMask = qa.bitwiseAnd(cloudShadowBit).eq(0);
  var snowMask = qa.bitwiseAnd(snowBit).eq(0);
  var cirrusMask = qa.bitwiseAnd(cirrusBit).eq(0);
  
  var scaled = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  
  var green = scaled.select('SR_B3');
  var swir = scaled.select('SR_B6');
  var mndwi = green.subtract(swir).divide(green.add(swir));
  var waterMask = mndwi.gt(0);
  
  var nir = scaled.select('SR_B5');
  var brightMask = nir.lt(0.15);
  
  var finalMask = cloudMask.and(cloudShadowMask).and(snowMask)
                           .and(cirrusMask).and(waterMask).and(brightMask);
  
  return image.updateMask(finalMask);
}

function maskVIIRSclouds(image) {
  var spectralBands = image.select(['M1','M2','M3','M4','M5','M7','M8','M10','M11']);
  var bandNames = image.bandNames();
  var hasQF1 = bandNames.contains('QF1');
  
  var masked = ee.Algorithms.If(
    hasQF1,
    function() {
      var qa = image.select('QF1');
      var cloudMask = qa.bitwiseAnd(1 << 2).eq(0);
      var shadowMask = qa.bitwiseAnd(1 << 3).eq(0);
      var combinedMask = cloudMask.or(shadowMask);
      return spectralBands.updateMask(combinedMask);
    }(),
    spectralBands
  );
  
  return ee.Image(masked);
}

function fillGaps(image, dataType) {
  if (dataType === 'VIIRS' || dataType === 'MODIS') {
    var filled1 = image.focal_mean({radius: 3, kernelType: 'square', units: 'pixels', iterations: 1});
    var result1 = image.unmask(filled1);
    var filled2 = result1.focal_mean({radius: 5, kernelType: 'square', units: 'pixels', iterations: 1});
    return result1.unmask(filled2);
  }
  
  var filled1 = image.focal_mean({radius: 1.5, kernelType: 'square', units: 'pixels', iterations: 1});
  var result1 = image.unmask(filled1);
  
  var filled2 = result1.focal_mean({radius: 2.5, kernelType: 'square', units: 'pixels', iterations: 1});
  var result2 = result1.unmask(filled2);
  
  var filled3 = result2.focal_mean({radius: 3.5, kernelType: 'square', units: 'pixels', iterations: 1});
  
  return result2.unmask(filled3);
}


function applyShoreBuffer(image, dataType, bufferMeters) {
  bufferMeters = bufferMeters || 90;
  
  var scaled;
  if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') {
    scaled = image.select('SR_B.').multiply(0.0000275).add(-0.2);
    var green = scaled.select('SR_B3');
    var swir = scaled.select('SR_B6');
    var mndwi = green.subtract(swir).divide(green.add(swir));
    var waterMask = mndwi.gt(0.1);
    var waterBuffered = waterMask.focal_min({radius: bufferMeters, units: 'meters'});
    return image.updateMask(waterBuffered);
  } else if (dataType === 'S2') {
    scaled = image.select(['B.*']).multiply(0.0001);
    var greenS2 = scaled.select('B3');
    var swirS2 = scaled.select('B11');
    var mndwiS2 = greenS2.subtract(swirS2).divide(greenS2.add(swirS2));
    var waterMaskS2 = mndwiS2.gt(0.1);
    var waterBufferedS2 = waterMaskS2.focal_min({radius: bufferMeters, units: 'meters'});
    return image.updateMask(waterBufferedS2);
  }
  
  return image;
}

function checkS2MixedType(collection) {
  var c_size = collection.size(), imgList = collection.toList(c_size);
  var type_1st_img = ee.String(ee.String(collection.first().get('PRODUCT_ID')).split('_').get(0));
  for (var i = 1; i < c_size.getInfo(); i++) {
    var current_type = ee.String(ee.String(ee.Image(imgList.get(i)).get('PRODUCT_ID')).split('_').get(0));
    if (ee.String(ee.Algorithms.If(type_1st_img.equals(current_type), 'true', 'false')).getInfo() === 'false') return 'false';
  }
  return type_1st_img.getInfo();
}

function computeIndex(img, dataType, S2img_type, indexType) {
  var image;
  if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') {
    image = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  } else if (dataType === 'S2') {
    image = img.select(['B.*']).multiply(0.0001);
  } else if (dataType === 'VIIRS') {
    image = img.select(['M1','M2','M3','M4','M5','M7','M8','M10','M11']).multiply(0.0001);
  } else if (dataType === 'MODIS') {
    image = img.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07']).multiply(0.0001);
  }
  
  var computedIndex;
  if (indexType === 'FAI') {
    computedIndex = computeFAI(image, dataType, S2img_type);
    computedIndex = computedIndex.where(computedIndex.lt(-0.5).or(computedIndex.gt(0.5)), -999);
  } else if (indexType === 'NDCI') {
    computedIndex = computeNDCI(image, dataType);
  } else if (indexType === 'FLH') {
    computedIndex = computeFLH(image, dataType, S2img_type);
    computedIndex = computedIndex.where(computedIndex.lt(-0.2).or(computedIndex.gt(0.2)), -999);
  } else if (indexType === '2BDA') {
    computedIndex = compute2BDA(image, dataType);
  } else if (indexType === '3BDA') {
    computedIndex = compute3BDA(image, dataType);
  } else if (indexType === 'CI') {
    computedIndex = computeCI(image, dataType);
  }
  return computedIndex.updateMask(computedIndex.neq(-999));
}

function computeFAI(image, dataType, S2img_type) {
  var exp = 'NIR - R - (SWIR - R) * (lambda_NIR - lambda_R)/(lambda_SWIR - lambda_R)';
  if (dataType === 'L8') {
    return image.expression(exp, {
      'NIR': image.select('SR_B5'), 
      'R': image.select('SR_B4'), 
      'SWIR': image.select('SR_B6'), 
      'lambda_SWIR': 1608.5, 
      'lambda_R': 654.59, 
      'lambda_NIR': 864.67
    }).rename('FAI');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression(exp, {
      'NIR': image.select('SR_B4'), 
      'R': image.select('SR_B3'), 
      'SWIR': image.select('SR_B5'), 
      'lambda_SWIR': 1650, 
      'lambda_R': 660, 
      'lambda_NIR': 835
    }).rename('FAI');
  } else if (dataType === 'S2') {
    var lS = 1613.7, lR = 664.6, lN = 832.8;
    if (S2img_type === 'S2B') { 
      lS = 1610.4; lR = 664.9; lN = 832.9; 
    } else if (S2img_type === 'false') { 
      lS = 1612.05; lR = 664.75; lN = 832.85; 
    }
    return image.expression(exp, {
      'NIR': image.select('B8'), 
      'R': image.select('B4'), 
      'SWIR': image.select('B11'), 
      'lambda_SWIR': lS, 
      'lambda_R': lR, 
      'lambda_NIR': lN
    }).rename('FAI');
  } else if (dataType === 'VIIRS') {
    return image.expression(exp, {
      'NIR': image.select('M7'), 
      'R': image.select('M5'), 
      'SWIR': image.select('M11'), 
      'lambda_SWIR': 2250, 
      'lambda_R': 672, 
      'lambda_NIR': 865
    }).rename('FAI');
  } else if (dataType === 'MODIS') {
    return image.expression(exp, {
      'NIR': image.select('sur_refl_b02'), 
      'R': image.select('sur_refl_b01'), 
      'SWIR': image.select('sur_refl_b06'), 
      'lambda_SWIR': 1640, 
      'lambda_R': 645, 
      'lambda_NIR': 858
    }).rename('FAI');
  }
}

function computeNDCI(image, dataType) {
  if (dataType === 'L8') {
    return image.expression('(NIR - RED) / (NIR + RED)', {
      'RED': image.select('SR_B4'),
      'NIR': image.select('SR_B5')
    }).rename('NDCI');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression('(NIR - RED) / (NIR + RED)', {
      'RED': image.select('SR_B3'),
      'NIR': image.select('SR_B4')
    }).rename('NDCI');
  } else if (dataType === 'S2') {
    return image.expression('(RE - R) / (RE + R)', {
      'RE': image.select('B5'), 
      'R': image.select('B4')
    }).rename('NDCI');
  } else if (dataType === 'VIIRS') {
    return image.expression('(NIR - RED) / (NIR + RED)', {
      'RED': image.select('M5'),
      'NIR': image.select('M7')
    }).rename('NDCI');
  } else if (dataType === 'MODIS') {
    return image.expression('(NIR - RED) / (NIR + RED)', {
      'RED': image.select('sur_refl_b01'),
      'NIR': image.select('sur_refl_b02')
    }).rename('NDCI');
  }
}

function computeFLH(image, dataType, S2img_type) {
  if (dataType === 'L8') {
    return image.expression('B5 - (B4 + B6) / 2', {
      'B4': image.select('SR_B4'), 
      'B5': image.select('SR_B5'), 
      'B6': image.select('SR_B6')
    }).rename('FLH');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression('B4 - (B3 + B5) / 2', {
      'B3': image.select('SR_B3'), 
      'B4': image.select('SR_B4'), 
      'B5': image.select('SR_B5')
    }).rename('FLH');
  } else if (dataType === 'S2') {
    return image.expression('RE2 - RE1 - (RE3 - RE1) * (l2 - l1) / (l3 - l1)', {
      'RE1': image.select('B5'), 
      'RE2': image.select('B6'), 
      'RE3': image.select('B7'), 
      'l1': 704.1, 
      'l2': 740.5, 
      'l3': 782.8
    }).rename('FLH');
  } else if (dataType === 'VIIRS') {
    return image.expression('M7 - (M5 + M8) / 2', {
      'M5': image.select('M5'), 
      'M7': image.select('M7'), 
      'M8': image.select('M8')
    }).rename('FLH');
  } else if (dataType === 'MODIS') {
    return image.expression('NIR - (RED + SWIR) / 2', {
      'RED': image.select('sur_refl_b01'), 
      'NIR': image.select('sur_refl_b02'), 
      'SWIR': image.select('sur_refl_b06')
    }).rename('FLH');
  }
}

function compute2BDA(image, dataType) {
  if (dataType === 'L8') {
    return image.expression('NIR / RED', {
      'NIR': image.select('SR_B5').add(0.001), 
      'RED': image.select('SR_B4').add(0.001)
    }).rename('2BDA');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression('NIR / RED', {
      'NIR': image.select('SR_B4').add(0.001), 
      'RED': image.select('SR_B3').add(0.001)
    }).rename('2BDA');
  } else if (dataType === 'S2') {
    return image.expression('RE / RED', {
      'RE': image.select('B5').add(0.001), 
      'RED': image.select('B4').add(0.001)
    }).rename('2BDA');
  } else if (dataType === 'VIIRS') {
    return image.expression('M7 / M5', {
      'M7': image.select('M7').add(0.001), 
      'M5': image.select('M5').add(0.001)
    }).rename('2BDA');
  } else if (dataType === 'MODIS') {
    return image.expression('NIR / RED', {
      'NIR': image.select('sur_refl_b02').add(0.001), 
      'RED': image.select('sur_refl_b01').add(0.001)
    }).rename('2BDA');
  }
}

function compute3BDA(image, dataType) {
  if (dataType === 'L8') {
    return image.expression('(1/GREEN - 1/RED) * NIR', {
      'GREEN': image.select('SR_B3').add(0.001), 
      'RED': image.select('SR_B4').add(0.001), 
      'NIR': image.select('SR_B5')
    }).rename('3BDA');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression('(1/GREEN - 1/RED) * NIR', {
      'GREEN': image.select('SR_B2').add(0.001), 
      'RED': image.select('SR_B3').add(0.001), 
      'NIR': image.select('SR_B4')
    }).rename('3BDA');
  } else if (dataType === 'S2') {
    return image.expression('(1/GREEN - 1/RED) * RE', {
      'GREEN': image.select('B3').add(0.001), 
      'RED': image.select('B4').add(0.001), 
      'RE': image.select('B5')
    }).rename('3BDA');
  } else if (dataType === 'VIIRS') {
    return image.expression('(1/M4 - 1/M5) * M7', {
      'M4': image.select('M4').add(0.001), 
      'M5': image.select('M5').add(0.001), 
      'M7': image.select('M7')
    }).rename('3BDA');
  } else if (dataType === 'MODIS') {
    return image.expression('(1/GREEN - 1/RED) * NIR', {
      'GREEN': image.select('sur_refl_b04').add(0.001), 
      'RED': image.select('sur_refl_b01').add(0.001), 
      'NIR': image.select('sur_refl_b02')
    }).rename('3BDA');
  }
}

function computeCI(image, dataType) {
  if (dataType === 'L8') {
    return image.expression('NIR - RED - (SWIR - RED) * (lN - lR) / (lS - lR)', {
      'NIR': image.select('SR_B5'), 
      'RED': image.select('SR_B4'), 
      'SWIR': image.select('SR_B6'), 
      'lN': 864.67, 
      'lR': 654.59, 
      'lS': 1608.5
    }).rename('CI');
  } else if (dataType === 'L7' || dataType === 'L5') {
    return image.expression('NIR - RED - (SWIR - RED) * (lN - lR) / (lS - lR)', {
      'NIR': image.select('SR_B4'), 
      'RED': image.select('SR_B3'), 
      'SWIR': image.select('SR_B5'), 
      'lN': 835, 
      'lR': 660, 
      'lS': 1650
    }).rename('CI');
  } else if (dataType === 'S2') {
    return image.expression('Rrs704 - Rrs665 - (Rrs740 - Rrs665) * (704 - 665) / (740 - 665)', {
      'Rrs704': image.select('B5'), 
      'Rrs665': image.select('B4'), 
      'Rrs740': image.select('B6')
    }).rename('CI');
  } else if (dataType === 'VIIRS') {
    return image.expression('M7 - M5 - (M10 - M5) * (l7 - l5) / (l10 - l5)', {
      'M7': image.select('M7'), 
      'M5': image.select('M5'), 
      'M10': image.select('M10'), 
      'l7': 865, 
      'l5': 672, 
      'l10': 1610
    }).rename('CI');
  } else if (dataType === 'MODIS') {
    return image.expression('NIR - RED - (SWIR - RED) * (lN - lR) / (lS - lR)', {
      'NIR': image.select('sur_refl_b02'), 
      'RED': image.select('sur_refl_b01'), 
      'SWIR': image.select('sur_refl_b06'), 
      'lN': 858, 
      'lR': 645, 
      'lS': 1640
    }).rename('CI');
  }
}

function computeChlorophyllA(img, dataType) {
  var image;
  if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') image = img.select('SR_B.').multiply(0.0000275).add(-0.2);
  else if (dataType === 'S2') image = img.select(['B.*']).multiply(0.0001);
  else if (dataType === 'VIIRS') image = img.select(['M1','M2','M3','M4','M5','M7','M8','M10','M11']).multiply(0.0001);
  else if (dataType === 'MODIS') image = img.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07']).multiply(0.0001);
  
  var ratio, chla;
  if (dataType === 'L8') {
    ratio = image.expression('log10(max(B, 0.0001) / max(G, 0.0001))', {B: image.select('SR_B2'), G: image.select('SR_B3')});
    chla = image.expression('pow(10, (a0 + a1*R + a2*pow(R,2) + a3*pow(R,3) + a4*pow(R,4)))', {R: ratio, a0: 0.2424, a1: -2.7423, a2: 1.8017, a3: 0.0015, a4: -1.2280}).rename('Chla');
    return chla.where(chla.lt(0), 0).where(chla.gt(100), -999).updateMask(chla.neq(-999));
  } else if (dataType === 'L7' || dataType === 'L5') {
    ratio = image.expression('log10(max(B, 0.0001) / max(G, 0.0001))', {B: image.select('SR_B1'), G: image.select('SR_B2')});
    chla = image.expression('pow(10, (a0 + a1*R + a2*pow(R,2) + a3*pow(R,3)))', {R: ratio, a0: 0.2830, a1: -2.7753, a2: 1.8272, a3: -0.5340}).rename('Chla');
    return chla.where(chla.lt(0), 0).where(chla.gt(100), -999).updateMask(chla.neq(-999));
  } else if (dataType === 'S2') {
    var reci = image.select('B5').divide(image.select('B4')).subtract(1);
    chla = reci.multiply(50).add(50).rename('Chla');
    return chla.where(chla.lt(0), 0).where(chla.gt(100), -999).updateMask(chla.neq(-999));
  } else if (dataType === 'VIIRS' || dataType === 'MODIS') {
    var b = (dataType === 'VIIRS') ? 'M3' : 'sur_refl_b03', g = (dataType === 'VIIRS') ? 'M4' : 'sur_refl_b04';
    ratio = image.expression('log10(max(B, 0.0001) / max(G, 0.0001))', {B: image.select(b), G: image.select(g)});
    chla = image.expression('pow(10, (a0 + a1*R + a2*pow(R,2) + a3*pow(R,3)))', {R: ratio, a0: 0.2830, a1: -2.7753, a2: 1.8272, a3: -0.5340}).rename('Chla');
    return chla.where(chla.lt(0), 0).where(chla.gt(100), -999).updateMask(chla.neq(-999));
  }
}


function toggleItem(side, shown) {
  if (side === 'left') {
    left.dataDate.chartLabel.style().set({shown: shown});
    left.dataDate.panel2.style().set({shown: shown});
    left.computeButton.style().set({shown: shown});
  } else {
    right.dataDate.chartLabel.style().set({shown: shown});
    right.dataDate.panel2.style().set({shown: shown});
    right.computeButton.style().set({shown: shown});
  }
}

function draw(featureColl, remove, side) {
  var chart = ui.Chart.feature.byFeature({features: featureColl.select('Number', 'label'), xProperty: 'label'}).setChartType('ColumnChart').setOptions({
    width: 100, height: 40, title: 'Image Availability',
    hAxis: {title: 'Dates', titleTextStyle: {italic: false, bold: true}},
    vAxis: {title: 'Number of images', titleTextStyle: {italic: false, bold: true}},
    colors: ['#e0440e'], isStacked: 'absolute'
  });
  if (remove === true) {
    if (side === 'left') {
      if (left.dataDate.chart) left.controlPanel.widgets().remove(left.dataDate.chart);
      left.controlPanel.widgets().insert(5, chart);
    } else {
      if (right.dataDate.chart) right.controlPanel.widgets().remove(right.dataDate.chart);
      right.controlPanel.widgets().insert(2, chart);
    }
    toggleItem(side, true);
  }
  return chart;
}

function getImg(collection, date, remove, side) {
  if (!collection || collection.size().getInfo() === 0) {
    toggleItem(side, false);
    return null;
  }
  var range = collection.reduceColumns(ee.Reducer.toList(), ["system:time_start"]).values().get(0);
  range = ee.List(range).map(function(n){ return ee.Date(n).format("YYYY-MM-dd"); });
  var dates = range.distinct().sort();
  var merged_collection = ee.FeatureCollection(dates.map(function(date){
    var nfc = collection.filter(ee.Filter.date(date, ee.Date(date).advance(1, 'day')));
    return nfc.union().first().set({label: date, Number: nfc.size()});
  }));
  return draw(merged_collection, remove, side);
}

function updateCollection(side) {
  if (side === 'left') {
    var start = ee.Date(x_left[0]), end = ee.Date(x_left[1]);
    if (end.difference(start, 'day').getInfo() > 0) {
      collection = getData('left');
      if (collection && collection.size().getInfo() > 0) {
        left.dataDate.chart = getImg(collection, x_left, true, 'left');
      } else {
        toggleItem('left', true);
      }
    }
  } else {
    var start = ee.Date(x_right[0]), end = ee.Date(x_right[1]);
    if (end.difference(start, 'day').getInfo() > 0) {
      collection_right = getData('right');
      if (collection_right && collection_right.size().getInfo() > 0) {
        right.dataDate.chart = getImg(collection_right, x_right, true, 'right');
      } else {
        toggleItem('right', true);
      }
    }
  }
}

function keepOneGeometry() {
  if (drawingTools.layers().get(0).geometries().length() > 1) {
    drawingTools.layers().get(0).geometries().remove(drawingTools.layers().get(0).geometries().get(0));
  }
  aoi = drawingTools.layers().get(0).getEeObject();
  updateCollection('left');
  updateCollection('right');
}

function drawRectangle() {
  clearGeometry();
  drawingTools.setShape('rectangle');
  drawingTools.draw();
  drawingTools.onDraw(keepOneGeometry);
}

function drawPolygon() {
  clearGeometry();
  drawingTools.setShape('polygon');
  drawingTools.draw();
  drawingTools.onDraw(keepOneGeometry);
}

function addPoint() {
  clearGeometry();
  drawingTools.setShape('point');
  drawingTools.draw();
  drawingTools.onDraw(keepOneGeometry);
}

function updateDownloadableItem() {
  var itemList = [];
  if (binaryMap) {
    itemList.push({label: 'Algal_Map (L)', value: '1'});
    itemList.push({label: 'Index_Map (L)', value: '2'});
    if (chlaMap) itemList.push({label: 'Chla_Map (L)', value: '3'});
  }
  if (binaryMap_right) {
    itemList.push({label: 'Algal_Map (R)', value: '4'});
    itemList.push({label: 'Index_Map (R)', value: '5'});
    if (chlaMap_right) itemList.push({label: 'Chla_Map (R)', value: '6'});
  }
  right.download.item.items().reset(itemList);
}

var left = {}, right = {};
left.controlPanel = ui.Panel();
left.divider = ui.Panel();
left.info = {};
left.info.titleLabel = ui.Label('TerraBloom');
left.info.aboutLabel = ui.Label("TerraBloom TerraBloom is a Google Earth Engine–powered web app designed to detect and analyze harmful algal blooms (HABs) worldwide using satellite imagery. Users can select an area of interest, choose a time range and satellite dataset, and calculate bloom indicators using indices such as NDVI, NDWI, and the Floating Algae Index. By combining real-time environmental data with an intuitive interface, TerraBloom empowers researchers, students, and policymakers to monitor water quality and better understand the global impact of harmful algal blooms.");
left.info.paper = ui.Label({value: 'References:', style: {width: '100px', margin: '0px 0px 0px 8px', height: '20px', fontWeight: 'bold'}});
left.info.paperLabel1 = ui.Label({value: 'FAI', targetUrl: 'https://doi.org/10.1016/j.rse.2009.05.012'});
left.info.paperLabel2 = ui.Label({value: 'NDCI', targetUrl: 'https://doi.org/10.1016/j.rse.2010.04.002'});
left.info.lab = ui.Label({value: 'Our Research Group', targetUrl: 'https://hydroinformatics.uiowa.edu/'});
left.info.horizontalPanel = ui.Panel({widgets: [left.info.paper, left.info.paperLabel1, left.info.paperLabel2], layout: ui.Panel.Layout.flow('horizontal')});
left.info.panel = ui.Panel([left.info.titleLabel, left.info.aboutLabel, left.info.horizontalPanel, left.info.lab]);

left.selectAOI = {};
left.selectAOI.label = ui.Label('Define Area of Interest');
left.selectAOI.aoi = ui.Select({
  items: [{label: 'Rectangle', value: 'Rectangle'}, {label: 'Polygon', value: 'Polygon'}, {label: 'Point', value: 'Point'}],
  placeholder: 'Choose & Draw AOI',
  onChange: function (value) {
    if (value === 'Rectangle') drawRectangle();
    else if (value === 'Polygon') drawPolygon();
    else addPoint();
  }
});
left.selectAOI.panel = ui.Panel([left.selectAOI.label].concat(left.selectAOI.aoi));

left.dataDate = {};
left.dataDate.label1 = ui.Label('Select Start Date:');
left.dataDate.slider1 = ui.DateSlider({
  start: ee.Date('1984-01-01'), value: '2017-08-01',
  onChange: function(dr) {
    x_left[0] = dr.start().format("YYYY-MM-dd");
    if (ee.Date(x_left[1]).difference(ee.Date(x_left[0]), 'day').getInfo() > 0) {
      collection = getData('left');
      if (collection && collection.size().getInfo() > 0) {
        left.dataDate.chart = getImg(collection, x_left, true, 'left');
      } else {
        toggleItem('left', true);
      }
    } else {
      toggleItem('left', true);
    }
  }
});

left.dataDate.label2 = ui.Label('Select End Date:');
left.dataDate.slider2 = ui.DateSlider({
  start: ee.Date('1984-01-01'), value: '2017-08-31', period: 1,
  onChange: function(dr) {
    x_left[1] = dr.start().format("YYYY-MM-dd");
    if (ee.Date(x_left[1]).difference(ee.Date(x_left[0]), 'day').getInfo() > 0) {
      collection = getData('left');
      if (collection && collection.size().getInfo() > 0) {
        left.dataDate.chart = getImg(collection, x_left, true, 'left');
      } else {
        toggleItem('left', true);
      }
    } else {
      toggleItem('left', true);
    }
  }
});

left.dataDate.dataLabel = ui.Label('Select Satellite:');
left.dataDate.data = ui.Select({
  items: [
    {label: 'Sentinel-2 (2015-present)', value: 'S2'},
    {label: 'Landsat 8 (2013-present)', value: 'L8'},
    {label: 'Landsat 7 (1999-present)', value: 'L7'},
    {label: 'Landsat 5 (1984-2012)', value: 'L5'},
    {label: 'VIIRS (2012-present)', value: 'VIIRS'},
    {label: 'MODIS (2000-present)', value: 'MODIS'}
  ],
  placeholder: 'Select the imagery',
  onChange: function() {
    collection = getData('left');
    if (collection && collection.size().getInfo() > 0) {
      left.dataDate.chart = getImg(collection, x_left, true, 'left');
    } else {
      toggleItem('left', true);
    }
  }
});

left.dataDate.indexLabel = ui.Label('Select Index:');
left.dataDate.indexSelect = ui.Select({
  items: [
    {label: 'FAI - Floating Algae Index', value: 'FAI'},
    {label: 'NDCI - Normalized Difference Chlorophyll Index', value: 'NDCI'},
    {label: 'FLH - Fluorescence Line Height', value: 'FLH'},
    {label: '2BDA - Two-Band Difference Algorithm', value: '2BDA'},
    {label: '3BDA - Three-Band Difference Algorithm', value: '3BDA'},
    {label: 'MCI - Maximum Chlorophyll Index', value: 'CI'}
  ],
  value: 'FAI',
  onChange: function() {
    var idx = left.dataDate.indexSelect.getValue();
    if (idx === 'FAI') { left.dataDate.trTextbox.setValue('0.005'); }
    else if (idx === 'NDCI') { left.dataDate.trTextbox.setValue('0.1'); }
    else if (idx === 'FLH') { left.dataDate.trTextbox.setValue('0.002'); }
    else if (idx === '2BDA') { left.dataDate.trTextbox.setValue('1.0'); }
    else if (idx === '3BDA') { left.dataDate.trTextbox.setValue('0.005'); }
    else if (idx === 'CI') { left.dataDate.trTextbox.setValue('0.002'); }
  }
});

left.dataDate.chartLabel = ui.Label('Image Availability:');
left.dataDate.mosaicLabel = ui.Label('Select Mosaic Approach:');
left.dataDate.mosaicc = ui.Select({
  items: [
    {label: 'Take the First', value: '1'},
    {label: 'Maximum (could be slow)', value: '2'},
    {label: 'Average (recommended)', value: '3'},
    {label: 'Minimum (could be slow)', value: '4'},
    {label: 'Spatial Mosaic', value: '5'},
    {label: 'Clear Composite (S2/MODIS only)', value: '6'}
  ],
  placeholder: 'Mosaic approaches',
  value: '3'
});

left.dataDate.thresholdLabel = ui.Label('Threshold for Classification:');
left.dataDate.trTextbox = ui.Textbox({
  value: '0.005',
  placeholder: 'Enter threshold value',
  style: {fontWeight: 'bold', textAlign: 'center', backgroundColor: 'EBEAEA'}
});

left.dataDate.bufferLabel = ui.Label('Shore Buffer (meters):');
left.dataDate.bufferTextbox = ui.Textbox({
  value: '90',
  placeholder: 'Buffer distance',
  style: {textAlign: 'center', backgroundColor: 'EBEAEA'}
});

left.dataDate.chlaCheckbox = ui.Checkbox({label: 'Calculate Chlorophyll-a', value: false});

left.dataDate.panel = ui.Panel([
  left.dataDate.label1, left.dataDate.slider1, left.dataDate.label2, left.dataDate.slider2,
  left.dataDate.dataLabel, left.dataDate.data, left.dataDate.indexLabel, left.dataDate.indexSelect,
  left.dataDate.mosaicLabel, left.dataDate.mosaicc, left.dataDate.chartLabel
]);

left.dataDate.panel2 = ui.Panel([
  left.dataDate.thresholdLabel, left.dataDate.trTextbox,
  left.dataDate.bufferLabel, left.dataDate.bufferTextbox,
  left.dataDate.chlaCheckbox
]);

left.computeButton = ui.Button({
  label: 'Get Algal Bloom Map',
  onClick: function() {
    var mosaic = left.dataDate.mosaicc.getValue() || '3';
    var dataType = (!left.dataDate.data.getValue()) ? 'S2' : left.dataDate.data.getValue();
    var indexType = left.dataDate.indexSelect.getValue() || 'FAI';
    var calcChla = left.dataDate.chlaCheckbox.getValue();
    var img, S2img_type;
    
    var Tr = parseFloat(left.dataDate.trTextbox.getValue());
    if (isNaN(Tr)) {
      print('Error: Invalid threshold value. Using default 0.005');
      Tr = 0.005;
      left.dataDate.trTextbox.setValue('0.005');
    }
    
    var bufferDist = parseFloat(left.dataDate.bufferTextbox.getValue());
    if (isNaN(bufferDist) || bufferDist < 0) {
      print('Error: Invalid buffer distance. Using default 90m');
      bufferDist = 90;
      left.dataDate.bufferTextbox.setValue('90');
    }
    
    if (mosaic === '6') {
      if (dataType === 'S2') {
        var aoiToUse = aoi || ee.Geometry.Point([-83.5379, 41.6528]).buffer(10000);
        img = getClearSentinel2Composite(aoiToUse, x_left[0], x_left[1]);
        S2img_type = 'S2A';
      } else if (dataType === 'MODIS') {
        var aoiToUse = aoi || ee.Geometry.Point([-83.5379, 41.6528]).buffer(10000);
        img = getClearMODISComposite(aoiToUse, x_left[0], x_left[1]);
      } else {
        print('Clear Composite only available for Sentinel-2 and MODIS. Using average mosaic.');
        mosaic = '3';
      }
    }
    
    if (mosaic !== '6') {
      function clp(img) {
        var newImg = img.clip(aoi);
        if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') return maskLandsatClouds(newImg);
        else if (dataType === 'S2') {
          var scl = newImg.select('SCL');
          var qa60 = newImg.select('QA60');
          var cloudBitMask = ee.Number(2).pow(10).int();
          var cirrusBitMask = ee.Number(2).pow(11).int();
          var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
          var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
          var water = scl.eq(6);
          return newImg.updateMask(qaMask.and(sclMask).and(water));
        }
        else if (dataType === 'VIIRS') return maskVIIRSclouds(newImg);
        else if (dataType === 'MODIS') {
          var qa = newImg.select('state_1km');
          return newImg.updateMask(qa.bitwiseAnd(3).eq(0));
        }
        return newImg;
      }

      if (mosaic === '1' || collection.size().getInfo() === 1) {
        img = collection.first();
        if (aoi && aoi.type().getInfo() === 'Polygon') img = img.clip(aoi);
        if (dataType === 'S2') {
          var scl = img.select('SCL');
          var qa60 = img.select('QA60');
          var cloudBitMask = ee.Number(2).pow(10).int();
          var cirrusBitMask = ee.Number(2).pow(11).int();
          var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
          var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
          var water = scl.eq(6);
          img = img.updateMask(qaMask.and(sclMask).and(water));
          S2img_type = (ee.String(img.get('PRODUCT_ID')).split('_').get(0)).getInfo();
        } else if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') {
          img = maskLandsatClouds(img);
        } else if (dataType === 'VIIRS') {
          img = maskVIIRSclouds(img);
        } else if (dataType === 'MODIS') {
          var qa = img.select('state_1km');
          img = img.updateMask(qa.bitwiseAnd(3).eq(0));
        }
      } else {
        var newColl;
        if (dataType === 'S2') S2img_type = checkS2MixedType(collection);
        if (aoi && aoi.type().getInfo() === 'Polygon') newColl = collection.map(clp);
        else {
          if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') newColl = collection.map(maskLandsatClouds);
          else if (dataType === 'S2') newColl = collection.map(function(image) {
            var scl = image.select('SCL');
            var qa60 = image.select('QA60');
            var cloudBitMask = ee.Number(2).pow(10).int();
            var cirrusBitMask = ee.Number(2).pow(11).int();
            var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
            var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
            var water = scl.eq(6);
            return image.updateMask(qaMask.and(sclMask).and(water));
          });
          else if (dataType === 'VIIRS') newColl = collection.map(maskVIIRSclouds);
          else if (dataType === 'MODIS') newColl = collection.map(function(image) {
            var qa = image.select('state_1km');
            return image.updateMask(qa.bitwiseAnd(3).eq(0));
          });
        }
        
        if (mosaic === '2') img = newColl.max();
        else if (mosaic === '3') img = newColl.mean();
        else if (mosaic === '4') img = newColl.min();
        else img = newColl.mosaic();
        
        var scope_left = newColl.geometry().dissolve({'maxError': 1});
        img = (aoi && aoi.type().getInfo() === 'Polygon') ? img.clip(aoi) : img.clip(scope_left);
      }
      
      img = fillGaps(img, dataType);
    }
    
    var waterMask = createWaterMask(img, dataType);
    waterMask = waterMask.and(permanentWater);
    img = img.updateMask(waterMask);
    
    if (bufferDist > 0) {
      img = applyShoreBuffer(img, dataType, bufferDist);
    }

    leftMap.layers().reset();
    reset('left');
    drawingTools.layers().get(0).setShown(false);

    indexMap = computeIndex(img, dataType, S2img_type, indexType);
    binaryMap = indexMap.gt(Tr).selfMask();
    
    if (calcChla) chlaMap = computeChlorophyllA(img, dataType);

    if (aoi) leftMap.layers().add(ui.Map.Layer(aoi, {color: 'A760D8'}, 'AOI', false));
    
    leftMap.layers().add(ui.Map.Layer(img, eval("visualization".concat("_", dataType)), 'Original Image'));

    var indexPalette, indexMin, indexMax;
    if (indexType === 'FAI') { indexPalette = index_palettes.colorbrewer.Greens[7]; indexMin = -0.05; indexMax = 0.1; }
    else if (indexType === 'NDCI') { indexPalette = index_palettes.colorbrewer.Blues[7]; indexMin = -0.5; indexMax = 0.5; }
    else if (indexType === 'FLH') { indexPalette = index_palettes.colorbrewer.Reds[7]; indexMin = -0.02; indexMax = 0.02; }
    else if (indexType === '2BDA') { indexPalette = index_palettes.colorbrewer.Oranges[7]; indexMin = 0.5; indexMax = 2.0; }
    else if (indexType === '3BDA') { indexPalette = index_palettes.colorbrewer.Purples[7]; indexMin = -0.02; indexMax = 0.1; }
    else if (indexType === 'CI') { indexPalette = index_palettes.colorbrewer.YlGn[7]; indexMin = -0.02; indexMax = 0.02; }

    leftMap.layers().add(ui.Map.Layer(indexMap, {palette: indexPalette, min: indexMin, max: indexMax}, indexType));
    leftMap.layers().add(ui.Map.Layer(binaryMap, {palette: "00FF00, FF0000"}, 'Algal Bloom'));
    
    if (calcChla) leftMap.layers().add(ui.Map.Layer(chlaMap, {palette: ['blue', 'cyan', 'yellow', 'red'], min: 0, max: 100}, 'Chlorophyll-a (μg/L)'));

    print('Images used: ' + collection.size().getInfo());
    updateDownloadableItem();
  }
});

right.controlPanel = ui.Panel();
right.divider = ui.Panel();
right.titleLabel = ui.Label('Comparison');

right.dataDate = {};
right.dataDate.label1 = ui.Label('Select Start Date:');
right.dataDate.slider1 = ui.DateSlider({
  start: ee.Date('1984-01-01'), value: '2017-08-01',
  onChange: function(dr) {
    x_right[0] = dr.start().format("YYYY-MM-dd");
    if (ee.Date(x_right[1]).difference(ee.Date(x_right[0]), 'day').getInfo() > 0) {
      collection_right = getData('right');
      if (collection_right && collection_right.size().getInfo() > 0) {
        right.dataDate.chart = getImg(collection_right, x_right, true, 'right');
      } else {
        toggleItem('right', true);
      }
    } else {
      toggleItem('right', true);
    }
  }
});

right.dataDate.label2 = ui.Label('Select End Date:');
right.dataDate.slider2 = ui.DateSlider({
  start: ee.Date('1984-01-01'), value: '2017-08-31', period: 1,
  onChange: function(dr) {
    x_right[1] = dr.start().format("YYYY-MM-dd");
    if (ee.Date(x_right[1]).difference(ee.Date(x_right[0]), 'day').getInfo() > 0) {
      collection_right = getData('right');
      if (collection_right && collection_right.size().getInfo() > 0) {
        right.dataDate.chart = getImg(collection_right, x_right, true, 'right');
      } else {
        toggleItem('right', true);
      }
    } else {
      toggleItem('right', true);
    }
  }
});

right.dataDate.dataLabel = ui.Label('Select Satellite:');
right.dataDate.data = ui.Select({
  items: [
    {label: 'Sentinel-2 (2015-present)', value: 'S2'},
    {label: 'Landsat 8 (2013-present)', value: 'L8'},
    {label: 'Landsat 7 (1999-present)', value: 'L7'},
    {label: 'Landsat 5 (1984-2012)', value: 'L5'},
    {label: 'VIIRS (2012-present)', value: 'VIIRS'},
    {label: 'MODIS (2000-present)', value: 'MODIS'}
  ],
  placeholder: 'Select the imagery',
  onChange: function() {
    collection_right = getData('right');
    if (collection_right && collection_right.size().getInfo() > 0) {
      right.dataDate.chart = getImg(collection_right, x_right, true, 'right');
    } else {
      toggleItem('right', true);
    }
  }
});

right.dataDate.indexLabel = ui.Label('Select Index:');
right.dataDate.indexSelect = ui.Select({
  items: [
    {label: 'FAI - Floating Algae Index', value: 'FAI'},
    {label: 'NDCI - Normalized Difference Chlorophyll Index', value: 'NDCI'},
    {label: 'FLH - Fluorescence Line Height', value: 'FLH'},
    {label: '2BDA - Two-Band Difference Algorithm', value: '2BDA'},
    {label: '3BDA - Three-Band Difference Algorithm', value: '3BDA'},
    {label: 'MCI - Maximum Chlorophyll Index', value: 'CI'}
  ],
  value: 'FAI',
  onChange: function() {
    var idx = right.dataDate.indexSelect.getValue();
    if (idx === 'FAI') { right.dataDate.trTextbox.setValue('0.005'); }
    else if (idx === 'NDCI') { right.dataDate.trTextbox.setValue('0.1'); }
    else if (idx === 'FLH') { right.dataDate.trTextbox.setValue('0.002'); }
    else if (idx === '2BDA') { right.dataDate.trTextbox.setValue('1.0'); }
    else if (idx === '3BDA') { right.dataDate.trTextbox.setValue('0.005'); }
    else if (idx === 'CI') { right.dataDate.trTextbox.setValue('0.002'); }
  }
});

right.dataDate.chartLabel = ui.Label('Image Availability:');
right.dataDate.mosaicLabel = ui.Label('Select Mosaic Approach:');
right.dataDate.mosaicc = ui.Select({
  items: [
    {label: 'Take the First', value: '1'},
    {label: 'Maximum (could be slow)', value: '2'},
    {label: 'Average (recommended)', value: '3'},
    {label: 'Minimum (could be slow)', value: '4'},
    {label: 'Spatial Mosaic', value: '5'},
    {label: 'Clear Composite (S2/MODIS only)', value: '6'}
  ],
  placeholder: 'Mosaic approaches',
  value: '3'
});

right.dataDate.thresholdLabel = ui.Label('Threshold for Classification:');
right.dataDate.trTextbox = ui.Textbox({
  value: '0.005',
  placeholder: 'Enter threshold value',
  style: {fontWeight: 'bold', textAlign: 'center', backgroundColor: 'EBEAEA'}
});

right.dataDate.bufferLabel = ui.Label('Shore Buffer (meters):');
right.dataDate.bufferTextbox = ui.Textbox({
  value: '90',
  placeholder: 'Buffer distance',
  style: {textAlign: 'center', backgroundColor: 'EBEAEA'}
});

right.dataDate.chlaCheckbox = ui.Checkbox({label: 'Calculate Chlorophyll-a', value: false});

right.dataDate.panel = ui.Panel([
  right.titleLabel, right.dataDate.label1, right.dataDate.slider1, right.dataDate.label2, right.dataDate.slider2,
  right.dataDate.dataLabel, right.dataDate.data, right.dataDate.indexLabel, right.dataDate.indexSelect,
  right.dataDate.mosaicLabel, right.dataDate.mosaicc, right.dataDate.chartLabel
]);

right.dataDate.panel2 = ui.Panel([
  right.dataDate.thresholdLabel, right.dataDate.trTextbox,
  right.dataDate.bufferLabel, right.dataDate.bufferTextbox,
  right.dataDate.chlaCheckbox
]);

right.computeButton = ui.Button({
  label: 'Get Algal Bloom Map',
  onClick: function() {
    var mosaic = right.dataDate.mosaicc.getValue() || '3';
    var dataType = (!right.dataDate.data.getValue()) ? 'S2' : right.dataDate.data.getValue();
    var indexType = right.dataDate.indexSelect.getValue() || 'FAI';
    var calcChla = right.dataDate.chlaCheckbox.getValue();
    var img, S2img_type;
    
    var Tr = parseFloat(right.dataDate.trTextbox.getValue());
    if (isNaN(Tr)) {
      print('Error: Invalid threshold value. Using default 0.005');
      Tr = 0.005;
      right.dataDate.trTextbox.setValue('0.005');
    }
    
    var bufferDist = parseFloat(right.dataDate.bufferTextbox.getValue());
    if (isNaN(bufferDist) || bufferDist < 0) {
      print('Error: Invalid buffer distance. Using default 90m');
      bufferDist = 90;
      right.dataDate.bufferTextbox.setValue('90');
    }
    
    if (mosaic === '6') {
      if (dataType === 'S2') {
        var aoiToUse = aoi || ee.Geometry.Point([-83.5379, 41.6528]).buffer(10000);
        img = getClearSentinel2Composite(aoiToUse, x_right[0], x_right[1]);
        S2img_type = 'S2A';
      } else if (dataType === 'MODIS') {
        var aoiToUse = aoi || ee.Geometry.Point([-83.5379, 41.6528]).buffer(10000);
        img = getClearMODISComposite(aoiToUse, x_right[0], x_right[1]);
      } else {
        print('Clear Composite only available for Sentinel-2 and MODIS. Using average mosaic.');
        mosaic = '3';
      }
    }
    
    if (mosaic !== '6') {
      function clp(img) {
        var newImg = img.clip(aoi);
        if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') return maskLandsatClouds(newImg);
        else if (dataType === 'S2') {
          var scl = newImg.select('SCL');
          var qa60 = newImg.select('QA60');
          var cloudBitMask = ee.Number(2).pow(10).int();
          var cirrusBitMask = ee.Number(2).pow(11).int();
          var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
          var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
          var water = scl.eq(6);
          return newImg.updateMask(qaMask.and(sclMask).and(water));
        }
        else if (dataType === 'VIIRS') return maskVIIRSclouds(newImg);
        else if (dataType === 'MODIS') {
          var qa = newImg.select('state_1km');
          return newImg.updateMask(qa.bitwiseAnd(3).eq(0));
        }
        return newImg;
      }

      if (mosaic === '1' || collection_right.size().getInfo() === 1) {
        img = collection_right.first();
        if (aoi && aoi.type().getInfo() === 'Polygon') img = img.clip(aoi);
        if (dataType === 'S2') {
          var scl = img.select('SCL');
          var qa60 = img.select('QA60');
          var cloudBitMask = ee.Number(2).pow(10).int();
          var cirrusBitMask = ee.Number(2).pow(11).int();
          var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
          var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
          var water = scl.eq(6);
          img = img.updateMask(qaMask.and(sclMask).and(water));
          S2img_type = (ee.String(img.get('PRODUCT_ID')).split('_').get(0)).getInfo();
        } else if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') {
          img = maskLandsatClouds(img);
        } else if (dataType === 'VIIRS') {
          img = maskVIIRSclouds(img);
        } else if (dataType === 'MODIS') {
          var qa = img.select('state_1km');
          img = img.updateMask(qa.bitwiseAnd(3).eq(0));
        }
      } else {
        var newColl;
        if (dataType === 'S2') S2img_type = checkS2MixedType(collection_right);
        if (aoi && aoi.type().getInfo() === 'Polygon') newColl = collection_right.map(clp);
        else {
          if (dataType === 'L8' || dataType === 'L7' || dataType === 'L5') newColl = collection_right.map(maskLandsatClouds);
          else if (dataType === 'S2') newColl = collection_right.map(function(image) {
            var scl = image.select('SCL');
            var qa60 = image.select('QA60');
            var cloudBitMask = ee.Number(2).pow(10).int();
            var cirrusBitMask = ee.Number(2).pow(11).int();
            var qaMask = qa60.bitwiseAnd(cloudBitMask).eq(0).and(qa60.bitwiseAnd(cirrusBitMask).eq(0));
            var sclMask = scl.neq(3).and(scl.neq(7)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
            var water = scl.eq(6);
            return image.updateMask(qaMask.and(sclMask).and(water));
          });
          else if (dataType === 'VIIRS') newColl = collection_right.map(maskVIIRSclouds);
          else if (dataType === 'MODIS') newColl = collection_right.map(function(image) {
            var qa = image.select('state_1km');
            return image.updateMask(qa.bitwiseAnd(3).eq(0));
          });
        }
        
        if (mosaic === '2') img = newColl.max();
        else if (mosaic === '3') img = newColl.mean();
        else if (mosaic === '4') img = newColl.min();
        else img = newColl.mosaic();
        
        var scope_right = newColl.geometry().dissolve({'maxError': 1});
        img = (aoi && aoi.type().getInfo() === 'Polygon') ? img.clip(aoi) : img.clip(scope_right);
      }
      
      img = fillGaps(img, dataType);
    }
    
    var waterMask = createWaterMask(img, dataType);
    img = img.updateMask(waterMask);
    
    if (bufferDist > 0) {
      img = applyShoreBuffer(img, dataType, bufferDist);
    }

    rightMap.layers().reset();
    reset('right');

    indexMap_right = computeIndex(img, dataType, S2img_type, indexType);
    binaryMap_right = indexMap_right.gt(Tr).selfMask();
    
    if (calcChla) chlaMap_right = computeChlorophyllA(img, dataType);

    if (aoi) rightMap.layers().add(ui.Map.Layer(aoi, {color: 'A760D8'}, 'AOI', false));
    
    rightMap.layers().add(ui.Map.Layer(img, eval("visualization".concat("_", dataType)), 'Original Image'));

    var indexPalette, indexMin, indexMax;
    if (indexType === 'FAI') { indexPalette = index_palettes.colorbrewer.Greens[7]; indexMin = -0.05; indexMax = 0.1; }
    else if (indexType === 'NDCI') { indexPalette = index_palettes.colorbrewer.Blues[7]; indexMin = -0.5; indexMax = 0.5; }
    else if (indexType === 'FLH') { indexPalette = index_palettes.colorbrewer.Reds[7]; indexMin = -0.02; indexMax = 0.02; }
    else if (indexType === '2BDA') { indexPalette = index_palettes.colorbrewer.Oranges[7]; indexMin = 0.5; indexMax = 2.0; }
    else if (indexType === '3BDA') { indexPalette = index_palettes.colorbrewer.Purples[7]; indexMin = -0.02; indexMax = 0.1; }
    else if (indexType === 'CI') { indexPalette = index_palettes.colorbrewer.YlGn[7]; indexMin = -0.02; indexMax = 0.02; }

    rightMap.layers().add(ui.Map.Layer(indexMap_right, {palette: indexPalette, min: indexMin, max: indexMax}, indexType));
    rightMap.layers().add(ui.Map.Layer(binaryMap_right, {palette: "00FF00, FF0000"}, 'Algal Bloom'));
    
    if (calcChla) rightMap.layers().add(ui.Map.Layer(chlaMap_right, {palette: ['blue', 'cyan', 'yellow', 'red'], min: 0, max: 100}, 'Chlorophyll-a (μg/L)'));

    print('Images used: ' + collection_right.size().getInfo());
    updateDownloadableItem();
  }
});

right.download = {};
right.download.label = ui.Label('Download Image');
right.download.KMLlink = ui.Label('KML File', {shown: false});
right.download.item = ui.Select({
  placeholder: 'Select Image',
  onChange: function(value) {
    var outputImg, outScope, fileName;
    if (value === '1') { outputImg = binaryMap; fileName = 'Algal_Map_Left'; }
    else if (value === '2') { outputImg = indexMap; fileName = 'Index_Map_Left'; }
    else if (value === '3') { outputImg = chlaMap; fileName = 'Chla_Map_Left'; }
    else if (value === '4') { outputImg = binaryMap_right; fileName = 'Algal_Map_Right'; }
    else if (value === '5') { outputImg = indexMap_right; fileName = 'Index_Map_Right'; }
    else { outputImg = chlaMap_right; fileName = 'Chla_Map_Right'; }

    outScope = (value === '1' || value === '2' || value === '3') ? binaryMap.geometry() : binaryMap_right.geometry();
    var scope = outputImg.reduceToVectors({geometry: outScope, crs: 'EPSG:4326', scale: 30, eightConnected: true, bestEffort: true, maxPixels: 1e13});
    var vector = ee.FeatureCollection(scope);
    var kml_url = vector.getDownloadURL({format: 'kml', filename: fileName});
    right.download.KMLlink.setUrl(kml_url);
    right.download.KMLlink.style().set({shown: true});
  }
});

clearGeometry();
var dummyGeometry = ui.Map.GeometryLayer({geometries: null, name: 'geometry', color: '23cba7'});
drawingTools.layers().add(dummyGeometry);

collection = getData('left');
collection_right = getData('right');
if (collection) left.dataDate.chart = getImg(collection, x_left, false, 'left');
if (collection_right) right.dataDate.chart = getImg(collection_right, x_right, false, 'right');

left.controlPanel.add(left.info.panel);
left.controlPanel.add(left.divider);
left.controlPanel.add(left.selectAOI.panel);
left.controlPanel.add(left.dataDate.panel);
if (left.dataDate.chart) left.controlPanel.add(left.dataDate.chart);
left.controlPanel.add(left.dataDate.panel2);
left.controlPanel.add(left.computeButton);

right.controlPanel.add(right.dataDate.panel);
if (right.dataDate.chart) right.controlPanel.add(right.dataDate.chart);
right.controlPanel.add(right.dataDate.panel2);
right.controlPanel.add(right.computeButton);
right.controlPanel.add(right.divider);
right.controlPanel.add(right.download.label);
right.controlPanel.add(right.download.item);
right.controlPanel.add(right.download.KMLlink);

var s = {};
s.aboutText = {fontSize: '13px', color: '505050', width: '180px'};
s.titleLabel = {fontSize: '20px', fontWeight: 'bold'};
s.widgetTitle = {fontWeight: 'bold', margin: '8px 8px 0px 8px', color: '383838'};
s.paperLabel = {margin: '0px 0px', height: '20px'};
s.divider = {backgroundColor: 'B4B3B3', height: '3px', margin: '5px 0px'};
s.panel = {width: '220px', padding: '0px', margin: '0px 0px 8px 0px'};
s.stretchHorizontal = {stretch: 'horizontal'};
s.normalTitle = {fontWeight: 'bold'};

left.controlPanel.style().set(s.panel);
left.info.titleLabel.style().set(s.titleLabel);
left.info.aboutLabel.style().set(s.aboutText);
left.info.paperLabel1.style().set(s.paperLabel);
left.info.paperLabel2.style().set(s.paperLabel);
left.divider.style().set(s.divider);
left.selectAOI.label.style().set(s.widgetTitle);
left.selectAOI.aoi.style().set(s.stretchHorizontal);
left.dataDate.label1.style().set(s.normalTitle);
left.dataDate.slider1.style().set(s.stretchHorizontal);
left.dataDate.label2.style().set(s.normalTitle);
left.dataDate.slider2.style().set(s.stretchHorizontal);
left.dataDate.dataLabel.style().set(s.normalTitle);
left.dataDate.data.style().set(s.stretchHorizontal);
left.dataDate.indexLabel.style().set(s.normalTitle);
left.dataDate.indexSelect.style().set(s.stretchHorizontal);
left.dataDate.mosaicLabel.style().set(s.normalTitle);
left.dataDate.mosaicc.style().set(s.stretchHorizontal);
left.dataDate.thresholdLabel.style().set(s.normalTitle);
left.dataDate.chartLabel.style().set(s.normalTitle);
left.dataDate.trTextbox.style().set(s.stretchHorizontal);
left.dataDate.bufferLabel.style().set(s.normalTitle);
left.dataDate.bufferTextbox.style().set(s.stretchHorizontal);
left.computeButton.style().set(s.stretchHorizontal);

right.controlPanel.style().set(s.panel);
right.titleLabel.style().set(s.titleLabel);
right.dataDate.label1.style().set(s.normalTitle);
right.dataDate.slider1.style().set(s.stretchHorizontal);
right.dataDate.label2.style().set(s.normalTitle);
right.dataDate.slider2.style().set(s.stretchHorizontal);
right.dataDate.dataLabel.style().set(s.normalTitle);
right.dataDate.data.style().set(s.stretchHorizontal);
right.dataDate.indexLabel.style().set(s.normalTitle);
right.dataDate.indexSelect.style().set(s.stretchHorizontal);
right.dataDate.mosaicLabel.style().set(s.normalTitle);
right.dataDate.mosaicc.style().set(s.stretchHorizontal);
right.dataDate.thresholdLabel.style().set(s.normalTitle);
right.dataDate.chartLabel.style().set(s.normalTitle);
right.dataDate.trTextbox.style().set(s.stretchHorizontal);
right.dataDate.bufferLabel.style().set(s.normalTitle);
right.dataDate.bufferTextbox.style().set(s.stretchHorizontal);
right.computeButton.style().set(s.stretchHorizontal);
right.divider.style().set(s.divider);
right.download.label.style().set(s.normalTitle);
right.download.item.style().set(s.stretchHorizontal);

ui.root.clear();
ui.root.add(left.controlPanel);
ui.root.add(splitPanel);
ui.root.add(right.controlPanel);

leftMap.setCenter(-83.5379, 41.6528, 9);
rightMap.setCenter(-83.5379, 41.6528, 9);
