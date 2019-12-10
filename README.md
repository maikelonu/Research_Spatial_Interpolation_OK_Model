![alt test](/sshot-1.png)
# Spatial_Interpolation_OK

# https://doi.org/10.1016/j.proeng.2016.07.595

Geostatistical Ordinary Kriging Interpolation Method

Instituto Tecnologico de Costa Rica (www.tec.ac.cr)

Maikel Mendez-M (mamendez@itcr.ac.cr);(maikel.mendez@gmail.com)

Luis Alexander Calvo-V (lcalvo@itcr.ac.cr);(lualcava.sa@gmail.com)

Precipitation Interpolation Methods

This script is structured in R (www.r-project.org)

General purpose: Generate temporal series of average precipitation for a waterbasin using the geostatistical interpolation method of Ordinary Kriging (OK).

Input files: "interpolation.dat","calibration.dat", "validation.dat", "blank.asc"

Output files: "output_interpolation.csv"

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT FILE: attri.txt (File Name must be respected). TAB delimited.

DESCRIPTION: it contains spatial attributes needed to generate a georeference plus an interpolation threshold

Following is a brief description of the variables:

	Column1		Column2
Row0	ATTRIBUTE	VALUE
Row1	spa.res		value
Row2	xmin		value
Row3	xmax		value
Row4	ymin		value
Row5	ymax		value
Row6	goef.col	value
Row7	geof.row	value
Row8	threshold	value
Row9	nugget.all	value				 
Row10	range.all 	value
Row11   sill.all	value
Row12	nugget.cal	value				 
Row13	range.cal 	value
Row14   sill.cal	value
Row15   cores		value

Explanation:

Row0-Column1:	label "ATTRIBUTE" must be respected.
Row1-Column1:	label "spa.res" must be respected. Metric geographical projection. e.g.: UTM.
Row2-Column1:	label "xmin" must be respected. Metric geographical projection. e.g.: UTM.
Row3-Column1:	label "xmax" must be respected. Metric geographical projection. e.g.: UTM.
Row4-Column1:	label "ymin" must be respected. Metric geographical projection. e.g.: UTM.
Row5-Column1:	label "ymax" must be respected. Metric geographical projection. e.g.: UTM.
Row6-Column1:	label "goef.col" must be respected. Number of columns in georeference.
Row7-Column1:	label "geof.row" must be respected. Number of rows in georeference.
Row8-Column1:	label "threshold" must be respected. Precipitation threshold to execute interpolation (mm/Time). Values below threshold are not interpolated.
Row9-Column1:	label "nugget.all" must be respected. Nugget value for AutoKrige method when using all stations available. If the value is not specified (NA), it will be assigned by the method.
Row10-Column1:	label "range.all" must be respected. Range value for AutoKrige method when using all stations available. If the value is not specified (NA), it will be assigned by the method.
Row11-Column1:	label "sill.all" must be respected. Sill value for AutoKrige method when using all stations available. If the value is not specified (NA), it will be assigned by the method.
Row12-Column1:	label "nugget.cal" must be respected. Nugget value for AutoKrige method when using a calibration subset of the all stations. If the value is not specified (NA), it will be assigned by the method.
Row13-Column1:	label "range.cal" must be respected. Range value for AutoKrige method when using a calibration subset of the all stations. If the value is not specified (NA), it will be assigned by the method.
Row14-Column1:	label "sill.cal" must be respected. Sill value for AutoKrige method when using a calibration subset of the all stations. If the value is not specified (NA), it will be assigned by the method.
Row15-Column1:	label "cores" must be respected. Number of hardware cores for cluster.


Row0-Column2:	Integer or Double. 
Row1-Column2:	Integer or Double.
Row2-Column2:	Integer or Double.
Row3-Column2:	Integer or Double.
Row4-Column2:	Integer or Double.
Row5-Column2:	Integer or Double.
Row6-Column2:	Integer or Double.
Row7-Column2:	Integer or Double.
Row8-Column2:	Integer or Double. 
Row9-Column2:	Integer or Double.
Row10-Column2:	Integer or Double.
Row11-Column2:	Integer or Double.
Row12-Column2:	Integer or Double.
Row13-Column2:	Integer or Double.
Row14-Column2:	Integer or Double.
Row15-Column2:	Integer.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT FILE: interpolation.txt (File Name must be respected). TAB delimited.

DESCRIPTION: it contains daily or sub-daily precipitation (or any other meteorological variable) observations (records) at specific X, Y and Z points (meteorological stations). 
             Observations must be consistent through the modeled period and no missing values are allowed. A consistent georeference must be defined defined.
             This file contains ALL STATIONS available (or selected) to generate a RasterLayer interpolated surface. 
             Once interpolated, the average precipitation value is calculated and the RasterLayer itself is disregarded (not saved) within the main loop.

Following is a brief description of the variables:

	Column1	Column2		Column3		Columnn
Row0	Name	Name_EST1	Name_EST2	Name_EST3
Row1	X	UTM		UTM		UTM
Row2	Y	UTM		UTM		UTM
Row3	Z	MASL		MASL		MASL
Row4	Date1	prec.value	prec.value	prec.value
Row5	Date2	prec.value	prec.value	prec.value
Rown	Daten	prec.value	prec.value	prec.value

Explanation:

Row0-Column1:	label "Name" must be respected.
Row1-Column1:	label "X" must be respected. Metric geographical projection. e.g.: UTM.
Row2-Column1:	label "Y" must be respected. Metric geographical projection. e.g.: UTM.
Row3-Column1:	label "Z" must be respected. Metres above sea level. WGS84 geoid.
Row4-Column1:	julian date (day-month-year) e.g.: 01-01-90. So on up to Rown.

Row0-Column2:	station name. No blank spaces allowed. Use only "-" or "_". e.g.: Agua_Caliente. So on up to Columnn.
Row1-Column2:	precipitation value (or any other variable). Integer or Double. So on up to Columnn.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT FILE: calibration.txt (File Name must be respected). TAB delimited.

DESCRIPTION: it contains daily or sub-daily precipitation (or any other meteorological variable) observations (records) at specific X, Y and Z points (meteorological stations). 
	     Observations must be consistent through the modeled period and no missing values are allowed. A consistent georeference must be defined defined.
	     This file contains ONLY a calibration SUBSET of the ALL STATIONS available (or selected) to generate a RasterLayer interpolated surface. 
	     Usually 2/3 of the available station, say 7 out of 10. The same georefence must be used for interpolation.
	     Once interpolated, the average precipitation value is calculated and saved (using ONLY the calibration SUBSET) AND the RasterLayer is KEPT for 
	     further comparison within the main loop (say validation).

Following is a brief description of the variables: SAME as interpolation.txt. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT FILE: validation.txt (File Name must be respected). TAB delimited.

DESCRIPTION: it contains daily or sub-daily precipitation (or any other meteorological variable) observations (records) at specific X, Y and Z points (meteorological stations). 
	     Observations must be consistent through the modeled period and no missing values are allowed. A consistent georeference must be defined defined.
	     This file contains ONLY a validation SUBSET of the ALL STATIONS available (or selected) to generate a RasterLayer interpolated surface. 
	     Usually 1/3 of the available station, say 3 out of 10. The same georefence must be used for rasterization.
	     In this case, prec.values are NOT interpolated but ONLY rasterized. Once rasterized, these prec.values are confronted against the calibration RasterLayer interpolated 
	     surface (as same georeference has been used) and RMSE + MAE are calculated.

Following is a brief description of the variables: SAME as interpolation.txt. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
INPUT FILE: blank.asc (File Name must be respected). SPACE delimited.

DESCRIPTION: universal ESRI ASCII Grid format used to "blank" the interpolation georeference. It may represent boundary condition of a waterbasin or a specific sub-basin.
	     It can be generated by R libraries such as "sp" + "raster" or externally, using for example: ILWIS, QGIS, SAGA, GRASS, etc.

Following is a brief description of the variables: SAME as interpolation.txt. 

The file format is:xxx (Where xxx is a number.)
<NCOLS xxx>
<NROWS xxx>
<XLLCENTER xxx | XLLCORNER xxx>
<YLLCENTER xxx | YLLCORNER xxx>
<CELLSIZE xxx>
{NODATA_VALUE xxx}
row 1
row 2
row n

e.g.

ncols 138
nrows 95
xllcenter 461815
yllcenter 1124046
cellsize 100.00
nodata_value -32767
-32767 -32767 -32767 -32767 -32767 -3276
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
OUTPUT FILE: output_general.csv (File Name must be respected). COMMA delimited.

DESCRIPTION: Average interpolated precipitation value for the blanked georeference (mm/Time) according to Ordinary Kriging method.

NOTE: By using libraries; "doParallel" and "foreach" among others, all output_general.csv outputs are exported as CHARACTER.
For instance, before processing externally (e.g. MS-EXCEL) character (") must be removed in order to be recognized as NUMERIC.

Following is a brief description of the variables: 

       	Column1      Column2	   Column3	 Column4       Column5	       Column6
Row0	  ""	     "DATE"	   "OK.ALL"	  "OK.CAL"	"MAE.OK"      "RMSE.OK"	
Row1  result.1      numeric	   numeric	  numeric       numeric        numeric	
Row2  result.2      numeric	   numeric	  numeric       numeric        numeric
Rown  result.n      numeric	   numeric	  numeric       numeric        numeric	

Explanation:

Row0-Column1:	Sequence number of results.
Row0-Column2:	"DATE".    Julian date (day-month-year). Data class: character.
Row0-Column3:	"OK.ALL".  Rounded Ordinary Kriging value. Calculated using ALL stations available. Data class: numeric. 
Row0-Column4:	"OK.CAL".  Rounded Ordinary Kriging value. Calculated using ONLY a calibration SUBSET of the ALL STATIONS available. Data class: numeric. 
Row0-Column5:	"MAE.OK".  Mean Absolute Error value (mm) according to Ordinary Kriging method. Data class: numeric.
Row0-Column6:	"RMSE.OK". Root Mean Square Error value (mm) according to Ordinary Kriging method. Data class: numeric.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
OUTPUT FILE: descri_output.csv. COMMA delimited.

DESCRIPTION: Various descriptive statistics ABOUT the output_general.csv data.frame

        Column1 Column2 Column3 Column4
Row0    numeric numeric numeric numeric
Row1    numeric numeric numeric numeric
Row2    numeric numeric numeric numeric
Row3    numeric numeric numeric numeric
Row4    numeric numeric numeric numeric
Row5    numeric numeric numeric numeric
Row6    numeric numeric numeric numeric
Row7    numeric numeric numeric numeric
Row8    numeric numeric numeric numeric
Row9    numeric numeric numeric numeric
Row10   numeric numeric numeric numeric
Row11   numeric numeric numeric numeric
Row12   numeric numeric numeric numeric
Row13   numeric numeric numeric numeric


Following is a brief description of the variables:

Row0-Column1:   nbr.val. Number of values. So on up to Columnn.
Row1-Column1:   nbr.null. Number of null values. So on up to Columnn.
Row2-Column1:   nbr.na. Number of missing values. So on up to Columnn.
Row3-Column1:   min. Minimal value. So on up to Columnn.
Row4-Column1:   max. Maximal value. So on up to Columnn.
Row5-Column1:   range. Range (range, that is, max-min). So on up to Columnn.
Row6-Column1:   sum. Sum of all non-missing values. So on up to Columnn.
Row7-Column1:   median. Median. So on up to Columnn.
Row8-Column1:   mean. Mean. So on up to Columnn.
Row9-Column1:   SE.mean. Standard error on the mean. So on up to Columnn.
Row10-Column1:  CI.mean.0.95. Confidence interval of the mean at 0.95  p-Value level. So on up to Columnn.
Row11-Column1:  var. Variance. So on up to Columnn.
Row12-Column1:  std.dev. Standard deviation. So on up to Columnn.
Row13-Column1:  coef.var. Variation coefficient. So on up to Columnn.
