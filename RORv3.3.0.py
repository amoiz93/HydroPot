#Author: Abdul Moiz
#River & Environmental Engineering Laboratory (REEL)
#Department of Civil Engineering
#The University of Tokyo
#last Modified: 2017/03/01 9:06 PM


#Change Log
#RORv1.0.0 - Basic Framework (Only Main River), Head Optimization
#RORv2.0.0 - Basic Framework (Including Tributaries), Performs watershed delineation within this code
#RORv2.0.1 - Optimization of Catchment_Area*Head/Length_of_Waterway
#RORv2.0.2 - Introduction to Interference Tests (Only Main River), Exclusion of Areas in WDPA
#RORv2.0.3 - Introduction to Interference Tests (Extension to Tributaries)
#RORv2.0.4 - CreatePointsLines tool has been incorporated into the script as a function. There is no need to use Supplementary Tools as an input
#RORv3.0.0 - Linked to ArcGIS Tool interface (First version in the Tool series, all the successive versions shall be based on this)
#RORv3.1.0 - Beta testing version
#RORv3.2.0 - 3 modes introduced (Note: Mode 3 hasn't been tested yet)
#RORv3.3.0 - Promethee Method II introduced for site selection optimimization (The algorithm had to be majorly restructured accordingly)(Note: Slope has been excluded as the minimum criterion for site selection)
#RORv3.3.0 - Minor modification to inlcude lvl 2 pfafstetter delineation


#Inputs
# 1 - DEM (Preferrably Clipped to the extent of the watershed)
# 2 - Basin Outlet (Point shapefile)
# 3 - Pfafstetter Subbasin (Raster, GRID file from WEBDHM Data Preprocessing)
# 4 - Flow Interval (Raster, GRID file from WEBDHM Data Preprocessing, from flow.asc)
# 5 - World Database for Protected Areas (Polygon shapefile, Preferrably Clipped to the extent of the watershed)
# 6 - Output river_folder from WEB-DHM (Optional)


#Criteria
# 1 - Threshold for stream delineation - Number of Cells
# 2 - Minimum Interval for Point Generation (Also, the minimum distance between consecutive schemes) - Meters
# 3 - Outer Radius of Annular Search Zone - Meters
# 4 - Inner Radius of Annular Search Zone - Meters
# 5 - Minimum Gross Head - Meters
# 6 - Minimum Slope - Percentage
# 7 - Minimum FLow Accumulation - Number of Cells
# 8 - Buffer zone around protected Areas - Meters
# 9 - Number of Validation Points (used in WEBDHM)

#Program Runtime
start_time = time.time()
#============================Importing Modules=================================
#print "Importing modules..."
#import archook
#archook.get_arcpy()
import arcpy
import math
import pandas as pd
import numpy as np
import os
import shutil
#==============================================================================

#=========================Defining Custom Functions============================
#------------------------------------------------------------------------------
# Purpose:     Creates points on lines at a specified distance, percentage, or
#              interval using a fixed or field-based value. All polyline fields
#              are included in the output Point feature class. Points can be
#              created starting from the beginning, or end of the line.
#
# Author:      Ian Broad
# Website:     www.ianbroad.com
#
# Created:     05/04/2016
#------------------------------------------------------------------------------

def CreatePointsLines(polyline, choice, start_from, use_field_for_value, field_with_value, distance, end_points, output):

    spatial_ref = arcpy.Describe(polyline).spatialReference
    
    mem_point = arcpy.CreateFeatureclass_management("in_memory", "mem_point", "POINT", "", "DISABLED", "DISABLED", polyline)
    arcpy.AddField_management(mem_point, "LineOID", "LONG")
    arcpy.AddField_management(mem_point, "Value", "FLOAT")
    
    result = arcpy.GetCount_management(polyline)
    features = int(result.getOutput(0))
    
    arcpy.SetProgressor("step", "Creating Points on Lines...", 0, features, 1)
    
    search_fields = ["SHAPE@", "OID@"]
    insert_fields = ["SHAPE@", "LineOID", "Value"]
    
    if use_field_for_value == "YES":
        search_fields.append(field_with_value)
    
    reverse_line = False
    if start_from == "END":
        reverse_line = True
    
    with arcpy.da.SearchCursor(polyline, (search_fields)) as search:
        with arcpy.da.InsertCursor(mem_point, (insert_fields)) as insert:
            for row in search:
                try:
                    line_geom = row[0]
                    length = float(line_geom.length)
                    count = distance
                    oid = str(row[1])
                    start = arcpy.PointGeometry(line_geom.firstPoint)
                    end = arcpy.PointGeometry(line_geom.lastPoint)
    
                    if reverse_line == True:
                       reversed_points = []
                       for part in line_geom:
                           for p in part:
                               reversed_points.append(p)
    
                       reversed_points.reverse()
                       array = arcpy.Array([reversed_points])
                       line_geom = arcpy.Polyline(array, spatial_ref)
    
                    if use_field_for_value == "YES":
                        count = float(row[2])
                        distance = float(row[2])
    
                    ################################################################
    
                    if choice == "DISTANCE":
                        point = line_geom.positionAlongLine(count, False)
                        insert.insertRow((point, oid, count))
    
                    elif choice == "PERCENTAGE":
                        point = line_geom.positionAlongLine(count, True)
                        insert.insertRow((point, oid, count))
    
                    elif choice == "INTERVAL BY DISTANCE":
                        while count <= length:
                            point = line_geom.positionAlongLine(count, False)
                            insert.insertRow((point, oid, count))
    
                            count += distance
    
                    elif choice == "INTERVAL BY PERCENTAGE":
                        percentage = float(count * 100.0)
                        total_runs = int(100.0 / percentage)
    
                        run = 1
                        while run <= total_runs:
                            current_percentage = float((percentage * run) / 100.0)
                            point = line_geom.positionAlongLine(current_percentage, True)
    
                            insert.insertRow((point, oid, current_percentage))
    
                            run += 1
    
                    elif choice == "START/END POINTS":
                        insert.insertRow((start, oid, 0))
                        insert.insertRow((end, oid, length))
    
                    ################################################################
    
                    if end_points == "START":
                        insert.insertRow((start, oid, 0))
    
                    elif end_points == "END":
                        insert.insertRow((end, oid, length))
    
                    elif end_points == "BOTH":
                        insert.insertRow((start, oid, 0)) #edited (Interchanged 0<->length)
                        insert.insertRow((end, oid, length)) #edited (Interchanged 0<->length)
    
                    arcpy.SetProgressorPosition()
    
                except Exception as e:
                    arcpy.AddMessage(str(e.message))
    
    
                    ################################################################
    
    
    line_keyfield = str(arcpy.ListFields(polyline, "", "OID")[0].name)
    
    mem_point_fl = arcpy.MakeFeatureLayer_management(mem_point, "Points_memory")
    
    arcpy.AddJoin_management(mem_point_fl, "LineOID", polyline, line_keyfield)
    
    if "in_memory" in output:
        arcpy.SetParameter(8, mem_point_fl)
    
    else:
        arcpy.CopyFeatures_management(mem_point_fl, output)
    
        arcpy.Delete_management(mem_point)
        arcpy.Delete_management(mem_point_fl)
    
    arcpy.ResetProgressor()
    arcpy.GetMessages()
    return str(output)

#------------------------------------------------------------------------------
# Purpose:          Calls the discharge value/range for a specific flow interval 
#                   at a specific time/time_range
# Author:           Abdul Moiz
# Organization:     River & Environment Engineering Laboratory (REEL)
#
# Created:          2017/05/15
#------------------------------------------------------------------------------

def flowint_to_Q(flowint_folder,subbasin,flowint,year,month,day,hr,time_range,year1,month1,day1,hr1):
    subbasin_str = str(subbasin)
    if len(subbasin_str) == 1:
        subbasin = subbasin*100
    elif len(subbasin_str) == 2:
        subbasin = subbasin*10
    elif len(subbasin_str) == 3:
        subbasin = subbasin*1
    flowint_df = pd.read_pickle(os.path.join(flowint_folder,'df'+str(subbasin)+'.pkl'))
    if time_range == False:
        time = pd.datetime(year,month,day,hr)
        flowint_q = flowint_df[flowint][time]
    if time_range ==  True:
        time = pd.datetime(year,month,day,hr)
        time1 = pd.datetime(year1,month1,day1,hr1)
        flowint_q = flowint_df[flowint][time:time1]
    return flowint_q


#------------------------------------------------------------------------------
#Purpose:   Calculates the discharge at a particular percentage flow exceedance 
#           and optionally plots the flow duration curves for a particular
#           Sub-Basin and Flow Interval
#
#Author:    Abdul Moiz
#------------------------------------------------------------------------------
def FlowDurationCurve(flowint_folder, subbasin, flowint, input_start_date, input_end_date, flow_duration):
    input_start_date_split = input_start_date.split('/')
    start_day = int(input_start_date_split[1])
    start_month = int(input_start_date_split[0])
    start_year = int(input_start_date_split[2])
    print start_day,start_month,start_year
    
    input_end_date_split = input_end_date.split('/')
    end_day = int(input_end_date_split[1])
    end_month = int(input_end_date_split[0])
    end_year = int(input_end_date_split[2])
    print end_day,end_month,end_year
    
    sim_df = flowint_to_Q(flowint_folder, subbasin, flowint ,start_year, start_month, start_day, 0, True, end_year, end_month, end_day, 23)
    S = sim_df.resample('D').mean()
    S = S[input_start_date:input_end_date]
    S = S.groupby(by=[S.index.day,S.index.month]).mean()
    
    S.sort_values(ascending=False, inplace=True, kind='mergesort')
    S_rank = S.rank(ascending=False)
    n = len(S_rank)
    P = 100*(S_rank/(n+1))
    database = pd.concat([P,S], axis=1, ignore_index=True)
    database.rename(index = None, columns = {0:"Flow Exceedance (%)", 1:"Discharge (cumecs)"}, inplace=True)
    
    database = database.reset_index()
    flow_exceedances = [10,20,30,40,50,60,70,80,90,95]
    Qx_list=[]
    for flow_exceedance in flow_exceedances:
        if flow_exceedance!=50:
    
            rank = float(flow_exceedance)*(float(n)+1.0)/100.0
            upperrank = int(math.ceil(rank))-1
            lowerrank = int(math.floor(rank))-1
        
            x = flow_exceedance
            
            x1 = database.at[lowerrank,"Flow Exceedance (%)"]
            x2 = database.at[upperrank,"Flow Exceedance (%)"]
            y1 = database.at[lowerrank,"Discharge (cumecs)"]
            y2 = database.at[upperrank,"Discharge (cumecs)"]
            Qx=((y2-y1)/(x2-x1))*(x-x1)+y1
    
        else:
            Qx = database.at[int(n/2), "Discharge (cumecs)"]
        Qx=round(Qx,5)
        Qx_list.append(Qx)
    
    if flow_duration=='ALL':
        return Qx_list
    elif int(flow_duration)==10:
        return [Qx_list[0]]
    elif int(flow_duration)==20:
        return [Qx_list[1]]
    elif int(flow_duration)==30:
        return [Qx_list[2]]
    elif int(flow_duration)==40:
        return [Qx_list[3]]
    elif int(flow_duration)==50:
        return [Qx_list[4]]
    elif int(flow_duration)==60:
        return [Qx_list[5]]
    elif int(flow_duration)==70:
        return [Qx_list[6]]
    elif int(flow_duration)==80:
        return [Qx_list[7]]
    elif int(flow_duration)==90:
        return [Qx_list[8]]
    elif int(flow_duration)==95:
        return [Qx_list[9]]
    

#------------------------------------------------------------------------------
def createfolder(folderName):
    if os.path.isdir(folderName)==True:
        shutil.rmtree(folderName)
        os.mkdir(folderName)
    else:
        os.mkdir(folderName)
#------------------------------------------------------------------------------
def rootpath(filepath):
    string=filepath.split('\\')
    folderpath=''
    for x in range(0,len(string)-1):
        folderpath=folderpath+str(string[x])+"\\"
    return folderpath
    del string
    del x
    del folderpath
#------------------------------------------------------------------------------
def AddFeatureLayer(shp_file_path):
    mxd = arcpy.mapping.MapDocument("CURRENT")
    df = arcpy.mapping.ListDataFrames(mxd, '*')[0]
    newlayer = arcpy.mapping.Layer(shp_file_path)
    arcpy.mapping.AddLayer(df, newlayer, 'TOP')
    arcpy.RefreshActiveView()
    arcpy.RefreshTOC()
    del mxd, df, newlayer
    
#------------------------------------------------------------------------------
#Purpose: Promethee II Outranking 
#Calcuates the net_flow to rank the alternatives according to the specified weights.
#Uses Promethee Method II
#
#Author: Abdul Moiz
#------------------------------------------------------------------------------
def promethee(alternative_matrix, max_column, min_column, p, q, weights):
    for column in alternative_matrix:
        if column in max_column:
            max_Aij = alternative_matrix[column].max()
            min_Aij = alternative_matrix[column].min()
            Aij = alternative_matrix[column]
            if max_Aij-min_Aij == 0:
                Rij = 1.0
            else:
                Rij = (Aij-min_Aij)/(max_Aij-min_Aij)
            alternative_matrix[column]=Rij
        elif column in min_column:
            max_Aij = alternative_matrix[column].max()
            min_Aij = alternative_matrix[column].min()
            Aij = alternative_matrix[column]
            if max_Aij-min_Aij == 0:
                Rij = 1.0
            else:
                Rij = (max_Aij-Aij)/(max_Aij-min_Aij)
            alternative_matrix[column]=Rij
                              
    normalized_matrix = alternative_matrix
    index_sites = range(0,len(normalized_matrix.index))
    aggregated_preference_matrix = pd.DataFrame(index=index_sites)
    
    #h,g are counter assingers for preference function table
    for a in range(0,len(normalized_matrix.index)):
        index1 = range(0,len(normalized_matrix.index)-1)
        column1 = range(0,len(normalized_matrix.columns))
        preference_function_matrix = pd.DataFrame(index=index1, columns=column1)
        h=0
        alternatives = range(0,len(normalized_matrix.index))
        alternatives = alternatives[:a] + alternatives[a+1:]
        for j in normalized_matrix:
            g=0
            for b in alternatives:
                A = normalized_matrix[j][a]
                B = normalized_matrix[j][b]
                D=A-B
                if D<=0:
                    Pref=0
                elif q<D and D<=p:
                    Pref = (D-q)/(p-q)
                elif D>p:
                    Pref = 1    #Add Preference Function Here
                
                preference_function_matrix[h][g] = Pref
                g=g+1
            h=h+1
        weight_counter = 0
        for weight in weights:
            preference_function_matrix[weight_counter] = preference_function_matrix[weight_counter]*(weights[weight_counter]/float(sum(weights)))
            weight_counter=weight_counter+1
        pie = preference_function_matrix.sum(1)
        pie.index = alternatives
        aggregated_preference_matrix = pd.concat([aggregated_preference_matrix, pie], axis=1)
    aggregated_preference_matrix.columns = index_sites
    aggregated_preference_matrix = aggregated_preference_matrix.transpose()
    final_outranking = pd.DataFrame(index = index_sites, columns=['Positive_Flow', 'Negative_Flow', 'Net_Flow'])
    for a in range(0,len(normalized_matrix.index)):
        final_outranking['Positive_Flow'][a] = (aggregated_preference_matrix.sum(1))[a]/(len(normalized_matrix.index)-1)
        final_outranking['Negative_Flow'][a] = (aggregated_preference_matrix.sum(0))[a]/(len(normalized_matrix.index)-1)
        final_outranking['Net_Flow'][a] = final_outranking['Positive_Flow'][a] - final_outranking['Negative_Flow'][a]
    return final_outranking
    del column, max_Aij, min_Aij, Aij, Rij, normalized_matrix, index_sites, aggregated_preference_matrix, a, index1, column1, preference_function_matrix, h, g, alternatives, j, b, A, B, D, Pref, weight_counter, weight, weights, pie, final_outranking
#==============================================================================

#====================Setting Environment Variables=============================
arcpy.AddMessage("Setting enviornment variables...")
mxd=arcpy.mapping.MapDocument("CURRENT")
wkspace = rootpath(mxd.filePath)
os.chdir(wkspace)
arcpy.env.workspace = wkspace
arcpy.env.scratchWorkspace = wkspace
arcpy.env.overwriteOutput = True
#==============================================================================


#==============================Output Paths====================================
output_folder = os.path.join(wkspace, 'Output')
preprocessing_folder = os.path.join(output_folder, 'preprocessing')
raster_folder = os.path.join(output_folder, 'raster')
shape_folder = os.path.join(output_folder, 'shape')
txt_folder = os.path.join(output_folder, 'text')
temp_folder = os.path.join(wkspace, 'temp')

txt_output = os.path.join(txt_folder, 'Results.txt')                        

hydropower_plants_output = os.path.join(shape_folder, 'hydropower_plant.shp')
dam_site_output = os.path.join(shape_folder, 'diversion_site.shp')          
schemes_output = os.path.join(shape_folder, 'schemes.shp')
#==============================================================================


#========================Checking out extensions===============================
arcpy.AddMessage("Checking out extensions...")
extensions = ["Spatial", "3D"]
for extension in extensions:
    if arcpy.CheckExtension(extension) == "Available":
        arcpy.AddMessage("Checking out " + extension + " Analyst...")
        arcpy.CheckOutExtension(extension)
    else:
        arcpy.AddError(extension + " Analyst license is unavailable")
#==============================================================================    

#==================================Inputs======================================
inDEM = arcpy.Raster(arcpy.GetParameterAsText(0))
outlet = arcpy.GetParameterAsText(1)

threshold = int(arcpy.GetParameterAsText(2))
min_interval = int(arcpy.GetParameterAsText(3))     
out_rad = int(arcpy.GetParameterAsText(4))   
in_rad = int(arcpy.GetParameterAsText(5))
min_head = float(arcpy.GetParameterAsText(6))       
min_slope = float(arcpy.GetParameterAsText(7))
min_FAC = float(arcpy.GetParameterAsText(8))


Protected_Area_Switch = arcpy.GetParameterAsText(9)
inWDPA = arcpy.GetParameterAsText(10)                                     
Protected_Area_Buffer = arcpy.GetParameterAsText(11)                                 

mode = arcpy.GetParameterAsText(12)

inSubbasin = arcpy.GetParameterAsText(13)                               
inFlowInt = arcpy.GetParameterAsText(14)                                  
river_folder = arcpy.GetParameterAsText(15)                                 

input_start_date = arcpy.GetParameterAsText(16)
input_end_date = arcpy.GetParameterAsText(17)

arcpy.AddMessage(input_start_date)
arcpy.AddMessage(input_end_date)

system_efficiency = arcpy.GetParameterAsText(18)
system_efficiency = float(int(system_efficiency)/100.0)

flow_exceedance = arcpy.GetParameterAsText(19)

weight_head = arcpy.GetParameterAsText(20)
weight_FAC = arcpy.GetParameterAsText(21)
weight_waterway = arcpy.GetParameterAsText(22)
weight_hydropower = arcpy.GetParameterAsText(23)

weights_mode_0 = [float(weight_head), float(weight_FAC), float(weight_waterway)]
weights_mode_2 = [float(weight_hydropower), float(weight_waterway)]

minimum_hydropower_potential = arcpy.GetParameterAsText(24)
maximum_hydropower_potential = arcpy.GetParameterAsText(25)
#==============================================================================


#================================Properties====================================
resolution = inDEM.meanCellHeight
cellsize = inDEM.meanCellHeight*inDEM.meanCellWidth
arcpy.AddMessage("DEM Resolution is " + str(resolution) + " meters")
arcpy.AddMessage("Cell area is " + str(cellsize/1000000) + " sq. km")
spatial_ref = arcpy.Describe(inDEM).spatialReference
#==============================================================================


#================================Counters======================================
WDPA_intersection_counter = 0
#==============================================================================


#========================Creating Necessary Folders============================
createfolder(output_folder)
createfolder(preprocessing_folder)
createfolder(raster_folder)
createfolder(shape_folder)
createfolder(txt_folder)
createfolder(temp_folder)
#==============================================================================

#=============================DEM Preprocessing================================
arcpy.AddMessage("Filling Sinks...")
outFill = arcpy.sa.Fill(inDEM)

arcpy.AddMessage("Generating Flow Direction...")
outFDR = arcpy.sa.FlowDirection(outFill)

arcpy.AddMessage("Generating Flow Accumulation...")
outFAC = arcpy.sa.FlowAccumulation(outFDR)

arcpy.AddMessage("Delineating Watershed...")
outbasin = arcpy.sa.Watershed(outFDR, outlet)
outbasin.save(os.path.join(preprocessing_folder, 'watershed.tif'))
watershed = arcpy.RasterToPolygon_conversion(os.path.join(preprocessing_folder, 'watershed.tif'), os.path.join(preprocessing_folder, 'watershed.shp'), "NO_SIMPLIFY")

outDEM = arcpy.sa.ExtractByMask(inDEM, outbasin)
outDEM.save(os.path.join(preprocessing_folder, 'DEM.tif'))

outFill = arcpy.sa.ExtractByMask(outFill, outbasin)
outFill.save(os.path.join(preprocessing_folder, 'Fill.tif'))

outFDR = arcpy.sa.ExtractByMask(outFDR, outbasin)
outFDR.save(os.path.join(preprocessing_folder, 'FDR.tif'))

outFAC = arcpy.sa.ExtractByMask(outFAC, outbasin)
outFAC.save(os.path.join(preprocessing_folder, 'FAC.tif'))

arcpy.AddMessage("Generating Stream Network...")

stream = arcpy.sa.Con(outFAC > int(threshold), 1)

arcpy.sa.StreamToFeature(stream, outFDR, \
                         os.path.join(preprocessing_folder, 'stream.shp'), "NO_SIMPLIFY")
#==============================================================================


#===========================WDPA Preprocessing=================================

arcpy.AddMessage("Buffering World Database for Protected Areas...")

if Protected_Area_Switch == 'true':
    WDPA_Buffer = arcpy.Buffer_analysis(inWDPA, os.path.join(preprocessing_folder, 'WDPA_Buffer.shp'), 
                                        int(Protected_Area_Buffer), '#', '#', 'ALL', '#', '#')

#==============================================================================


#==========================Strahler Stream Order===============================
arcpy.AddMessage("Ordering Streams...")
strahler = arcpy.sa.StreamOrder(stream, outFDR)
strahler.save(os.path.join(preprocessing_folder, 'strahler.tif'))

streamorders = []
cursorA = arcpy.da.SearchCursor(os.path.join(preprocessing_folder, 'strahler.tif'), "VALUE")
for rowA in cursorA:
    streamorders.append(rowA[0])
#==============================================================================

#======================ROR Site Selection (MODE==0)============================
#Note: i, j and k are variables to keep track of the optimization procedures

if mode == "Approach (i): Using only Topographic Factors (Only Site Selection)":
    arcpy.AddMessage("Initializing Run-of-River Site Selection Module (" + mode +")")
    
    #Creating Empty Results DataFrame
    columns = ['D/S_ID', 'U/S_ID', 'D/S_X', 'D/S_Y', 'D/S_Z', 'U/S_X', 'U/S_Y', 'U/S_Z', 
               'Headrace_Length_3D(m)', 'Headrace_Length_2D(m)', 'Head(m)', 'Slope(%)', 
                'FAC_at_U/S', 'Stream_Segment', 'Stream_Order', 'Chain_DS', 'Chain_US']
                
    results_df = pd.DataFrame([], columns = columns)
    
    for streamorder in streamorders:
        order_select = arcpy.sa.ExtractByAttributes(strahler, "VALUE = " + str(streamorder))
        order_select.save(os.path.join(temp_folder, 'strahler_select.tif'))
        order = arcpy.sa.StreamToFeature(order_select, outFDR, os.path.join(temp_folder, 'strahler_select.shp'), "NO_SIMPLIFY")
        streamsegs = []
        cursorB = arcpy.da.SearchCursor(order, "FID")
        for rowB in cursorB:
            streamsegs.append(rowB[0])
        
        for streamseg in streamsegs:
            streamseg_select = arcpy.Select_analysis(order, os.path.join(temp_folder, 'streamseg_select.shp'), "FID = " + str(streamseg))
            inPoints = CreatePointsLines(streamseg_select, "INTERVAL BY DISTANCE", 
                                                           "END", "NO", '', float(min_interval), "NO", os.path.join(temp_folder, 'densified_stream.shp'))
            
            total_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            if Protected_Area_Switch == 'true':
                inPoints_WDPA = arcpy.Erase_analysis(inPoints, WDPA_Buffer, os.path.join(temp_folder, 'wdpa_stream.shp'))
                arcpy.Delete_management(inPoints)
                inPoints = inPoints_WDPA
                
            remaining_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            excluded_searching_points = total_searching_points - remaining_searching_points
            WDPA_intersection_counter += excluded_searching_points
            
            arcpy.AddXY_management(inPoints)
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outFAC, "FAC"]])
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outDEM, "DEM"]], 'BILINEAR')
            inPoints = arcpy.Sort_management(inPoints, os.path.join(temp_folder, 'densified_stream_sorted.shp'), [["mem_point1", "ASCENDING"]])
    
            streampoints = []
            cursorC = arcpy.da.SearchCursor(inPoints, "FID")
            for rowC in cursorC:
                streampoints.append(rowC[0])
            
            count = len(streampoints)
            
            i=0
            
            while i < count:
                arcpy.AddMessage("Processing SO " + str(streamorder) + 
                      " of " + str(streamorders[-1]) + " SS " +  str(streamseg + 1) + 
                      " of " + str(streamsegs[-1] + 1) + " Point " + str(i + 1) + " of " + str(count))
    
    
                #Accessing Hydropower Plant Data Base
                with arcpy.da.SearchCursor(inPoints, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM"], "FID = " + str(i)) as cursor1:
                    for row1 in cursor1:
                        Orig_ID_1 = row1[0]
                        Chainage_1 = row1[1]
                        x1 = row1[2]
                        y1 = row1[3]
                        fac1 = row1[4]*(cellsize)/1000000.0
                        z1 = row1[5]
                        segment = int(streamseg + 1)
                        
                        
                        #Selecting Hydropower Plant Point       
                        focus_point = arcpy.Select_analysis(inPoints, os.path.join(temp_folder, 'focus_point_' + str(i) +'.shp'), "FID = " + str(i))
                        
                        #Creating Search Radius
                        out_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'out_buf_' + str(i) + '.shp'), out_rad, method = "PLANAR")
                        in_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'in_buf_' + str(i) + '.shp'), in_rad, method = "PLANAR")
                        search_rad = arcpy.Erase_analysis(out_buf, in_buf, os.path.join(temp_folder, 'search_rad_' + str(i) + '.shp'))
                
                        #Creating Potential Diversion Points Feature Class
                        pot_div_points = arcpy.Clip_analysis(inPoints, search_rad, os.path.join(temp_folder, 'pot_div_' + str(i) + '.shp'))
                        
                        #Counting Total Alternatives that satisfy Minimum Criteria
                        total_alternatives = 0
                        with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM"]) as cursor2:
                            for row2 in cursor2:
                                Orig_ID_2 = row2[0]
                                Chainage_2 = row2[1]
                                x2 = row2[2]
                                y2 = row2[3]
                                fac2 = row2[4]*(cellsize)/1000000.0
                                z2 = row2[5]
            
                                head = z2-z1
                                distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                slope = head/distance2D
                                slope_percentage = slope*100
                                head_slope = head*slope_percentage
                                
                                if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC):
                                    total_alternatives = total_alternatives + 1
                        
                        #Accessing Diversion Site Database
                        if total_alternatives == 0:
                            arcpy.AddMessage('No Suitable Diversion Site')
                            
                        elif total_alternatives == 1:
                            Chainage_3 = Chainage_2
                            Orig_ID_3 = Orig_ID_2
                            x3 = x2
                            y3 = y2
                            fac3 = fac2
                            z3 = z2
                            
                            head_fin = head
                            distance3D_fin = distance3D
                            distance2D_fin = distance2D
                            slope_fin = slope
                            slope_percentage_fin = slope_percentage
                            head_slope_fin = head_slope
                            j = int((Chainage_3/min_interval)-1)
    
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " sq. km.")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " m")
        
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                            'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                            'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                            'FAC_at_U/S':[fac3], 'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                            'Chain_US':[Chainage_3]}, columns=columns)
                        
                            results_df = pd.concat([results_df, df1])
                        
                        else:
                            alternatives = pd.DataFrame(index = range(0,total_alternatives), columns=["Chainage", "ID", "X", "Y", "Z", "Slope", "Distance2D", "Head", "Flow_Accumulation", "Waterway"])
                            alternatives_counter = 0
                            with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM"]) as cursor2:
                                for row2 in cursor2:
                                    Orig_ID_2 = row2[0]
                                    Chainage_2 = row2[1]
                                    x2 = row2[2]
                                    y2 = row2[3]
                                    fac2 = row2[4]*(cellsize)/1000000.0
                                    z2 = row2[5]
            
                                    head = z2-z1
                                    distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                    distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                    slope = head/distance2D
                                    slope_percentage = slope*100
                                    head_slope = head*slope_percentage
                                
                                    if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC):
                                        Orig_ID_3 = row2[0]
                                        Chainage_3 = row2[1]
                                        x3 = row2[2]
                                        y3 = row2[3]
                                        fac3 = row2[4]*(cellsize)/1000000.0
                                        z3 = row2[5]
                                
                                        head_fin = z3-z1
                                        distance3D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2)+((z1-z3)**2))
                                        distance2D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2))
                                        slope_fin = head_fin/distance2D_fin
                                        slope_percentage_fin = slope_fin*100
                                        head_slope_fin = head_fin*slope_percentage_fin
                                        
                                        #Creating the input for Promethee Module
                                        alternatives["Chainage"][alternatives_counter] = Chainage_3
                                        alternatives["ID"][alternatives_counter] = Orig_ID_3
                                        alternatives["X"][alternatives_counter] = x3
                                        alternatives["Y"][alternatives_counter] = y3
                                        alternatives["Z"][alternatives_counter] = z3
                                        alternatives["Slope"][alternatives_counter] = slope_percentage_fin
                                        alternatives["Distance2D"][alternatives_counter] = distance2D_fin
                                        alternatives["Head"][alternatives_counter] = head_fin
                                        alternatives["Flow_Accumulation"][alternatives_counter] = fac3
                                        alternatives["Waterway"][alternatives_counter] = distance3D_fin
                                        
                                        alternatives_counter = alternatives_counter + 1
                        
                            promethee_alternatives = alternatives[["Head","Flow_Accumulation","Waterway"]]
                            outranking_flows = promethee(promethee_alternatives, ["Head", "Flow_Accumulation"], ["Waterway"], 1, 0, weights_mode_0)
                            outranking_flows.index = alternatives.index
                            final_alternatives = pd.concat([alternatives,outranking_flows], axis=1)
                            final_alternatives.sort(["Net_Flow"], ascending=False, inplace=True)
                            final_alternatives.reset_index(inplace=True, drop=True)
                                
                            Chainage_3 = final_alternatives["Chainage"][0]
                            Orig_ID_3 = final_alternatives["ID"][0]
                            x3 = final_alternatives["X"][0]
                            y3 = final_alternatives["Y"][0]
                            fac3 = final_alternatives["Flow_Accumulation"][0]
                            z3 = final_alternatives["Z"][0]
                                
                            head_fin = final_alternatives["Head"][0]
                            distance3D_fin = final_alternatives["Waterway"][0]
                            distance2D_fin = final_alternatives["Distance2D"][0]
                            slope_fin = head_fin/distance2D_fin
                            slope_percentage_fin = slope_fin*100
                            head_slope_fin = head_fin*slope_percentage_fin
                            j = int((Chainage_3/min_interval)-1)
                                
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " sq. km.")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " m")
    
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                        'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                        'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                        'FAC_at_U/S':[fac3], 'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                        'Chain_US':[Chainage_3]}, columns=columns)
                            results_df = pd.concat([results_df, df1])
                            
                        arcpy.Delete_management(focus_point)
                        arcpy.Delete_management(out_buf)
                        arcpy.Delete_management(in_buf)
                        arcpy.Delete_management(search_rad)
                        arcpy.Delete_management(pot_div_points)
                        
                        if arcpy.Exists(os.path.join(temp_folder, 'pot_div.shp'))==True:
                            arcpy.Delete_management(os.path.join(temp_folder, 'pot_div.shp'))
                
                i+=1
                arcpy.AddMessage('\n')
            arcpy.Delete_management(inPoints)
            arcpy.Delete_management(os.path.join(temp_folder, 'densified_stream.shp'))
            arcpy.Delete_management(streamseg_select)
            
        arcpy.Delete_management(order_select)
        arcpy.Delete_management(order)
  
    results_df.reset_index(drop = True, inplace = True)
    results_df = results_df[columns]
    results_df.index = np.arange(1, len(results_df)+1)
    
    alternative_schemes_counter = int(len(results_df.index))
#==============================================================================

#======================ROR Site Selection (MODE==1)============================
#Note: i, j and k are variables to keep track of the optimization procedures

elif mode == "Approach (i): Using only Topographic Factors (Hydropower Potential Evaluation)":
    arcpy.AddMessage("Initializing Run-of-River Site Selection Module (" + mode + ")")
    
    #Creating Empty Results DataFrame
    columns = ['D/S_ID', 'U/S_ID', 'D/S_X', 'D/S_Y', 'D/S_Z', 'U/S_X', 'U/S_Y', 'U/S_Z', 
               'Headrace_Length_3D(m)', 'Headrace_Length_2D(m)', 'Head(m)', 'Slope(%)', 
                'FAC_at_U/S', 'Flow_Interval_U/S', 'Subbasin_U/S', 'Stream_Segment', 
                'Stream_Order', 'Chain_DS', 'Chain_US', 'Q10', 'Q20', 'Q30', 'Q40', 'Q50', 'Q60', 'Q70', 'Q80', 'Q90', 'Q95',
                'Hydro10', 'Hydro20', 'Hydro30', 'Hydro40', 'Hydro50', 'Hydro60', 'Hydro70', 'Hydro80', 'Hydro90', 'Hydro95']

    results_df = pd.DataFrame([], columns = columns)
    
    for streamorder in streamorders:
        order_select = arcpy.sa.ExtractByAttributes(strahler, "VALUE = " + str(streamorder))
        order_select.save(os.path.join(temp_folder, 'strahler_select.tif'))
        order = arcpy.sa.StreamToFeature(order_select, outFDR, os.path.join(temp_folder, 'strahler_select.shp'), "NO_SIMPLIFY")
        streamsegs = []
        cursorB = arcpy.da.SearchCursor(order, "FID")
        for rowB in cursorB:
            streamsegs.append(rowB[0])
        
        for streamseg in streamsegs:
            streamseg_select = arcpy.Select_analysis(order, os.path.join(temp_folder, 'streamseg_select.shp'), "FID = " + str(streamseg))
            inPoints = CreatePointsLines(streamseg_select, "INTERVAL BY DISTANCE", 
                                                           "END", "NO", '', float(min_interval), "NO", os.path.join(temp_folder, 'densified_stream.shp'))
            
            total_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            if Protected_Area_Switch == 'true':
                inPoints_WDPA = arcpy.Erase_analysis(inPoints, WDPA_Buffer, os.path.join(temp_folder, 'wdpa_stream.shp'))
                arcpy.Delete_management(inPoints)
                inPoints = inPoints_WDPA
                
            remaining_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            excluded_searching_points = total_searching_points - remaining_searching_points
            WDPA_intersection_counter += excluded_searching_points
            
            arcpy.AddXY_management(inPoints)
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outFAC, "FAC"]])
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outDEM, "DEM"]], 'BILINEAR')
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[inSubbasin, "SUBBASIN"]])
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[inFlowInt, "FLOW_INT"]])            
            inPoints = arcpy.Sort_management(inPoints, os.path.join(temp_folder, 'densified_stream_sorted.shp'), [["mem_point1", "ASCENDING"]])
    
            streampoints = []
            cursorC = arcpy.da.SearchCursor(inPoints, "FID")
            for rowC in cursorC:
                streampoints.append(rowC[0])
            
            count = len(streampoints)
            
            i=0
            
            while i < count:
                arcpy.AddMessage("Processing SO " + str(streamorder) + 
                      " of " + str(streamorders[-1]) + " SS " +  str(streamseg + 1) + 
                      " of " + str(streamsegs[-1] + 1) + " Point " + str(i + 1) + " of " + str(count))
    
    
                #Accessing Hydropower Plant Database
                with arcpy.da.SearchCursor(inPoints, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM"], "FID = " + str(i)) as cursor1:
                    for row1 in cursor1:
                        Orig_ID_1 = row1[0]
                        Chainage_1 = row1[1]
                        x1 = row1[2]
                        y1 = row1[3]
                        fac1 = row1[4]*(cellsize)/1000000.0
                        z1 = row1[5]
                        segment = int(streamseg + 1)
                        
                        
                        #Selecting Hydropower Plant Point       
                        focus_point = arcpy.Select_analysis(inPoints, os.path.join(temp_folder, 'focus_point.shp'), "FID = " + str(i))
                        
                        #Creating Search Radius
                        out_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'out_buf.shp'), out_rad, method = "PLANAR")
                        in_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'in_buf.shp'), in_rad, method = "PLANAR")
                        search_rad = arcpy.Erase_analysis(out_buf, in_buf, os.path.join(temp_folder, 'search_rad.shp'))
                
                        #Creating Potential Diversion Points Feature Class
                        pot_div_points = arcpy.Clip_analysis(inPoints, search_rad, os.path.join(temp_folder, 'pot_div.shp'))
                        
                        #Counting Total Alternatives that satisfy Minimum Criteria
                        with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM", "SUBBASIN", "FLOW_INT"]) as cursor2:
                            total_alternatives = 0
                            for row2 in cursor2:
                                Orig_ID_2 = row2[0]
                                Chainage_2 = row2[1]
                                x2 = row2[2]
                                y2 = row2[3]
                                fac2 = row2[4]*(cellsize)/1000000.0
                                z2 = row2[5]

                                flow_int2 = int(row2[7]) 
                                subbasin2 = int(row2[6])
                                
            
                                head = z2-z1
                                distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                slope = head/distance2D
                                slope_percentage = slope*100
                                head_slope = head*slope_percentage
                                
                                if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC):
                                    total_alternatives = total_alternatives + 1
                        
                        #Accessing Diversion Site Database
                        if total_alternatives == 0:
                            arcpy.AddMessage('No Suitable Diversion Site')
                        elif total_alternatives == 1:
                            Chainage_3 = Chainage_2
                            Orig_ID_3 = Orig_ID_2
                            x3 = x2
                            y3 = y2
                            fac3 = fac2
                            z3 = z2
                            
                            flow_int_fin = flow_int2
                            subbasin_fin = subbasin2
                            
                            head_fin = head
                            distance3D_fin = distance3D
                            distance2D_fin = distance2D
                            slope_fin = slope
                            slope_percentage_fin = slope_percentage
                            head_slope_fin = head_slope
                            
                            FDC = FlowDurationCurve(river_folder, subbasin_fin, flow_int_fin, input_start_date, input_end_date, 'ALL')
                            Q10 = FDC[0]
                            Q20 = FDC[1]
                            Q30 = FDC[2]
                            Q40 = FDC[3]
                            Q50 = FDC[4]
                            Q60 = FDC[5]
                            Q70 = FDC[6]
                            Q80 = FDC[7]
                            Q90 = FDC[8]
                            Q95 = FDC[9]

                            Hydro10 = 0.00981*system_efficiency*head_fin*Q10
                            Hydro20 = 0.00981*system_efficiency*head_fin*Q20
                            Hydro30 = 0.00981*system_efficiency*head_fin*Q30
                            Hydro40 = 0.00981*system_efficiency*head_fin*Q40
                            Hydro50 = 0.00981*system_efficiency*head_fin*Q50
                            Hydro60 = 0.00981*system_efficiency*head_fin*Q60
                            Hydro70 = 0.00981*system_efficiency*head_fin*Q70
                            Hydro80 = 0.00981*system_efficiency*head_fin*Q80
                            Hydro90 = 0.00981*system_efficiency*head_fin*Q90
                            Hydro95 = 0.00981*system_efficiency*head_fin*Q95
                    
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " sq. km.")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " m")
                            arcpy.AddMessage("Q70 (cumecs) = " + str(Q70))
                            arcpy.AddMessage("Hydropower Potential (MW) at Q70 = " + str(Hydro70) + " MW")
                            
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                                'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                                'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                                'FAC_at_U/S':[fac3], 'Flow_Interval_U/S':[flow_int_fin], 'Subbasin_U/S':[subbasin_fin], 
                                                'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                                'Chain_US':[Chainage_3], 'Q10':[Q10], 'Q20':[Q20], 'Q30':[Q30], 'Q40':[Q40], 'Q50':[Q50],
                                                'Q60':[Q60], 'Q70':[Q70], 'Q80':[Q80], 'Q90':[Q90], 'Q95':[Q95],
                                                'Hydro10':[Hydro10], 'Hydro20':[Hydro20], 'Hydro30':[Hydro30], 'Hydro40':[Hydro40], 'Hydro50':[Hydro50],
                                                'Hydro60':[Hydro60], 'Hydro70':[Hydro70], 'Hydro80':[Hydro80], 'Hydro90':[Hydro90], 'Hydro95':[Hydro95]}, columns=columns)
                                            
                            results_df = pd.concat([results_df, df1])
                        
                        else:
                            alternatives = pd.DataFrame(index = range(0,total_alternatives), columns=["Chainage", "ID", "X", "Y", "Z", "Slope", "Distance2D", "Head", "Flow_Accumulation", "Waterway", "Flow_int", "Subbasin"])
                            alternatives_counter = 0
                            with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM", "SUBBASIN", "FLOW_INT"]) as cursor2:
                                for row2 in cursor2:
                                    Orig_ID_2 = row2[0]
                                    Chainage_2 = row2[1]
                                    x2 = row2[2]
                                    y2 = row2[3]
                                    fac2 = row2[4]*(cellsize)/1000000.0
                                    z2 = row2[5]
            
                                    head = z2-z1
                                    distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                    distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                    slope = head/distance2D
                                    slope_percentage = slope*100
                                    head_slope = head*slope_percentage
                                
                                    if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC):
                                        Orig_ID_3 = row2[0]
                                        Chainage_3 = row2[1]
                                        x3 = row2[2]
                                        y3 = row2[3]
                                        fac3 = row2[4]*(cellsize)/1000000.0
                                        z3 = row2[5]
    
                                        flow_int_fin = int(row2[7]) 
                                        subbasin_fin = int(row2[6])
                                
                                        head_fin = z3-z1
                                        distance3D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2)+((z1-z3)**2))
                                        distance2D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2))
                                        slope_fin = head_fin/distance2D_fin
                                        slope_percentage_fin = slope_fin*100
                                        head_slope_fin = head_fin*slope_percentage_fin
                                        
                                        
                                        
                                        #Creating the input for Promethee Module
                                        alternatives["Chainage"][alternatives_counter] = Chainage_3
                                        alternatives["ID"][alternatives_counter] = Orig_ID_3
                                        alternatives["X"][alternatives_counter] = x3
                                        alternatives["Y"][alternatives_counter] = y3
                                        alternatives["Z"][alternatives_counter] = z3
                                        alternatives["Slope"][alternatives_counter] = slope_percentage_fin
                                        alternatives["Distance2D"][alternatives_counter] = distance2D_fin
                                        alternatives["Head"][alternatives_counter] = head_fin
                                        alternatives["Flow_Accumulation"][alternatives_counter] = fac3
                                        alternatives["Waterway"][alternatives_counter] = distance3D_fin
                                        alternatives["Flow_int"][alternatives_counter] = flow_int_fin
                                        alternatives["Subbasin"][alternatives_counter] = subbasin_fin

                                        alternatives_counter = alternatives_counter + 1
                        
                            promethee_alternatives = alternatives[["Head","Flow_Accumulation","Waterway"]]
                            outranking_flows = promethee(promethee_alternatives, ["Head", "Flow_Accumulation"], ["Waterway"], 1, 0, weights_mode_0)
                            outranking_flows.index = alternatives.index
                            final_alternatives = pd.concat([alternatives,outranking_flows], axis=1)
                            final_alternatives.sort(["Net_Flow"], ascending=False, inplace=True)
                            final_alternatives.reset_index(inplace=True, drop=True)
                                
                            Chainage_3 = final_alternatives["Chainage"][0]
                            Orig_ID_3 = final_alternatives["ID"][0]
                            x3 = final_alternatives["X"][0]
                            y3 = final_alternatives["Y"][0]
                            fac3 = final_alternatives["Flow_Accumulation"][0]
                            z3 = final_alternatives["Z"][0]

                            flow_int_fin = final_alternatives["Flow_int"][0]
                            subbasin_fin = final_alternatives["Subbasin"][0]
                                
                            head_fin = final_alternatives["Head"][0]
                            distance3D_fin = final_alternatives["Waterway"][0]
                            distance2D_fin = final_alternatives["Distance2D"][0]
                            slope_fin = head_fin/distance2D_fin
                            slope_percentage_fin = slope_fin*100
                            head_slope_fin = head_fin*slope_percentage_fin
                            
                            j = int((Chainage_3/min_interval)-1)
                            
                            arcpy.AddMessage('subbasin_fin = '+ str(subbasin_fin))
                            arcpy.AddMessage('flow_int_fin = '+ str(flow_int_fin))
                            FDC = FlowDurationCurve(river_folder, subbasin_fin, flow_int_fin, input_start_date, input_end_date, 'ALL')
                            Q10 = FDC[0]
                            Q20 = FDC[1]
                            Q30 = FDC[2]
                            Q40 = FDC[3]
                            Q50 = FDC[4]
                            Q60 = FDC[5]
                            Q70 = FDC[6]
                            Q80 = FDC[7]
                            Q90 = FDC[8]
                            Q95 = FDC[9]

                            Hydro10 = 0.00981*system_efficiency*head_fin*Q10
                            Hydro20 = 0.00981*system_efficiency*head_fin*Q20
                            Hydro30 = 0.00981*system_efficiency*head_fin*Q30
                            Hydro40 = 0.00981*system_efficiency*head_fin*Q40
                            Hydro50 = 0.00981*system_efficiency*head_fin*Q50
                            Hydro60 = 0.00981*system_efficiency*head_fin*Q60
                            Hydro70 = 0.00981*system_efficiency*head_fin*Q70
                            Hydro80 = 0.00981*system_efficiency*head_fin*Q80
                            Hydro90 = 0.00981*system_efficiency*head_fin*Q90
                            Hydro95 = 0.00981*system_efficiency*head_fin*Q95
                    
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " sq. km.")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " m")
                            arcpy.AddMessage("Q70 (cumecs) = " + str(Q70))
                            arcpy.AddMessage("Hydropower Potential (MW) at Q70 = " + str(Hydro70) + " MW")
                            
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                                'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                                'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                                'FAC_at_U/S':[fac3], 'Flow_Interval_U/S':[flow_int_fin], 'Subbasin_U/S':[subbasin_fin], 
                                                'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                                'Chain_US':[Chainage_3], 'Q10':[Q10], 'Q20':[Q20], 'Q30':[Q30], 'Q40':[Q40], 'Q50':[Q50],
                                                'Q60':[Q60], 'Q70':[Q70], 'Q80':[Q80], 'Q90':[Q90], 'Q95':[Q95],
                                                'Hydro10':[Hydro10], 'Hydro20':[Hydro20], 'Hydro30':[Hydro30], 'Hydro40':[Hydro40], 'Hydro50':[Hydro50],
                                                'Hydro60':[Hydro60], 'Hydro70':[Hydro70], 'Hydro80':[Hydro80], 'Hydro90':[Hydro90], 'Hydro95':[Hydro95]}, columns=columns)
                            results_df = pd.concat([results_df, df1])

                        arcpy.Delete_management(focus_point)
                        arcpy.Delete_management(out_buf)
                        arcpy.Delete_management(in_buf)
                        arcpy.Delete_management(search_rad)
                        arcpy.Delete_management(pot_div_points)
                        
                        if arcpy.Exists(os.path.join(temp_folder, 'pot_div.shp'))==True:
                            arcpy.Delete_management(os.path.join(temp_folder, 'pot_div.shp'))
                
                i+=1
                arcpy.AddMessage('\n')
            arcpy.Delete_management(inPoints)
            arcpy.Delete_management(os.path.join(temp_folder, 'densified_stream.shp'))
            arcpy.Delete_management(streamseg_select)
            
        arcpy.Delete_management(order_select)
        arcpy.Delete_management(order)
  
    results_df.reset_index(drop = True, inplace = True)
    results_df = results_df[columns]
    results_df.index = np.arange(1, len(results_df)+1)
    
    alternative_schemes_counter = int(len(results_df.index))
#==============================================================================


#======================ROR Site Selection (MODE==2)============================
#Note: i, j and k are variables to keep track of the optimization procedures

elif mode == "Appraoch (ii): Using both Topographic and Hydrological Factors (For use with WEB-DHM)":
    arcpy.AddMessage("Initializing Run-of-River Site Selection Module (" + mode +")")

    
    k1 = float(minimum_hydropower_potential)
    k2 = float(maximum_hydropower_potential)
    flow_exceedance = int(flow_exceedance)
    
    #Creating Empty Results DataFrame
    columns = ['D/S_ID', 'U/S_ID', 'D/S_X', 'D/S_Y', 'D/S_Z', 'U/S_X', 'U/S_Y', 'U/S_Z', 
               'Headrace_Length_3D(m)', 'Headrace_Length_2D(m)', 'Head(m)', 'Slope(%)', 
                'FAC_at_U/S', 'Flow_Interval_U/S', 'Subbasin_U/S', 'Stream_Segment', 
                'Stream_Order', 'Chain_DS', 'Chain_US', 'Q10', 'Q20', 'Q30', 'Q40', 'Q50', 'Q60', 'Q70', 'Q80', 'Q90', 'Q95',
                'Hydro10', 'Hydro20', 'Hydro30', 'Hydro40', 'Hydro50', 'Hydro60', 'Hydro70', 'Hydro80', 'Hydro90', 'Hydro95']

    results_df = pd.DataFrame([], columns = columns)
    
    for streamorder in streamorders:
        order_select = arcpy.sa.ExtractByAttributes(strahler, "VALUE = " + str(streamorder))
        order_select.save(os.path.join(temp_folder, 'strahler_select.tif'))
        order = arcpy.sa.StreamToFeature(order_select, outFDR, os.path.join(temp_folder, 'strahler_select.shp'), "NO_SIMPLIFY")
        streamsegs = []
        cursorB = arcpy.da.SearchCursor(order, "FID")
        for rowB in cursorB:
            streamsegs.append(rowB[0])
        
        for streamseg in streamsegs:
            streamseg_select = arcpy.Select_analysis(order, os.path.join(temp_folder, 'streamseg_select.shp'), "FID = " + str(streamseg))
            inPoints = CreatePointsLines(streamseg_select, "INTERVAL BY DISTANCE", 
                                                           "END", "NO", '', float(min_interval), "NO", os.path.join(temp_folder, 'densified_stream.shp'))

            total_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            if Protected_Area_Switch == 'true':
                inPoints_WDPA = arcpy.Erase_analysis(inPoints, WDPA_Buffer, os.path.join(temp_folder, 'wdpa_stream.shp'))
                arcpy.Delete_management(inPoints)
                inPoints = inPoints_WDPA
                
            remaining_searching_points = int(arcpy.GetCount_management(inPoints).getOutput(0))
            excluded_searching_points = total_searching_points - remaining_searching_points
            WDPA_intersection_counter += excluded_searching_points
            
            arcpy.AddXY_management(inPoints)
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outFAC, "FAC"]])
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[outDEM, "DEM"]], 'BILINEAR')
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[inSubbasin, "SUBBASIN"]])
            arcpy.sa.ExtractMultiValuesToPoints(inPoints, [[inFlowInt, "FLOW_INT"]])
            inPoints = arcpy.Sort_management(inPoints, os.path.join(temp_folder, 'densified_stream_sorted.shp'), [["mem_point1", "ASCENDING"]])
    
            streampoints = []
            cursorC = arcpy.da.SearchCursor(inPoints, "FID")
            for rowC in cursorC:
                streampoints.append(rowC[0])
            
            count = len(streampoints)
            
            i=0
            
            while i < count:
                arcpy.AddMessage("Processing SO " + str(streamorder) + 
                      " of " + str(streamorders[-1]) + " SS " +  str(streamseg + 1) + 
                      " of " + str(streamsegs[-1] + 1) + " Point " + str(i + 1) + " of " + str(count))
    
    
                #Accessing Hydropower Plant Data Base
                with arcpy.da.SearchCursor(inPoints, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM"], "FID = " + str(i)) as cursor1:
                    for row1 in cursor1:
                        Orig_ID_1 = row1[0]
                        Chainage_1 = row1[1]
                        x1 = row1[2]
                        y1 = row1[3]
                        fac1 = row1[4]*(cellsize)/1000000.0
                        z1 = row1[5]
                        segment = int(streamseg + 1)
                        
                        
                        #Selecting Hydropower Plant Point       
                        focus_point = arcpy.Select_analysis(inPoints, os.path.join(temp_folder, 'focus_point.shp'), "FID = " + str(i))
                        
                        
                        #Creating Search Radius
                        out_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'out_buf.shp'), out_rad, method = "PLANAR")
                        in_buf = arcpy.Buffer_analysis(focus_point, os.path.join(temp_folder, 'in_buf.shp'), in_rad, method = "PLANAR")
                        search_rad = arcpy.Erase_analysis(out_buf, in_buf, os.path.join(temp_folder, 'search_rad.shp'))
                
                        #Creating Potential Diversion Points Feature Class
                        pot_div_points = arcpy.Clip_analysis(inPoints, search_rad, os.path.join(temp_folder, 'pot_div.shp'))
                        
                        #Counting Total Alternatives that satisfy Minimum Criteria
                        with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM", "SUBBASIN", "FLOW_INT"]) as cursor2:
                            total_alternatives = 0
                            for row2 in cursor2:
                                Orig_ID_2 = row2[0]
                                Chainage_2 = row2[1]
                                x2 = row2[2]
                                y2 = row2[3]
                                fac2 = row2[4]*(cellsize)/1000000.0
                                z2 = row2[5]

                                flow_int2 = int(row2[7]) 
                                subbasin2 = int(row2[6])
                                
                                if (flow_int2 > 0 and subbasin2 > 0):
                                    head = z2-z1
                                    distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                    distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                    slope = head/distance2D
                                    slope_percentage = slope*100
    #                                    head_slope = head*slope_percentage
                                    
                                    FDC = FlowDurationCurve(river_folder, subbasin2, flow_int2, input_start_date, input_end_date, 'ALL')
                                    Q10 = FDC[0]
                                    Q20 = FDC[1]
                                    Q30 = FDC[2]
                                    Q40 = FDC[3]
                                    Q50 = FDC[4]
                                    Q60 = FDC[5]
                                    Q70 = FDC[6]
                                    Q80 = FDC[7]
                                    Q90 = FDC[8]
                                    Q95 = FDC[9]
    
                                    Hydro10 = 0.00981*system_efficiency*head*Q10
                                    Hydro20 = 0.00981*system_efficiency*head*Q20
                                    Hydro30 = 0.00981*system_efficiency*head*Q30
                                    Hydro40 = 0.00981*system_efficiency*head*Q40
                                    Hydro50 = 0.00981*system_efficiency*head*Q50
                                    Hydro60 = 0.00981*system_efficiency*head*Q60
                                    Hydro70 = 0.00981*system_efficiency*head*Q70
                                    Hydro80 = 0.00981*system_efficiency*head*Q80
                                    Hydro90 = 0.00981*system_efficiency*head*Q90
                                    Hydro95 = 0.00981*system_efficiency*head*Q95
                                    
                                    if flow_exceedance == 10:
                                        hydropower = Hydro10
                                    elif flow_exceedance == 20:
                                        hydropower = Hydro20
                                    elif flow_exceedance == 30:
                                        hydropower_criteria = Hydro30
                                    elif flow_exceedance == 40:
                                        hydropower_criteria = Hydro40
                                    elif flow_exceedance == 50:
                                        hydropower_criteria = Hydro50
                                    elif flow_exceedance == 60:
                                        hydropower_criteria = Hydro60
                                    elif flow_exceedance == 70:
                                        hydropower_criteria = Hydro70
                                    elif flow_exceedance == 80:
                                        hydropower_criteria = Hydro80
                                    elif flow_exceedance == 90:
                                        hydropower_criteria = Hydro90
                                    elif flow_exceedance == 95:
                                        hydropower_criteria = Hydro95
                                    
                                    if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC and float(hydropower_criteria) >= k1 and float(hydropower_criteria) <= k2):
                                        total_alternatives = total_alternatives + 1
                        
                        #Accessing Diversion Site Database
                        if total_alternatives == 0:
                            arcpy.AddMessage('No Suitable Diversion Site')
                            
                        elif total_alternatives == 1:
                            Chainage_3 = Chainage_2
                            Orig_ID_3 = Orig_ID_2
                            x3 = x2
                            y3 = y2
                            fac3 = fac2
                            z3 = z2
                            
                            flow_int_fin = flow_int2
                            subbasin_fin = subbasin2
                            
                            head_fin = head
                            distance3D_fin = distance3D
                            distance2D_fin = distance2D
                            slope_fin = slope
                            slope_percentage_fin = slope_percentage
#                                head_slope_fin = head_slope
                            
                            FDC = FlowDurationCurve(river_folder, subbasin_fin, flow_int_fin, input_start_date, input_end_date, 'ALL')
                            Q10 = FDC[0]
                            Q20 = FDC[1]
                            Q30 = FDC[2]
                            Q40 = FDC[3]
                            Q50 = FDC[4]
                            Q60 = FDC[5]
                            Q70 = FDC[6]
                            Q80 = FDC[7]
                            Q90 = FDC[8]
                            Q95 = FDC[9]

                            Hydro10 = 0.00981*system_efficiency*head_fin*Q10
                            Hydro20 = 0.00981*system_efficiency*head_fin*Q20
                            Hydro30 = 0.00981*system_efficiency*head_fin*Q30
                            Hydro40 = 0.00981*system_efficiency*head_fin*Q40
                            Hydro50 = 0.00981*system_efficiency*head_fin*Q50
                            Hydro60 = 0.00981*system_efficiency*head_fin*Q60
                            Hydro70 = 0.00981*system_efficiency*head_fin*Q70
                            Hydro80 = 0.00981*system_efficiency*head_fin*Q80
                            Hydro90 = 0.00981*system_efficiency*head_fin*Q90
                            Hydro95 = 0.00981*system_efficiency*head_fin*Q95
                    
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " cells")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " sq. km.")
                            arcpy.AddMessage("Q70 (cumecs) = " + str(Q70))
                            arcpy.AddMessage("Hydropower Potential (MW) at Q70 = " + str(Hydro70) + " MW")
            
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                                'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                                'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                                'FAC_at_U/S':[fac3], 'Flow_Interval_U/S':[flow_int_fin], 'Subbasin_U/S':[subbasin_fin], 
                                                'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                                'Chain_US':[Chainage_3], 'Q10':[Q10], 'Q20':[Q20], 'Q30':[Q30], 'Q40':[Q40], 'Q50':[Q50],
                                                'Q60':[Q60], 'Q70':[Q70], 'Q80':[Q80], 'Q90':[Q90], 'Q95':[Q95],
                                                'Hydro10':[Hydro10], 'Hydro20':[Hydro20], 'Hydro30':[Hydro30], 'Hydro40':[Hydro40], 'Hydro50':[Hydro50],
                                                'Hydro60':[Hydro60], 'Hydro70':[Hydro70], 'Hydro80':[Hydro80], 'Hydro90':[Hydro90], 'Hydro95':[Hydro95]}, columns=columns)
                                            
                            results_df = pd.concat([results_df, df1])
                        
                        else:
                            alternatives = pd.DataFrame(index = range(0,total_alternatives), columns=["Chainage", "ID", "X", "Y", "Z", "Slope", "Distance2D", "Head", "Flow_Accumulation", "Waterway", "Flow_int", "Subbasin", "Hydropower"])
                            alternatives_counter = 0
                            with arcpy.da.SearchCursor(pot_div_points, ["FID", "mem_point1", "Point_X", "Point_Y", "FAC", "DEM", "SUBBASIN", "FLOW_INT"]) as cursor2:
                                for row2 in cursor2:
                                    Orig_ID_2 = row2[0]
                                    Chainage_2 = row2[1]
                                    x2 = row2[2]
                                    y2 = row2[3]
                                    fac2 = row2[4]*(cellsize)/1000000.0
                                    z2 = row2[5]

                                    flow_int2 = int(row2[7]) 
                                    subbasin2 = int(row2[6])
                                    
                                    if (flow_int2 > 0  and subbasin2 > 0):
                
                                        head = z2-z1
                                        distance3D = math.sqrt(((x1-x2)**2) +((y1-y2)**2)+((z1-z2)**2))
                                        distance2D = math.sqrt(((x1-x2)**2) +((y1-y2)**2))
                                        slope = head/distance2D
                                        slope_percentage = slope*100
    #                                        head_slope = head*slope_percentage
                                    
                                        FDC = FlowDurationCurve(river_folder, subbasin2, flow_int2, input_start_date, input_end_date, 'ALL')
                                        Q10 = FDC[0]
                                        Q20 = FDC[1]
                                        Q30 = FDC[2]
                                        Q40 = FDC[3]
                                        Q50 = FDC[4]
                                        Q60 = FDC[5]
                                        Q70 = FDC[6]
                                        Q80 = FDC[7]
                                        Q90 = FDC[8]
                                        Q95 = FDC[9]
        
                                        Hydro10 = 0.00981*system_efficiency*head*Q10
                                        Hydro20 = 0.00981*system_efficiency*head*Q20
                                        Hydro30 = 0.00981*system_efficiency*head*Q30
                                        Hydro40 = 0.00981*system_efficiency*head*Q40
                                        Hydro50 = 0.00981*system_efficiency*head*Q50
                                        Hydro60 = 0.00981*system_efficiency*head*Q60
                                        Hydro70 = 0.00981*system_efficiency*head*Q70
                                        Hydro80 = 0.00981*system_efficiency*head*Q80
                                        Hydro90 = 0.00981*system_efficiency*head*Q90
                                        Hydro95 = 0.00981*system_efficiency*head*Q95
                                        
                                        if flow_exceedance == 10:
                                            hydropower = Hydro10
                                        elif flow_exceedance == 20:
                                            hydropower = Hydro20
                                        elif flow_exceedance == 30:
                                            hydropower_criteria = Hydro30
                                        elif flow_exceedance == 40:
                                            hydropower_criteria = Hydro40
                                        elif flow_exceedance == 50:
                                            hydropower_criteria = Hydro50
                                        elif flow_exceedance == 60:
                                            hydropower_criteria = Hydro60
                                        elif flow_exceedance == 70:
                                            hydropower_criteria = Hydro70
                                        elif flow_exceedance == 80:
                                            hydropower_criteria = Hydro80
                                        elif flow_exceedance == 90:
                                            hydropower_criteria = Hydro90
                                        elif flow_exceedance == 95:
                                            hydropower_criteria = Hydro95
                                        
                                        if (head > min_head and slope_percentage > min_slope and fac2 > min_FAC and float(hydropower_criteria) >= k1 and float(hydropower_criteria) <= k2):
                                            Orig_ID_3 = row2[0]
                                            Chainage_3 = row2[1]
                                            x3 = row2[2]
                                            y3 = row2[3]
                                            fac3 = row2[4]*(cellsize)/1000000.0
                                            z3 = row2[5]
        
                                            flow_int_fin = int(row2[7]) 
                                            subbasin_fin = int(row2[6])
                                    
                                            head_fin = z3-z1
                                            distance3D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2)+((z1-z3)**2))
                                            distance2D_fin = math.sqrt(((x1-x3)**2) +((y1-y3)**2))
                                            slope_fin = head_fin/distance2D_fin
                                            slope_percentage_fin = slope_fin*100
    #                                            head_slope_fin = head_fin*slope_percentage_fin
                                            
                                            FDC = FlowDurationCurve(river_folder, subbasin_fin, flow_int_fin, input_start_date, input_end_date, 'ALL')
                                            Q10 = FDC[0]
                                            Q20 = FDC[1]
                                            Q30 = FDC[2]
                                            Q40 = FDC[3]
                                            Q50 = FDC[4]
                                            Q60 = FDC[5]
                                            Q70 = FDC[6]
                                            Q80 = FDC[7]
                                            Q90 = FDC[8]
                                            Q95 = FDC[9]
            
                                            Hydro10 = 0.00981*system_efficiency*head_fin*Q10
                                            Hydro20 = 0.00981*system_efficiency*head_fin*Q20
                                            Hydro30 = 0.00981*system_efficiency*head_fin*Q30
                                            Hydro40 = 0.00981*system_efficiency*head_fin*Q40
                                            Hydro50 = 0.00981*system_efficiency*head_fin*Q50
                                            Hydro60 = 0.00981*system_efficiency*head_fin*Q60
                                            Hydro70 = 0.00981*system_efficiency*head_fin*Q70
                                            Hydro80 = 0.00981*system_efficiency*head_fin*Q80
                                            Hydro90 = 0.00981*system_efficiency*head_fin*Q90
                                            Hydro95 = 0.00981*system_efficiency*head_fin*Q95
                                            
                                            if flow_exceedance == 10:
                                                hydropower = Hydro10
                                            elif flow_exceedance == 20:
                                                hydropower = Hydro20
                                            elif flow_exceedance == 30:
                                                hydropower_criteria = Hydro30
                                            elif flow_exceedance == 40:
                                                hydropower_criteria = Hydro40
                                            elif flow_exceedance == 50:
                                                hydropower_criteria = Hydro50
                                            elif flow_exceedance == 60:
                                                hydropower_criteria = Hydro60
                                            elif flow_exceedance == 70:
                                                hydropower_criteria = Hydro70
                                            elif flow_exceedance == 80:
                                                hydropower_criteria = Hydro80
                                            elif flow_exceedance == 90:
                                                hydropower_criteria = Hydro90
                                            elif flow_exceedance == 95:
                                                hydropower_criteria = Hydro95
                                              
                                            #Creating the input for Promethee Module
                                            alternatives["Chainage"][alternatives_counter] = Chainage_3
                                            alternatives["ID"][alternatives_counter] = Orig_ID_3
                                            alternatives["X"][alternatives_counter] = x3
                                            alternatives["Y"][alternatives_counter] = y3
                                            alternatives["Z"][alternatives_counter] = z3
                                            alternatives["Slope"][alternatives_counter] = slope_percentage_fin
                                            alternatives["Distance2D"][alternatives_counter] = distance2D_fin
                                            alternatives["Head"][alternatives_counter] = head_fin
                                            alternatives["Flow_Accumulation"][alternatives_counter] = fac3
                                            alternatives["Waterway"][alternatives_counter] = distance3D_fin
                                            alternatives["Flow_int"][alternatives_counter] = flow_int_fin
                                            alternatives["Subbasin"][alternatives_counter] = subbasin_fin
                                            alternatives["Hydropower"][alternatives_counter] = hydropower_criteria
    
                                            alternatives_counter = alternatives_counter + 1
                        
                            promethee_alternatives = alternatives[["Hydropower","Waterway"]]
                            outranking_flows = promethee(promethee_alternatives, ["Hydropower"], ["Waterway"], 1, 0, weights_mode_2)
                            outranking_flows.index = alternatives.index
                            final_alternatives = pd.concat([alternatives,outranking_flows], axis=1)
                            final_alternatives.sort(["Net_Flow"], ascending=False, inplace=True)
                            final_alternatives.reset_index(inplace=True, drop=True)
                                
                            Chainage_3 = final_alternatives["Chainage"][0]
                            Orig_ID_3 = final_alternatives["ID"][0]
                            x3 = final_alternatives["X"][0]
                            y3 = final_alternatives["Y"][0]
                            fac3 = final_alternatives["Flow_Accumulation"][0]
                            z3 = final_alternatives["Z"][0]

                            flow_int_fin = final_alternatives["Flow_int"][0]
                            subbasin_fin = final_alternatives["Subbasin"][0]
                                
                            head_fin = final_alternatives["Head"][0]
                            distance3D_fin = final_alternatives["Waterway"][0]
                            distance2D_fin = final_alternatives["Distance2D"][0]
                            slope_fin = head_fin/distance2D_fin
                            slope_percentage_fin = slope_fin*100
#                                head_slope_fin = head_fin*slope_percentage_fin
                            
                            j = int((Chainage_3/min_interval)-1)
                                
                            FDC = FlowDurationCurve(river_folder, subbasin_fin, flow_int_fin, input_start_date, input_end_date, 'ALL')
                            Q10 = FDC[0]
                            Q20 = FDC[1]
                            Q30 = FDC[2]
                            Q40 = FDC[3]
                            Q50 = FDC[4]
                            Q60 = FDC[5]
                            Q70 = FDC[6]
                            Q80 = FDC[7]
                            Q90 = FDC[8]
                            Q95 = FDC[9]

                            Hydro10 = 0.00981*system_efficiency*head_fin*Q10
                            Hydro20 = 0.00981*system_efficiency*head_fin*Q20
                            Hydro30 = 0.00981*system_efficiency*head_fin*Q30
                            Hydro40 = 0.00981*system_efficiency*head_fin*Q40
                            Hydro50 = 0.00981*system_efficiency*head_fin*Q50
                            Hydro60 = 0.00981*system_efficiency*head_fin*Q60
                            Hydro70 = 0.00981*system_efficiency*head_fin*Q70
                            Hydro80 = 0.00981*system_efficiency*head_fin*Q80
                            Hydro90 = 0.00981*system_efficiency*head_fin*Q90
                            Hydro95 = 0.00981*system_efficiency*head_fin*Q95
                            
                            arcpy.AddMessage("Head = " + str(head_fin) + " m")
                            arcpy.AddMessage("Flow Accumulation = " + str(fac3) + " sq. km.")
                            arcpy.AddMessage("Length of Waterway = " + str(distance3D_fin) + " m")
                            arcpy.AddMessage("Q70 (cumecs) = " + str(Q70))
                            arcpy.AddMessage("Hydropower Potential (MW) at Q70 = " + str(Hydro70) + " MW")
            
                            df1 = pd.DataFrame({'D/S_ID':[Orig_ID_1], 'U/S_ID':[j], 'D/S_X':[x1], 'D/S_Y':[y1], 'D/S_Z':[z1], 
                                                'U/S_X':[x3], 'U/S_Y':[y3], 'U/S_Z':[z3],  'Headrace_Length_3D(m)':[distance3D_fin], 
                                                'Headrace_Length_2D(m)':[distance2D_fin], 'Head(m)':[head_fin], 'Slope(%)':[slope_percentage_fin], 
                                                'FAC_at_U/S':[fac3], 'Flow_Interval_U/S':[flow_int_fin], 'Subbasin_U/S':[subbasin_fin], 
                                                'Stream_Segment':[segment], 'Stream_Order':[streamorder], 'Chain_DS':[Chainage_1], 
                                                'Chain_US':[Chainage_3], 'Q10':[Q10], 'Q20':[Q20], 'Q30':[Q30], 'Q40':[Q40], 'Q50':[Q50],
                                                'Q60':[Q60], 'Q70':[Q70], 'Q80':[Q80], 'Q90':[Q90], 'Q95':[Q95],
                                                'Hydro10':[Hydro10], 'Hydro20':[Hydro20], 'Hydro30':[Hydro30], 'Hydro40':[Hydro40], 'Hydro50':[Hydro50],
                                                'Hydro60':[Hydro60], 'Hydro70':[Hydro70], 'Hydro80':[Hydro80], 'Hydro90':[Hydro90], 'Hydro95':[Hydro95]}, columns=columns)
                            results_df = pd.concat([results_df, df1])

                        arcpy.Delete_management(focus_point)
                        arcpy.Delete_management(out_buf)
                        arcpy.Delete_management(in_buf)
                        arcpy.Delete_management(search_rad)
                        arcpy.Delete_management(pot_div_points)
                        
                        if arcpy.Exists(os.path.join(temp_folder, 'pot_div.shp'))==True:
                            arcpy.Delete_management(os.path.join(temp_folder, 'pot_div.shp'))

                
                i+=1
                arcpy.AddMessage('\n')
            arcpy.Delete_management(inPoints)
            arcpy.Delete_management(os.path.join(temp_folder, 'densified_stream.shp'))
            arcpy.Delete_management(streamseg_select)
            
        arcpy.Delete_management(order_select)
        arcpy.Delete_management(order)
  
    results_df.reset_index(drop = True, inplace = True)
    results_df = results_df[columns]
    results_df.index = np.arange(1, len(results_df)+1)
    
    alternative_schemes_counter = int(len(results_df.index))
#==============================================================================


#======================Removing Interfering Schemes============================
results_df_final = pd.DataFrame([], columns = columns)
results_df.reset_index(drop = True, inplace = True)

if mode == "SITE_SELECTION_USING_POWER_CRITERIA":
    if flow_exceedance == 10:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro10'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 20:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro20'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 30:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro30'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 40:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro40'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 50:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro50'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 60:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro60'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 70:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro70'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 80:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro80'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 90:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro90'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)
    elif flow_exceedance == 95:
        outranking_flows = promethee(results_df[['Hydro10', 'Headrace_Length_3D(m)']],['Hydro95'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_2)

else:
    outranking_flows = promethee(results_df[['Head(m)', 'FAC_at_U/S', 'Headrace_Length_3D(m)']],['Head(m)', 'FAC_at_U/S'], ['Headrace_Length_3D(m)'], 1, 0, weights_mode_0)        

outranking_flows.index = results_df.index
results_df = pd.concat([results_df,outranking_flows], axis=1)
results_df.index = np.arange(1, len(results_df)+1)
results_df = results_df.sort(['Net_Flow'], ascending = [False])

results_df['Interference']=False
results_df['Chain_DS'] = results_df['Chain_DS'].astype(int)
results_df['Chain_US'] = results_df['Chain_US'].astype(int)
results_df['Stream_Segment'] = results_df['Stream_Segment'].astype(int)
results_df['Stream_Order'] = results_df['Stream_Order'].astype(int)

for streamorder in streamorders:
    results_df_streamorder = results_df[results_df['Stream_Order']==streamorder]
    streamsegs = []
    streamsegs = map(int, np.array(results_df_streamorder['Stream_Segment'].unique()))
    for streamseg in streamsegs:
        results_df_streamseg = results_df_streamorder[results_df_streamorder['Stream_Segment']==streamseg]
        results_df_streamseg.reset_index(drop = True, inplace = True)
        for i in range(0, len(results_df_streamseg.index)):
            
            if results_df_streamseg.loc[i, 'Interference'] == False:
                final_DS_ID = results_df_streamseg.loc[i, 'Chain_DS']
                final_US_ID = results_df_streamseg.loc[i, 'Chain_US']
            
                for j in range(i+1, len(results_df_streamseg.index)):

                    if results_df_streamseg.loc[j, 'Interference'] == False:
                        check_DS_ID = results_df_streamseg.loc[j, 'Chain_DS']
                        check_US_ID = results_df_streamseg.loc[j, 'Chain_US']

                        if check_DS_ID in range(final_DS_ID, final_US_ID + 1):
                            results_df_streamseg.loc[j,'Interference'] = True

                        elif check_US_ID in range(final_DS_ID, final_US_ID + 1):
                            results_df_streamseg.loc[j, 'Interference'] = True

                        elif final_DS_ID in range(check_DS_ID, check_US_ID + 1) and final_US_ID in range(check_DS_ID, check_US_ID + 1):
                            results_df_streamseg.loc[j, 'Interference'] = True

        results_df_streamseg = results_df_streamseg[results_df_streamseg['Interference']==False]
        results_df_final = pd.concat([results_df_streamseg, results_df_final])

#results_df_final = results_df #Remove comment to see intersections #Temporary

results_df_final = results_df_final[columns+['Net_Flow']]
results_df_final.reset_index(drop = True, inplace = True)
arcpy.AddMessage(results_df_final)
results_df_final.to_csv(txt_output, sep = " ", mode = 'w')

final_schemes_counter = int(len(results_df_final.index))
#==============================================================================


#===================Generating ShapeFile of Points & Schemes===================
Hydropower_Plants_Event = arcpy.MakeXYEventLayer_management(txt_output, 'D/S_X', 'D/S_Y', 'Hydropower_Plants', spatial_ref)
Hydropower_Plants_Shape = arcpy.CopyFeatures_management(Hydropower_Plants_Event, hydropower_plants_output)

Dam_Site_Event = arcpy.MakeXYEventLayer_management(txt_output, 'U/S_X', 'U/S_Y', 'Dam_Sites', spatial_ref)
Dam_Site_Shape = arcpy.CopyFeatures_management(Dam_Site_Event, dam_site_output)

Schemes = arcpy.XYToLine_management(txt_output, schemes_output, 'D/S_X', 'D/S_Y', 'U/S_X', 'U/S_Y', '#', '#', spatial_ref)
Schemes = arcpy.DeleteField_management(Schemes, ['D/S_Y', 'U/S_X', 'U/S_Y'])

Join_Fields = (arcpy.ListFields(Hydropower_Plants_Shape))[3:]                                                                   #May need to change
Join_Fields1 = []
for Join_Field in Join_Fields:
    Join_Fields1.append(Join_Field.name)

Schemes = arcpy.JoinField_management(Schemes, 'FID', Hydropower_Plants_Shape, 'FID', Join_Fields1)
Schemes = arcpy.DeleteField_management(Schemes, ['D/S_X'])

arcpy.DeleteField_management(Hydropower_Plants_Shape, ['Field1'])
arcpy.DeleteField_management(Dam_Site_Shape, ['Field1'])
#==============================================================================


#========================Deleting Temporary Files==============================
temp_rasters = arcpy.ListRasters()
for temp_raster in temp_rasters:
    arcpy.Delete_management(temp_raster)
#==============================================================================

#=======================Adding Layers to Current DataFrame=====================
AddFeatureLayer(os.path.join(preprocessing_folder, 'watershed.shp'))
AddFeatureLayer(os.path.join(preprocessing_folder, 'stream.shp'))
AddFeatureLayer(schemes_output)
AddFeatureLayer(hydropower_plants_output)
AddFeatureLayer(dam_site_output)    
#==========================Checking in extensions==============================
arcpy.AddMessage("Checking in extensions...")
if arcpy.CheckExtension("spatial") == "Available":
    arcpy.CheckInExtension("sptial")
    
if arcpy.CheckExtension("3D") == "Available":
    arcpy.CheckInExtension("3D")
#==============================================================================


#===================Printing Counters for Rejected Sites=======================
arcpy.AddMessage("WDPA Intersections = " + str(WDPA_intersection_counter))
arcpy.AddMessage("Rejected by Interference Tests = " + str(alternative_schemes_counter-final_schemes_counter))
#==============================================================================

arcpy.AddMessage("Successfully Finished!")

end_time = time.time()
program_runtime = end_time - start_time
program_runtime = 'Program Execution Time = ' + str(program_runtime) + ' seconds'
arcpy.AddMessage(program_runtime)