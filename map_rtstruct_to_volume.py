'''
Created on Jan 30, 2020

@author: davisr28
'''

import pydicom
import os
import pandas as pd
import numpy as np

# mapToOriginalVolume takes a image volume & RTStruct output by slicer and maps that RTStruct back to 
#     the original volume that slicer originally read in.
#
# Input:
#     original_volume_loc: path to folder containing the original volume that we want the RTStruct to reference
#     slicer_volume_loc: path to folder containing the image volume ant RTStruct that was output by slicer
#     rtss_save_loc: path to folder where the RTStruct (which references the original volume) should be stored.
#
# Caveats:
#    1) This function throws an error when there is not a 1-to-1 mapping of z positions between the original and Slicer
#    output volumes. This situation (no 1-to-1 mapping) occurs when Slicer loads an irregular geometry, such as
#    jumps in z position.
#
#    2) In the situation where there is irregular geometry but still a 1-to-1 mapping of z positions between volumes
#    this code maps the volume SOPInstanceUIDs by sorting z positions of both volumes, then pairing by the sorted order.
#    
#    3) Take caution in general when changing referenced UIDs, and especially when using this function on volumes
#    with irregular geometries.

def mapRTStructToVolume(original_volume_loc, slicer_volume_loc, rtss_save_loc):

    # read the z position to SOPInstanceUID mapping into dataframes
    original_uids_df,original_vol_attr=getUIDAndZPos(original_volume_loc)
    original_uids_df=original_uids_df.rename(mapper={"SOPInstanceUID":"SOPInstanceUID_Original"},axis=1)
    new_uids_df=getUIDAndZPos(slicer_volume_loc)[0].rename(mapper={"SOPInstanceUID":"SOPInstanceUID_New"},axis=1)
    
    # double check that the mapping based on Z position will be 1-to-1 as expected
    original_z=set(original_uids_df["ImageZPosMm"])
    new_z=set(new_uids_df["ImageZPosMm"])
    union=new_z.union(original_z)
    
    # generate the mapping from slicer_output_volume SOPInstanceUIDs to original volume SOPInstanceUIDs
    if len(union) == len(original_z): # if this is false, then the new and original volumes are not matched
    
        # Find the SOPInstanceUID mapping and drop the Z position
        new_uids_df=new_uids_df.rename(axis=1,mapper={"ImageZPosMm":"ImageZPosMm_New"})
        original_uids_df=original_uids_df.rename(axis=1,mapper={"ImageZPosMm":"ImageZPosMm_Original"})
    
        new_to_original_UID_mapping_df=new_uids_df[["ImageZPosMm_New","SOPInstanceUID_New"]].merge(original_uids_df[["ImageZPosMm_Original","SOPInstanceUID_Original"]],left_on="ImageZPosMm_New",right_on="ImageZPosMm_Original",how="left")
    
    # often this happens when there is irregular slice geometry
    elif len(original_z) == len(new_z):
        new_to_original_UID_mapping_df = mapBySliceSorting(new_uids_df,original_uids_df)
        if len(new_to_original_UID_mapping_df)==0:
            return False
    else:
        return False
    
    # get the rtss location and make sure there is only one RTSS in that directory as expected:
    new_vol_files=os.listdir(slicer_volume_loc)
    rtss_files=[file for file in new_vol_files if "rtss" in file]
    assert len(rtss_files) == 1
    rtss_file=rtss_files[0]
    
    # read the rtss
    rtss=pydicom.dcmread(os.path.join(slicer_volume_loc,rtss_file))
    
    # Set top level attributes of the RTStruct to output
    rtss[0x0020,0x000d].value=original_vol_attr["StudyInstanceUID"]
    
    # Set attributes within ReferencedFrameOfReferenceSequence
    for i_RefFrameOfRef in range(len(list(rtss[0x3006,0x0010]))):
        for i_RTRefStudy in range(len(list(rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012]))): # (3006, 0012)  RT Referenced Study Sequence
            
            # make sure for this simple code that RTStruct only references one series
            if len(list(rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014])) != 1:
                raise RuntimeError("This function does not handle RTStructs referencing more than one series.")
            
            for i_RTRefSeries in range(len(list(rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014]))):
                for i_ContourImage in range(len(list(rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014][i_RTRefSeries][0x3006,0x0016]))):
    #             for RTReferencedSeries in RTReferencedStudy[0x3006,0x0014]:
                    new_vol_SOPInstanceUID=rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014][i_RTRefSeries][0x3006,0x0016][i_ContourImage][0x0008,0x1155].value
                    original_vol_SOPInstanceUID=new_to_original_UID_mapping_df[new_to_original_UID_mapping_df["SOPInstanceUID_New"]==new_vol_SOPInstanceUID]["SOPInstanceUID_Original"].iloc[0]
                    rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014][i_RTRefSeries][0x3006,0x0016][i_ContourImage][0x0008,0x1155].value=original_vol_SOPInstanceUID
                    
                    # ReferenceSOPClassUID                
                    rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014][i_RTRefSeries][0x3006,0x0016][i_ContourImage][0x0008,0x1150].value=original_vol_attr["SOPClassUID"]
                    
                # Referenced SeriesInstanceUID
                original_series=original_uids_df[original_uids_df["SOPInstanceUID_Original"] == original_vol_SOPInstanceUID]["SeriesInstanceUID"].iloc[0]
                rtss[0x3006,0x0010][i_RefFrameOfRef][0x3006,0x0012][i_RTRefStudy][0x3006,0x0014][i_RTRefSeries][0x0020,0x000e].value=original_series
  
    # Set Attributes within ROIContourSequence
    for i_ROIContourSequence in range(len(list(rtss[0x3006,0x0039]))):
        for i_ContourSequence in range(len(list(rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040]))):
            
            # check if slicer altered the z values for this contour (as happens when there is an irregular geometry
            contour_seq_z_val=round(float(rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0050].value[2]),1)
            df_row=new_to_original_UID_mapping_df[new_to_original_UID_mapping_df["ImageZPosMm_New"]==contour_seq_z_val]

            # if slicer has altered the z values, then change the z values in the RTStruct to match those of the original image volume
            if not (df_row["ImageZPosMm_Original"].iloc[0] == df_row["ImageZPosMm_New"].iloc[0]):
                original_z_pos=df_row["ImageZPosMm_Original"].iloc[0]
                old_contour=[float(el) for el in rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0050].value]
                new_contour_float=[old_contour[ii] if (ii-2) % 3 != 0 else original_z_pos for ii in range(len(old_contour))]
                new_contour_str=[str(el) for el in new_contour_float]
                rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0050].value = new_contour_str
                raise UserWarning("Have not verified that sorting-based SOPInstanceUID mapping works. Compare the contour locations on the slicer and original volumes to make sure.")
                
            for i_ContourImageSequence in range(len(list(rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0016]))):

                # Referenced SOPInstanceUID
                new_vol_ReferencedSOPInstanceUID=rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0016][i_ContourImageSequence][0x0008,0x1155].value
                original_vol_ReferencedSOPInstanceUID=new_to_original_UID_mapping_df[new_to_original_UID_mapping_df["SOPInstanceUID_New"]==new_vol_ReferencedSOPInstanceUID]["SOPInstanceUID_Original"].iloc[0]
                rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0016][i_ContourImageSequence][0x0008,0x1155].value=original_vol_ReferencedSOPInstanceUID
                
                # ReferenceSOPClassUID
                rtss[0x3006,0x0039][i_ROIContourSequence][0x3006,0x0040][i_ContourSequence][0x3006,0x0016][i_ContourImageSequence][0x0008,0x1150].value=original_vol_attr["SOPClassUID"]

    pydicom.filewriter.dcmwrite(os.path.join(rtss_save_loc,rtss_file),rtss)
            
    return True

def mapBySliceSorting(new_uids_df,original_uids_df):
    
    # order the z position of the slices in both volumes
    slicer_vol_z_list=list(new_uids_df["ImageZPosMm"])
    slicer_vol_z_list_sort=slicer_vol_z_list.copy()
    slicer_vol_z_list_sort.sort(reverse=slicer_vol_z_list[0]>slicer_vol_z_list[-1])
    
    orig_vol_z_list=list(original_uids_df["ImageZPosMm"])
    orig_vol_z_list_sort=orig_vol_z_list.copy()
    orig_vol_z_list_sort.sort(reverse=orig_vol_z_list[0]>orig_vol_z_list[-1])
    
    # check if the slices are in order based on z position:
    if all(np.array(orig_vol_z_list) == np.array(orig_vol_z_list_sort)) and all(np.array(slicer_vol_z_list) == np.array(slicer_vol_z_list_sort)):
        
        original_uids_df=original_uids_df.rename(axis=1,mapper={"ImageZPosMm":"ImageZPosMm_Original"})
        new_uids_df=new_uids_df.rename(axis=1,mapper={"ImageZPosMm":"ImageZPosMm_New"})
        
        # if the slice orders are reversed, reverse the new volume order
        if np.sign(original_uids_df["ImageZPosMm_Original"].iloc[0] - original_uids_df["ImageZPosMm_Original"].iloc[len(original_uids_df)-1]) != np.sign(new_uids_df["ImageZPosMm_New"].iloc[0] - new_uids_df["ImageZPosMm_New"].iloc[len(new_uids_df)-1]):
            new_uids_df=new_uids_df.sort_index(axis=0, ascending=False)
            new_uids_df.index=range(len(new_uids_df))
            
        mapping_df=new_uids_df[["SOPInstanceUID_New","ImageZPosMm_New"]].join(original_uids_df[["SOPInstanceUID_Original","ImageZPosMm_Original"]])
        return mapping_df
    else:
        return pd.DataFrame([])
    print("t")


def getUIDAndZPos(directory):
    
    # get mapping of SOPInstanceUIDs to z position
    uids_df=pd.DataFrame()
    image_files=os.listdir(directory)
    image_files=[im_f for im_f in image_files if (".dcm" in im_f and "rtss" not in im_f)]
#     count=0
    for image_file in image_files:
        try:
            dcm=pydicom.dcmread(directory+"\\"+image_file)
        except:
            print("t")
        to_append=pd.DataFrame({"SOPInstanceUID":[dcm[0x0008,0x0018].value],"ImageZPosMm":round(float(dcm[0x0020,0x0032].value[2]),1),"SOPClassUID":[dcm[0x0008,0x0016].value],"SeriesInstanceUID":[dcm[0x0020,0x000e].value],"StudyInstanceUID":[dcm[0x0020,0x000d].value],"FileName":[image_file]})
        uids_df=uids_df.append(to_append,ignore_index=True,sort=False)
    
    attr={}
    
    # get SOPClassUID_
    SOPClassUID_set=set(uids_df["SOPClassUID"])
    assert len(SOPClassUID_set) == 1
    attr["SOPClassUID"]=SOPClassUID_set.pop()
    
    #get SeriesInstanceUID
    SeriesInstanceUID_set=set(uids_df["SeriesInstanceUID"])
    if len(SeriesInstanceUID_set) > 1:
        raise RuntimeError("!"*100 + "\nWARNING!!!!!! More than one series in dir: " + directory + "\n" + "!"*100)
    attr["SeriesInstanceUID"]=SeriesInstanceUID_set.pop()    

    #get StudyInstanceUID
    StudyInstanceUID_set=set(uids_df["StudyInstanceUID"])
    if len(StudyInstanceUID_set) > 1:
        raise RuntimeError("!"*100 + "\nWARNING!!!!!! More than one study in dir: " + directory + "\n" + "!"*100)
    attr["StudyInstanceUID"]=StudyInstanceUID_set.pop()    
    
    return uids_df,attr
