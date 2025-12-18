import pandas as pd
import numpy as np
import datetime, os
import yaml
import sdmc_tools.process as sdmc
import sdmc_tools.constants as constants
import sdmc_tools.access_ldms as access_ldms


def main():
    input_data_path = "/trials/vaccine/p137/s001/qdata/LabData/IHC_Pass-through/Cell Count Analysis from Youyi  with cell densities per um2 with CD3 test v3 with SHARP columns.xlsx"
    data = pd.read_excel(input_data_path)

    rename = {
        'SpecID':'guspec',
        'spectype':'spectype_submitted',
        'PTID':'ptid_submitted',
        'Visit':'visit_submitted',
        'Collection Date':'drawdt_submitted',
        'protocol':'protocol_submitted',    
        'network': 'network',
        'specrole': 'specrole',
        
        'lab ID': 'lab_id',
        'assay': 'assay',
        'assaysub': 'assaysub',
        'method': 'method',
        'runbatchno': 'runbatchno',
        'rundate': 'rundate',
        'imagebatchno': 'imagebatchno',
        'imagedate': 'imagedate',
        'analysisbatchno': 'analysisbatchno',
        'analysisfilename': 'analysisfilename',
        'reliable': 'reliable',
        'studyimageid': 'studyimageid',
        'control_sample_name': 'control_sample_name',
        'runnum': 'runnum',
        'comments': 'comments',
        
        'Lamina_Propria_CV2  Annotations files': 'lamina_propria_cv2_annotations_files',
        'Epithelium_CV2 Annotations files': 'epithelium_cv2_annotations_files',
        'Epi  Manual Revisions Comments': 'epi_manual_revisions_comments',
        'Epithelium_CV2_and_Manual Annotationss File': 'epithelium_cv2_and_manual_annotations_file',
        'Lamina_Propria_CV2 Manual Revisions Comments': 'lamina_propria_cv2_manual_revisions_comments',
        'Lamina_Propria_CV2 +Manual Revisions File': 'lamina_propria_cv2_+manual_revisions_file',
        'Nuclear Detection: Minimal Nuclear Intensity': 'nuclear_detection_minimal_nuclear_intensity',
        
        'CD3 Threshold': 'cd3_threshold',
        'CCR5 Threshhold': 'ccr5_threshhold',
        'CD4 Threshold': 'cd4_threshold',
        'T-Cell Algorithm': 't-cell_algorithm',
        'T-Cell algorithm Comments': 't-cell_algorithm_comments',
        'Image Tag of Cell Count Analysis from Youyi': 'image_tag_of_cell_count_analysis_youyi',
        'tissue_area_type': 'tissue_area_type',
        'area_size_px': 'area_size_px',
        
        'CD3+': 'result_cd3+_cell_count',
        'CD3+CD4+': 'result_cd3+cd4+_cell_count',
        'CD3+CCR5+CD4+': 'result_cd3+ccr5+cd4+_cell_count',
        
        'Python Density CD3+/um2': 'result_cd3+_density_python',
        'Python Density CD3+ CD4+/um2': 'result_cd3+cd4+_density_python',
        'Python Density CD3+ CD4+ CCR5+/um2': 'result_cd3+cd4+ccr5+_density_python',
        
        'HALO Image Tag of T Cell Density analysis': 'image_tag_of_t_cell_density_analysis_halo',
        'HALO Area Analyzed (μm²)': 'halo_area_analyzed',
        
        'HALO CD3 Density cells/um2': 'result_cd3+_density_halo',
        'HALO CD4+ Density cells/um2': 'result_cd4+_density_halo',
        'HALO CD3+ CD4+ Density cells/um2': 'result_cd3+cd4+_density_halo',
        'HALO CCR5+ Density cells/um2': 'result_ccr5+_density_halo',
        'HALO CCR5+ CD3+ Density cells/um2': 'result_ccr5+cd3+_density_halo',
        'HALO CD3+ CD4+ CCR5+ Density cells/um2': 'result_cd3+cd4+ccr5+_density_halo',
    }

    data = data.rename(columns=rename)

    md = {
        'halo_area_analyzed_units':'micrometers squared',
        'result_units_density':'cells per micrometer squared',
        'lab_software':'HALO',
        'upload_lab_id':'FH',
        'assay_precision':'Quantitative',
        'assay_lab_name':'McElrath Lab (Fred Hutch)',
        'assay_type':'IHC (Immunohistochemistry)',
        'assay_subtype':'FFPE (Formalin-Fixed, Paraffin-Embedded)',
        'instrument':'PhenoImager',
    }

    ldms = access_ldms.pull_one_protocol('hvtn', 137)

    # set(data.guspec).difference(ldms.guspec)
    # data.loc[data.guspec=="J9T00CJC-01",['guspec','specrole']]

    ldms = ldms.loc[ldms.guspec.isin(data.guspec)]

    check_cols = ['protocol','spectype','specrole','PTID','Visit','Collection Date', 'SpecID']
    # tissue_area_type	area_size_px	

    outputs = sdmc.standard_processing(
            input_data=data, #.loc[data.guspec!="J9T00CJC-01"],
            input_data_path=input_data_path,
            guspec_col='guspec',
            network='hvtn',
            metadata_dict=md,
            ldms=ldms
    )


    assert (outputs.loc[outputs.specrole=='Sample'].spectype_submitted!=outputs.loc[outputs.specrole=='Sample'].spectype).sum() == 0
    assert (outputs.loc[outputs.specrole=='Sample'].ptid_submitted.astype(int)!=outputs.loc[outputs.specrole=='Sample'].ptid.astype(int)).sum() == 0
    assert (outputs.loc[outputs.specrole=='Sample'].visit_submitted.astype(int)!=outputs.loc[outputs.specrole=='Sample'].visitno.astype(int)).sum() == 0
    assert (outputs.loc[outputs.specrole=='Sample'].drawdt_submitted!=outputs.loc[outputs.specrole=='Sample'].drawdt).sum() == 0
    assert (outputs.loc[outputs.specrole=='Sample'].protocol_submitted.astype(float)!=outputs.loc[outputs.specrole=='Sample'].protocol.astype(float)).sum() == 0
    assert (outputs.image_tag_of_t_cell_density_analysis_halo!=outputs.image_tag_of_cell_count_analysis_youyi).sum() == 0

    outputs = outputs.drop(columns=[
        'spectype_submitted',
        'ptid_submitted',
        'visit_submitted',
        'drawdt_submitted',
        'protocol_submitted',
        'analysisfilename',
        'assay',
        'assaysub',
        'lab_id',
        'rundate', # Maria wasnt able to add the run dates
    ])

    outputs.control_sample_name.unique()

    reorder = [
        'network',
        'protocol',
        'specrole',
        'guspec',
        'control_sample_name',
        'ptid',
        'visitno',
        'drawdt',
        'spectype',
        'spec_primary',
        'spec_additive',
        'spec_derivative',
        'upload_lab_id',
        'assay_lab_name',
        'assay_type',
        'assay_subtype',
        'assay_precision',
        'instrument',
        'lab_software',
        
        'method',
        
        'cd3_threshold',
        'ccr5_threshhold',
        'cd4_threshold',
        't-cell_algorithm',
        't-cell_algorithm_comments',
        'image_tag_of_cell_count_analysis_youyi',
        'image_tag_of_t_cell_density_analysis_halo',
        'tissue_area_type',
        'area_size_px',
        
        'result_cd3+_cell_count',
        'result_cd3+cd4+_cell_count',
        'result_cd3+ccr5+cd4+_cell_count',
        
        'result_cd3+_density_python',
        'result_cd3+cd4+_density_python',
        'result_cd3+cd4+ccr5+_density_python',

        'halo_area_analyzed',
        'halo_area_analyzed_units',

        'result_cd3+_density_halo',
        'result_cd4+_density_halo',
        'result_cd3+cd4+_density_halo',
        'result_ccr5+_density_halo',
        'result_ccr5+cd3+_density_halo',
        'result_cd3+cd4+ccr5+_density_halo',
        'result_units_density',
        
        'runbatchno',
        'runnum',
        'imagebatchno',
        'imagedate',
        'analysisbatchno',
        'studyimageid',
        'reliable',
        'comments',
        
        'lamina_propria_cv2_annotations_files',
        'epithelium_cv2_annotations_files',
        'epi_manual_revisions_comments',
        'epithelium_cv2_and_manual_annotations_file',
        'lamina_propria_cv2_manual_revisions_comments',
        'lamina_propria_cv2_+manual_revisions_file',
        'nuclear_detection_minimal_nuclear_intensity',

        'sdmc_processing_datetime',
        'sdmc_data_receipt_datetime',
        'input_file_name',
    ]

    set(reorder).symmetric_difference(outputs.columns)
    outputs = outputs[reorder]

    outputs.loc[outputs.specrole!='Sample', ['specrole','guspec','control_sample_name']]
    outputs.loc[outputs.specrole!='Sample', 'guspec'] = np.nan

    outputs.control_sample_name.iloc[0]

    outputs.visitno = outputs.visitno.astype(float)
    outputs.ptid = outputs.ptid.astype(float)

    savedir = '/networks/vtn/lab/SDMC_labscience/studies/HVTN/HVTN137/assays/IHC/misc_files/data_processing/'
    today = datetime.date.today().isoformat()
    outputs.to_csv(savedir + f"HVTN137_McElrath_IHC_processed_{today}.txt", sep="\t", index=False)



    summary = pd.pivot_table(
        outputs,
        index=['specrole','ptid'],
        columns='visitno',
        aggfunc='count',
        values='input_file_name',
        dropna=False
    ).dropna(how='all').fillna(0)

    summary.to_excel(savedir + "HVTN137_IHC_sample_summary.xlsx")

    # outputs[[i for i in outputs.columns if 'unit' in i]]

    # outputs.spectype.unique()
    # outputs.control_sample_name.unique()
    # outputs.visitno.unique()
    # outputs[[i for i in outputs.columns if 'unit' in i]]

if __name__=="__main__":
    main()