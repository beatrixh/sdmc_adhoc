Beatrix Haddock
2026-05-28

The code/outputs in this dir are for the purpose of generating Fig S5 for the AMP manuscript resubmission. The fig depicts 16 representative sensitive vs resistant curves from 703/704.

See AMP_nAb_curve_Fig_S5_final.ipynb for all inputs used.

## File notes

* code
	+ AMP_nAb_curve_Fig_S5_final.ipynb: generated final outputs
	+ develop_AMP_representative_curve_plots_2026-05-27.ipynb: drafting before outputs were finalized
* outputs 
	+ AMP_nAb_representative_curves_2026-05-28.pdf PDF version of final fig
	+ AMP_nAb_representative_curves_2026-05-28.svg: SVG version of final fig
	+ HVTN703_704_nAb_pnmle_p0/: files from first run of monolix pnlme model fitting all 6 replicates simultaneously
	+ HVTN703_704_nAb_pnmle_p1/: files from second (final) run of monolix pnlme model fitting all 6 replicates simultaneously. This one fixes the population lower limit to zero but allows each isolate to wander.
	+ curves_from_atlas.pdf: earlier draft of fig. This one depicts the curves that labkey fit that I could find -- never finalized, one ppt still missing curves, and some look wrong.
	+ curves_from_atlas_force_lower_limit_0.pdf: earlier draft of fig. This one takes curves_from_atlas and forces the lower limit to be zero. Still has issues.
	+ monolix_fivepl.txt: the model file used for the monolix model. Note I renamed it since running so to rerun this will have to resource this file.
	+ refit_data_simultaneously.pdf: earlier draft of fig. This uses the published IC50/80s but has my fits from HVTN703_704_nAb_pnmle_p0. I accidentally didn't correct the lower limit from prop space to percentage space in this one.
