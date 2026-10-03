# Method dictionary

Scan the public Python API by name. Each entry jumps to the heading on
the hand-written module page. Examples live in `tests/<module>_examples/`
and are included on those pages.

Regenerate with `python scripts/gen_api_dictionary.py`.

<div class="api-dictionary" markdown="1">

## Structure Manager

- [`assign_secondary_structure_map`](structure_manager.md#assign_secondary_structure_map)
- [`parse_pdb`](structure_manager.md#parse_pdb)
- [`StructureManager.__init__`](structure_manager.md#constructor)
- [`StructureManager.align_to_axis`](structure_manager.md#align_to_axis)
- [`StructureManager.apply_mempro_orientation`](structure_manager.md#apply_mempro_orientation)
- [`StructureManager.assign_secondary_structure`](structure_manager.md#assign_secondary_structure)
- [`StructureManager.auto_detect_molecules`](structure_manager.md#auto_detect_molecules)
- [`StructureManager.center_atoms`](structure_manager.md#center_atoms)
- [`StructureManager.delete_atoms`](structure_manager.md#delete_atoms)
- [`StructureManager.get_chains`](structure_manager.md#get_chains)
- [`StructureManager.get_residues`](structure_manager.md#get_residues)
- [`StructureManager.get_secondary_structure_summary`](structure_manager.md#get_secondary_structure_summary)
- [`StructureManager.get_structure_info`](structure_manager.md#get_structure_info)
- [`StructureManager.load_from_pdb_id`](structure_manager.md#load_from_pdb_id)
- [`StructureManager.load_structure`](structure_manager.md#load_structure)
- [`StructureManager.rename_chain`](structure_manager.md#rename_chain)
- [`StructureManager.rename_chain_by_indices`](structure_manager.md#rename_chain_by_indices)
- [`StructureManager.rename_residues`](structure_manager.md#rename_residues)
- [`StructureManager.rename_residues_by_indices`](structure_manager.md#rename_residues_by_indices)
- [`StructureManager.renumber_residues`](structure_manager.md#renumber_residues)
- [`StructureManager.renumber_residues_by_indices`](structure_manager.md#renumber_residues_by_indices)
- [`StructureManager.rotate_atoms`](structure_manager.md#rotate_atoms)
- [`StructureManager.save_pdb`](structure_manager.md#save_pdb)
- [`StructureManager.select_atoms`](structure_manager.md#select_atoms)
- [`StructureManager.select_by_criteria`](structure_manager.md#select_by_criteria)
- [`StructureManager.translate_atoms`](structure_manager.md#translate_atoms)

## Preparation

- [`complete_missing_heavy_atoms`](preparation.md#complete_missing_heavy_atoms)
- [`count_protein_hydrogens`](preparation.md#count_protein_hydrogens)
- [`extract_summary_section`](preparation.md#extract_summary_section)
- [`modify_pdb_based_on_summary`](preparation.md#modify_pdb_based_on_summary)
- [`parse_summary_section`](preparation.md#parse_summary_section)
- [`PreparationManager.__init__`](preparation.md#constructor)
- [`PreparationManager.apply_disulfide_bonds`](preparation.md#apply_disulfide_bonds)
- [`PreparationManager.apply_protonation_states`](preparation.md#apply_protonation_states)
- [`PreparationManager.build_titration_figure`](preparation.md#build_titration_figure)
- [`PreparationManager.complete_missing_heavy_atoms`](preparation.md#complete_missing_heavy_atoms)
- [`PreparationManager.detect_disulfide_bonds`](preparation.md#detect_disulfide_bonds)
- [`PreparationManager.extract_summary`](preparation.md#extract_summary)
- [`PreparationManager.get_available_states`](preparation.md#get_available_states)
- [`PreparationManager.get_default_protonation_state`](preparation.md#get_default_protonation_state)
- [`PreparationManager.get_ph_titration_curve`](preparation.md#get_ph_titration_curve)
- [`PreparationManager.get_residue_statistics`](preparation.md#get_residue_statistics)
- [`PreparationManager.parse_summary`](preparation.md#parse_summary)
- [`PreparationManager.restore_original_residue_numbers`](preparation.md#restore_original_residue_numbers)
- [`PreparationManager.run_analysis`](preparation.md#run_analysis)
- [`PreparationManager.run_pdb4amber_with_cap_fix`](preparation.md#run_pdb4amber_with_cap_fix)
- [`run_propka`](preparation.md#run_propka)
- [`strip_protein_hydrogens`](preparation.md#strip_protein_hydrogens)

## Residue mapping

- [`parse_residue_mapping_by_chain`](preparation.md#parse_residue_mapping_by_chain)
- [`parse_residue_mapping_file`](preparation.md#parse_residue_mapping_file)

## Protein capping

- [`cap_protein`](preparation.md#cap_protein)
- [`detect_terminal_caps`](preparation.md#detect_terminal_caps)
- [`is_already_capped`](preparation.md#is_already_capped)
- [`ProteinCapper.__init__`](preparation.md#constructor)
- [`ProteinCapper.remove_hydrogens_and_cap`](preparation.md#remove_hydrogens_and_cap)

## Builder

- [`Builder.__init__`](builder.md#constructor)
- [`Builder.cancel_preparation`](builder.md#cancel_preparation)
- [`Builder.generate_preparation_inputs`](builder.md)
- [`Builder.modify_residue_names_for_propka`](builder.md#modify_residue_names_for_propka)
- [`Builder.prepare_system`](builder.md#prepare_system)
- [`Builder.prepare_system_stage1_for_propka`](builder.md#prepare_system_stage1_for_propka)
- [`Builder.prepare_system_stage2_for_propka`](builder.md#prepare_system_stage2_for_propka)
- [`Builder.run_preparation`](builder.md#run_preparation)
- [`Builder.set_configuration`](builder.md#set_configuration)
- [`Builder.validate_system_inputs`](builder.md#validate_system_inputs)
- [`Builder.wait_for_completion`](builder.md#wait_for_completion)

## Job monitor

- [`JobMonitor.__init__`](builder.md#constructor)
- [`JobMonitor.cleanup_stale_jobs`](builder.md#cleanup_stale_jobs)
- [`JobMonitor.export_job_summary`](builder.md)
- [`JobMonitor.get_active_jobs`](builder.md#get_active_jobs)
- [`JobMonitor.get_completed_jobs`](builder.md#get_completed_jobs)
- [`JobMonitor.get_job`](builder.md#get_job)
- [`JobMonitor.get_job_statistics`](builder.md#get_job_statistics)
- [`JobMonitor.refresh_job`](builder.md#refresh_job)
- [`JobMonitor.remove_job`](builder.md#remove_job)
- [`JobMonitor.scan_for_jobs`](builder.md#scan_for_jobs)
- [`JobMonitor.set_working_directory`](builder.md)

## Force fields

- [`ForceFieldManager.__init__`](builder.md#constructor)
- [`ForceFieldManager.get_available_anions`](builder.md#get_available_anions)
- [`ForceFieldManager.get_available_cations`](builder.md#get_available_cations)
- [`ForceFieldManager.get_available_lipids`](builder.md#get_available_lipids)
- [`ForceFieldManager.get_force_field_info`](builder.md#get_force_field_info)
- [`ForceFieldManager.get_lipid_force_fields`](builder.md#get_lipid_force_fields)
- [`ForceFieldManager.get_protein_force_fields`](builder.md#get_protein_force_fields)
- [`ForceFieldManager.get_recommendations`](builder.md#get_recommendations)
- [`ForceFieldManager.get_water_models`](builder.md#get_water_models)
- [`ForceFieldManager.validate_combination`](builder.md#validate_combination)
- [`ForceFieldManager.validate_ion`](builder.md#validate_ion)
- [`ForceFieldManager.validate_lipid`](builder.md#validate_lipid)

## MemPrO

- [`apply_orientation_transform`](mempro.md#apply_orientation_transform)
- [`compute_orientation_transform`](mempro.md#compute_orientation_transform)
- [`MemPrO.__init__`](mempro.md#constructor)
- [`MemPrO.build_command`](mempro.md#build_command)
- [`MemPrO.get_oriented_pdb`](mempro.md#get_oriented_pdb)
- [`MemPrO.is_available`](mempro.md#is_available)
- [`MemPrO.parse_results`](mempro.md#parse_results)
- [`MemPrO.run`](mempro.md#run)

## Hydration

- [`build_hydrate_inp_text`](hydration.md#build_hydrate_inp_text)
- [`check_packmol_available`](hydration.md#check_packmol_available)
- [`detect_hydrogen_status`](hydration.md#detect_hydrogen_status)
- [`estimate_cavity_volume`](hydration.md#estimate_cavity_volume)
- [`hydrate_cavity`](hydration.md#hydrate_cavity)
- [`prepare_hydration_job`](hydration.md#prepare_hydration_job)
- [`preview_hydrate_inp`](hydration.md#preview_hydrate_inp)
- [`run_custom_packmol`](hydration.md#run_custom_packmol)

## Equilibration

- [`adjust_protocol_for_missing_protein`](equilibration.md#adjust_protocol_for_missing_protein)
- [`AmberEquilibrationManager.__init__`](equilibration.md#constructor)
- [`AmberEquilibrationManager.build_group_restraint_block`](equilibration.md)
- [`AmberEquilibrationManager.ensure_inpcrd_box`](equilibration.md)
- [`AmberEquilibrationManager.ensure_prmtop_box`](equilibration.md)
- [`AmberEquilibrationManager.find_system_files`](equilibration.md#find_system_files)
- [`AmberEquilibrationManager.generate_input_file`](equilibration.md)
- [`AmberEquilibrationManager.generate_mdin_file`](equilibration.md)
- [`AmberEquilibrationManager.generate_run_script`](equilibration.md)
- [`AmberEquilibrationManager.get_default_selections`](equilibration.md)
- [`AmberEquilibrationManager.get_default_stage_params`](equilibration.md#get_default_stage_params)
- [`AmberEquilibrationManager.setup_amber_equilibration`](equilibration.md)
- [`EquilibrationStage.replace`](equilibration.md)
- [`EquilibrationStage.to_dict`](equilibration.md)
- [`GROMACSEquilibrationManager.__init__`](equilibration.md#constructor)
- [`GROMACSEquilibrationManager.convert_from_amber`](equilibration.md#convert_from_amber)
- [`GROMACSEquilibrationManager.find_system_files`](equilibration.md#find_system_files)
- [`GROMACSEquilibrationManager.generate_com_colvars_config`](equilibration.md#generate_com_colvars_config)
- [`GROMACSEquilibrationManager.generate_index_ndx`](equilibration.md)
- [`GROMACSEquilibrationManager.generate_mdp_file`](equilibration.md#generate_mdp_file)
- [`GROMACSEquilibrationManager.generate_posres_itp`](equilibration.md)
- [`GROMACSEquilibrationManager.generate_run_script`](equilibration.md)
- [`GROMACSEquilibrationManager.get_default_selections`](equilibration.md)
- [`GROMACSEquilibrationManager.get_default_stage_params`](equilibration.md#get_default_stage_params)
- [`GROMACSEquilibrationManager.setup_gromacs_equilibration`](equilibration.md#setup_gromacs_equilibration)
- [`NAMDEquilibrationManager.__init__`](equilibration.md#constructor)
- [`NAMDEquilibrationManager.count_all_selections`](equilibration.md)
- [`NAMDEquilibrationManager.count_selection_atoms`](equilibration.md)
- [`NAMDEquilibrationManager.find_system_files`](equilibration.md#find_system_files)
- [`NAMDEquilibrationManager.generate_bilayer_thickness_colvar`](equilibration.md)
- [`NAMDEquilibrationManager.generate_charmm_gui_config_file`](equilibration.md#generate_charmm_gui_config_file)
- [`NAMDEquilibrationManager.generate_colvars_file`](equilibration.md)
- [`NAMDEquilibrationManager.generate_com_colvars_config`](equilibration.md#generate_com_colvars_config)
- [`NAMDEquilibrationManager.generate_config_file`](equilibration.md)
- [`NAMDEquilibrationManager.generate_restraints_file`](equilibration.md#generate_restraints_file)
- [`NAMDEquilibrationManager.generate_restraints_file_mda`](equilibration.md#generate_restraints_file_mda)
- [`NAMDEquilibrationManager.generate_run_script`](equilibration.md)
- [`NAMDEquilibrationManager.get_default_selections`](equilibration.md)
- [`NAMDEquilibrationManager.get_default_stage_params`](equilibration.md#get_default_stage_params)
- [`NAMDEquilibrationManager.load_charmm_gui_template`](equilibration.md)
- [`NAMDEquilibrationManager.run_equilibration`](equilibration.md)
- [`NAMDEquilibrationManager.setup_namd_equilibration`](equilibration.md#setup_namd_equilibration)
- [`OpenMMEquilibrationManager.__init__`](equilibration.md#constructor)
- [`OpenMMEquilibrationManager.find_system_files`](equilibration.md#find_system_files)
- [`OpenMMEquilibrationManager.generate_openmm_config`](equilibration.md)
- [`OpenMMEquilibrationManager.generate_openmm_restraint_files`](equilibration.md)
- [`OpenMMEquilibrationManager.generate_run_script`](equilibration.md)
- [`OpenMMEquilibrationManager.get_default_selections`](equilibration.md)
- [`OpenMMEquilibrationManager.get_default_stage_params`](equilibration.md#get_default_stage_params)
- [`OpenMMEquilibrationManager.setup_openmm_equilibration`](equilibration.md#setup_openmm_equilibration)
- [`pdb_has_protein`](equilibration.md#pdb_has_protein)

## Equilibration extras

- [`cluster_engine_executable`](equilibration.md#cluster_engine_executable)
- [`resolve_cluster_launch_script`](equilibration.md#resolve_cluster_launch_script)
- [`write_cluster_run_script`](equilibration.md#write_cluster_run_script)

## NAMD / energy analysis

- [`EnergyAnalyzer.__init__`](analysis.md#constructor)
- [`EnergyAnalyzer.get_available_properties`](analysis.md#get_available_properties)
- [`EnergyAnalyzer.get_statistics`](analysis.md#get_statistics)
- [`EnergyAnalyzer.plot_energy`](analysis.md#plot_energy)
- [`EnergyAnalyzer.plot_properties`](analysis.md#plot_properties)
- [`parse_namd_log`](analysis.md#parse_namd_log)
- [`run_energetic_analysis`](analysis.md#run_energetic_analysis)

## Trajectory analysis

- [`prepare_structural_inputs`](analysis.md#prepare_structural_inputs)
- [`run_structural_analysis`](analysis.md#run_structural_analysis)
- [`TrajectoryAnalyzer.__init__`](analysis.md#constructor)
- [`TrajectoryAnalyzer.calculate_distances`](analysis.md#calculate_distances)
- [`TrajectoryAnalyzer.calculate_radius_of_gyration`](analysis.md#calculate_radius_of_gyration)
- [`TrajectoryAnalyzer.calculate_rmsd`](analysis.md#calculate_rmsd)
- [`TrajectoryAnalyzer.calculate_rmsf`](analysis.md#calculate_rmsf)
- [`TrajectoryAnalyzer.clear_analysis_cache`](analysis.md#clear_analysis_cache)
- [`TrajectoryAnalyzer.plot_distances`](analysis.md#plot_distances)
- [`TrajectoryAnalyzer.plot_radius_of_gyration`](analysis.md#plot_radius_of_gyration)
- [`TrajectoryAnalyzer.plot_rmsd`](analysis.md#plot_rmsd)
- [`TrajectoryAnalyzer.plot_rmsf`](analysis.md#plot_rmsf)
- [`TrajectoryAnalyzer.plot_summary`](analysis.md#plot_summary)
- [`TrajectoryAnalyzer.time_array_for_analysis`](analysis.md)

## Bilayer analysis

- [`BilayerTrajectoryAnalyzer.__init__`](analysis.md#constructor)
- [`BilayerTrajectoryAnalyzer.calculate_area_per_lipid`](analysis.md#calculate_area_per_lipid)
- [`BilayerTrajectoryAnalyzer.calculate_membrane_thickness`](analysis.md#calculate_membrane_thickness)
- [`BilayerTrajectoryAnalyzer.plot_area_per_lipid`](analysis.md)
- [`BilayerTrajectoryAnalyzer.plot_membrane_thickness`](analysis.md)
- [`BilayerTrajectoryAnalyzer.universe`](analysis.md)
- [`run_bilayer_analysis`](analysis.md#run_bilayer_analysis)

## OpenMM analysis

- [`OpenMMLogAnalyzer.__init__`](analysis.md#constructor)
- [`OpenMMLogAnalyzer.get_statistics`](analysis.md#get_statistics)
- [`OpenMMLogAnalyzer.plot_energy`](analysis.md#plot_energy)
- [`OpenMMLogAnalyzer.plot_properties`](analysis.md#plot_properties)
- [`run_openmm_energetic_analysis`](analysis.md#run_openmm_energetic_analysis)

## GROMACS analysis

- [`GROMACSLogEnergyAnalyzer.__init__`](analysis.md#constructor)
- [`GROMACSLogEnergyAnalyzer.get_available_properties`](analysis.md#get_available_properties)
- [`GROMACSLogEnergyAnalyzer.get_statistics`](analysis.md#get_statistics)

## Amber analysis

- [`AmberLogEnergyAnalyzer.__init__`](analysis.md#constructor)
- [`AmberLogEnergyAnalyzer.get_available_properties`](analysis.md#get_available_properties)
- [`run_amber_energetic_analysis`](analysis.md#run_amber_energetic_analysis)

## FATSLiM / GridMAT

- [`run_fatslim_apl`](analysis.md#run_fatslim_apl)
- [`run_gridmat_md_apl`](analysis.md#run_gridmat_md_apl)

## Cluster helpers

- [`connect_ssh`](equilibration.md#connect_ssh)
- [`rsync_from_remote`](equilibration.md#rsync_from_remote)
- [`rsync_to_remote`](equilibration.md#rsync_to_remote)
- [`run_remote`](equilibration.md#run_remote)

</div>
