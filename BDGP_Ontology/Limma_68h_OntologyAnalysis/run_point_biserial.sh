RESULT_DIR='test_results'
RESULT_YAML_DIR='test_result_yamls'

#generate links to the ontology categories
ln -s /g/furlong/garfield/categorical_enrichment_software/annotation_files/GO/V_1.2/flies GO
ln -s /g/furlong/garfield/categorical_enrichment_software/annotation_files/PANTHER_8/newOntologyFiles/GO new_go
ln -s /g/furlong/garfield/categorical_enrichment_software/annotation_files/PANTHER_8/newOntologyFiles/PANTHER new_panther
ln -s /g/furlong/garfield/categorical_enrichment_software/annotation_files/BDGP BDGP
#and a link to the script directory
ln -s /g/furlong/garfield/categorical_enrichment_software/software/point_biserial scripts
#and a directory to allow us to convert names
ln -s /g/furlong/garfield/categorical_enrichment_software/annotation_files/name_conversion name_conversion
#the first step is to restructure the data
./scripts/restr.rb > restr_stdout
#add, as needed, the names of the genes that can be linked to the various databases
./scripts/add_ids_to_data.rb > add_ids_stdout
#and now we do some calculations, first with the standard GO packages
./scripts/calc_rpbs.rb go_biol_proc flyBase > calc_rpbs_new_go_biol_proc_stdout
./scripts/calc_rpbs.rb go_mol_funct flyBase > calc_rpbs_go_mol_funct_stdout
./scripts/calc_rpbs.rb go_cell_comp flyBase > calc_rpbs_go_cell_comp_stdout
./scripts/anal_rpbs.rb go_biol_proc 10
./scripts/anal_rpbs.rb go_mol_funct 10
./scripts/anal_rpbs.rb go_cell_comp 10
#and now for the PANTHER functions
./scripts/calc_rpbs.rb new_panther_prot_class flyBase > calc_rpbs_new_panther_prot_class_stdout
./scripts/calc_rpbs.rb new_panther_pathw flyBase > calc_new_panther_pathw_stdout
./scripts/calc_rpbs.rb new_go_biol_proc flyBase > calc_rpbs_new_go_biol_proc_stdout
./scripts/calc_rpbs.rb new_go_mol_funct flyBase > calc_rpbs_new_go_mol_funct_stdout
./scripts/calc_rpbs.rb new_go_cell_comp flyBase > calc_rpbs_new_go_cell_comp_stdout
./scripts/anal_rpbs.rb new_panther_prot_class 5
./scripts/anal_rpbs.rb new_panther_pathw 5
./scripts/anal_rpbs.rb new_go_biol_proc 5
./scripts/anal_rpbs.rb new_go_mol_funct 5
./scripts/anal_rpbs.rb new_go_cell_comp 5
#and, finally, the BDGP....
./scripts/calc_rpbs.rb bdgp flyBase > calc_rpbs_bdgp_stdout
./scripts/anal_rpbs.rb bdgp 5

#and collect our results in one place
mkdir $RESULT_DIR
mv *results $RESULT_DIR
mkdir $RESULT_YAML_DIR
mv *.yml $RESULT_YAML_DIR

#and do something useful with the names of the categories..this will create files that you can view if you want to see which genes are in a given category
#warning, the script is a little on the slow side
cd $RESULT_YAML_DIR
python ../scripts/add_category_names.py go_biol_proc 
python ../scripts/add_category_names.py go_mol_funct 
python ../scripts/add_category_names.py go_cell_comp
python ../scripts/add_category_names.py new_panther_prot_class 
python ../scripts/add_category_names.py new_panther_pathw 
python ../scripts/add_category_names.py new_go_biol_proc 
python ../scripts/add_category_names.py new_go_mol_funct 
python ../scripts/add_category_names.py new_go_cell_comp 
python ../scripts/add_category_names.py bdgp
cd ..


##next, visualization. This section is optional
ln -s /g/furlong/garfield/categorical_enrichment_software/software/pyEnrichment pyEnrichment
python scripts/create_clusters.py go_biol_proc 0.01 $RESULT_DIR $RESULT_YAML_DIR 10
python scripts/create_clusters.py go_mol_funct 0.01 $RESULT_DIR $RESULT_YAML_DIR 10
python scripts/create_clusters.py go_cell_comp 0.01 $RESULT_DIR $RESULT_YAML_DIR 10
python scripts/create_clusters.py new_panther_prot_class 0.01 $RESULT_DIR $RESULT_YAML_DIR 5
python scripts/create_clusters.py new_panther_pathw 0.01 $RESULT_DIR $RESULT_YAML_DIR 5
python scripts/create_clusters.py new_go_biol_proc 0.01 $RESULT_DIR $RESULT_YAML_DIR 5
python scripts/create_clusters.py new_go_mol_funct 0.01 $RESULT_DIR $RESULT_YAML_DIR 5
python scripts/create_clusters.py new_go_cell_comp 0.01 $RESULT_DIR $RESULT_YAML_DIR 5
python scripts/create_clusters.py bdgp 0.01 $RESULT_DIR $RESULT_YAML_DIR 5


