module GenomeAnnotation
{
    typedef string md5;
    typedef list<md5> md5s;
    typedef string genome_id;
    typedef string feature_id;
    typedef string contig_id;
    typedef string feature_type;

    /* A region of DNA is maintained as a tuple of four components:

		the contig
		the beginning position (from 1)
		the strand
		the length

	   We often speak of "a region".  By "location", we mean a sequence
	   of regions from the same genome (perhaps from distinct contigs).
        */
    typedef tuple<contig_id, int begin, string strand,int length> region_of_dna;

    /*
	a "location" refers to a sequence of regions
    */
    typedef list<region_of_dna> location;
    
    typedef tuple<string comment, string annotator, int annotation_time> annotation;

    typedef structure {
	feature_id id;
	location location;
	feature_type type;
	string function;
	string protein_translation;
	list<string> aliases;
	list<annotation> annotations;
    } feature;

    typedef structure {
	contig_id id;
	string dna;
    } contig;

    typedef structure {
	genome_id id;
	string scientific_name;
	string domain;
	int genetic_code;
	string source;
	string source_id;
	
	list<contig> contigs;
	list<feature> features;
    } genomeTO;

    typedef string subsystem;
    typedef string variant;
    typedef tuple<subsystem,variant> variant_of_subsystem;
    typedef list<variant_of_subsystem> variant_subsystem_pairs;
    typedef string fid;
    typedef string role;
    typedef string function;
    typedef tuple<fid,role> fid_role_pair;
    typedef list<fid_role_pair> fid_role_pairs;
    typedef tuple<fid,function> fid_function_pair;
    typedef list<fid_function_pair> fid_function_pairs;

    typedef structure {
	variant_subsystem_pairs subsystems;
	fid_role_pairs bindings;
	fid_function_pairs assignments;
    } reconstructionTO;

    typedef tuple<fid,md5,location,function> fid_data_tuple;
    typedef list<fid_data_tuple> fid_data_tuples;

    funcdef genomeTO_to_reconstructionTO (genomeTO) returns (reconstructionTO);
    funcdef genomeTO_to_feature_data (genomeTO) returns (fid_data_tuples);
    funcdef reconstructionTO_to_roles (reconstructionTO) returns (list<role>);
    funcdef reconstructionTO_to_subsystems(reconstructionTO) returns (variant_subsystem_pairs);

    /*
     * Given a genome object populated with contig data, perform gene calling
     * and functional annotation and return the annotated genome.
     *
     *  NOTE: Many of these "transformations" modify the input hash and
     *        copy the pointer.  Be warned.
     */
    funcdef annotate_genome(genomeTO) returns (genomeTO);
    funcdef call_RNAs(genomeTO) returns (genomeTO);
    funcdef call_CDSs(genomeTO) returns (genomeTO);
    funcdef find_close_neighbors(genomeTO) returns (genomeTO);
    funcdef assign_functions_to_CDSs(genomeTO) returns (genomeTO);

    /*
     * Given a genome object populated with feature data, reannotate
     * the features that have protein translations. Return the updated
     * genome object.
     */
    funcdef annotate_proteins(genomeTO) returns (genomeTO);

};
