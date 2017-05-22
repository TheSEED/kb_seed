# -*- perl -*-
#       This is a SAS Component.
########################################################################
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
#
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License.
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
########################################################################

# usage:  svr_genbank_to_seed  OrgDir genbank.file
my $usage = "svr_genbank_to_seed  OrgDir genbank.file";

use strict;
use warnings;

use FIGV;
use gjogenbank;
use gjoseqlib;

use File::Path;
use Data::Dumper;

my $orgID;
my $org_dir;
my $genbank_file;
my $extension = 0;
use Getopt::Long;
my $rc    = GetOptions("orgID=s"        => \$orgID,
		       "orgDir=s"       => \$org_dir,
                       "gbk"            => \$genbank_file,
		       "extension=i"    => \$extension,
		       );
if (! $rc) { print STDERR "$usage\n"; exit }


if (!defined($org_dir) && !defined($genbank_file) && (@ARGV == 2)) {
    ($org_dir, $genbank_file) = @ARGV;
}

if (!defined($orgID) && ($org_dir =~ m/(\d+\.\d+)$/o)) {
    $orgID = $1;
}


if (-d $org_dir) {
    warn "WARNING: orgDir='$org_dir' already exists\n";
}
else {
    ($rc = mkpath($org_dir)) || die qq(Could not create orgDir='$org_dir', rc=$rc);
}

my $accession = gjogenbank::parse_next_genbank($genbank_file);
my $db_xref   = $accession->{FEATURES}->{source}->[0]->[1]->{db_xref}->[0];
if (!defined($orgID)) {
    if (defined($db_xref)) {
	if ($db_xref =~ m/^taxon:(\d+)$/) {
	    $orgID = "$1.$extension";
	    warn "WARNING: Using orgID='$orgID' based on GenBank db_xref='$db_xref'";
	}
	else {
	    die "ERROR: No orgID and could not parse db_xref='$db_xref' --- aborting";
	}
    }
    else {
	die "ERROR: No orgID and no db_xref --- aborting";
    }
}

if (defined($orgID)) {
    open(GENOME_ID, q(>), "$org_dir/GENOME_ID")
	|| die "Could not write-open file '$org_dir/GENOME_ID'";
    print GENOME_ID ($orgID, "\n");
    close(GENOME_ID);
}
else {
    die "ERROR: no orgID found on command-line, at end of orgDir path, or within GenBank file --- aborting";
}


my $figV = FIGV->new($org_dir);

my $contigs_file = $org_dir.q(/contigs);
open(my $contigs_fh, q(>), $contigs_file)
    || die qq(Could not write-open contigs file \'$contigs_file\');

do {
    my $contig_id   = $accession->{ACCESSION}->[0];
    my $contig_dna  = $accession->{SEQUENCE};
    
    my $this_xref   = $accession->{FEATURES}->{source}->[0]->[1]->{db_xref}->[0];
    if (defined($this_xref) && ($this_xref ne $db_xref)) {
	warn "WARNING: db_xref mismatch for contig='$contig_id': '$this_xref' differs from '$db_xref'\n";
    }
    
    $figV->display_id_and_seq( $contig_id, \$contig_dna, $contigs_fh);
    
    foreach my $cds (@ { $accession->{FEATURES}->{CDS} }) {
#	die Dumper($cds);
	my $gb_loc      = gjogenbank::location( $cds, $accession );
	my $locus       = gjogenbank::genbank_loc_2_seed($contig_id, $gb_loc);
	my $func        = gjogenbank::product( $cds ) || q();
	my $translation = gjogenbank::CDS_translation($cds);
	my $pseudo      = defined( $cds->[1]->{pseudo}->[0] );
	
	my $gene_name   = defined($cds->[1]->{gene}->[0])       ? $cds->[1]->{gene}->[0]       : q();
	my $locus_tag   = defined($cds->[1]->{locus_tag}->[0])  ? $cds->[1]->{locus_tag}->[0]  : q();
	my $protein_id  = defined($cds->[1]->{protein_id}->[0]) ? $cds->[1]->{protein_id}->[0] : q();
	
	my @db_xrefs    = defined($cds->[1]->{db_xref}->[0])    ? @ { $cds->[1]->{db_xref} }   : ();
	
	my @gi_nums     = map { m/GI\:(\d+)/o     ? (q(gi|).$1)     : () } @db_xrefs;
	my @gene_nums   = map { m/GeneID\:(\d+)/o ? (q(GeneID|).$1) : () } @db_xrefs;
	
	my @aliases     = grep { $_ } ($gene_name, $locus_tag, $protein_id, @gi_nums, @db_xrefs, @gene_nums);
	my $aliases     = join(q(,), @aliases);
	
	if ($locus) {
	    my $fid;
	    if ($translation) {
		if ($fid = $figV->add_feature(q(Initial Import), $orgID, q(peg), $locus, $aliases, $translation)) {
		    if ($func) {
			$figV->assign_function($fid, q(master:Initial Import), $func);
		    }		    
		}
	    }
	    elsif ($pseudo) {
		my $sequence = gjogenbank::ftr_seq( $cds, $contig_dna);
		if ($fid = $figV->add_feature(q(Initial Import), $orgID, q(pseudo), $locus, $aliases, $sequence)) {
		    if ($func) {
			if ($func !~ m/pseudogene/i) {
			    $func .= " \# pseudogene";
			}
			
			$figV->assign_function($fid, q(master:Initial Import), $func);
		    }
		}
	    }
	    else {
		die (qq(Could not add feature\n), Dumper($cds));
	    }
	}
	else {
	    warn (qq(Could not parse CDS feature in accession '$contig_id':\n),
		  Dumper($cds),
		  qq(\n));
	}
    }
    
    foreach my $rna (map { $_ ? @$_ : () }
		     map { my $x = $accession->{FEATURES}->{$_}
		       } qw( rRNA tRNA misc_RNA ncRNA )
		     ) {
#	die Dumper($rna);
	
	my $gb_loc      = gjogenbank::location( $rna, $accession );
	my $locus       = gjogenbank::genbank_loc_2_seed($contig_id, $gb_loc);
	my $func        = gjogenbank::product( $rna ) || q();
	my $sequence    = gjogenbank::ftr_seq( $rna, $contig_dna);
	
	my $gene_name   = defined($rna->[1]->{gene}->[0])       ? $rna->[1]->{gene}->[0]       : q();
	my $locus_tag   = defined($rna->[1]->{locus_tag}->[0])  ? $rna->[1]->{locus_tag}->[0]  : q();
	
	my @db_xrefs    = defined($rna->[1]->{db_xref}->[0])    ? @ { $rna->[1]->{db_xref} }   : ();
	
	my @gi_nums     = map { m/GI\:(\d+)/o     ? (q(gi|).$1)     : () } @db_xrefs;
	my @gene_nums   = map { m/GeneID\:(\d+)/o ? (q(GeneID|).$1) : () } @db_xrefs;
	
	my @aliases     = grep { $_ } ($gene_name, $locus_tag, @gi_nums, @db_xrefs, @gene_nums);
	my $aliases     = join(q(,), @aliases);
	
	if ($locus && defined($func) && $sequence) {
	    if (my $fid = $figV->add_feature(q(Initial Import), $orgID, q(rna), $locus, $aliases, $sequence)) {
		if ($func) {
		    if ($func =~ m/23S\s+(ribosomal)?\s+RNA/i) { 
			$func =  "LSU rRNA \#\# 23S rRNA, large subunit ribosomal RNA";
		    }
		    elsif ($func =~ m/16S\s+(ribosomal)?\s+RNA/i) { 
			$func = "SSU rRNA \#\# 16S rRNA, small subunit ribosomal RNA";
		    }
		    elsif ($func =~ m/5S\s+(ribosomal)?\s+RNA/i) {
			$func = "5S rRNA \#\# 5S ribosomal RNA";
		    }
		    
		    $figV->assign_function($fid, q(master:Initial Import), $func);
		}
	    }
	    else {
		die (qq(Could not add feature\n), Dumper($rna));
	    }
	}
	else {
	    warn (qq(Could not parse RNA feature in accession '$contig_id': locus='$locus', func='$func', sequence='$sequence'\n),
		  Dumper($rna),
		  qq(\n));
	}
    }
} until (not defined($accession = gjogenbank::parse_next_genbank($genbank_file)));
exit(0);
