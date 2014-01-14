#!/usr/bin/perl -w

#
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
#

    use strict;
    use Stats;
    use SeedUtils;
    use Data::Dumper;
    use Bio::KBase::CDMI::CDMI;
    use Bio::KBase::CDMI::CDMILoader;


=head1 CDMI Test Script

    CDMITest [options]

The command-line options are as specified in L<Bio::KBase::CDMI::CDMI/new_for_script>. There are
no positional parameters.

=cut

$| = 1;
# Connect to the database.
my $cdmi = Bio::KBase::CDMI::CDMI->new_for_script();
if (! $cdmi) {
    print "usage: CDMITest [options]\n";
} else {
    my @genomes = $cdmi->GetAll("Genome WasSubmittedBy", "Genome(dna-size) < ? AND Genome(prokaryotic) = 1 AND Genome(complete) = 1 ORDER BY Genome(dna-size) DESC",
            [500000], 'Genome(id) Genome(dna-size) Genome(pegs) WasSubmittedBy(to-link) Genome(source-id)');
    print "ID\tdna size\tpegs\tsource\tsource-id\ttaxonomy\n";
    for my $genome (@genomes) {
        my @taxonomy = $cdmi->Taxonomy($genome->[0]);
        print join("\t", @$genome, @taxonomy) . "\n";
    }
}