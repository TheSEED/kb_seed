use strict;
use Data::Dumper;
use Bio::KBase::Utilities::ScriptThing;
use Carp;

#
# This is a SAS Component
#


=head1 get_entity_Family

The Kbase will support the maintenance of protein families (as sets of Features
with associated translations).  We are initially only supporting the notion of a family
as composed of a set of isofunctional homologs.  That is, the families we
initially support should be thought of as containing protein-encoding genes whose
associated sequences all implement the same function
(we do understand that the notion of "function" is somewhat ambiguous, so let
us sweep this under the rug by calling a functional role a "primitive concept").
We currently support families in which the members are
translations of features, and we think of Features as
having an associated function. Identical protein sequences
as products of translating distinct genes may or may not
have identical functions, and we allow multiple members of
the same Family to share identical protein sequences.  This
may be justified, since in a very, very, very few cases
identical proteins do, in fact, have distinct functions.
We would prefer to reach the point where our Families are
sets of protein sequence, rather than sets of
protein-encoding Features.


Example:

    get_entity_Family -a < ids > table.with.fields.added

would read in a file of ids and add a column for each filed in the entity.

The standard input should be a tab-separated table (i.e., each line
is a tab-separated set of fields).  Normally, the last field in each
line would contain the id. If some other column contains the id,
use

    -c N

where N is the column (from 1) that contains the id.

This is a pipe command. The input is taken from the standard input, and the
output is to the standard output.

=head2 Related entities

The Family entity has the following relationship links:

=over 4
    
=item HasMember Feature

=item HasProteinMember ProteinSequence

=item IsCoupledTo Family

=item IsCoupledWith Family

=item IsFamilyFor Role

=item IsRepresentedIn Genome


=back

=head2 Command-Line Options

=over 4

=item -c Column

Use the specified column to define the id of the entity to retrieve.

=item -h

Display a list of the fields available for use.

=item -fields field-list

Choose a set of fields to return. Field-list is a comma-separated list of 
strings. The following fields are available:

=over 4

=item type

=item release

=item family_function

=item alignment

=back    

=back

=head2 Output Format

The standard output is a tab-delimited file. It consists of the input
file with an extra column added for each requested field.  Input lines that cannot
be extended are written to stderr.  

=cut

use Bio::KBase::CDMI::CDMIClient;
use Getopt::Long;

#Default fields

my @all_fields = ( 'type', 'release', 'family_function', 'alignment' );
my %all_fields = map { $_ => 1 } @all_fields;

my $usage = "usage: get_entity_Family [-h] [-c column] [-a | -f field list] < ids > extended.by.a.column(s)";

my $column;
my $a;
my $f;
my $i = "-";
my @fields;
my $show_fields;
my $geO = Bio::KBase::CDMI::CDMIClient->new_get_entity_for_script('c=i'	   	=> \$column,
								  "a"	   	=> \$a,
								  "h"	   	=> \$show_fields,
								  "show-fields"	=> \$show_fields,
								  "fields=s" 	=> \$f,
								  'i=s'	   	=> \$i);
if ($show_fields)
{
    print STDERR "Available fields: @all_fields\n";
    exit 0;
}
if ($a && $f) { print STDERR $usage; exit 1 }
if ($a)
{
    @fields = @all_fields;
}
elsif ($f) {
    my @err;
    for my $field (split(",", $f))
    {
	if (!$all_fields{$field})
	{
	    push(@err, $field);
	}
	else
	{
	    push(@fields, $field);
	}
    }
    if (@err)
    {
	print STDERR "get_entity_Family: unknown fields @err. Valid fields are: @all_fields\n";
	exit 1;
    }
} else {
    print STDERR $usage;
    exit 1;
}

my $ih;
if ($i eq '-')
{
    $ih = \*STDIN;
}
else
{
    open($ih, "<", $i) or die "Cannot open input file $i: $!\n";
}

while (my @tuples = Bio::KBase::Utilities::ScriptThing::GetBatch($ih, undef, $column)) {
    my @h = map { $_->[0] } @tuples;
    my $h = $geO->get_entity_Family(\@h, \@fields);
    for my $tuple (@tuples) {
        my @values;
        my ($id, $line) = @$tuple;
        my $v = $h->{$id};
	if (! defined($v))
	{
	    #nothing found for this id
	    print STDERR $line,"\n";
     	} else {
	    foreach $_ (@fields) {
		my $val = $v->{$_};
		push (@values, ref($val) eq 'ARRAY' ? join(",", @$val) : $val);
	    }
	    my $tail = join("\t", @values);
	    print "$line\t$tail\n";
        }
    }
}
