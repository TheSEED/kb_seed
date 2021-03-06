use strict;
use Data::Dumper;
use Carp;

#
# This is a SAS Component
#

=head1 NAME

external_ids_to_fids

=head1 SYNOPSIS

external_ids_to_fids [arguments] < input > output

=head1 DESCRIPTION


external_ids_to_fids is used to search for the feature ids that a set of given
external_ids maps to.


Example:

    external_ids_to_fids [arguments] < input > output

The standard input should be a tab-separated table (i.e., each line
is a tab-separated set of fields).  Normally, the last field in each
line would contain the identifer. If another column contains the identifier
use

    -c N

where N is the column (from 1) that contains the subsystem.

This is a pipe command. The input is taken from the standard input,
and the output is to the standard output. For each input line, there
can be many output lines, one per feature. The feature id is added to
the end of the line.

=head1 COMMAND-LINE OPTIONS

Usage: external_ids_to_fids [arguments] < input > output


    -c num        Select the identifier from column num
    -i filename   Use filename rather than stdin for input

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut


our $usage = "usage: external_ids_to_fids [-c column] [-prefix] < input > output";

use Bio::KBase::CDMI::CDMIClient;
use Bio::KBase::Utilities::ScriptThing;

my $column;

my $input_file;
my $prefix_flag = 0;

my $kbO = Bio::KBase::CDMI::CDMIClient->new_for_script('c=i' => \$column,
						       'prefix' => \$prefix_flag,
						       'i=s' => \$input_file);
if (! $kbO) { print STDERR $usage; exit }

my $ih;
if ($input_file)
{
    open $ih, "<", $input_file or die "Cannot open input file $input_file: $!";
}
else
{
    $ih = \*STDIN;
}

while (my @tuples = Bio::KBase::Utilities::ScriptThing::GetBatch($ih, 10, $column)) {
    my @h = map { $_->[0] } @tuples;
    my $h = $kbO->external_ids_to_fids(\@h, $prefix_flag);
    for my $tuple (@tuples) {
        #
        # Process output here and print.
        #
        my ($id, $line) = @$tuple;
        my $v = $h->{$id};

        if (! defined($v))
        {
            print STDERR $line,"\n";
        }
        elsif (ref($v) eq 'ARRAY')
        {
            foreach $_ (@$v)
            {
                print "$line\t$_\n";
            }
        }
        else
        {
            print "$line\t$v\n";
        }
    }
}

__DATA__
