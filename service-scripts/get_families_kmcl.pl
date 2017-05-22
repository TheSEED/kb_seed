#
# This is a SAS Component
#


=head1 get_families

Generate protein families (of isofunctional homologs) using kmer technology.

------

Example:

    get_families -d Data.kmers -f Families/families -s Seqs.Fasta < genomes > families

This uses a Data.kmer directory built to support kmer_guts processing.
We suggest using the one in pubSEED (Global/Data.kmers).  The invocation causes a
set of "families" files to be generated in the existing Families directory.  They will all
be prefixed with the word "families.".  The final set of protein
families is written to STDOUT.

Seqs.Fasta is a directory that contains protein fasta files.  The file names
must be genome IDs.  Thus, it is assumed that 

    Seqs.Fasta/83333.1

would be the peg translations for E.coli (assuming that you wished E.coli
to be one of the genomes from which families get produced).

The files in Seqs.Fasta used in constructing the families is determined
by the contents of STDIN (each input line contains just a genome ID to be included).
Each included genome must have a corresponding fasta file in Seqs.Fasta.

------

The standard input should be a file of genome IDs (an input line must end in
a genome ID (/\d+\.\d+$/)).  These will be the genomes from which families are constructed.

------

Now, let us summarize the steps used to generate the families.  We go through the following steps:

    1. the script get_families_1.pl runs all of the PEG translations from all of
       the genomes (specified in STDIN) through kmer_search, which uses kmers to attempt
       assignment of function.  Successfully called PEGs are written to tmp.$$.calls.  Those
       that were not assigned a function are written to tmp.$$.missed.

    2. The PEGs in tmp.$$.missed are thn processed using svr_representative_sequences, which
       generates sets based on blast for those kmers not handled by kmer_search.  They are
       all assigned the function "hypothetical protein", and the sets are written to
       families.missed.

    3. Then, we go through the PEGs that were called by kmers (and recorded in tmp.$$.calls).
       This is done by get_families_3.
       We form potential sets as all PEGs assigned the same function.  For each "function-based set"
       we count the number of PEGs from each genome.  If 90% (i.e., the cutoff parameter 
       defines this value, which defaults to 0.9) of the genomes represented
       in the set have only one PEG in the set, the set is considered "good" and written to
       "families.good".  Otherwise, the set is written to tmp.$$.bad.

    4. Now, get_families-4 is used to  process the families written to tmp.$$.bad.
       Note that kmer assignment of function may "group" disparate
       sequences into a single function.  If the manual assignments of
       function upon which the kmers were derived correctly assigned
       one set of sequences to a function F and incorrectly assigned a
       second set to function F, then sequences that get assigned a
       sequence F by kmers may have gotten the assignment due to signature
       kmers from either of 2 distinct sets.  For that matter, if two
       non-homologous classes of sequences both have proteins with a
       common function (due to non-orthologous replacement), you will
       have distinct kmers that produce a common call, and there may
       legitimately be multiple instances of the function in a single
       genome. 

       Anyway, for each set we wish to split, we compute the kmers
       that are associated with each peg in the set.  Then, we sort
       the pegs in the set based on the number of kmers that hit each
       peg.  Then, we make passes through the sorted set of pegs,
       seeding a new set and adding pegs that share at least MatchN
       common kmers (we set MatchN to 3, usually).  Each pass induces
       a new set written to families.bad.fixed (possibly singletons).

Finally, the sets from families.good, families.bad.fixed, and
families.missed are all gathered and renumbered and written to STDOUT..

=head2 Command-Line Options

=over 4

=item -d Data

This is a Data directory usable by kmer_guts.  I suggest using the one in
the Global directory (FIGfisk/FIG/Data/Global/Data.kmers).

=item -m MatchN

During one step, families that may need to be split use an algorithm in which
two PEGs are kept in the same family iff they share at least MatchN kmers (see above)

=item -i IdentityFraction

This is the fraction used by Gary's representative_sequences when forming families
of the sequences left uncalled by kmers (see above)

=item -f FamilyFilesPrefix

The prefix used when writing files recording subfamilies.

=item -c cutoff used to differntiate between "good" and "bad" "called families"

if a fraction more than "cutoff" genomes in a family have just one PEG,
the family is "good"; else it is "bad", and an attempt will be made to split it.

=item -s Seqs.Fasta

The directory from which the translations of PEGs from each genome are 
used.

=back

=head2 Output Format

Output is written to STDOUT and constitutes the derived protein families (which
include singletons).  An 8-column, tab-separated table is written:

    FamilyID - an integer
    Function - function assigned to family
    SubFunction - the Function and an integer (SubFunction) together uniquely
                  determine the FamilyID.  Another way to look at it is

                    a) each family is assigned a unique ID and a function
                    b) multiple families can have the same function (consider
                       "hypothetical protein")
                    c) the Function+SubFunction uniquely determine the FamilyID
    PEG
    LengthProt - the length of the translated PEG
    Mean       - the mean length of PEGs in the family
    StdDev     - standard deviation of lengths for family
    Z-sc       - the Z-score associated with the length of this PEG

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use SeedEnv;
use gjoseqlib;
use File::Slurp;
use File::Basename;
use File::Temp 'tempdir';
use POSIX ':sys_wait_h';
use IPC::Run 'run';
use IO::Handle;
use Time::HiRes 'gettimeofday';

my $usage = "usage: get_families -d Data -s Seqs [-l logfile]  < genomes\n";
my $dataD;
my $seqsD;
my $origSeqsD;
my $matchN = 3;
my $iden = 0.5;
my $families;
my $tmpdir;
my $inflation = 3.0;
my $cutoff = 0.9;  # fraction of members with uniq genomes to be "good"
my $logfile;
my $parallel = 1;
 
my $rc  = GetOptions('d=s' => \$dataD,
		     'l=s' => \$logfile,
		     'm=i' => \$matchN,
		     'i=f' => \$iden,
		     'I=f' => \$inflation,
		     'f=s' => \$families,
		     'c=f' => \$cutoff,
		     't=s' => \$tmpdir,
		     'p=i' => \$parallel,
                     's=s' => \$seqsD);

if ((! $rc) || (! $dataD) || (! $seqsD) || (! $families))
{ 
    print STDERR $usage; exit ;
}

$origSeqsD = $seqsD;

if ($logfile)
{
    open(LOG, ">>", $logfile) or die "Cannot open $logfile for append: $!";
}
else
{
    open(LOG, ">&STDERR");
}
LOG->autoflush(1);

my $overall_tstart = gettimeofday;
print LOG "start overall_run $overall_tstart\n";

#
# Start up our kmer servers.
#

#my $skip_compute = defined($tmpdir);
my $skip_compute;
if ($tmpdir)
{
    if (! -d $tmpdir)
    {
	mkdir($tmpdir);
    }
}
else
{
    $tmpdir = tempdir();
}

#goto x;

if (0)
{
    
    my $guts_port_file = "$tmpdir/guts_port";
    my @guts_cmd = ("kmer_guts", "-D", $dataD, "-l", 0, "-L", $guts_port_file, "-P", $$);
    
    unlink($guts_port_file);
    
    my $guts_pid = fork;
    defined($guts_pid) or die "fork failed: $!";
    if ($guts_pid == 0)
    {
	print "Starting guts_cmd in $$: @guts_cmd\n";
	exec(@guts_cmd);
	die "guts cmd failed with $?: @guts_cmd\n";
    }
    
    #
    # Wait for server to start
    #
    
    my $guts_port;
    while (1)
    {
	my $kid = waitpid($guts_pid, WNOHANG);
	if ($kid)
	{
	    die "kmer_guts server did not start\n";
	}
	if (-s $guts_port_file)
	{
	    $guts_port = read_file($guts_port_file);
	    print "Read '$guts_port' from $guts_port_file\n";
	    chomp $guts_port;
	    if ($guts_port !~ /^\d+$/)
	    {
		kill(1, $guts_pid);
		kill(9, $guts_pid);
		die "Invalid guts port '$guts_port'\n";
	    }
	    last;
	}
	sleep(0.1);
    }

    print "Running guts $guts_pid on port $guts_port\n";
}


#
# Start the stateful kmer server
#

my $kser_port_file = "$tmpdir/kser_port";
#my @kser_cmd = ("kser", "--listen-port-file", $kser_port_file, 0, 'localhost', $guts_port, $dataD);
my @kser_cmd = ("/scratch/olson/close_kmers/kser",
		"--reserve-mapping", 10000000,
		"--n-kmer-threads", $parallel,
		"--listen-port-file", $kser_port_file, 0, $dataD);

unlink($kser_port_file);

my $kser_pid = fork;
if (!defined($kser_pid))
{
#    kill(1, $guts_pid);
#    kill(9, $guts_pid);
    die "fork failed: $!";
}
if ($kser_pid == 0)
{
    print "Starting kser_cmd in $$: @kser_cmd\n";
    exec(@kser_cmd);
    die "Kser did not start $?: @kser_cmd\n";
}

my $kser_port;
while (1)
{
    my $kid = waitpid($kser_pid, WNOHANG);
    if ($kid)
    {
#	kill(1, $guts_pid);
#	kill(9, $guts_pid);	
	die "kser server did not start\n";
    }
    if (-s $kser_port_file)
    {
	$kser_port = read_file($kser_port_file);
	print "Read '$kser_port' from $kser_port_file\n";
	chomp $kser_port;
	if ($kser_port !~ /^\d+$/)
	{
#	    kill(1, $guts_pid);
#	    kill(9, $guts_pid);
	    kill(1, $kser_pid);
	    kill(9, $kser_pid);
	    die "Invalid kser port '$kser_port'\n";
	}
	last;
    }
    sleep(0.1);
}

print "Kser running pid $kser_pid on port $kser_port\n";

my $add_url = "http://localhost:$kser_port/add";
my $matrix_url = "http://localhost:$kser_port/matrix";

$SIG{__DIE__} = sub {
    kill 1, $kser_pid;
#    kill 1, $guts_pid;
};

#
# If we have an nr directory, construct a mapped sequence file and use its sequence directory
# for our sequences.
#
x:

my $fam_dir = dirname($families);
my %pegs_to_inflate;
my %pegs_for_singleton_fams;

if (-s "$fam_dir/nr/peg.synonyms")
{
    if (1 || ! -d "$fam_dir/nr-seqs")
    {
	#
	# In the event files already exist in nr-seqs, remove them to eliminate duplicates.
	#

	my $sdir = "$fam_dir/nr-seqs";
	-d $sdir or mkdir($sdir);

	my $open_genome;
	my $open_fh;

	opendir(DH, $sdir) or die "Cannot opendir $sdir: $!";
	while (my $f = readdir(DH))
	{
	    next if $f =~ /^\./;
	    my $p = "$sdir/$f";
	    if (-s $p)
	    {
		print STDERR "Unlinking existing NR file $p\n";
		unlink($p) or die "Cannot unlink $p: $!";
	    }
	}
	closedir(DH);

	my $tstart = gettimeofday;
	print LOG "start pegsyn processing $tstart\n";

	my %f;
	my %g;
	open(N, "<", "$fam_dir/nr/peg.synonyms") or die "Cannot open $fam_dir/nr/peg.synonyms: $!";

	# gnl|md5|ffdb15eded2cbc4a89b6e9ece42b30e9,108  fig|1005475.3.peg.3229,108;fig|1182692.3.peg.3088,108;fig|1182722.3.peg.1758,108
	while (<N>)
	{
	    if (/^([^,]+),\d+\t(fig\|(\d+\.\d+)\.peg\.\d+),\d+(.*)/)
	    {
		my $md5 = $1;
		my $ref = $2;
		my $genome = $3;
		my $rest = $4;
		    
		$f{$md5} = $ref;
		$g{$md5} = $genome;

		$pegs_for_singleton_fams{$ref} = 1;
		my $lst = [];
		while ($rest =~ /(fig\|\d+\.\d+\.peg\.\d+)/mg)
		{
		    push(@$lst, $1);
		}
		$pegs_to_inflate{$ref} = $lst;
	    }
	}

	open(N, "<", "$fam_dir/nr/nr") or die "Cannot open $fam_dir/nr/nr: $!";
	while (<N>)
	{
	    if (/^>(\S+)(.*)/)
	    {
		my $fid = $f{$1};
		my $genome = $g{$1};

		if (!$genome)
		{
		    warn "Cannot map $1 to genome\n";
		    next;
		}

		if ($genome ne $open_genome)
		{
		    close($open_fh) if $open_fh;
		    open($open_fh, ">>", "$sdir/$genome") or die "Cannot append $sdir/$genome: $!";
		    $open_genome = $genome;
		}		    

		print $open_fh ">$f{$1}$2\n";
	    }
	    else
	    {
		print $open_fh $_;
	    }

	}
	close(N);
	close($open_fh);
	my $tend = gettimeofday;
	my $elap = $tstart - $tend;
	print LOG "finish pegsyn processing $tend $elap\n";
    }

    $seqsD = "$fam_dir/nr-seqs";
}
# &SeedUtils::run("get_families_1 -u $add_url -d $dataD -s $seqsD > $tmpdir/calls 2> $tmpdir/missed");
my $tstart = gettimeofday;
print LOG "start get_families_1 $tstart\n";
&SeedUtils::ipc_run(["get_families_1", "-u", $add_url, "-d", $dataD, "-s", $seqsD],
			     '>', "$tmpdir/calls",
			     '2>', "$tmpdir/missed");
my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish get_families_1 $tend $elap\n";

mkdir "$families.fasta";
mkdir "$families.dist";
mkdir "$families.mcl";

if (!$skip_compute)
{
#    &SeedUtils::run("get_families_5 -I $matrix_url $seqsD $tmpdir/calls $tmpdir/missed $families.map $families.fasta $families.dist $families.mcl > $families.good $tmpdir/missed2");
#    &SeedUtils::run("get_families_2 -i $iden -s $seqsD < $tmpdir/missed2 > $families.missed");

    my $tstart = gettimeofday;
    print LOG "start get_families_5 $tstart\n";
    my $cmd = ["get_families_5",
			 "-o", "$families.good",
			 "-I", $inflation,
			 "-p", $parallel,
			 $matrix_url, $seqsD, "$tmpdir/calls",
			 "$tmpdir/missed", "$families.map", "$families.fasta",
			 "$families.dist", "$families.mcl", "$tmpdir/missed2"];
    print STDERR "Run: @$cmd\n";
    &SeedUtils::ipc_run($cmd);
    my $tend = gettimeofday;
    my $elap = $tend - $tstart;
    print LOG "finish get_families_5 $tend $elap\n";
    
    my $tstart = gettimeofday;
    print LOG "start get_families_2 $tstart\n";
    &SeedUtils::ipc_run(["get_families_2", "-i", $iden, "-s", $seqsD],
			'<', "$tmpdir/missed2",
			'>', "$families.missed");
    
    my $tend = gettimeofday;
    my $elap = $tend - $tstart;
    print LOG "finish get_families_2 $tend $elap\n";
    open(N, ">", "$families.bad.fixed");
    close(N);
}

#
# Inflate generated families if we ran from NR.
#

my $tstart = gettimeofday;
print LOG "start inflate $tstart\n";

if (%pegs_to_inflate)
{
    inflate("$families.good", "$families.good.refs");
    inflate("$families.missed", "$families.missed.refs");

    #
    # Use families.bad.fixed for our singletons.
    #

    my %subfam;
    open(CALLS, "<", "$tmpdir/calls") or die "Cannot open $tmpdir/calls: $!";
    open(B, ">", "$families.bad.fixed") or die "Cannot write $families.bad.fixed: $!";
    while (<CALLS>)
    {
	chomp;
	my($peg, $fun) = split(/\t/);
	if ($pegs_for_singleton_fams{$peg})
	{
	    my $sub;
	    if (exists($subfam{$fun}))
	    {
		$subfam{$fun}++;
		$sub = $subfam{$fun};
	    }
	    else
	    {
		$sub = 1;
		$subfam{$fun} = 1;
	    }
		
	    print B "$fun\t$sub\t$peg\n";
	    delete $pegs_for_singleton_fams{$peg};
	    if (!%pegs_for_singleton_fams)
	    {
		last;
	    }
	}
    }
    close(CALLS);
    my $fun = "hypothetical protein";
    my $sub;
    if (exists($subfam{$fun}))
    {
	$subfam{$fun}++;
	$sub = $subfam{$fun};
    }
    else
    {
	$sub = 1;
	$subfam{$fun} = 1;
    }
    for my $peg (sort keys %pegs_for_singleton_fams)
    {
	print B "$fun\t$sub\t$peg\n";
	$sub++;
    }
    close(B);

    inflate("$families.bad.fixed", "$families.bad.fixed.refs");

}
my $tend = gettimeofday;
my $elap = $tend - $tstart;
print LOG "finish inflate $tend $elap\n";

my $tstart = gettimeofday;
print LOG "start get_families_final $tstart\n";

#&SeedUtils::run("get_families_final -f $families -s $origSeqsD > $families.all");
&SeedUtils::ipc_run(["get_families_final", "-f", $families, "-s", $origSeqsD],
		    '>', "$families.all");

my $tend = gettimeofday;
my $elap = $tend - $tstart;
my $overall_elap = $tend - $overall_tstart;
print LOG "finish get_families_final $tend $elap\n";
print LOG "finish overall_run $tend $overall_elap\n";

kill 1, $kser_pid;
# kill 1, $guts_pid;

exit;

sub inflate
{
    my($file, $bak) = @_;
    rename($file, $bak);

    open(I, "<", $bak) or die "Cannot open $bak: $!";
    open(O, ">", $file) or die "Cannot write $file: $!";

    while (<I>)
    {
	print O $_;
	chomp;
	my($fam, $subfam, $peg) = split(/\t/);

	for my $exp (@{$pegs_to_inflate{$peg}})
	{
	    print O "$fam\t$subfam\t$exp\n";
	}
	delete $pegs_for_singleton_fams{$peg};
    }
    close(I);
    close(O);
}



#&SeedUtils::run("get_families_3 -c $cutoff < $tmpdir/calls > $families.good 2> $tmpdir/bad");
#&SeedUtils::run("get_families_4 -d $dataD -s $seqsD -m $matchN < $tmpdir/bad > $families.bad.fixed");

#unlink("$tmpdir/tmp.$$.missed","$tmpdir/tmp.$$.calls","$tmpdir/tmp.$$.bad");
#system("rm", "-r", $tmpdir);

