# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.160",
	package_date => 1508872672,
	package_date_str => "Oct 24, 2017 14:17:52",
    };
    return bless $self, $class;
}
1;
