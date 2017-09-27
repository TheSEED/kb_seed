# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.156",
	package_date => 1506541023,
	package_date_str => "Sep 27, 2017 14:37:03",
    };
    return bless $self, $class;
}
1;
