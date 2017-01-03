# This is a SAS component.
package myRASTVersion;
use base 'Class::Accessor';
__PACKAGE__->mk_accessors(qw(release));
sub new
{
    my($class) = @_;
    my $self = {
	release => "1.148",
	package_date => 1483470218,
	package_date_str => "Jan 03, 2017 13:03:38",
    };
    return bless $self, $class;
}
1;
